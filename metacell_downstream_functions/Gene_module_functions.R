require("WGCNA")
require("flashClust")
require("genefilter")
require("zoo")
require("viridis")
require("tgstat")
require("pheatmap")
#disableWGCNAThreads()
options(stringsAsFactors = FALSE)
require("topGO")
require("RColorBrewer")
require("scales")
require("matrixStats")
require("flashClust")

### Example typical usage:
# var_genes = names(which(apply(mc@mc_fp, 1, max)>1.8))
# gmod_determineSoftPowerWGCNA(mc@mc_fp[var_genes, ])   #Look at the plot and decide power
# x = gmod_runWGCNA(mc@mc_fp[var_genes, ], propGenes = 1, softPower = 8, cor_method = "pearson", signedNetwork = TRUE)
# gmod_plotModulesCut(x, "WGCNA_mods_pearson_signed_8")
# me = gmod_calculateModuleEigengenes(x, split = 4)
# gmods = gmod_moduleMembership(x, me, kME_threshold = 0.5)
# gmod_plotMCheatmap_annotate_modules(mc@mc_fp, gmods, me, tf_file = "../Annot_gene_lists/Sros_TF_list", gene_annot_file = "../Annot_gene_lists/Sros_FINAL_ANNOTATION")
# gmod_moduleGeneOntology(gmods, n_cpu = 10, p.thrs = 0.01, p.adj = FALSE)

#' @param data data.frame of normalised gene expression values with samples as columns and rows as genes.
#' @param output_file name of output file
#' @param propGenes integer, default is 1
#' @return power estimate
#' @description A wrapper function for constructing a WGCNA network
#'
gmod_determineSoftPowerWGCNA = function(
	data,
	propGenes = 1,
	output_file = NULL, width = 8, height = 4, res = NA
) {
	
	options(stringsAsFactors = FALSE)
	
	# Remove bad genes (missing values, ...)
	nGenesInput = dim(data)[1]
	data = data[ WGCNA::goodSamplesGenes(datExpr = t(data))$goodGenes, ]
	nGenesGood = dim(data)[1]
	nGenesRemoved = nGenesInput - nGenesGood
	message(paste(nGenesRemoved, " genes filtered from dataset"))
	propGenes = round(propGenes * dim(data)[1])
	
	# Filter genes based in variance
	keepGenesExpr1 = rank( - matrixStats::rowVars(data) ) <= propGenes
	data = data[keepGenesExpr1, ]
	genesRetained = dim(data)[1]
	message(paste(genesRetained, " genes retained in dataset"))
	
	# plot powers to work out soft-power
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
	
	# Call the network topology analysis function
	sft = WGCNA::pickSoftThreshold(t(data), powerVector = powers, verbose = 5)
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		# side by side
		par(mfrow = c(1, 2))
		
		# Scale-free topology fit index as a function of the soft-thresholding power
		plot(
			sft$fitIndices[, 1],
			-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
			xlab = "Soft Threshold (power)",
			ylab = "Scale Free Topology Model Fit, signed R^2",
			col = "gray",
			main = paste("Scale independence"))
		
		text(
			sft$fitIndices[, 1],
			-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
			labels = powers, cex = 0.8, col = "blue")
		
		# this line corresponds to using an R^2 cut-off of h
		abline(h = 0.80, col = "red", lty = 2)
		
		# Mean connectivity as a function of the soft-thresholding power
		plot(
			sft$fitIndices[, 1], sft$fitIndices[, 5],
			xlab = "Soft Threshold (power)",
			ylab = "Mean Connectivity",
			col = "gray",
			main = paste("Mean connectivity"))
		
		text(
			sft$fitIndices[, 1], sft$fitIndices[, 5],
			labels = powers, cex = 0.8, col = "blue")
		
	})
	return(sft$powerEstimate)
}

#' @param data data.frame of normalised gene expression values with samples as columns and rows as genes.
#' @param propGenes A numeric value indicating the proportion of most variable genes to be retained in the analysis between 0 and 1.
#' @param softPower Integer. The soft thresholding power used to construct the WGCNA network. This can be determined using the determineSoftPowerWGCNA function.
#' @param cor_method Character. One of "spearman", "pearson", or "bicor", defines which correlation metric to calculate.
#' @param hclust_method Character. Passed to `?flashClust()`, defines which clustering method to use.
#' @param signedNetwork Logical. Keep track of the sign of the correlation in network construction?
#' @return A list of 3. data is a data.frame of the input data after filtering. geneTreeA1 is the gene tree constructed by WGCNA. dissTOMA1 is the distance matirx.
#' @description A wrapper function for constructing a WGCNA network
gmod_runWGCNA = function(
    data, 
    propGenes = 1, 
    softPower = 10, 
    cor_method = "pearson", 
    hclust_method = "average",
    signedNetwork = TRUE
) {
	
	# Try to always use unsigned network and average hclust for gene modules.
	options(stringsAsFactors = FALSE)
	
	type = ifelse(test = signedNetwork == TRUE, yes = "signed", no = "unsigned")
	# Remove bad genes (missing values, ...)
	nGenesInput = dim(data)[1]
	data = data[goodSamplesGenes(datExpr = t(data))$goodGenes, ]
	nGenesGood = dim(data)[1]
	nGenesRemoved = nGenesInput - nGenesGood
	message(paste(nGenesRemoved, " genes filtered from dataset"))
	propGenes = round(propGenes * dim(data)[1])
	
	# Filter genes based in variance
	keepGenesExpr1 = rank( -matrixStats::rowVars(data) ) <= propGenes
	data = data[keepGenesExpr1, ]
	genesRetained = dim(data)[1]
	message(paste(genesRetained, " genes retained in dataset"))
	
	# Run WGCNA on the datasets
	datExprA1g = data
	if (cor_method == "pearson") {
		message("Using pearson correlation...")
		adjacencyA1 = adjacency(t(datExprA1g), power = softPower, type = type)
	} else if (cor_method == "spearman") {
		message("Using spearman correlation...")
		adjacencyA1 = adjacency(t(datExprA1g), power = softPower, type = type, corOptions = list(use = "p", method = "spearman"))
	} else if (cor_method == "bicor") {
		message("Using weighted bicorrelation...")
		adjacencyA1 = adjacency(t(datExprA1g), power = softPower, type = type, corFnc = "bicor", corOptions = list(maxPOutliers = 0.5))
	}else{
		message("Unkown correlation function!")
		break
	}
	
	diag(adjacencyA1)=0
	dissTOMA1 = 1 - TOMsimilarity(adjacencyA1, TOMType = type)
	
	geneTreeA1 = flashClust::flashClust(as.dist(dissTOMA1), method = hclust_method)
	
	# Return the relevant objects
	return(list(data = data, geneTreeA1 = geneTreeA1, dissTOMA1 = dissTOMA1))
}

#' Plot dendogram cuts
#' @param minClusterSize Integer. Target minimum size of a gene module, see `?cutreeHybrid`. Default 10.
#' @param cutHeight Numeric. Cut dendogram at this height, see `?cutreeHybrid`. Default 0.99.
#' @param maxds Numeric. Maximum dynamic splits to plot. Default 4.
#' @param output_file
#' @param width,height,res 
gmod_plotModulesCut = function(
	referenceDataset, minClusterSize = 10, cutHeight = 0.99, maxds = 4,
	output_file = NULL, width = 12, height = 6, res = NA
) {
	
	# Seperate list objects
	mColorh = NULL
	for (ds in 0:maxds) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE, minClusterSize = minClusterSize,
			cutHeight = cutHeight, deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh = cbind(mColorh, labels2colors(tree$labels))
	}
	
	for (ds in 0:maxds) {
		tree = cutreeDynamic(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE, minClusterSize = minClusterSize,
			cutHeight = cutHeight, deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh = cbind(mColorh, labels2colors(tree))
	}
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		WGCNA::plotDendroAndColors(
			referenceDataset[["geneTreeA1"]], mColorh,
			c(paste("Hybrid_dpSplt =", 0:maxds), paste("dynamic_dpSplt =", 0:maxds)), main = "",
			dendroLabels = FALSE)
			
	})
}

gmod_plotTOMheatmap = function(
	referenceDataset,
	minClusterSize = minClusterSize, splitdepth = c(2, 3),
	output_file = NULL, width = 6, height = 6, res = NA
) {
	
	mColorh = list()
	for (ds in splitdepth) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE,
			minClusterSize = minClusterSize,
			cutHeight = cutHeight,
			deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh[[as.character(ds)]]=labels2colors(tree$labels)
	}
	
	# plot
	plotting_function(output_file, width, height, res, EXP = {
		
		dissim = referenceDataset[["dissTOMA1"]]
		diag(dissim)=NA
		TOMplot(dissim = dissim, dendro = referenceDataset[["geneTreeA1"]], Colors = as.character(unlist(mColorh[[1]])), ColorsLeft = as.character(unlist(mColorh[[2]])))
		
	})
}


gmod_calculateModuleEigengenes = function(referenceDataset, split, minClusterSize = 10, cutHeight = 0.99) {
	#split is to be selected from examining plotModulesCut output. TBA: why default I use hybrid not dynamic?
	# Seperate list objects
	mColorh = NULL
	for (ds in 0:4) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE,
			minClusterSize = minClusterSize,
			cutHeight = cutHeight,
			deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]])
		
		mColorh = cbind(mColorh, labels2colors(tree$labels))
	}
	modulesA1 = mColorh[ , split] # (Chosen based on plot)
	PCs = moduleEigengenes(t(referenceDataset$data), colors = modulesA1)
	ME = PCs$eigengenes
	rownames(ME)=colnames(referenceDataset$data)
	
	return(ME)
}

gmod_moduleHubGenes = function(referenceDataset, MEs, nGenes, split = 1) {
	
	message("ranking genes")
	kMEs = WGCNA::signedKME(datExpr = t(referenceDataset$data), datME = MEs)
	
	# rank the genes for each module on kMEs
	rankGenes= function(x) {
		kMErank = rank(-kMEs[ , x])
		genes = rownames(kMEs)
		genes = genes[order(kMErank)]
		genes[1:nGenes]
	}
	
	topGenes = lapply(1:ncol(kMEs), rankGenes)
	
	# Get the top results in a data.frame
	topGenes = do.call(cbind, topGenes)
	colnames(topGenes)=substr(colnames(kMEs), start = 4, stop = 30)
	return(topGenes)
}

gmod_moduleMembership = function(referenceDataset, MEs, kME_threshold = 0.7) {
	
	# KEY FUNCTION. Increase threshold for tighter, conservative modules.
	message("calculating kME's")
	
	kMEs = WGCNA::signedKME(datExpr = t(referenceDataset$data), datME = MEs)
	# kMEs <<- kMEs
	
	# Get gene names with kME > kME_threshold
	modGenes= function(x) {
		rownames(kMEs[kMEs[ , x] > kME_threshold, ])  # IF you want unsigned, use abs(kMEs[, x]) thresholding
	}
	gmods = lapply(X = 1:ncol(kMEs), FUN = modGenes)
	names(gmods)=substr(colnames(kMEs), start = 4, stop = 20)
	
	gene_module_membership = apply(kMEs, 1, max)
	gmods = lapply(gmods, function(x) x[order(gene_module_membership[x], decreasing = TRUE)])
	return(gmods)
	
}

gmod_moduleNonoverlapingMembership = function(referenceDataset, MEs, kME_threshold = 0.7) {
	
	# Force genes to be included ONLY in 1 module (the one with best correlation)
	message("calculating kME's")
	
	kMEs = WGCNA::signedKME(datExpr = t(referenceDataset$data), datME = MEs)
	# kMEs <<- kMEs
	
	g_to_keep = rownames(kMEs)[which(apply(kMEs, 1, max) > kME_threshold)]  # we you want unsigned, use abs(kMEs) thresholding
	kMEs = kMEs[g_to_keep, ]
	g_to_gmod = apply(kMEs, 1, which.max)  # we you want unsigned, use abs(kMEs) thresholding
	gmods = split(names(g_to_gmod), f = g_to_gmod)
	names(gmods)=substr(colnames(kMEs), start = 4, stop = 20)
	
	gene_module_membership = apply(kMEs, 1, max)
	gmods = lapply(gmods, function(x) x[order(gene_module_membership[x], decreasing = TRUE)])
	return(gmods)
	
}


gmod_moduleOverlapingMembership = function(referenceDataset, MEs, kME_threshold = 0.7) {
	
	mc_wgcna_kmes = WGCNA::signedKME(datExpr = t(referenceDataset$data), datME = MEs)
	mc_wgcna_kmes_top = apply(mc_wgcna_kmes, 2, function(v) names(which(v > kME_threshold)))
	names(mc_wgcna_kmes_top) = gsub("^kME", "", names(mc_wgcna_kmes_top))
	return(mc_wgcna_kmes_top)
	
}

gmod_moduleGeneOntology = function(gmods, go_annotation = "../Annot_gene_lists/Sros_Gene_Ontology", n_cpu = 32, p.thrs = 0.05, min_genes_in_GO = 5, p.adj = TRUE, out_folder = "./module_GO_analysis/") {
	print("I AM SLOW AS HELL, DAMN YOU TopGO!!")
	thread_num = n_cpu
	doMC::registerDoMC(thread_num)
	
	go = readMappings(go_annotation)
	dir.create(out_folder, showWarnings = FALSE)
	
	go_analysis= function(i) {
		
		fg_genes = as.factor(as.integer(names(go) %in% gmods[[i]]))
		names(fg_genes)=names(go)
		
		GOdataBP = new("topGOdata", ontology = c("BP"), allGenes = fg_genes, annot = annFUN.gene2GO, gene2GO = go)
		fisher = runTest(GOdataBP, "classic", "fisher")
		parent_child = runTest(GOdataBP, "parentchild", "fisher")
		#suppressMessages(ks = runTest(GOdata, "classic", "ks"))
		r1 = GenTable(GOdataBP, fisher = fisher, parent_child = parent_child, topNodes= length(usedGO(GOdataBP)), orderBy = "fisher")
		r1 = r1[which(r1$Annotated > min_genes_in_GO), ]
		r1$fisher = as.numeric(r1$fisher)
		r1$parent_child = as.numeric(r1$parent_child)  #retarded TopGO reports values as characters...
		if (p.adj) {
			r1$fisher = p.adjust(r1$fisher, method = "BH")
			r1$parent_child = p.adjust(r1$parent_child, method = "BH")
		}
		r1[is.na(r1)]=1e-30 # TopGO doesn't report smaller evals
		r1 = r1[which(r1$parent_child < p.thrs | r1$fisher < p.thrs), ]
		r1 = cbind.data.frame(r1[, 1:2], rep("BP", nrow(r1)), r1[, 3:ncol(r1)])
		colnames(r1)[3]="GO_cat"
		
		GOdataMF = new("topGOdata", ontology = c("MF"), allGenes = fg_genes, annot = annFUN.gene2GO, gene2GO = go)
		fisher = runTest(GOdataMF, "classic", "fisher")
		parent_child = runTest(GOdataMF, "parentchild", "fisher")
		#suppressMessages(ks = runTest(GOdata, "classic", "ks"))
		r2 = GenTable(GOdataMF, fisher = fisher, parent_child = parent_child, topNodes= length(usedGO(GOdataMF)), orderBy = "fisher")
		r2 = r2[which(r2$Annotated > min_genes_in_GO), ]
		r2$fisher = as.numeric(r2$fisher)
		r2$parent_child = as.numeric(r2$parent_child)
		if (p.adj) {
			r2$fisher = p.adjust(r2$fisher, method = "BH")
			#r2$parent_child = p.adjust(r2$parent_child, method = "BH")
		}
		r2[is.na(r2)]=1e-30 # TopGO doesn't report smaller evals.
		r2 = r2[which(r2$parent_child < p.thrs | r2$fisher < p.thrs), ]
		r2 = cbind.data.frame(r2[, 1:2], rep("MF", nrow(r2)), r2[, 3:ncol(r2)])
		colnames(r2)[3]="GO_cat"
		
		r_final = rbind.data.frame(r1, r2)
		r_final = r_final[order(r_final$fisher), ]
		r_final$fisher = formatC(r_final$fisher, format = "e", digits = 2) #let's reports few dec positions
		r_final$parent_child = formatC(r_final$parent_child, format = "e", digits = 2)
		write.table(r_final, file = paste0(out_folder, "module_", i, "_GO.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	}
	
	plyr::alply(names(gmods), 1, function(x) go_analysis(x), .parallel = TRUE)
	
}

#' function to visualize WGCNA modules on metacells and to export annotated gene module table
#' 
#' @param expr_matrix expression matrix (e.g. `mc@mc_fp`)
#' @param gmods expression matrix
#' @param me eigenvalue matrix
#' @param expr_matrix_colors vector of metacell colors (same order as mc@mc_fp)
#' @param an_output_file,me_output_file,ex_output_file output files for the gene module annotations, the eigen values heatmap, and the expression heatmap
#' @param me_width,me_height size of the eigenvalues heatmap
#' @param do_expression boolean, whether to output the huge expression map (it's usually unadvisable, it's too big)
#' @param ex_width,ex_height size of the expression heatmap
#' @param resolution_rate downsample matrix to this rate, by average matrix subsampling
#' @param highlight_genes vector of genes to highlight in the expression heatmap
#' @param highlight_genes_annot vector of annotations for the genes to highlight in the expression heatmap (same order)
#' @param cor_cutoff_min,cor_cutoff_max low and high values for the correlation heatmap (default: 0 and 0.5). NULL sets max to quantile 0.99
#' @param eigen_min,eigen_ax low and high values for the eigenvalue heatmap (default: -0.1 to 0.2). NULL sets max to quantile 0.98
#' @param heatmap_colors,heatmap_colors_cor vector of colors to map to the heatmaps
#' 
#' @return nothing, only output plots
#' 
gmod_plotMCheatmap_annotate_modules = function(
	expr_matrix,
	gmods,
	me,
	expr_matrix_colors = NULL,
	an_output_file = "wgcna.gmod_annotation.csv", 
	me_output_file = "wgcna.gmod_eigenvalues.pdf", 
	ex_output_file = "wgcna.gmod_expression.pdf", 
	do_expression = TRUE,
	me_width = 10,  me_height = 5,
	ex_width = 20, ex_height = 10,
	res = NA,
	resolution_rate = 0.1,
	annotation = NULL,
	highlight_genes = NULL, 
	highlight_genes_annot = NULL,
	cor_cutoff_min = 0, 
	cor_cutoff_max = 0.5,
	eigen_min = -0.1,
	eigen_max = 0.2,
	heatmap_colors = c("white","orange","orangered2","#520c52"),
	heatmap_colors_cor = c("white","#d6e72e","#6fb600","#003f4d")
	
) {

	# me derive from calculateModuleEigengenes
	message("gmod plot | define input...")
	me = t(me) 
	
	# plot eigenvalues and reorder gmods based on them
	me = me[order(apply(me, 1, function(x) which.max(rollmean(x, 1)))), ]
	
	# do we have colors for the metacells?
	if (!is.null(expr_matrix_colors)) {
		mc_colors = expr_matrix_colors
	} else {
		mc_colors = FALSE
	}
	
	# plot
	# browser()
	message("gmod plot | ME plot...")
	plotting_function(me_output_file, width = me_width, height = me_height, res, EXP = {
		
		# color vector
		me_plot = me [ rev(rownames(me)) , ]
		rownames(me_plot) = gsub("^ME", "", rownames(me_plot))
		gmod_color_me = rownames(me_plot)
		names(gmod_color_me) = rownames(me_plot)
		
		# min and max values for the eigen heatmap
		if (is.null(eigen_max)) {
			eigen_max = quantile(me_plot, 0.98)
		}
		if (is.null(eigen_min)) {
			eigen_min = min(me_plot)
		}
		
		# plot
		hm = gmod_plot_complex_heatmap(
			mat = me_plot,
			name = "eigenvalue",
			color_mat = heatmap_colors_cor,
			color_min = eigen_min, 
			color_max = eigen_max,
			cluster_row = FALSE, 
			cluster_col = FALSE, 
			colors_row = gmod_color_me,
			colors_col = mc_colors,
			fontsize = 5, 
			use_raster = FALSE,
			title_row = "gene modules")
			
		print(hm)
		
	})
	gmods = gmods[gsub("ME", "", rownames(me))]
	gmods = sapply(gmods, function(x) x[x %in% rownames(expr_matrix)], USE.NAMES = TRUE, simplify = FALSE)
	unlist_gmods = unlist(gmods)
	
	# sort genes inside gmod by membership score?
	message("gmod plot | ME membership scores...")
	gene_module_membership = apply(WGCNA::signedKME(t(expr_matrix[unlist_gmods, ]), t(me)), 1, max)
	gmods = lapply(gmods, function(x) x[order(gene_module_membership[x], decreasing = TRUE)])

	# save gene module annotation table	
	message("gmod plot | gene module annotations...")
	tab = data.frame(
		gene = unlist_gmods, 
		gene_module = rep(names(gmods), lengths(gmods)), 
		membership_score = round(gene_module_membership[unlist_gmods], 3)
	)
	# add gene annotations, if available
	if (!is.null(annotation)) {
		tab = cbind(tab, annotation[unlist(gmods), ])
		colnames(tab) [ ( ncol(tab) - ncol(annotation) + 1 ) : ncol(tab) ] = colnames(annotation)
	}
	write.table(tab, file = an_output_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	
	# compute pearson correlation matrix (^ 5) and scale it
	if (do_expression) {
	
		message("gmod plot | expression/correlation plots, calculate correlation...")
		x_cor = tgs_cor(t(expr_matrix[unlist(gmods), ]), spearman = FALSE)
		x_cor[is.na(x_cor)]=0
		diag(x_cor)=0
		x_cor = (( 1 + x_cor ) / 2) ^ 5
	
		message("gmod plot | expression/correlation plots, set up plots...")
		plotting_function(ex_output_file, width = ex_width, height = ex_height, res, EXP = {
			
			# subsample x_cor matrix
			if (resolution_rate < 1) {
				x_cor_s = matrix_subsample(x_cor, subsample_rate = resolution_rate)
			} else {
				x_cor_s = x_cor
			}
			# rownames(x_cor_s) = 1:nrow(x_cor_s)
			# colnames(x_cor_s) = 1:ncol(x_cor_s)
			
			# subsample gmod
			g_cor_v = lengths(gmods)
			g_cor_s = g_cor_v * resolution_rate
			g_cor_s = round_smart(g_cor_s)
			gmod_color = names(g_cor_s)
			names(gmod_color) = names(g_cor_s)
			gmod_gene_color_category = rep(names(g_cor_s), g_cor_s)
			names(gmod_gene_color_category) = unlist(gmods)
			
			# colnames are names of gmods (located at the start position of each module)
			v = as.numeric(as.factor(gmod_color))
			ixs_names = c(1, 1 + which(diff( v ) != 0))
			ixs_other = c(1 + which(diff( v ) == 0))
			ixs_names = pmin(ixs_names, dim(x_cor_s)[1])
			ixs_names = pmax(ixs_names, 1)
			ixs_other = pmin(ixs_other, dim(x_cor_s)[1])
			ixs_other = pmax(ixs_other, 1)
			colnames(x_cor_s) [ ixs_names ] = names(g_cor_s)
			colnames(x_cor_s) [ ixs_other ] = ""
			
			# min and max values for the correlation map
			if (is.null(cor_cutoff_max)) {
				cor_cutoff_max = quantile(x_cor_s, 0.99)
			}
			if (is.null(cor_cutoff_min)) {
				cor_cutoff_min = 0
			}
			
			# first heatmap: gene-gene correlation
			message("gmod plot | expression/correlation plots, gene-gene correlation...")
			hm1 = gmod_plot_complex_heatmap(
				mat = x_cor_s,
				name = "cor^5",
				color_mat = heatmap_colors_cor,
				color_min = cor_cutoff_min,
				color_max = cor_cutoff_max,
				cluster_row = FALSE,
				cluster_col = FALSE,
				name_row_show = FALSE,
				name_col_show = TRUE,
				colors_row = gmod_color,
				colors_col = gmod_color,
				fontsize = 5,
				use_raster = TRUE,
				raster_quality = 1,
				title_col = "genes sorted by gene module",
				title_row = sprintf("%i genes in %i modules", nrow(x_cor), length(gmods)))
			
			# second heatmap: footprints
			if (resolution_rate < 1) {
				expr_s = matrix_subsample(expr_matrix[rownames(x_cor), ], subsample_rate_row = resolution_rate, subsample_rate_col = 1)
			} else {
				expr_s = expr_matrix[rownames(x_cor), ]
			}
			rownames(expr_s) = rep("", nrow(expr_s))
			colnames(expr_s) = colnames(expr_matrix)
			
			# rownames are names of selected genes (if needed)
			if (!is.null(highlight_genes)) {
				ixs_tfs_o = which(rownames(x_cor) %in% highlight_genes)
				vec_tfs_o = rownames(x_cor) [ ixs_tfs_o ]
				ixs_tfs_s = round(ixs_tfs_o * resolution_rate)
				ixs_tfs_s = pmin(ixs_tfs_s, dim(expr_s)[1])
				ixs_tfs_s = pmax(ixs_tfs_s, 1)
				ixs_oth_s = which(!as.character(1:nrow(expr_s)) %in% as.character(ixs_tfs_s))
				ixs_oth_s = pmin(ixs_oth_s, dim(expr_s)[1])
				ixs_oth_s = pmax(ixs_oth_s, 1)
				rownames(expr_s) [ ixs_tfs_s ] = vec_tfs_o
				rownames(expr_s) [ is.na(rownames(expr_s)) ] = ""
				
				if (!is.null(highlight_genes_annot)) {
					names(highlight_genes_annot) = highlight_genes
					vec_ann_o = highlight_genes_annot [ rownames(expr_s) [ ixs_tfs_s ] ]
					rownames(expr_s) [ ixs_tfs_s ] = paste(vec_tfs_o, stringr::str_trunc(vec_ann_o, width = 40), sep = " | ")
				}
				
			}
			
			message("gmod plot | expression/correlation plots, gene expression...")
			hm2 = gmod_plot_complex_heatmap(
				mat = expr_s,
				name = "fp",
				color_mat = heatmap_colors,
				color_min = 0,
				color_max = 4,
				cluster_row = TRUE,
				cluster_col = FALSE,
				name_row_show = TRUE,
				name_col_show = TRUE,
				colors_row = gmod_color,
				colors_col = mc_colors,
				fontsize = 3.5,
				use_raster = TRUE,
				raster_quality = 1,
				title_row = "genes")
			
			# relative sizes of each heatmap
			hm1@matrix_param$width = unit(0.7, "npc")
			hm2@matrix_param$width = unit(0.5, "npc")

			message("gmod plot | expression/correlation plots, draw...")
			print(hm1 + hm2)

		})
	}
}


gmod_geneModuleColors = function(
	referenceDataset, 
	split,
	minClusterSize = 20) {
	# Seperate list objects
	mColorh = NULL
	for (ds in 0:3) {
		tree = cutreeHybrid(
			dendro = referenceDataset[["geneTreeA1"]],
			pamStage = FALSE, minClusterSize = minClusterSize,
			cutHeight = 0.99, deepSplit = ds,
			distM = referenceDataset[["dissTOMA1"]]
		)
		mColorh = cbind(mColorh, labels2colors(tree$labels))
	}
	
	modules = mColorh[ , split] # (Chosen based on plot)
	return(modules)
}

gmod_getModuleLabels = function(
	referenceDataset,
	split,
	minClusterSize = 20) {
	tree = cutreeHybrid(
		dendro = referenceDataset[["geneTreeA1"]],
		pamStage = FALSE, minClusterSize = minClusterSize,
		cutHeight = 0.99, deepSplit = split,
		distM = referenceDataset[["dissTOMA1"]]
	)
	return(tree)
}

#compares the preservation of WGCNA solution in another dataset.
gmod_preservationStatsWGCNA = function(
	referenceDataset, data2, split = 3, type = "unsigned",
	nPermutations = 100, minClusterSize = 20,
	greyName = "grey", qVal = FALSE
) {
	
	colors = geneModuleColors(referenceDataset, split, minClusterSize)
	# Remove bad genes (missing values, ...)
	data2 = data2[goodSamplesGenes(datExpr = t(data2))$goodGenes, ]
	
	geneModules = cbind(rownames(referenceDataset$data), colors)
	
	# Quantify module preservation
	data = referenceDataset[["data"]]
	multiExpr = list(A1 = list(data = t(data)), A2 = list(data = t(data2)))
	multiColor = list(A1 = colors)
	
	mp = modulePreservation(
		multiData = multiExpr, multiColor = multiColor,
		referenceNetworks = 1, verbose = 3,
		calculateQvalue = qVal,
		networkType = type,
		nPermutations = nPermutations,
		maxGoldModuleSize = 100,
		greyName = greyName
		#maxModuleSize = 400
	)
	
	stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
	stats = stats[order(-stats[, 2]), ]
	print(stats[ , 1:3])
	return(list(mp = mp, geneModules = geneModules))
}


# Convert module colors to module names (so thing look more professional)
gmod_convertColorToLabel = function(colors, n = 100, prefix = "mod") {
	
	# Get the WGCNA colors
	cols = standardColors(n)
	labels = 1:n
	labels = paste(prefix, labels, sep = "")
	colToLabel = data.frame(colors = cols, labels = labels)
	
	# Add grey as label 0
	grey = c("grey", paste(prefix, "0", sep = ""))
	colToLabel = rbind(grey, colToLabel)
	
	# Function to return corresponding label for a given color
	convertFunc= function(x) {
		colToLabel[colToLabel$colors == x, ]$labels
	}
	
	# Apply the function to a list of colors
	unlist(lapply(colors, convertFunc))
}

gmod_compare_module_solutions = function(mods1, mods2) {
	# input must be a vector of module membership in same order.
	f = overlapTable(mods1, mods2)
	return(f)
}


gmod_compare_gene_module_membership = function(list1, list2) {
	common = intersect(unlist(list1), unlist(list2))
	
	perc_list1 = round(length(common) / length(unlist(list1)), 3)
	perc_list2 = round(length(common) / length(unlist(list2)), 3)
	
	message("Fraction of list1 genes in common: ", perc_list1)
	message("Fraction of list2 genes in common: ", perc_list2)
	
	list1_red = lapply(list1, function(x) intersect(x, common))
	list2_red = lapply(list2, function(x) intersect(x, common))
	m1 = matrix(rep(0, length(common) * length(common)), ncol = length(common), nrow = length(common))
	colnames(m1)=rownames(m1)=unlist(list1_red)
	m2 = matrix(rep(0, length(common) * length(common)), ncol = length(common), nrow = length(common))
	colnames(m2)=rownames(m2)=unlist(list2_red)
	
	# sapply(names(list1_red), function(x) m1[unlist(list1_red[[x]]), unlist(list1_red[[x]])]=1)
	# sapply(names(list2_red), function(x) m2[unlist(list2_red[[x]]), unlist(list2_red[[x]])]=1)
	
	for (module in names(which(lengths(list1_red) > 0))) {
		m1[list1_red[[module]], list1_red[[module]]]=1
	}
	
	
	for (module in names(which(lengths(list2_red) > 0))) {
		m2[list2_red[[module]], list2_red[[module]]]=1
	}
	
	m1_focused = m1
	m2_reord = m2[unlist(list1_red), unlist(list1_red)]
	m1_focused[lower.tri(m1_focused)]=m2_reord[lower.tri(m2_reord)]
	
	m2_focused = m2
	m1_reord = m1[unlist(list2_red), unlist(list2_red)]
	m2_focused[upper.tri(m2_focused)]=m1_reord[lower.tri(m1_reord)]
	
	png("test.pdf", height = 4, width = 4)
	par(mar = c(0, 0, 0, 0))
	par(fig = c(0.1, 0.5, 0.1, 0.95))
	pheatmap(m1_focused, cluster_cols = FALSE, cluster_rows = FALSE)
	
	par(fig = c(0.55, 0.9, 0.1, 0.95))
	pheatmap(m2_focused, cluster_cols = FALSE, cluster_rows = FALSE)
	dev.off()
}






gmod_moduleTraitCorrelations = function(MEs, trait) {
	
	# Order traits based on clustering tree
	d = dist(x = t(trait))
	traitTree = hclust(d, method = "a")
	plot(traitTree, xlab = "", ylab = "", main = "", sub = "")
	traitOrder = traitTree$labels[traitTree$order]
	trait = trait[ , match(traitOrder, colnames(trait))]
	
	# Order modules based on clustering tree
	d = dist(x = t(MEs))
	modTree = hclust(d, method = "a")
	plot(modTree, xlab = "", ylab = "", main = "", sub = "")
	modOrder = modTree$labels[modTree$order]
	MEs = MEs[ , match(modOrder, colnames(MEs))]
	
	# Calculate module-trait correlations
	moduleTraitCor = cor(x = trait, y = MEs, use = "pairwise.complete")
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs))
	
	return(list(moduleTraitCor = moduleTraitCor, moduleTraitPvalue = moduleTraitPvalue))
}

gmod_plotModuleTraitCorrelations = function(
	eigenTraitCor, pThresh,
	output_file = NULL, width = 8, height = 4, res = NA
) {
	
	moduleTraitCor = eigenTraitCor[["moduleTraitCor"]]
	moduleTraitPvalue = eigenTraitCor[["moduleTraitPvalue"]]
	
	
	# *** Experimental: organse traits by clustering correlation values
	d = dist(x = moduleTraitCor)
	traitTree = hclust(d, method = "a")
	
	d = dist(x = t(moduleTraitCor))
	moduleTree = hclust(d, method = "a")
	
	moduleTraitCor = moduleTraitCor[traitTree$order, moduleTree$order]
	moduleTraitPvalue = moduleTraitPvalue[traitTree$order, moduleTree$order]
	
	# Put the p-values in a text matrix for plotting
	# textMatrix = paste( signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "" )
	textMatrix = sprintf("%.2f (%.2f)", moduleTraitCor, moduleTraitPvalue )
	
	# Change text matrix to only contain significant p-values
	textMatrix = signif(moduleTraitPvalue, 1)
	textMatrix[textMatrix >= pThresh]=""
	dim(textMatrix)=dim(moduleTraitCor)
	
	# Reorder the the plot by dendrogram of correlation values
	
	# Display the correlation values within a heatmap plot
	myColRamp = colorRampPalette(colors = c("#1F8A70", "#FFFFFF", "#FD7400"))
	
	# Graphical parameters
	par(mar = c(6, 8.5, 3, 3), cex.lab = 0.5)
	
	modSymbols = gmod_convertColorToLabel(substr(x = colnames(moduleTraitCor), 3, 30), prefix = "M")
	
	plotting_function(output_file, width, height, res, EXP = {
		
		labeledHeatmap(
			Matrix = moduleTraitCor,
			xColorLabels = TRUE,
			xLabels = colnames(moduleTraitCor),
			yLabels = rownames(moduleTraitCor),
			xSymbols = modSymbols,
			textMatrix = textMatrix,
			colors = myColRamp(50),
			main = paste("Gene Module - Clinical Variable Correlations"),
			cex.text = 0.5,
			cex.lab.y = 0.5,
			cex.lab.x = 0.7)
		
		dev.off()
	})
	
}



# Generic function to plot complex heatmaps
gmod_plot_complex_heatmap = function(
		mat,
		name = "heatmap",
		color_mat = c("white","#d6e72e","#6fb600","#003f4d"),
		color_min = 0,
		color_max = 1,
		fontsize = 10,
		categories_col = NULL,
		categories_row = NULL,
		separate_col = FALSE,
		separate_row = FALSE,
		colors_col = NULL,
		colors_row = NULL,
		title_row = NULL,
		title_col = NULL,
		name_row_show = TRUE,
		name_col_show = TRUE,
		cluster_row = TRUE,
		cluster_col = TRUE,
		use_raster = TRUE,
		raster_quality = 1,
		show_legend_row = FALSE,
		show_legend_col = FALSE,
		both_sides_row = TRUE,
		both_sides_col = TRUE,
		cell_border = gpar(col = NA, lwd = 1, lty = 1),
		heatmap_border = gpar(col = NA, lwd = 1, lty = 1),
		do_dotplot = FALSE,
		dot_size_mat = NULL,
		dot_size_min = NULL,
		dot_size_max = NULL,
		cex_dotplot = 0.02
) {
	
	require("ComplexHeatmap")
	require("circlize")
	ht_opt$message = FALSE
	
	# color function
	col_fun = circlize::colorRamp2(seq(color_min, color_max, length.out = length(color_mat)), color_mat)
	
	
	# # vector of clusters (row/col annotation)
	# cluster_vector = as.character(mot_merge_d$cluster)
	# # categorical colors for clusters (get recycled)
	# catcol_vec = rep(catcol_lis, length.out = length(unique(cluster_vector)))
	# names(catcol_vec) = unique(cluster_vector)
	# cluster_colors = catcol_vec [ cluster_vector ]
	# names(cluster_colors) = cluster_colors
	# # quantitative colorscale for heatmap
	# catcol_fun = circlize::colorRamp2(1:length(catcol_vec), catcol_vec)
	# names(cluster_colors) = cluster_vector
	
	if (is.null(title_row)) {
		title_row = sprintf("n = %i rows", nrow(mat))
	}
	if (is.null(title_col)) {
		title_col = sprintf("n = %i columns", ncol(mat))
	}
	
	# left row annotations
	if (is.null(categories_row)) {
		categories_row = rownames(mat)
	}
	if (is.null(colors_row)) {
		colors_row = rep(NA, nrow(mat))
		left_annotation = NULL
	} else {
		if (is.null(names(colors_row))) { names(colors_row) = categories_row }
		left_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_row, name = title_row, col = list(c = colors_row), which = "row", show_legend = show_legend_row)
	}
	
	# top col annotations
	if (is.null(categories_col)) {
		categories_col = colnames(mat)
	}
	if (is.null(colors_col)) {
		colors_col = rep(NA, ncol(mat))
		top_annotation = NULL
	} else {
		if (is.null(names(colors_col))) { names(colors_col) = categories_col }
		top_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_col, name = title_col, col = list(c = colors_col), which = "column", show_legend = show_legend_col)
	}
	
	# add annotations to both sides of the rows or columns?
	if (both_sides_col) {
		bottom_annotation = top_annotation
	} else {
		bottom_annotation = NULL
	}
	if (both_sides_row) {
		right_annotation = left_annotation
	} else {
		right_annotation = NULL
	}
	
	# split columns and rows?
	if (separate_row) {
		split_row = categories_row
	} else {
		split_row = NULL
	}
	if (separate_col) {
		split_col = categories_col
	} else {
		split_col = NULL
	}
	
	# should this be a dot plot?
	if (do_dotplot) {
		if (is.null(dot_size_mat)) {
			dot_size_mat = mat
		}
		if (!is.null(dot_size_max) & !is.null(dot_size_min)) {
			# dot_size_mat [ dot_size_mat > dot_size_max ] = dot_size_max
			# dot_size_min [ dot_size_mat < dot_size_min ] = dot_size_min
			rangesize = function(x) { (x - dot_size_min) / (dot_size_max - dot_size_min) }
		} else {
			rangesize = function(x) { (x - min(x)) / (max(x) - min(x)) }
		}
		cell_fun_dotplot = function(j, i, x, y, width, height, fill) {
			grid.circle(
				x = x, y = y, 
				r = sqrt(rangesize(dot_size_mat)[i, j]) * cex_dotplot, 
				gp = gpar(col = NA, fill = col_fun(mat[i, j]))
			)
		}
		cell_border = gpar(type = "none")
		# cell_border = gpar(col = "gray", lwd = 1, lty = 1, fill = NULL)
	} else {
		cell_fun_dotplot = NULL
	}
	
	
	# plot
	hm = ComplexHeatmap::Heatmap(
		mat,
		name = name,
		cell_fun = cell_fun_dotplot,
		use_raster = use_raster,
		raster_quality = raster_quality,
		cluster_rows = cluster_row,
		cluster_columns = cluster_col,
		row_title = title_row,
		column_title = title_col,
		show_row_names = name_row_show,
		show_column_names = name_col_show,
		row_names_gp = gpar(fontsize = fontsize),
		column_names_gp = gpar(fontsize = fontsize),
		top_annotation = top_annotation,
		left_annotation = left_annotation,
		right_annotation = right_annotation,
		bottom_annotation = bottom_annotation,
		column_split = split_col,
		row_split = split_row,
		row_gap = unit(0.5, "mm"),
		column_gap = unit(0.5, "mm"),
		rect_gp = cell_border,
		border_gp = heatmap_border,
		col = col_fun)
	
	# return heatmap
	return(hm)
	
}

# Helper function to round values without changing total sum
round_smart <- function(x, digits = 0) {
	up <- 10 ^ digits
	x <- x * up
	y <- floor(x)
	indices <- tail(order(x - y), round(sum(x)) - sum(y))
	y[indices] <- y[indices] + 1
	y / up
}


# subsample matrix at a certain rate (to reduce resolution of a matrix in a brute-force way)
matrix_subsample = function(mat, subsample_rate = 0.1, subsample_rate_row = NULL, subsample_rate_col = NULL) {
	
	# subsample rates
	if (is.null(subsample_rate_row)) { 
		subsample_rate_row = subsample_rate
	}
	if (is.null(subsample_rate_col)) { 
		subsample_rate_col = subsample_rate
	}
	
	# num bins and widths
	nbin_r = nrow(mat) * subsample_rate_row
	nbin_c = ncol(mat) * subsample_rate_col
	wbin_r = round(nrow(mat) / nbin_r)
	wbin_c = round(ncol(mat) / nbin_c)
	
	# init subsampled matrix
	smat = matrix(nrow = nbin_r, ncol = nbin_c)
	
	# log
	message(sprintf("subsample matrix at %.2f x %.2f rate: %i x %i to %i x %i", subsample_rate_row, subsample_rate_col, nrow(mat), ncol(mat), nrow(smat), ncol(smat)))
	
	# loop
	for (rs in 0:(nbin_r - 1)) {
		ro_s = (rs * wbin_r) + 1
		ro_e = rs * wbin_r + wbin_r
		for (cs in 0:(nbin_c - 1)) {
			co_s = (cs * wbin_c) + 1
			co_e = cs * wbin_c + wbin_c
			smat_i = mat [ ro_s : ro_e , co_s : co_e ]
			smat[rs + 1, cs + 1]   = mean(smat_i)
		}
	}
	
	# out
	return(smat)
	
}


#' Map modules to metacells based on eigenvalue matrix
#' 
#' @param me eigenvalue matrix
#' @param ct_vector a named vector: names are metacells (corresponding to the rownames of the `me` matrix) and contents are cell types (default is NULL, i.e. the function simply maps the top mc to each module, but doesn't report cell types)
#' @param clean_module_names whether to remove the `ME` prefix in the module names (default is TRUE)
#' 
#' @return dataframe with module to mc and cell type annotations
#' 
gmod_annotate_modules_to_mc_and_ct = function(me, ct_vector = NULL, clean_module_names = TRUE, output_fn = NULL) {
	
	# clean module names? (remove ME prefix)
	if (clean_module_names) {
		colnames(me) = gsub("^ME","",colnames(me))
	}
	
	# find best mc per module
	top_mc_per_module = rownames(me) [ apply(me, 2, which.max) ]
	names(top_mc_per_module) = colnames(me)

	# if cell type table is given, map mcs to cell types
	if (!is.null(ct_vector)) {
		top_ct_per_module = ct_vector [ top_mc_per_module ]
	} else {
		top_ct_per_module = rep(NA, length(top_mc_per_module))
	}

	# output
	module_annots = data.frame(
		module = names(top_mc_per_module),
		mc_top = top_mc_per_module,
		ct_top = top_ct_per_module
	)
	
	if (!is.null(output_fn)) {
		write.table(module_annots, output_fn, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
	
	return(module_annots)
	
}