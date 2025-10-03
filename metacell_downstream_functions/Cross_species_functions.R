# libraries
require("metacell")
require("tgconfig")
require("tgstat")
require("tglkmeans")
require("zoo")
require("scales")
require("RColorBrewer")
require("viridis")
require("data.table")
require("stringr")
require("circlize")
require("ComplexHeatmap")
require("LaplacesDemon")
require("preprocessCore")
require("ggpubr")
require("doParallel")
require("philentropy")
require("biotools")
require("igraph")
require("WGCNA")
require("tgstat")


# deactivate complexheatmap warnings
ht_opt$message = FALSE


# # # # # # # # # # # # # # # # #
#                               #
#       UTILITY FUNCTIONS       #
#                               #
# # # # # # # # # # # # # # # # #

#' Overalap between all cells in two matrices
overlap <- function(M1, M2) {
	# make matrices conformable
	missing_in_M2 <- setdiff(rownames(M1), rownames(M2))
	missing_in_M1 <- setdiff(rownames(M2), rownames(M1))
	add_M2 <- matrix(0, nrow = length(missing_in_M2), ncol = ncol(M2))
	add_M1 <- matrix(0, nrow = length(missing_in_M1), ncol = ncol(M1))
	rownames(add_M2) <- missing_in_M2
	rownames(add_M1) <- missing_in_M1
	M1 <- data.matrix(rbind(M1, add_M1))
	M2 <- data.matrix(rbind(M2, add_M2))
	# intersect: matrix multiplication t(M1) x M2
	outM <- t(M1) %*% M2
	rownames(outM) <- colnames(M1)
	colnames(outM) <- colnames(M2)
	outM
}

#' Intersect normalized by smaller group size
intersect_by_min <- function(M1, M2) {
	# make matrices conformable
	missing_in_M2 <- setdiff(rownames(M1), rownames(M2))
	missing_in_M1 <- setdiff(rownames(M2), rownames(M1))
	add_M2 <- matrix(0, nrow = length(missing_in_M2), ncol = ncol(M2))
	add_M1 <- matrix(0, nrow = length(missing_in_M1), ncol = ncol(M1))
	rownames(add_M2) <- missing_in_M2
	rownames(add_M1) <- missing_in_M1
	M1 <- data.matrix(rbind(M1, add_M1))
	M2 <- data.matrix(rbind(M2, add_M2))
	# intersect: matrix multiplication t(M1) x M2
	intM <- t(M1) %*% M2
	rownames(intM) <- colnames(M1)
	colnames(intM) <- colnames(M2)
	
	modlen1 <- colSums(M1)
	modlen2 <- colSums(M2)
	minM <- intM
	for (i in 1:nrow(minM)) {
		for (j in 1:ncol(minM)) {
			minM[i,j] <- min(modlen1[i], modlen2[j])
		}
	}
	
	outM <- intM / minM
	return(outM)
}

#' Jaccard distance between columns of two binary matrices
#' @param M1,M2 binary matrix
jaccard <- function(M1, M2) {
	# make matrices conformable
	missing_in_M2 <- setdiff(rownames(M1), rownames(M2))
	missing_in_M1 <- setdiff(rownames(M2), rownames(M1))
	add_M2 <- matrix(0, nrow = length(missing_in_M2), ncol = ncol(M2))
	add_M1 <- matrix(0, nrow = length(missing_in_M1), ncol = ncol(M1))
	rownames(add_M2) <- missing_in_M2
	rownames(add_M1) <- missing_in_M1
	M1 <- data.matrix(rbind(M1, add_M1))
	M2 <- data.matrix(rbind(M2, add_M2))
	# intersect: matrix multiplication t(M1) x M2
	intersectM <- t(M1) %*% M2
	# union: sum by rows and columns
	unionM <- matrix(nrow = ncol(M1), ncol = ncol(M2))
	for (i in 1:ncol(M1)) {
		unionM[i,] <- colSums((M1[,i] + M2) > 0)
	}
	# jaccard
	outM <- intersectM / unionM
	outM[is.na(outM)] <- 0
	rownames(outM) <- colnames(M1)
	colnames(outM) <- colnames(M2)
	outM
}

#' Shafer similarity index between columns of two binary matrices
#' @param M1,M2 binary matrix
shaferindex <- function(M1, M2) {
	# make matrices conformable
	missing_in_M2 <- setdiff(rownames(M1), rownames(M2))
	missing_in_M1 <- setdiff(rownames(M2), rownames(M1))
	add_M2 <- matrix(0, nrow = length(missing_in_M2), ncol = ncol(M2))
	add_M1 <- matrix(0, nrow = length(missing_in_M1), ncol = ncol(M1))
	rownames(add_M2) <- missing_in_M2
	rownames(add_M1) <- missing_in_M1
	M1 <- data.matrix(rbind(M1, add_M1))
	M2 <- data.matrix(rbind(M2, add_M2))
	# intersect: matrix multiplication t(M1) x M2
	intersectM <- t(M1) %*% M2
	# union: sum by rows and columns
	unionM <- matrix(nrow = ncol(M1), ncol = ncol(M2))
	for (i in 1:ncol(M1)) {
		unionM[i,] <- colSums((M1[,i] + M2) > 0)
	}
	# shafer index
	outM <- 1 - sqrt( ( 1 - intersectM / nrow(M1) ) * ( 1 - intersectM / nrow(M2) ) ) 
	outM[is.na(outM)] <- 0
	rownames(outM) <- colnames(M1)
	colnames(outM) <- colnames(M2)
	outM
}

#' KLD between columns (metacells, cell types) of matrix (footprint, genes in rows)
calcKLD = function(mc_fp) {
	nc=ncol(mc_fp)
	kld_mat=matrix(NA,nrow=nc,ncol=nc)
	kld_list=vector("list",nc)
	for (i in 1:nc) {
		# message(sprintf("%s of %s",i,nc))
		for (j in 1:nc) {
			if (!i > j) {
				kld = LaplacesDemon::KLD(px=mc_fp[,i],py=mc_fp[,j])
				kld_mat[i,j]=kld$sum.KLD.px.py
				kld_mat[j,i]=kld$sum.KLD.py.px
				kld_list[[i]][[j]]=kld$KLD.px.py
				kld_list[[j]][[i]]=kld$KLD.py.px
			}
		}
	}
	colnames(kld_mat) = colnames(mc_fp)
	rownames(kld_mat) = colnames(mc_fp)
	list(mat=kld_mat,probs=kld_list)
}

#' Quantile normalization
#' @param x matrix or data.frame
quantile_normalisation <- function(x){
	df_rank <- data.frame(apply(x,2,rank,ties.method="min"))
	df_sorted <- data.frame(apply(x, 2, sort))
	df_mean <- apply(df_sorted, 1, mean)
	
	index_to_mean <- function(my_index, my_mean){
		return(my_mean[my_index])
	}
	
	df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
	rownames(df_final) <- rownames(x)
	return(df_final)
}

# # # # # # # # # # # # # # # # #
#                               #
#     CSPS INIT FUNCTIONS       #
# create csps object as output  #
#                               #
# # # # # # # # # # # # # # # # #

#' Create cross-species comparison object.
#' 
#' @param sp1_fp_fn,sp2_fp_fn path to expression matrices for the analyzed species, 
#'   rows are genes and columns are cells, metacells or cell types (usually these 
#'   are metacell footprint matrices, `mc@mc_fp`); this can be either raw fold changes 
#'   UMI counts or UMI fraction matrices in either RDS or tab-separated text format.
#'   Alternatively, you can provide a matrix object (like `mc@mc_fp`) as input too.
#' @param OG_pairs_fn path to file with broccoli-style pairwise orthologies in
#'   tab-separated text format; make sure they are in the right order
#' @param sp_names character of length 2, optional species names/abbreviations;
#'   if not specified, will use `c("sp1","sp2")`
#' @param make_sp_colnames logical, append the species abbreviations at the 
#'   begining of column names (default: TRUE)
#' @param quant_norm logical (default: TRUE)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` 
#'   we require in each species (default: 1, set to NULL to skip filtering by fc)
#' @param cross_fc_thrs numeric, fold change threshold (default: 2, set to NULL 
#'   to skip filtering by fc)
#' @param one2one logical (default: TRUE)
#'
#' @return a list with the following elements:
#'   * merged, a matrix with selected orthologous genes in rows and metacells 
#'     (or cell types) from both species in columns
#'   * og_pairs, data.frame with orthologous pairs in columns
#'   * sp1, a matrix with selected orthologous genes in rows and metacells 
#'     of first species in columns
#'   * sp2, a matrix with selected orthologous genes in rows and metacells 
#'     of the second species in columns
#'   * top_cross_sp1, top orthologous genes in the first species
#'   * top_cross_sp2, top orthologous genes in the second species
#' This will be a standard object for cross-species analyses
#' 
csps_create_crossspecies_object <- function(
	sp1_fp_fn, sp2_fp_fn, OG_pairs_fn, sp_names=c("sp1","sp2"), 
	make_sp_colnames=TRUE, quant_norm=TRUE, cross_n=1, cross_fc_thrs=2, one2one=TRUE
){ 
	
	# Read and parse data  
	# species 1
	extension = rev(stringr::str_split(sp1_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp1_fp_fn)) {
		sp1_fp = sp1_fp_fn
	} else if ("character" %in% class(sp1_fp_fn) & extension == ".RDS") {
		sp1_fp=readRDS(sp1_fp_fn)
	} else if ("character" %in% class(sp1_fp_fn)) {
		sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp1_fp_fn is neither a matrix object nor a path to such an object, I quit.")    
	}
	# species 2
	extension = rev(stringr::str_split(sp2_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp2_fp_fn)) {
		sp2_fp = sp2_fp_fn
	} else if ("character" %in% class(sp2_fp_fn) & extension == ".RDS") {
		sp2_fp=readRDS(sp2_fp_fn)
	} else if ("character" %in% class(sp2_fp_fn)) {
		sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp2_fp_fn is neither a matrix object nor a path to such an object, I quit.")    
	}
	OG_pairs=fread(OG_pairs_fn,header=FALSE,sep="\t",stringsAsFactors=FALSE)
	
	setnames(OG_pairs, c(sp_names))
	sp1=sp_names[1]
	sp2=sp_names[2]
	sp1_fp=sp1_fp[intersect(rownames(sp1_fp),OG_pairs[[sp1]]),]
	sp2_fp=sp2_fp[intersect(rownames(sp2_fp),OG_pairs[[sp2]]),]
	if (make_sp_colnames) {
		colnames(sp1_fp) = paste(sp1, colnames(sp1_fp), sep = "_")
		colnames(sp2_fp) = paste(sp2, colnames(sp2_fp), sep = "_")
	}
	
	# Reduce orthology table to only expressed genes
	OG_pairs=OG_pairs[which(OG_pairs[[sp1]] %in% rownames(sp1_fp) & OG_pairs[[sp2]] %in% rownames(sp2_fp)),]
	
	# Eliminate non-one2one orthologs (relaxed non-expression criterion)  
	OG_pairs_one2one=as.data.frame(
		OG_pairs[!duplicated(OG_pairs[[sp1]]) & !duplicated(OG_pairs[[sp2]]),]
	)
	
	if (one2one) {
		
		OG_pairs = OG_pairs_one2one
		cross1 <- OG_pairs[[sp1]]
		
	} else { 
		
		# Allow for one2many (2-3 max) relationships, DUPLICATING entries for the paralogs.
		# one2one
		one2one <- OG_pairs_one2one[[sp1]]
		# one2many (2 or three)
		OG_pairs[[sp1]] <- make.names(OG_pairs[[sp1]],unique=TRUE)
		one2many <- grep("\\.[123]$",OG_pairs[[sp1]],value=TRUE)
		one2many <- unique(sort(c(one2many, str_remove(one2many,"\\.[123]$"))))
		
		OG_pairs = OG_pairs[OG_pairs[[sp1]] %in% c(one2one,one2many),]
		cross1 <- str_remove(OG_pairs[[sp1]],"\\.[123]$")
		
	}
	
	# Reorder matrices, quantile_norm? and merge matrices
	if (class(sp1_fp)[1] == "dgCMatrix" | class(sp2_fp)[1] == "dgCMatrix") {
		merged=Matrix::Matrix(cbind(sp1_fp[cross1,],sp2_fp[OG_pairs[[sp2]],]),sparse=TRUE)
	} else {
		merged=data.matrix(cbind(sp1_fp[cross1,],sp2_fp[OG_pairs[[sp2]],]))
	}
	rnm=OG_pairs[[sp1]]
	rownames(merged)=rnm
	cnm=colnames(merged)
	
	if (quant_norm) {
		if (any(class(merged) == "dgCMatrix")) {
			sparse=TRUE
			merged=as.matrix(merged)
			copymat=FALSE
		} else {
			sparse=FALSE
			copymat=TRUE
		}
		merged=tryCatch({
			preprocessCore::normalize.quantiles(merged,copy=copymat)
		}, error = function(e){
			warning(e)
			quantile_normalisation(merged)
		})
		if (sparse)
			merged=Matrix::Matrix(merged,sparse=TRUE)
		rownames(merged)=rnm
		colnames(merged)=cnm
	}
	
	# Select variable genes in BOTH species
	if (!is.null(cross_fc_thrs) & !is.null(cross_n)) {
		top_cross_list = csps_select_covariable_genes(sp1_fp = sp1_fp, sp2_fp = sp2_fp, merged =  merged, cross_fc_thrs = cross_fc_thrs, cross_n = cross_n)
		top_cross_sp1 = top_cross_list[[1]]
		top_cross_sp2 = top_cross_list[[2]]
		# old
		# top_cross=names(which(
		# 	apply(merged[,1:ncol(sp1_fp)],1,function(x) 
		# 		sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs & 
		# 		apply(merged[,(ncol(sp1_fp) + 1):ncol(merged)],1,function(x) 
		# 			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		# ))
	} else {
		top_cross_sp1 = rownames(sp1_fp)
		top_cross_sp2 = rownames(sp2_fp)
	}
	
	# Return a list with both matrices
	out_m1=merged[ ,1:ncol(sp1_fp) ]
	rownames(out_m1)=OG_pairs[[sp1]]
	out_m2=merged[ ,(ncol(sp1_fp) + 1):ncol(merged) ]
	rownames(out_m2)=make.names(OG_pairs[[sp2]],unique=TRUE)
	
	# commented-out when adding `csps_select_covariable_genes`
	# ids <- unlist(lapply(top_cross,function(x) which(OG_pairs[[sp1]] == x)))
	# top_cross_sp2=OG_pairs[ids,][[sp2]]
	#top_cross <- str_remove(top_cross,"\\.\\d+")
	#top_cross_sp2 <- str_remove(top_cross2,"\\.\\d+")
	
	csps <- list(
		merged=merged,
		og_pairs=OG_pairs,
		sp1=out_m1, sp2=out_m2,
		top_cross_sp1 = top_cross_sp1, top_cross_sp2 = top_cross_sp2
	)
	return(csps)
}

#' @inheritParams csps_create_crossspecies_object
#' @seealso [csps_create_crossspecies_object()]
csps_create_3way_crossspecies_object <- function(
	sp1_fp_fn, sp2_fp_fn, sp3_fp_fn, 
	OG_pairs1_2_fn, OG_pairs1_3_fn, OG_pairs2_3_fn, 
	sp_names=c("sp1","sp2","sp3"), make_sp_colnames=TRUE,
	one2one=TRUE, quant_norm=TRUE, cross_n=1, cross_fc_thrs=2
){
	# Read and parse data  
	if (length(sp_names) != 3) {
		stop("Length of sp_names should be 3!")
	}
	# species 1
	extension = rev(stringr::str_split(sp1_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp1_fp_fn)) {
		sp1_fp = sp1_fp_fn
	} else if ("character" %in% class(sp1_fp_fn) & extension == ".RDS") {
		sp1_fp=readRDS(sp1_fp_fn)
	} else if ("character" %in% class(sp1_fp_fn)) {
		sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp1_fp_fn is neither a mc@mc_fp matrix object nor a path to such an object, I quit.")    
	}
	# species 2
	extension = rev(stringr::str_split(sp2_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp2_fp_fn)) {
		sp2_fp = sp2_fp_fn
	} else if ("character" %in% class(sp2_fp_fn) & extension == ".RDS") {
		sp2_fp=readRDS(sp2_fp_fn)
	} else if ("character" %in% class(sp2_fp_fn)) {
		sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp2_fp_fn is neither a mc@mc_fp matrix object nor a path to such an object, I quit.")    
	}
	# species 3
	extension = rev(stringr::str_split(sp3_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp3_fp_fn)) {
		sp3_fp = sp3_fp_fn
	} else if ("character" %in% class(sp3_fp_fn) & extension == ".RDS") {
		sp3_fp=readRDS(sp3_fp_fn)
	} else if ("character" %in% class(sp3_fp_fn)) {
		sp3_fp=read.table(sp3_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp3_fp_fn is neither a mc@mc_fp matrix object nor a path to such an object, I quit.")    
	}
	
	OG_pairs1_2=fread(OG_pairs1_2_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	OG_pairs1_3=fread(OG_pairs1_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	OG_pairs2_3=fread(OG_pairs2_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	sp1=sp_names[1]
	sp2=sp_names[2]
	sp3=sp_names[3]
	setnames(OG_pairs1_2, c(sp1,sp2))
	setnames(OG_pairs1_3, c(sp1,sp3))
	setnames(OG_pairs2_3, c(sp2,sp3))
	intersect1=intersect(intersect(rownames(sp1_fp),OG_pairs1_2[[sp1]]),OG_pairs1_3[[sp1]])
	intersect2=intersect(intersect(rownames(sp2_fp),OG_pairs1_2[[sp2]]),OG_pairs2_3[[sp2]])
	intersect3=intersect(intersect(rownames(sp3_fp),OG_pairs1_3[[sp3]]),OG_pairs2_3[[sp3]])
	sp1_fp=sp1_fp[intersect1,]
	sp2_fp=sp2_fp[intersect2,]
	sp3_fp=sp3_fp[intersect3,]
	if (make_sp_colnames) {
		colnames(sp1_fp) = paste(sp1, colnames(sp1_fp), sep = "_")
		colnames(sp2_fp) = paste(sp2, colnames(sp2_fp), sep = "_")
		colnames(sp3_fp) = paste(sp3, colnames(sp3_fp), sep = "_")
	}
	
	# Reduce orthology table to only expressed genes
	OG_pairs1_2=OG_pairs1_2[which(OG_pairs1_2[[sp1]] %in% rownames(sp1_fp) & 
								  	OG_pairs1_2[[sp2]] %in% rownames(sp2_fp)),]
	OG_pairs1_3=OG_pairs1_3[which(OG_pairs1_3[[sp1]] %in% rownames(sp1_fp) & 
								  	OG_pairs1_3[[sp3]] %in% rownames(sp3_fp)),]
	OG_pairs2_3=OG_pairs2_3[which(OG_pairs2_3[[sp2]] %in% rownames(sp2_fp) & 
								  	OG_pairs2_3[[sp3]] %in% rownames(sp3_fp)),]
	
	# Eliminate non-one2one orthologs (relaxed non-expression criterion)  
	OG_pairs1_2_one2one=as.data.frame(
		OG_pairs1_2[!duplicated(OG_pairs1_2[[sp1]]) & !duplicated(OG_pairs1_2[[sp2]]),]
	)
	OG_pairs1_3_one2one=as.data.frame(
		OG_pairs1_3[!duplicated(OG_pairs1_3[[sp1]]) & !duplicated(OG_pairs1_3[[sp3]]),]
	)
	OG_pairs2_3_one2one=as.data.frame(
		OG_pairs2_3[!duplicated(OG_pairs2_3[[sp2]]) & !duplicated(OG_pairs2_3[[sp3]]),]
	)
	
	if (one2one) {
		
		OG_pairs1_2 = OG_pairs1_2_one2one
		OG_pairs1_3 = OG_pairs1_3_one2one
		OG_pairs2_3 = OG_pairs2_3_one2one
		mdt1=merge.data.table(OG_pairs1_2,OG_pairs1_3,by=sp1)
		mdt=merge.data.table(mdt1,OG_pairs2_3,by=c(sp2,sp3))
		cross1=mdt[[sp1]]
		cross2=mdt[[sp2]]
		cross3=mdt[[sp3]]
		
	} else { 
		
		stop("one2many not implemented yet, set one2one=TRUE")
		
	}
	
	# Reorder matrices, quantile_norm and merge matrices
	if (any(c(class(sp1_fp),class(sp2_fp),class(sp3_fp)) == "dgCMatrix")) {
		merged1=Matrix::Matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]),sparse=TRUE)
		merged=Matrix::Matrix(cbind(merged1,sp3_fp[cross3,]),sparse=TRUE)
	} else {
		merged1=data.matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]))
		merged=data.matrix(cbind(merged1,sp3_fp[cross3,]))
	}
	rnm=cross1
	rownames(merged)=rnm
	cnm=colnames(merged)
	
	if (quant_norm){
		if (any(class(merged) == "dgCMatrix")) {
			sparse=TRUE
			merged=as.matrix(merged)
			copymat=FALSE
		} else {
			sparse=FALSE
			copymat=TRUE
		}
		merged=tryCatch({
			preprocessCore::normalize.quantiles(merged,copy=copymat)
		}, error = function(e){
			quantile_normalisation(merged)
		})
		if (sparse)
			merged=Matrix::Matrix(merged,sparse=TRUE)
		rownames(merged)=rnm
		colnames(merged)=cnm
	}
	
	# Select variable genes in ALL species
	nc1=ncol(sp1_fp)
	nc2=ncol(sp1_fp) + ncol(sp2_fp)
	nc3=ncol(merged)
	if (!is.null(cross_fc_thrs) & !is.null(cross_n)) {
		var1= apply(merged[,1:nc1],1,function(x) 
			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		var2= apply(merged[,c(nc1 + 1):nc2],1,function(x) 
			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		var3= apply(merged[,c(nc2 + 1):nc3],1,function(x) 
			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		top_cross=names(which(var1 & var2 & var3))
	} else {
		top_cross=rownames(merged)
	}
	
	# Return a list with both matrices
	out_m1=merged[,1:nc1]
	rownames(out_m1)=cross1
	out_m2=merged[,(nc1 + 1):nc2]
	rownames(out_m2)=make.names(cross2,unique=TRUE)
	out_m3=merged[,(nc2 + 1):nc3]
	rownames(out_m3)=make.names(cross3,unique=TRUE)
	
	ids2 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_2[[sp1]] == x)))
	top_cross_sp2=OG_pairs1_2[ids2,][[sp2]]
	ids3 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_3[[sp1]] == x)))
	top_cross_sp3=OG_pairs1_3[ids3,][[sp3]]
	
	csps <- list(
		merged=merged,
		OG_pairs1_2=OG_pairs1_2, OG_pairs1_3=OG_pairs1_3, OG_pairs2_3=OG_pairs2_3,
		sp1=out_m1, sp2=out_m2, sp3=out_m3,
		top_cross_sp1=top_cross, top_cross_sp2=top_cross_sp2, top_cross_sp3=top_cross_sp3
	)
	return(csps)
	
}

#' @inheritParams csps_create_crossspecies_object
#' @param restrict_paralogs logical, restrict number of possible paralogs for a gene 
#'   (default: FALSE, if TRUE, genes with more than 4 paralogs will not be considered)
#' @seealso [csps_create_crossspecies_object()]
csps_create_4way_crossspecies_object <- function(
	sp1_fp_fn, sp2_fp_fn, sp3_fp_fn, sp4_fp_fn,  
	OG_pairs1_2_fn, OG_pairs1_3_fn,OG_pairs1_4_fn, OG_pairs2_3_fn, OG_pairs2_4_fn, OG_pairs3_4_fn,
	sp_names=c("sp1","sp2","sp3","sp4"), make_sp_colnames=TRUE,
	one2one=TRUE, restrict_paralogs=FALSE, quant_norm=TRUE, 
	cross_n=1, cross_fc_thrs=2
){
	# Read and parse data  
	if (length(sp_names) != 4)
		stop("Length of sp_names should be 4!")
	# species 1
	extension = rev(stringr::str_split(sp1_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp1_fp_fn)) {
		sp1_fp = sp1_fp_fn
	} else if ("character" %in% class(sp1_fp_fn) & extension == ".RDS") {
		sp1_fp=readRDS(sp1_fp_fn)
	} else if ("character" %in% class(sp1_fp_fn)) {
		sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp1_fp_fn is neither a mc@mc_fp matrix object nor a path to such an object, I quit.")    
	}
	# species 2
	extension = rev(stringr::str_split(sp2_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp2_fp_fn)) {
		sp2_fp = sp2_fp_fn
	} else if ("character" %in% class(sp2_fp_fn) & extension == ".RDS") {
		sp2_fp=readRDS(sp2_fp_fn)
	} else if ("character" %in% class(sp2_fp_fn)) {
		sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp2_fp_fn is neither a mc@mc_fp matrix object nor a path to such an object, I quit.")    
	}
	# species 3
	extension = rev(stringr::str_split(sp3_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp3_fp_fn)) {
		sp3_fp = sp3_fp_fn
	} else if ("character" %in% class(sp3_fp_fn) & extension == ".RDS") {
		sp3_fp=readRDS(sp3_fp_fn)
	} else if ("character" %in% class(sp3_fp_fn)) {
		sp3_fp=read.table(sp3_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp3_fp_fn is neither a mc@mc_fp matrix object nor a path to such an object, I quit.")    
	}
	# species 4
	extension = rev(stringr::str_split(sp4_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp4_fp_fn)) {
		sp4_fp = sp4_fp_fn
	} else if ("character" %in% class(sp4_fp_fn) & extension == ".RDS") {
		sp4_fp=readRDS(sp4_fp_fn)
	} else if ("character" %in% class(sp4_fp_fn)) {
		sp4_fp=read.table(sp4_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp4_fp_fn is neither a mc@mc_fp matrix object nor a path to such an object, I quit.")    
	}
	
	# load orthology
	OG_pairs1_2=fread(OG_pairs1_2_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	OG_pairs1_3=fread(OG_pairs1_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	OG_pairs1_4=fread(OG_pairs1_4_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	OG_pairs2_3=fread(OG_pairs2_3_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	OG_pairs2_4=fread(OG_pairs2_4_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	OG_pairs3_4=fread(OG_pairs3_4_fn,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
	sp1=sp_names[1]
	sp2=sp_names[2]
	sp3=sp_names[3]
	sp4=sp_names[4]
	setnames(OG_pairs1_2, c(sp1,sp2))
	setnames(OG_pairs1_3, c(sp1,sp3))
	setnames(OG_pairs1_4, c(sp1,sp4))
	setnames(OG_pairs2_3, c(sp2,sp3))
	setnames(OG_pairs2_4, c(sp2,sp4))
	setnames(OG_pairs3_4, c(sp3,sp4))
	intersect1=intersect(intersect(intersect(rownames(sp1_fp),OG_pairs1_2[[sp1]]),OG_pairs1_3[[sp1]]),OG_pairs1_4[[sp1]])
	intersect2=intersect(intersect(intersect(rownames(sp2_fp),OG_pairs1_2[[sp2]]),OG_pairs2_3[[sp2]]),OG_pairs2_4[[sp2]])
	intersect3=intersect(intersect(intersect(rownames(sp3_fp),OG_pairs1_3[[sp3]]),OG_pairs2_3[[sp3]]),OG_pairs3_4[[sp3]])
	intersect4=intersect(intersect(intersect(rownames(sp4_fp),OG_pairs1_4[[sp4]]),OG_pairs2_4[[sp4]]),OG_pairs3_4[[sp4]])
	sp1_fp=sp1_fp[intersect1,]
	sp2_fp=sp2_fp[intersect2,]
	sp3_fp=sp3_fp[intersect3,]
	sp4_fp=sp4_fp[intersect4,]
	if (make_sp_colnames) {
		colnames(sp1_fp) = paste(sp1, colnames(sp1_fp), sep = "_")
		colnames(sp2_fp) = paste(sp2, colnames(sp2_fp), sep = "_")
		colnames(sp3_fp) = paste(sp3, colnames(sp3_fp), sep = "_")
		colnames(sp4_fp) = paste(sp4, colnames(sp4_fp), sep = "_")
	}
	
	# Reduce orthology table to only expressed genes
	OG_pairs1_2=OG_pairs1_2[which(OG_pairs1_2[[sp1]] %in% rownames(sp1_fp) & 
								  	OG_pairs1_2[[sp2]] %in% rownames(sp2_fp)),]
	OG_pairs1_3=OG_pairs1_3[which(OG_pairs1_3[[sp1]] %in% rownames(sp1_fp) & 
								  	OG_pairs1_3[[sp3]] %in% rownames(sp3_fp)),]
	OG_pairs1_4=OG_pairs1_4[which(OG_pairs1_4[[sp1]] %in% rownames(sp1_fp) & 
								  	OG_pairs1_4[[sp4]] %in% rownames(sp4_fp)),]
	OG_pairs2_3=OG_pairs2_3[which(OG_pairs2_3[[sp2]] %in% rownames(sp2_fp) & 
								  	OG_pairs2_3[[sp3]] %in% rownames(sp3_fp)),]
	OG_pairs2_4=OG_pairs2_4[which(OG_pairs2_4[[sp2]] %in% rownames(sp2_fp) & 
								  	OG_pairs2_4[[sp4]] %in% rownames(sp4_fp)),]
	OG_pairs3_4=OG_pairs3_4[which(OG_pairs3_4[[sp3]] %in% rownames(sp3_fp) & 
								  	OG_pairs3_4[[sp4]] %in% rownames(sp4_fp)),]
	
	# Eliminate non-one2one orthologs (relaxed non-expression criterion)  
	OG_pairs1_2_one2one=as.data.frame(
		OG_pairs1_2[!duplicated(OG_pairs1_2[[sp1]]) & !duplicated(OG_pairs1_2[[sp2]]),]
	)
	OG_pairs1_3_one2one=as.data.frame(
		OG_pairs1_3[!duplicated(OG_pairs1_3[[sp1]]) & !duplicated(OG_pairs1_3[[sp3]]),]
	)
	OG_pairs1_4_one2one=as.data.frame(
		OG_pairs1_4[!duplicated(OG_pairs1_4[[sp1]]) & !duplicated(OG_pairs1_4[[sp4]]),]
	)
	OG_pairs2_3_one2one=as.data.frame(
		OG_pairs2_3[!duplicated(OG_pairs2_3[[sp2]]) & !duplicated(OG_pairs2_3[[sp3]]),]
	)
	OG_pairs2_4_one2one=as.data.frame(
		OG_pairs2_4[!duplicated(OG_pairs2_4[[sp2]]) & !duplicated(OG_pairs2_4[[sp4]]),]
	)
	OG_pairs3_4_one2one=as.data.frame(
		OG_pairs3_4[!duplicated(OG_pairs3_4[[sp3]]) & !duplicated(OG_pairs3_4[[sp4]]),]
	)
	
	if (one2one) {
		
		OG_pairs1_2 = OG_pairs1_2_one2one
		OG_pairs1_3 = OG_pairs1_3_one2one
		OG_pairs1_4 = OG_pairs1_4_one2one
		OG_pairs2_3 = OG_pairs2_3_one2one
		OG_pairs2_4 = OG_pairs2_4_one2one
		OG_pairs3_4 = OG_pairs3_4_one2one
		mdt1=merge.data.table(OG_pairs1_2,OG_pairs1_3,by=sp1)
		mdt2=merge.data.table(mdt1,OG_pairs1_4,by=sp1)
		mdt3=merge.data.table(mdt2,OG_pairs2_3,by=c(sp2,sp3))
		mdt4=merge.data.table(mdt3,OG_pairs2_4,by=c(sp2,sp4))
		mdt=merge.data.table(mdt4,OG_pairs3_4,by=c(sp3,sp4))
		cross1=mdt[[sp1]]
		cross2=mdt[[sp2]]
		cross3=mdt[[sp3]]
		cross4=mdt[[sp4]]
		
	} else { 
		
		#stop("one2many not implemented yet, set one2one=TRUE")
		
		# Allow for one2many (3 max) relationships, DUPLICATING entries for the paralogs.
		multiOGs <- function(OG_pairs) {
			# one2one
			one2one <- OG_pairs[[1]]
			# one2many (2 or three)
			OG_pairs[[1]] <- make.names(OG_pairs[[1]],unique=TRUE)
			one2many <- grep("\\.[12]$",OG_pairs[[1]],value=TRUE)
			one2many <- unique(sort(c(one2many, str_remove(one2many,"\\.[12]$"))))
			
			OG_pairs = OG_pairs[OG_pairs[[1]] %in% c(one2one,one2many),]
			cross <- str_remove(OG_pairs[[1]],"\\.[12]$")
			#list(OG_pairs, cross)
			OG_pairs[[1]] <- cross
			OG_pairs
		}
		
		ogs_dts <- list(
			OG_pairs1_2_one2one, OG_pairs1_3_one2one, OG_pairs1_4_one2one,
			OG_pairs2_3_one2one, OG_pairs2_4_one2one, OG_pairs3_4_one2one
		)
		og_dts <- lapply(list(
			OG_pairs1_2, OG_pairs1_3, OG_pairs1_4, 
			OG_pairs2_3, OG_pairs2_4, OG_pairs3_4
		), multiOGs)
		
		og_list <- lapply(sp_names, function(sp) {
			ogids <- grep(sp,lapply(og_dts, colnames))
			multis <- og_dts[ogids]
			ogsids <- setdiff(1:length(og_dts),ogids)
			singles <- ogs_dts[ogsids]
			
			all_multis <- na.omit(Reduce(function(...) merge(..., all=TRUE, by=sp, allow.cartesian=TRUE), multis))
			mdt1 <- merge.data.table(singles[[1]],singles[[2]])
			mdt2 <- merge.data.table(mdt1,singles[[3]], by=intersect(colnames(mdt1),colnames(singles[[3]])))
			dt <- merge.data.table(all_multis,mdt2,by=setdiff(sp_names,sp))
			setcolorder(dt,sp_names)
			dt
		})
		
		mdt <- unique(rbindlist(og_list))
		
		cross1=mdt[[sp1]]
		cross2=mdt[[sp2]]
		cross3=mdt[[sp3]]
		cross4=mdt[[sp4]]
		
		mdt[[sp1]] <- make.names(mdt[[sp1]], unique=TRUE)
		mdt[[sp2]] <- make.names(mdt[[sp2]], unique=TRUE)
		mdt[[sp3]] <- make.names(mdt[[sp3]], unique=TRUE)
		mdt[[sp4]] <- make.names(mdt[[sp4]], unique=TRUE)
		
		if (restrict_paralogs) {
			ll <- lapply(sp_names, function(sp) grepl("*\\.[3-9]",mdt[[sp]]))
			ol <- ll[[1]] | ll[[2]] | ll[[3]] | ll[[4]]
			remove_genes <- unique(str_remove(mdt[[sp1]][ol],"\\.\\d+"))
			remove_genes_length <- length(remove_genes)
			remove_genes_message <- paste(head(remove_genes,pmin(remove_genes_length,6)),collapse=", ")
			if (remove_genes_length > 6) remove_genes_message <- paste0(remove_genes_message, ", ...")
			message("Removed ", remove_genes_length, " genes with more than three paralogs.")
			message(remove_genes_message)
			mdt <- mdt[grep(paste(remove_genes,collapse="|"),mdt[[sp1]],invert=TRUE),]
			
			cid <- !(cross1 %in% remove_genes)
			cross1=cross1[cid]
			cross2=cross2[cid]
			cross3=cross3[cid]
			cross4=cross4[cid]
		}
		
	}
	
	# Reorder matrices, quantile_norm and merge matrices
	if (any(c(class(sp1_fp),class(sp2_fp),class(sp3_fp)) == "dgCMatrix")) {
		merged1=Matrix::Matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]),sparse=TRUE)
		merged2=Matrix::Matrix(cbind(merged1,sp3_fp[cross3,]),sparse=TRUE)
		merged=Matrix::Matrix(cbind(merged2,sp4_fp[cross4,]),sparse=TRUE)
	} else {
		merged1=data.matrix(cbind(sp1_fp[cross1,],sp2_fp[cross2,]))
		merged2=data.matrix(cbind(merged1,sp3_fp[cross3,]))
		merged=data.matrix(cbind(merged2,sp4_fp[cross4,]))
	}
	rnm=mdt[[sp1]]
	rownames(merged)=rnm
	cnm=colnames(merged)
	
	if (quant_norm) {
		if (any(class(merged) == "dgCMatrix")) {
			sparse=TRUE
			merged=as.matrix(merged)
			copymat=FALSE
		} else {
			sparse=FALSE
			copymat=TRUE
		}
		merged=tryCatch({
			preprocessCore::normalize.quantiles(merged,copy=copymat)
		}, error = function(e){
			quantile_normalisation(merged)
		})
		if (sparse)
			merged=Matrix::Matrix(merged,sparse=TRUE)
		rownames(merged)=rnm
		colnames(merged)=cnm
	}
	
	# Select variable genes in ALL species
	nc1=ncol(sp1_fp)
	nc2=ncol(sp1_fp) + ncol(sp2_fp)
	nc3=nc2 + ncol(sp3_fp)
	nc4=ncol(merged)
	if (!is.null(cross_fc_thrs) & !is.null(cross_n)) {
		var1= apply(merged[,1:nc1],1,function(x) 
			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		var2= apply(merged[,c(nc1 + 1):nc2],1,function(x) 
			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		var3= apply(merged[,c(nc2 + 1):nc3],1,function(x) 
			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		var4= apply(merged[,c(nc3 + 1):nc4],1,function(x) 
			sort(x,decreasing=TRUE)[cross_n]) > cross_fc_thrs
		top_cross=names(which(var1 & var2 & var3 & var4))
	} else {
		top_cross=rownames(merged)
	}
	
	# Return a list with both matrices
	out_m1=merged[,1:nc1]
	rownames(out_m1)=cross1
	out_m2=merged[,(nc1 + 1):nc2]
	rownames(out_m2)=make.names(cross2,unique=TRUE)
	out_m3=merged[,(nc2 + 1):nc3]
	rownames(out_m3)=make.names(cross3,unique=TRUE)
	out_m4=merged[,(nc3 + 1):nc4]
	rownames(out_m4)=make.names(cross4,unique=TRUE)
	
	ids2 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_2[[sp1]] == x)))
	top_cross_sp2=OG_pairs1_2[ids2,][[sp2]]
	ids3 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_3[[sp1]] == x)))
	top_cross_sp3=OG_pairs1_3[ids3,][[sp3]]
	ids4 <- unlist(lapply(top_cross,function(x) which(OG_pairs1_4[[sp1]] == x)))
	top_cross_sp4=OG_pairs1_4[ids4,][[sp4]]
	
	csps <- list(
		merged=merged,
		OG_pairs1_2=OG_pairs1_2, OG_pairs1_3=OG_pairs1_3, OG_pairs1_4=OG_pairs1_4, 
		OG_pairs2_3=OG_pairs2_3, OG_pairs2_4=OG_pairs2_4, OG_pairs3_4=OG_pairs3_4,
		sp1=out_m1, sp2=out_m2, sp3=out_m3, sp4=out_m4,
		top_cross_sp1=top_cross, top_cross_sp2=top_cross_sp2, top_cross_sp3=top_cross_sp3, top_cross_sp4=top_cross_sp4
	)
	return(csps)
	
}


# # # # # # # # # # # # # # # # #
#                               #
#    CSPS ANALYSIS FUNCTIONS    #
#   take csps object as input   #
#                               #
# # # # # # # # # # # # # # # # #

#' Plots cross-species correlation matrix
#' 
#' Takes a csps object, and optionally a path to metacell annotation table files.
#' Saves heatmap of correlation between two species and return correlation matrix
#' used for plotting, as well as overlapping genes supporting each pairwise relation 
#' (defined by binarizing gene expression with `fc_thrs`).
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param output_file output png image filename
#' @param cluster_together logical, whether to cluster together metacells from 
#'   both species (i.e. rows and columns)
#' @param var_genes either a character vector of genes for comparison, or logical, 
#'   if TRUE, using the csps-defined top-cross genes, if FALSE, using all genes
#' @param annotation_file_1,annotation_file_2 metacell annotation tsv files with 
#'   three columns: metacell, cell_type, color 
#' @param species_1_level,species_2_level character, column name in `annotation_file_1` 
#'   that is being used as colnames in the `csps$merged` object. Defaults to `metacell`.
#' @param cor_method character, corelation method to use, one of the following: 
#'   `c("pearson","spearman","kendall","jaccard")` (default: "jaccard")
#' @param cor_max numeric, color scaling max value (default: 1)
#' @param cor_min numeric, color scaling min value (default: 0)
#' @param cex_dot numeric, dot size scaling factor (default: 1)
#' @param pow numeric, for plotting, raise the correlation to the power of 
#'   `pow` (default: 1)
#' @param fc_thrs numeric, fold change threshold for binarizing gene expression 
#'   (default: 1.2)
#'   when calculating Jaccard distance (default: 2); only used when `cor_method=='jaccard'`
#' @param reorder_sp2 logical, whether to plot reordered metacells in second species
#'   (default: FALSE)
#' @param reorder_by_ann1,reorder_by_ann2 logical, whether to plot metacells in order 
#'   in which they appear in annotation files (default: FALSE)
#' @param width numeric, figure height (in px)
#' @param height numeric, figure width (in px)
#' @param res numeric (default: NA)
#' @param annotation_size numeric, height of the annotation color bar
#' @param label_font_size numeric, size of annotation labels
#' @param annot_length_chr numeric, max length of column and row annotations in heatmap 
#' @param show_legend_heatmap_annot logical, whether to show annotations for heatmap color tracks (default is FALSE)
#' 
#' @return a list with following elements: 
#'   1) `heatmap` complex heatmap object
#'   2) `cor_matrix` correlation matrix used for plotting
#'   3) `overlap_matrix` matrix with number of overlapping genes
#'   4) `overlapping_genes` list of overlapping genes, nested list where at 
#'   the first level are the columns from the first matrix, and at the second 
#'   level are the columns from the second matrix.
#' 
csps_plot_correlation_matrix=function(
	csps, cluster_together=FALSE, var_genes=TRUE, 
	annotation_file_1=NULL, annotation_file_2=NULL,
	species_1_level = "metacell", species_2_level = "metacell",
	cor_method="jaccard", fc_thrs=1.2, 
	heatmap_colors = c("white","gray99","orange","orangered2","#520c52"),
	dotplot_colors = c("white", "red"),
	output_file=NULL, width=3000, height=3000, res=NA,
	reorder_sps2=FALSE, reorder_by_ann1=FALSE, reorder_by_ann2=FALSE,
	pow = 1, plot_type="heatmap", 
	annotation_size = 3, label_font_size = 5, 
	cor_max=1, 
	cor_min=0, 
	cex_dot=1,
	grid = TRUE,  annotation_grid_1 = TRUE, annotation_grid_2 = TRUE,
	annot_length_chr = 40,
	show_legend_heatmap_annot = FALSE
){ 
	
	merged = cbind.data.frame(csps$sp1,csps$sp2)

	# get variable genes: if var_genes is a vector, keep it.
	# if it's TRUE or FALSE, retrieve them from csps object
	# - TRUE: variable in species 1
	# - FALSE: all
	if (class(var_genes) == "logical") {
		if (var_genes) {
			var_genes=csps$top_cross_sp1
		} else {
			var_genes=rownames(csps$merged)
		}
	}
	
	# ordering
	if (cluster_together) {
		cor_joint=cor(merged[var_genes,],method=cor_method) ^ pow
		diag(cor_joint)=NA
		hc=hclust(dist(cor_joint),method="ward.D2")
		cor_mat <- cor_joint[hc$order,hc$order]
	} else {
		m1 <- merged[var_genes,1:ncol(csps$sp1)]
		m2 <- merged[var_genes,(ncol(csps$sp1) + 1):ncol(merged)]
		# binarize for overlap
		m1j <- m1
		m1j[] <- 0 
		m2j <- m2
		m2j[] <- 0
		m1j[!(m1 < fc_thrs)] <- 1
		m2j[!(m2 < fc_thrs)] <- 1
		m1j <- data.matrix(m1j)
		m2j <- data.matrix(m2j)
		# genes in common
		out_list <- vector("list",length = ncol(m1j))
		for (i in 1:ncol(m1j)) {
			intglist <- lapply(1:ncol(m2j), function(j) {
				ints <- which(m1j[,i] > 0 & m2j[,j] > 0)
				if (length(ints) == 0) {
					NA 
				} else {
					pairs_ids <- match(rownames(m1j)[ints], rownames(csps$sp1))
					list(rownames(csps$sp1)[pairs_ids],rownames(csps$sp2)[pairs_ids])
				}
			})
			names(intglist) <- colnames(m2j)
			out_list[[i]] <- intglist
		}
		names(out_list) <- colnames(m1j)
		
		# correlation matrix calculation
		# message(sprintf("%s correlation…", cor_method))
		if (cor_method == "jaccard") {
			cor_m = jaccard(m1j,m2j) ^ pow
			# } else if (cor_method=="mahalanobis") {
			#   mt <- t(merged[var_genes,])
			#   mhdist <- biotools::D2.dist(data=mt,cov=cov(mt))
			#   hcl <- hclust(mhdist)
		} else if (cor_method == "kld") {
			mm=cbind(m1,m2)
			cor_m=calcKLD(mm)$mat ^ pow
			cor_m=cor_m[colnames(m1),colnames(m2)]
		} else {
			cor_m=cor(m1,m2,method=cor_method) ^ pow
		}

		# reorder sps 2?
		if (reorder_sps2) {
			order_cols=names(sort(apply(cor_m,2,function(x) which.max(x)))) 
		} else { 
			order_cols=colnames(cor_m) 
		}
		cor_mat <- cor_m[,order_cols]

	}

	# plotting matrix
	col_ord <- colnames(cor_mat)
	row_ord <- rownames(cor_mat)
	cor_mat_name <- cor_method
	
	
	# METACELL ANNOTATIONS
	# first, for species 1
	if (!is.null(annotation_file_1)) {
		clust_anno_size  <- unit(annotation_size,"mm")
		annr = fread(annotation_file_1,header=TRUE)
		setnames(annr,c("metacell","cell_type","color"))
		# if (!any(annr$metacell %in% colnames(csps$merged))) {
		# 	preffix1 = stringr::str_split(colnames(csps$sp1), "_")[[1]] [1]
		# 	annr$metacell = paste(preffix1,annr$metacell, sep = "_")
		# }
		if (reorder_by_ann1) {
			tryCatch({
				row_ord <- annr$metacell
				cor_mat <- cor_mat[row_ord,]
				
			}, error=function(e) stop(
				"Different annotation and matrix names: ", 
				setdiff(row_ord,rownames(cor_mat)), " vs ",
				setdiff(rownames(cor_mat),row_ord),
			)
			)
		}
		
		# get metacell order, color and cluster information
		rid <- unlist(lapply(row_ord, match, table=as.data.frame(annr)[,species_1_level]))
		rid_label = stringr::str_trunc(rid, annot_length_chr)
		rid_label_l = stringr::str_pad(rid, annot_length_chr, side="left")
		rid_label_r = stringr::str_pad(rid, annot_length_chr, side="left")
		row_clusts <- annr[rid,]$cell_type
		row_clust_col <- annr[rid,]$color
		names(row_clust_col) <- row_clusts
		
		# map it to complex heatmap object, as row-wise annotations
		right_row_col_ha <- HeatmapAnnotation(
			which = "row", MC = row_clusts, col = list(MC = row_clust_col),
			lab = anno_text(
				x = rid_label_l,
				which = "row",
				gp = gpar(fontsize = label_font_size)),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = show_legend_heatmap_annot, gap = unit(annotation_size / 3,"mm")
		)
		left_row_col_ha <- HeatmapAnnotation(
			which = "row",
			lab = anno_text(
				x = rid_label_r,
				which = "row",
				gp = gpar(fontsize = label_font_size)),
			MC = row_clusts,
			col = list(MC = row_clust_col), 
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = show_legend_heatmap_annot, gap = unit(annotation_size / 3,"mm")
		)
		
	} else {
		
		# if no mc annotation is available for species 1, do this
		right_row_col_ha <- HeatmapAnnotation(
			which = "row",
			lab = anno_text(which = "row", row_ord, gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = show_legend_heatmap_annot
		)
		left_row_col_ha <- right_row_col_ha
		
	}
	
	# same, for species 2
	if (!is.null(annotation_file_2)) {
		clust_anno_size  <- unit(annotation_size,"mm")
		annc <- fread(annotation_file_2,header=TRUE)
		setnames(annc,c("metacell","cell_type","color"))
		# if (!any(annc$metacell %in% colnames(csps$merged))) {
		# 	preffix2 = stringr::str_split(colnames(csps$sp2), "_")[[1]] [1]
		# 	annc$metacell <- paste(preffix2, annc$metacell, sep="_")
		# }
		if (reorder_by_ann2) {
			col_ord <- annc$metacell
			cor_mat <- cor_mat[,col_ord]
		} 
		
		# get metacell order, color and cluster information
		cid <- unlist(lapply(col_ord, match, table=as.data.frame(annc)[,species_2_level]))
		col_label = stringr::str_trunc(col_ord, annot_length_chr)
		col_label_t = stringr::str_pad(col_label, annot_length_chr, side="right")
		col_label_b = stringr::str_pad(col_label, annot_length_chr, side="left")
		col_clusts <- annc[cid,]$cell_type
		col_clust_col <- annc[cid,]$color
		names(col_clust_col) <- col_clusts
		# print(cid)
		# print(cid_label)
		# print(col_clusts)
		# print(col_ord)
		
		# map it to complex heatmap object, as col-wise annotations
		top_column_col_ha <- HeatmapAnnotation(
			which = "column", 
			lab = anno_text(
				x = col_label_t,
				which = "column",
				just = "left", 
				location = unit(0, "npc"), 
				gp = gpar(fontsize = label_font_size)), 
			MC = col_clusts, 
			col = list(MC = col_clust_col),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = show_legend_heatmap_annot, gap = unit(annotation_size / 3,"mm")
		)
		bottom_column_col_ha <- HeatmapAnnotation(
			which = "column", col = list(MC = col_clust_col), MC = col_clusts,
			lab = anno_text(
				x = col_label_b,
				which = "column", 
				gp = gpar(fontsize = label_font_size)),
			border = TRUE, simple_anno_size = clust_anno_size,
			show_annotation_name = FALSE, show_legend = show_legend_heatmap_annot, gap = unit(annotation_size / 3,"mm")
		)
		
	} else {
		
		# if no mc annotation is available for species 2, do this
		top_column_col_ha <- HeatmapAnnotation(
			which = "column", 
			lab = anno_text(
				which = "column", col_ord, 
				gp = gpar(fontsize = label_font_size)),
			border = FALSE, show_annotation_name = FALSE, show_legend = show_legend_heatmap_annot
		)
		bottom_column_col_ha <- top_column_col_ha
	}
	
	
	# intersect
	cl_int <- overlap(m1j,m2j)
	ovrl_mat <- cl_int[row_ord,col_ord]
	sf <- max(ovrl_mat,na.rm=TRUE)
	cor_mat_sc <- ovrl_mat / sf
	
	# heatmap
	cor_mat_plot = pmin(cor_mat,cor_max)
	cor_mat_plot = pmax(cor_mat_plot,cor_min)
	if (plot_type == "dotplot") {
		
		# dotplot colors
		dotplot_colors = c("white", "red")
		if (cor_method == "kld") dotplot_colors <- rev(dotplot_colors)
		cor_color = circlize::colorRamp2(c(0, cor_max), dotplot_colors)
		
		hm <- Heatmap(
			cor_mat_plot, col = cor_color, name = cor_mat_name, border = TRUE, 
			rect_gp = gpar(type = "none"),
			cell_fun = function(j, i, x, y, width, height, fill) {
				if (grid == TRUE)
					grid.rect(
						x = x, y = y, width = width, height = height, 
						gp = gpar(col = "gray50", fill = NA, lty = 1, lwd = 0.1)
					)
				grid.circle(
					x = x, y = y, r = abs(cor_mat_sc[i, j]) / 2 * max(unit.c(width, height)) * cex_dot, 
					gp = gpar(fill = cor_color(cor_mat[i, j]), col = "white")
				)
			},
			cluster_rows = FALSE, cluster_columns = FALSE, 
			show_row_names = FALSE, show_column_names = FALSE,
			right_annotation = right_row_col_ha, left_annotation = left_row_col_ha,  
			top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
			heatmap_legend_param = list(
				#col_fun = cor_color, at = cols_range,
				title = cor_mat_name, border = TRUE,
				legend_height = unit(2, "cm"), 
				grid_width = unit(annotation_size / 2,"mm"),
				title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
				labels_gp = gpar(fontsize = label_font_size)
			)
		)
		
		
	} else if (plot_type == "heatmap") {
		
		# Plot heatmap (default)
		# create color palette
		#col_fun = circlize::colorRamp2(c(-2,0,2), heatmap_colors)
		col_fun = colorRampPalette(colors = heatmap_colors)
		cor_color = col_fun( 40 )
		if (cor_method == "kld") {
			cor_color = rev(cor_color)
		}
		
		# This checks for an old ComplexHeatmap version to adapt the final plot 
		# but this looks like an unnecessary source of downstream problems. I'd drop it.
		cv <- as.numeric(stringr::str_extract(as.character(packageVersion("ComplexHeatmap")),"\\d\\.\\d")) 
		if (cv < 2) { 
			hm <- Heatmap(  
				cor_mat_plot, 
				name=cor_mat_name,
				rect_gp = gpar(col = "gray50", lwd = 0.2), 
				cluster_rows = FALSE, cluster_columns = FALSE,
				show_row_names = FALSE, show_column_names = FALSE,
				top_annotation = top_column_col_ha, bottom_annotation = bottom_column_col_ha,
				heatmap_legend_param = list(
					title = cor_mat_name, border = TRUE,
					legend_height = unit(1, "cm"), 
					grid_width = unit(4,"mm"),
					title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
					labels_gp = gpar(fontsize = label_font_size)
				)
			)
		} else {
			hm <- Heatmap(
				cor_mat_plot, 
				name=cor_mat_name,
				col=cor_color, 
				rect_gp = gpar(col = "black", lwd = 0.5), 
				border=TRUE, 
				cluster_rows = FALSE, 
				cluster_columns = FALSE,
				show_row_names = FALSE, 
				show_column_names = FALSE,
				right_annotation = right_row_col_ha, 
				left_annotation = left_row_col_ha,  
				top_annotation = top_column_col_ha, 
				bottom_annotation = bottom_column_col_ha,
				heatmap_legend_param = list(
					title = cor_mat_name, border = TRUE,
					legend_height = unit(6, "cm"), 
					grid_width = unit(annotation_size,"mm"),
					title_position = "leftcenter-rot", 
					title_gp = gpar(fontsize = label_font_size),
					labels_gp = gpar(fontsize = label_font_size)
				)
			)
		}
		
	}
	
	# if output_file is null, do not plot heatmap; return instead
	plotting_function(output_file, width, height, res, EXP = draw(hm) )
	
	ht_opt(RESET = TRUE)
	
	# output
	return(list(heatmap=hm, cor_matrix=cor_mat, overlap_matrix=ovrl_mat, overlapping_genes=out_list))
	
}


#' Create cross-species correlation matrix
#' 
#' Takes a csps object and returns a sp1 (rows) v. sp2 (columns) correlation matrix
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param use_var_genes either a character vector of genes for comparison, or logical;
#'   if TRUE, using the csps-defined top-cross genes; if FALSE, use all genes
#' @param cor_method character, corelation method to use, one of the following: 
#'   `c("pearson","spearman","kendall","jaccard","jsd","kld","wpearson","wspearman","shaferindex")` 
#'   (default: "jaccard")
#' @param gene_weights a numeric vector of non-negative weights used to calculated weighted
#'   correlations when method is "wpearson" o "wspearman". It has to match the order and length 
#    of the vector defined by the `use_var_genes` variable (i.e. either the supplied vector, 
#'   or the csps-defined vector of gene names).
#' @param pow numeric, for plotting, raise the correlation to the power of 
#'   `pow` (default: 1)
#' @param fc_thrs numeric, fold change threshold for binarising gene expression 
#'   (default: 1.2) when calculating Jaccard distance; only used when `cor_method="jaccard" or `"shaferindex"`
#' @param report_overlaps logical, whether to report shared genes between cell types (relies
# '   on the `fc_thrs` binarisation threshold
#' @param quantile_discretisation whether to convert input matrix to a discretised matrix
#'   where gene expression levels are discretised to quantiles (default: FALSE)
#' @param quantile_n if `quantile_discretisation=TRUE`, use this many quantiles (default = 10)
#' 
#' @return a list with following elements: 
#'   1) `cor_matrix` correlation matrix used for plotting
#'   2) `overlap_matrix` matrix with number of overlapping genes
#'   3) `overlapping_genes` list of overlapping genes, nested list where at 
#'   the first level are the columns from the first matrix, and at the second 
#'   level are the columns from the second matrix.
#' 
csps_correlation_matrix = function(
	csps, 
	use_var_genes = TRUE, 
	cor_method = "jaccard", 
	gene_weights = NULL,
	report_overlaps = TRUE,
	fc_thrs = 1.2, 
	pow = 1,
	add_sps_prefix = TRUE,
	prefix_sp1 = NULL,
	prefix_sp2 = NULL,
	quantile_discretisation = FALSE,
	quantile_n = 10) { 
	
	# get matrices
	mm = csps$merged
	m1 = csps$sp1
	m2 = csps$sp2

	# vector of gene names from sp1	
	gene_names = rownames(m1)
	
	if (quantile_discretisation) {
		mm = apply(mm, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
		m1 = apply(m1, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
		m2 = apply(m2, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
	}
	
	# use same rownames across all matrices (from sp1)
	rownames(mm) = gene_names
	rownames(m1) = gene_names
	rownames(m2) = gene_names

	# use unique colnames
	if (add_sps_prefix) {
		if (is.null(prefix_sp1)) prefix_sp1 = "sp1"
		if (is.null(prefix_sp2)) prefix_sp1 = "sp2"
		colnames(m1) = paste(prefix_sp1, colnames(m1), sep = "|")
		colnames(m2) = paste(prefix_sp2, colnames(m2), sep = "|")
		colnames(mm) = c(colnames(m1), colnames(m2))
	}

	# get variable genes: if use_var_genes is a vector, keep it.
	# if it's TRUE or FALSE, retrieve them from the csps object
	if (class(use_var_genes) == "logical" & any(use_var_genes == TRUE)) {
		var_genes = csps$top_cross_sp1
		message(sprintf("csps matrix | use %i genes (from `csps` object)", length(var_genes)))
	} else if (class(use_var_genes) == "logical" & any(use_var_genes == FALSE)) {
		var_genes = rownames(csps$merged)
		message(sprintf("csps matrix | use %i genes (complete matrix)", length(var_genes)))
	} else if (class(use_var_genes) == "character") {
		var_genes = use_var_genes
		message(sprintf("csps matrix | use %i variable genes (given)", length(var_genes)))
	} else {
		stop("`use_var_genes` has to be either logical, or a vector of genes from species 1")
	}
	
	# subset matrices to variable genes
	mm = mm [ var_genes, ]
	m1 = m1 [ var_genes, ]
	m2 = m2 [ var_genes, ]
	var_genes_sp2 = csps$og_pairs$sp2 [ csps$og_pairs$sp1 %in% var_genes ]
	
	# sanity check
	if (nrow(mm) == 0) {
		stop("No rows (genes) left in merged matrix!")
	}
	
	# list of genes for each species
	genes_sp1 = rownames(csps$sp1) [ var_genes %in% rownames(csps$sp1) ]
	genes_sp2 = rownames(csps$sp2) [ var_genes %in% rownames(csps$sp1) ]
	
	# correlation matrix calculation
	message(sprintf("csps matrix | method: %s", cor_method))
	if (cor_method == "jaccard") {
		
		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2, threshold = fc_thrs)
		
		# get jaccard distance
		com = jaccard(m1_bin, m2_bin) ^ pow
		
	} else if (cor_method == "jsd") {

		# calculate similarity based on Jensen-Shannon divergence, for a count-like matrix 
		# `est.prob = "empirical"` converts the count matrix to a probability matrix, like this:
		# counts -> 1:10
		# probs  -> 1:10 / sum(1:10)
		require("philentropy")
		com_d = philentropy::JSD(t(as.matrix(mm)), unit = "log2", est.prob = "empirical")
		# convert Jensen-Shannon distance (i.e. sqrt divergence) to similarity matrix
		com_s = as.matrix(1 - sqrt(com_d))
		colnames(com_s) = colnames(mm)
		rownames(com_s) = colnames(mm)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1) , colnames(com_s) %in% colnames(m2) ]

	} else if (cor_method == "jsdnp") {

		# calculate similarity based on Jensen-Shannon divergence, for a non-count matrix
		require("philentropy")
		com_s = philentropy::JSD(t(as.matrix(mm)), unit = "log2", est.prob = NULL)
		com_s = sqrt(com_s)
		colnames(com_s) = colnames(mm)
		rownames(com_s) = colnames(mm)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1) , colnames(com_s) %in% colnames(m2) ]

	} else if (cor_method == "kld") {

		# old approach: our own KLD function
		# com = calcKLD(mm)$mat ^ pow
		# com = com[ colnames(m1), colnames(m2) ]
		# com = 1 - com
		
		#  new approach: philentropy KLD (faster!)
		require("philentropy")
		com_d = philentropy::KL(t(as.matrix(mm)), unit = "log2", est.prob = "empirical")
		# convert KL distance (i.e. sqrt divergence) to similarity matrix
		com_s = as.matrix(1 - sqrt(com_d))
		colnames(com_s) = colnames(mm)
		rownames(com_s) = colnames(mm)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1) , colnames(com_s) %in% colnames(m2) ]
		
	} else if (cor_method == "shaferindex") {
		
		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2, threshold = fc_thrs)

		# get shaferindex distance
		com = shaferindex(m1_bin, m2_bin) ^ pow
		
	} else if (cor_method == "onls") {
		
		require("pracma")
		cod = matrix(NA, nrow = ncol(m1), ncol = ncol(m2))
		for (i in 1:ncol(m1)) {
			for (j in 1:ncol(m2)) {
				mi = log10(m1[,i])
				mj = log10(m2[,j])
				po = pracma::odregress(mi, mj)
				sq = po$ssq / nrow(m1)
				cod[i,j] = sq
			}
		}
		rownames(cod) = colnames(m1)
		colnames(cod) = colnames(m2)
		com = scale(1 / cod)

	} else if (cor_method %in% c("wpearson","wspearman")) {
		
		require("WGCNA")
		if (is.null(gene_weights)) {
			stop(sprintf("I can't use weighted correlation method %s because no weights were supplied by `gene_weights`!", cor_method))
		} else if (length(gene_weights) != length(var_genes)) {
			stop(sprintf("I can't use weighted correlation method %s because length of `gene_weights` (%i) does not match length of `var_genes` (%i)!", cor_method, length(gene_weights), length(var_genes)))
		}
		
		# create matrix of weights to apply to each footprint (based on gene-level weights applied across all tissues)
		gene_weights_m1 = matrix(rep(gene_weights, ncol(m1)), ncol = ncol(m1))
		gene_weights_m2 = matrix(rep(gene_weights, ncol(m2)), ncol = ncol(m2))
		
		# weighted correlation
		wgcna_cor_method = gsub("^w","", cor_method)
		com = WGCNA::cor(m1, m2, method = wgcna_cor_method, weights.x = gene_weights_m1, weights.y = gene_weights_m2)
		com = com ^ pow
		
	
	} else {
		
		# any other correlation metric
		com = cor(m1, m2, method = cor_method) ^ pow
		
	}

	# get shared genes between each pair of cell types
	# this is a nested named list: first level are cell types in sp1, second level are cell types in sp2
	if (report_overlaps) {

		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2, threshold = fc_thrs)
		
		# empty upper level lists (shi_l1 will contain gene names from sp1, shi_l2 for sp2)
		shi_l1 = vector("list", length = ncol(m1_bin))
		shi_l2 = vector("list", length = ncol(m1_bin))
		names(shi_l1) = colnames(m1_bin)
		names(shi_l2) = colnames(m1_bin)

		# loop through sp1 cell types
		for (i in 1:ncol(m1_bin)) {
			
			# get vector of gene presence in sp1
			v1 = m1_bin[,i]
			
			# empty lower level lists for sp2 cell types
			# create named list of cell types in sp2
			shj_l1 = vector("list", length = ncol(m2_bin))
			shj_l2 = vector("list", length = ncol(m2_bin))
			names(shj_l1) = colnames(m2_bin)
			names(shj_l2) = colnames(m2_bin)
			
			# loop through sp2 cell types
			# populate list with a vector of genes shared in cell type i and j
			for (j in 1:ncol(m2_bin)) {
				v2 = m2_bin[,j]
				shj_l1[[j]] = var_genes     [ which(v1 > 0 & v2 > 0) ]
				shj_l2[[j]] = var_genes_sp2 [ which(v1 > 0 & v2 > 0) ]
			}
			
			# populate upper level entry (sp1 cell type i) with lower level lists (sp2 cell types)
			shi_l1[[i]] = shj_l1
			shi_l2[[i]] = shj_l2
			
		}
		
		# get counts of shared genes between cell types
		sha_m = t(sapply(1:length(shi_l1), function (i) { 
			lengths(shi_l1[[i]])
		}))
		rownames(sha_m) = colnames(m1_bin)
		colnames(sha_m) = colnames(m2_bin)
		
		# shared lists
		sha_l = list(sp1 = shi_l1, sp2 = shi_l2)
		
	} else {
		
		sha_l = list(sp1 = NULL, sp2 = NULL)
		sha_m = NULL
		
	}
	
	# return
	return(list(cor_matrix = com, overlap_matrix = sha_m, overlap_genes = sha_l, method = cor_method, var_genes = var_genes))
	
}


#' Create cross-species correlation matrix
#' 
#' Takes a csps object and returns a sp1 (rows) v. sp2 (columns) correlation matrix
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param use_var_genes either a character vector of genes for comparison, or logical;
#'   if TRUE, using the csps-defined top-cross genes; if FALSE, use all genes
#' @param cor_method character, corelation method to use, one of the following: 
#'   `c("pearson","spearman","kendall","jaccard","jsd","kld","wpearson","wspearman","shaferindex")` 
#'   (default: "jaccard")
#' @param gene_weights a numeric vector of non-negative weights used to calculated weighted
#'   correlations when method is "wpearson" o "wspearman". It has to match the order and length 
#    of the vector defined by the `use_var_genes` variable (i.e. either the supplied vector, 
#'   or the csps-defined vector of gene names).
#' @param pow numeric, for plotting, raise the correlation to the power of 
#'   `pow` (default: 1)
#' @param fc_thrs numeric, fold change threshold for binarising gene expression 
#'   (default: 1.2) when calculating Jaccard distance; only used when `cor_method="jaccard" or `"shaferindex"`
#' @param report_overlaps logical, whether to report shared genes between cell types (relies
# '   on the `fc_thrs` binarisation threshold
#' @param quantile_discretisation whether to convert input matrix to a discretised matrix
#'   where gene expression levels are discretised to quantiles (default: FALSE)
#' @param quantile_n if `quantile_discretisation=TRUE`, use this many quantiles (default = 10)
#' 
#' @return a list with following elements: 
#'   1) `cor_matrix` correlation matrix used for plotting
#'   2) `overlap_matrix` matrix with number of overlapping genes
#'   3) `overlapping_genes` list of overlapping genes, nested list where at 
#'   the first level are the columns from the first matrix, and at the second 
#'   level are the columns from the second matrix.
#' 
csps_correlation_matrix_noobj = function(
	mm,
	m1,
	m2,
	use_var_genes = FALSE, 
	cor_method = "jaccard", 
	gene_weights = NULL,
	report_overlaps = TRUE,
	fc_thrs = 1.2, 
	pow = 1,
	add_sps_prefix = TRUE,
	prefix_sp1 = NULL,
	prefix_sp2 = NULL,
	quantile_discretisation = FALSE,
	quantile_n = 10) { 
	
	# vector of gene names from sp1	
	gene_names = rownames(m1)
	
	if (quantile_discretisation) {
		mm = apply(mm, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
		m1 = apply(m1, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
		m2 = apply(m2, 2, function(i) as.numeric(cut(i, quantile(i), prob = seq(0, 1, length = quantile_n + 1), include.lowest = TRUE)) )
	}
	
	# use same rownames across all matrices (from sp1)
	rownames(mm) = gene_names
	rownames(m1) = gene_names
	rownames(m2) = gene_names

	# use unique colnames
	if (add_sps_prefix) {
		if (is.null(prefix_sp1)) prefix_sp1 = "sp1"
		if (is.null(prefix_sp2)) prefix_sp1 = "sp2"
		colnames(m1) = paste(prefix_sp1, colnames(m1), sep = "|")
		colnames(m2) = paste(prefix_sp2, colnames(m2), sep = "|")
		colnames(mm) = c(colnames(m1), colnames(m2))
	}

	# get variable genes: if use_var_genes is a vector, keep it.
	# if it's FALSE, retrieve them from the merged matrix
	if (class(use_var_genes) == "logical" & any(use_var_genes == FALSE)) {
		var_genes = rownames(mm)
		message(sprintf("csps matrix | use %i genes (complete matrix)", length(var_genes)))
	} else if (class(use_var_genes) == "character") {
		var_genes = use_var_genes
		message(sprintf("csps matrix | use %i variable genes (given)", length(var_genes)))
	} else {
		stop("`use_var_genes` has to be either logical, or a vector of genes from species 1")
	}
	
	# subset matrices to variable genes
	mm_f = mm [ var_genes, ]
	m1_f = m1 [ var_genes, ]
	m2_f = m2 [ var_genes, ]
	
	# sanity check
	if (nrow(mm) == 0) {
		stop("No rows (genes) left in merged matrix!")
	}
	
	# list of genes for each species
	genes_sp1 = rownames(m1) [ var_genes %in% rownames(m1) ]
	genes_sp2 = rownames(m2) [ var_genes %in% rownames(m2) ]
	
	# correlation matrix calculation
	message(sprintf("csps matrix | method: %s", cor_method))
	if (cor_method == "jaccard") {
		
		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1_f, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2_f, threshold = fc_thrs)
		
		# get jaccard distance
		com = jaccard(m1_bin, m2_bin) ^ pow
		
	} else if (cor_method == "jsd") {

		# calculate similarity based on Jensen-Shannon divergence, for a count-like matrix 
		# `est.prob = "empirical"` converts the count matrix to a probability matrix, like this:
		# counts -> 1:10
		# probs  -> 1:10 / sum(1:10)
		require("philentropy")
		com_d = philentropy::JSD(t(as.matrix(mm_f)), unit = "log2", est.prob = "empirical")
		# convert Jensen-Shannon distance (i.e. sqrt divergence) to similarity matrix
		com_s = as.matrix(1 - sqrt(com_d))
		colnames(com_s) = colnames(mm_f)
		rownames(com_s) = colnames(mm_f)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1_f) , colnames(com_s) %in% colnames(m2_f) ]

	} else if (cor_method == "jsdnp") {

		# calculate similarity based on Jensen-Shannon divergence, for a non-count matrix
		require("philentropy")
		com_s = philentropy::JSD(t(as.matrix(mm_f)), unit = "log2", est.prob = NULL)
		com_s = sqrt(com_s)
		colnames(com_s) = colnames(mm_f)
		rownames(com_s) = colnames(mm_f)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1_f) , colnames(com_s) %in% colnames(m2_f) ]

	} else if (cor_method == "kld") {

		# old approach: our own KLD function
		# com = calcKLD(mm)$mat ^ pow
		# com = com[ colnames(m1), colnames(m2) ]
		# com = 1 - com
		
		#  new approach: philentropy KLD (faster!)
		require("philentropy")
		com_d = philentropy::KL(t(as.matrix(mm_f)), unit = "log2", est.prob = "empirical")
		# convert KL distance (i.e. sqrt divergence) to similarity matrix
		com_s = as.matrix(1 - sqrt(com_d))
		colnames(com_s) = colnames(mm_f)
		rownames(com_s) = colnames(mm_f)
		# subset to keep only sp1 as rows and sp2 as columns
		com = com_s [ rownames(com_s) %in% colnames(m1_f) , colnames(com_s) %in% colnames(m2_f) ]
		
	} else if (cor_method == "shaferindex") {
		
		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1_f, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2_f, threshold = fc_thrs)

		# get shaferindex distance
		com = shaferindex(m1_bin, m2_bin) ^ pow
		
	} else if (cor_method == "onls") {
		
		require("pracma")
		cod = matrix(NA, nrow = ncol(m1_f), ncol = ncol(m2_f))
		for (i in 1:ncol(m1_f)) {
			for (j in 1:ncol(m2_f)) {
				mi = log10(m1_f[,i])
				mj = log10(m2_f[,j])
				po = pracma::odregress(mi, mj)
				sq = po$ssq / nrow(m1_f)
				cod[i,j] = sq
			}
		}
		rownames(cod) = colnames(m1_f)
		colnames(cod) = colnames(m2_f)
		com = scale(1 / cod)

	} else if (cor_method %in% c("wpearson","wspearman")) {
		
		require("WGCNA")
		if (is.null(gene_weights)) {
			stop(sprintf("I can't use weighted correlation method %s because no weights were supplied by `gene_weights`!", cor_method))
		} else if (length(gene_weights) != length(var_genes)) {
			stop(sprintf("I can't use weighted correlation method %s because length of `gene_weights` (%i) does not match length of `var_genes` (%i)!", cor_method, length(gene_weights), length(var_genes)))
		}
		
		# create matrix of weights to apply to each footprint (based on gene-level weights applied across all tissues)
		gene_weights_m1 = matrix(rep(gene_weights, ncol(m1_f)), ncol = ncol(m1_f))
		gene_weights_m2 = matrix(rep(gene_weights, ncol(m2_f)), ncol = ncol(m2_f))
		
		# weighted correlation
		wgcna_cor_method = gsub("^w","", cor_method)
		com = WGCNA::cor(m1_f, m2_f, method = wgcna_cor_method, weights.x = gene_weights_m1, weights.y = gene_weights_m2)
		com = com ^ pow
		
	
	} else {
		
		# any other correlation metric
		com = cor(m1_f, m2_f, method = cor_method) ^ pow
		
	}

	# get shared genes between each pair of cell types
	# this is a nested named list: first level are cell types in sp1, second level are cell types in sp2
	if (report_overlaps) {

		# binarise footprints
		m1_bin = sca_binarise_expression(counts = m1_f, threshold = fc_thrs)
		m2_bin = sca_binarise_expression(counts = m2_f, threshold = fc_thrs)
		
		# empty upper level lists (shi_l1 will contain gene names from sp1, shi_l2 for sp2)
		shi_l1 = vector("list", length = ncol(m1_bin))
		shi_l2 = vector("list", length = ncol(m1_bin))
		names(shi_l1) = colnames(m1_bin)
		names(shi_l2) = colnames(m1_bin)

		# loop through sp1 cell types
		for (i in 1:ncol(m1_bin)) {
			
			# get vector of gene presence in sp1
			v1 = m1_bin[,i]
			
			# empty lower level lists for sp2 cell types
			# create named list of cell types in sp2
			shj_l1 = vector("list", length = ncol(m2_bin))
			shj_l2 = vector("list", length = ncol(m2_bin))
			names(shj_l1) = colnames(m2_bin)
			names(shj_l2) = colnames(m2_bin)
			
			# loop through sp2 cell types
			# populate list with a vector of genes shared in cell type i and j
			for (j in 1:ncol(m2_bin)) {
				v2 = m2_bin[,j]
				shj_l1[[j]] = var_genes     [ which(v1 > 0 & v2 > 0) ]
				shj_l2[[j]] = var_genes     [ which(v1 > 0 & v2 > 0) ] ### this should be genes from sp2?
			}
			
			# populate upper level entry (sp1 cell type i) with lower level lists (sp2 cell types)
			shi_l1[[i]] = shj_l1
			shi_l2[[i]] = shj_l2
			
		}
		
		# get counts of shared genes between cell types
		sha_m = t(sapply(1:length(shi_l1), function (i) { 
			lengths(shi_l1[[i]])
		}))
		rownames(sha_m) = colnames(m1_bin)
		colnames(sha_m) = colnames(m2_bin)
		
		# shared lists
		sha_l = list(sp1 = shi_l1, sp2 = shi_l2)
		
	} else {
		
		sha_l = list(sp1 = NULL, sp2 = NULL)
		sha_m = NULL
		
	}
	
	# return
	return(list(cor_matrix = com, overlap_matrix = sha_m, overlap_genes = sha_l, method = cor_method, var_genes = var_genes))
	
}



#' Plot annotated matrix
#' 
#' @param mat any type of data matrix
#' @param name name of the type of data in the matrix (default: "data")
#' @param heatmap_colors vector of colors to map to the data
#' @param min_val,max_val min and max values of the colorscale
#' @param use_raster whether to rasterise
#' @param row_title,col_title titles for rows and columns
#' @param row_labels,col_labels ad-hoc labels for rows and columns (if set to NULL, they are taken from `mat` object)
#' @param max_length_labels truncate `row_labels` and `col_labels` to this maximum length, in characters (default 40)
#' @param fontsize size of labels (default 5 pts)
#' @param row_annot,col_annot either dataframes where the 1st column is a vector of categories for each row/column (same order is assumed) and 2nd is a vector of colors, or simply a vector of categories. Default is NULL, i.e. no annotations.
#' @param row_annot_cols,col_annot_cols named vector of colors, where names are categories that match the vector in `row_annot`/`col_annot` (not necessary if `row_annot`/`col_annot` are dataframes).
#' @param row_annot_legend,col_annot_legend whether to plot row/col annotation legends
#' @param row_cluster,col_cluster if TRUE/FALSE, whether to cluster rows/columns with default parameters. Other built-in options are "pearson", "euclidean", or a precomputed `hclust` object.
#' @param cex_dotplot transformation factor for dot size, if `do_dotplot=TRUE`
#' @param do_dotplot draw a dot plot instead of a heatmap
#' 
#' @return a ComplexHeatmap object
#' 
csps_plot_annotated_matrix = function(
	mat,
	name = "data",
	heatmap_colors = c("white","orange","orangered2","#520c52"),
	min_val = 0,
	max_val = 1,
	use_raster = TRUE,
	row_title = NULL,
	col_title = NULL,
	row_labels = NULL,
	col_labels = NULL,
	max_length_labels = 40,
	fontsize = 5,
	row_annot = NULL,
	row_annot_cols = NULL,
	row_annot_legend = FALSE,
	row_cluster = FALSE,
	col_annot = NULL,
	col_annot_cols = NULL,
	col_annot_legend = FALSE,
	col_cluster = FALSE,
	cex_dotplot = 0.02,
	do_dotplot = FALSE) { 
	
	# libraries
	require("ComplexHeatmap")
	require("circlize")
	
	# get colors for heatmap
	col_fun = circlize::colorRamp2(
		breaks = seq(from = min_val, to = max_val, length.out = length(heatmap_colors)), 
		colors = heatmap_colors)
	
	# Titles	
	# get row title
	if (is.null(row_title)) { 
		row_title = sprintf("n = %i", nrow(mat))
	} else {
		row_title = sprintf("%s, n = %i", row_title, nrow(mat))
	}
	# get col title
	if (is.null(col_title)) { 
		col_title = sprintf("n = %i", ncol(mat))
	} else {
		col_title = sprintf("%s, n = %i", col_title, ncol(mat))
	}
	
	# Row and column names
	if (is.null(row_labels)) {
		row_labels = rownames(mat)
	}
	if (is.null(col_labels)) {
		col_labels = colnames(mat)
	}
	# truncate if necessary
	if (!is.null(max_length_labels)) {
		row_labels = stringr::str_trunc(row_labels, max_length_labels)
		col_labels = stringr::str_trunc(col_labels, max_length_labels)
		row_labels = stringr::str_pad(row_labels, max_length_labels, side = "right")
		col_labels = stringr::str_pad(col_labels, max_length_labels)
	}

	# get row annotations
	# by default, get labels
	ha_row_base = ComplexHeatmap::HeatmapAnnotation(lab = anno_text(row_labels, which = "row", gp = gpar(fontsize = fontsize), just = "left"), which = "row")
	# add coloring info if available
	if (!is.null(row_annot)) {
		# ensure that row_annot is a list of dataframes
		if ("data.frame" %in% class(row_annot)) { row_annot = list(row_annot) }
		# loop through list of dataframes, adding colors
		for (ni in 1:length(row_annot)) {
			# get unique named list of colors
			row_annot_u = unique(row_annot[[ni]])
			row_annot_cols = row_annot_u[,2]
			names(row_annot_cols) = row_annot_u[,1]
			# get vector of categories
			row_annot_cats = row_annot[[ni]][,1]
			# add new annotation track
			ha_row_left  = c(ha_row_base, ComplexHeatmap::HeatmapAnnotation(clusters = row_annot_cats, col = list(clusters = row_annot_cols), which = "row", show_annotation_name = FALSE, show_legend = row_annot_legend))
			ha_row_right = c(ComplexHeatmap::HeatmapAnnotation(clusters = row_annot_cats, col = list(clusters = row_annot_cols), which = "row", show_annotation_name = FALSE, show_legend = row_annot_legend), ha_row_base)
		}
	} else {
		ha_row_left  = ha_row_base
		ha_row_right = ha_row_base
	}
		
	# get column annotations
	# by default, get labels
	ha_col_base = ComplexHeatmap::HeatmapAnnotation(lab = anno_text(col_labels, which = "column", gp = gpar(fontsize = fontsize), just = "right"), which = "column")
	# add coloring info if available
	if (!is.null(col_annot)) {
		# ensure that col_annot is a list of dataframes
		if ("data.frame" %in% class(col_annot)) { col_annot = list(col_annot) }
		# loop through list of dataframes, adding colors
		for (ni in 1:length(col_annot)) {
			# get unique named list of colors 
			col_annot_u = unique(col_annot[[ni]])
			col_annot_cols = col_annot_u[,2]
			names(col_annot_cols) = col_annot_u[,1]
			# get vector of categories
			col_annot_cats = col_annot[[ni]][,1]
			# add new annotation track
			ha_col_top = c(ha_col_base, ComplexHeatmap::HeatmapAnnotation(clusters = col_annot_cats, col = list(clusters = col_annot_cols), which = "column", show_annotation_name = FALSE, show_legend = col_annot_legend))
			ha_col_bot = c(ComplexHeatmap::HeatmapAnnotation(clusters = col_annot_cats, col = list(clusters = col_annot_cols), which = "column", show_annotation_name = FALSE, show_legend = col_annot_legend), ha_col_base)
		}
	} else {
		ha_col_top = ha_col_base
		ha_col_bot = ha_col_base
	}

	
	# should this be a dot plot?
	if (do_dotplot) {
		cell_fun_dotplot = function(j, i, x, y, width, height, fill) {
			range01 = function(x) { (x - min(x)) / (max(x) - min(x)) }
			grid.circle(
				x = x, y = y, 
				r = sqrt(range01(mat)[i, j]) * cex_dotplot, 
				gp = gpar(col = col_fun(mat[i, j]), fill = col_fun(mat[i, j]))
			)
		}
		rect_gp = gpar(type = "none")
	} else {
		cell_fun_dotplot = NULL
		rect_gp = gpar(col = "white", lwd = 0.3)
	}
	
	# how to perform row-wise clustering
	if (row_cluster %in% c("pearson","spearman","kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
		row_cluster_method = row_cluster
		row_cluster = TRUE
	} else if (row_cluster == TRUE) {
		row_cluster_method = "euclidean"
	} else {
		row_cluster = FALSE
		row_cluster_method = NULL
	}
	# how to perform column-wise clustering
	if (col_cluster %in% c("pearson","spearman","kendall", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
		col_cluster_method = col_cluster
		col_cluster = TRUE
	} else if (col_cluster == TRUE) {
		col_cluster_method = "euclidean"
	} else {
		col_cluster = FALSE
		col_cluster_method = NULL
	}			
	
	# heatmap object
	hm = ComplexHeatmap::Heatmap(
		mat, 
		name = name,
		cell_fun = cell_fun_dotplot,
		rect_gp = rect_gp,
		use_raster = use_raster,
		col = col_fun,
		border = TRUE,
		cluster_rows = row_cluster,
		clustering_distance_rows = row_cluster_method,
		cluster_columns = col_cluster,
		clustering_distance_columns = col_cluster_method,
		row_title = row_title,
		column_title = col_title,
		show_row_names = FALSE,
		show_column_names = FALSE,
		# column annotations
		top_annotation = ha_col_top,
		bottom_annotation = ha_col_bot,
		# row annotations
		left_annotation = ha_row_left,
		right_annotation = ha_row_right
	)

	# return heatmap object
	return(hm)		
	
}


#' Identify co-expressed genes in a defined set of metacells/cell types in two compared species
#' Optionally can be subseted by gene list (e.g. showing only TFs)
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param sp1_focus,sp2_focus character, one or more metacell (or cell type) names 
#'   (i.e. columns in `csps$merged`)
#' @param glist character, optional gene list for subsetting the genes (default: NULL)
#' @param fc_thrs numeric, gene fold change threshold (default: 1.5)
#' @param gene_annot_sp1_fn,gene_annot_sp1_fn character, path to gene annotation file for 
#'   the first and the second species in `csps`
#' @param out_fn character, path to png file where expression heatmap will be saved, same filename
#'   will be used to generate other output files (barplot figure and txt file)
#'   (default: NULL, saves files to current working directory)
#' 
csps_identify_coexpressed_genes=function(
	csps, sp1_focus, sp2_focus,
	glist=NULL, fc_thrs=1.5,
	gene_annot_sp1_fn, gene_annot_sp2_fn,
	out_fn=NULL, width_cex=20, height_cex=20, res=NA,
	fc_max=3, annotation_size = 10, label_font_size=12
){ 
	
	if (is.null(out_fn)) {
		out_bfn=paste0("Coexpressed_genes_",paste(c(sp1_focus,sp2_focus),collapse="_"))
	} else {
		out_bfn=str_remove(out_fn,".png")
	}
	message(out_bfn)
	
	annot_sp1=read.table(gene_annot_sp1_fn,h=TRUE,row.names=1,sep="\t",quote="",stringsAsFactors=FALSE,fill=TRUE)
	annot_sp2=read.table(gene_annot_sp2_fn,h=TRUE,row.names=1,sep="\t",quote="",stringsAsFactors=FALSE,fill=TRUE)
	bckg_sp1=setdiff(colnames(csps$sp1),sp1_focus)
	bckg_sp2=setdiff(colnames(csps$sp2),sp2_focus)
	
	f_sp1=apply(csps$merged[,sp1_focus,drop=FALSE],1,function(x) sort(x,decreasing=TRUE)[1]) >= fc_thrs & apply(csps$merged[,bckg_sp1],1,median) < fc_thrs
	f_sp2=apply(csps$merged[,sp2_focus,drop=FALSE],1, function(x) sort(x,decreasing=TRUE)[1]) >= fc_thrs & apply(csps$merged[,bckg_sp2],1,median) < fc_thrs
	
	sps1_selected_ids=names(which(f_sp1 & f_sp2))
	hc=hclust(dist(cor(t(csps$merged[sps1_selected_ids,]))),method="ward.D2")
	gene_order=sps1_selected_ids[hc$order]
	sps2_selected_ids=csps$og_pairs[match(gene_order,csps$og_pairs[[1]]),2]
	
	output_table=cbind.data.frame(annot_sp1[str_remove(gene_order,"\\.\\d"),],annot_sp2[sps2_selected_ids,])
	if (!is.null(out_fn))
		write.table(output_table,file=paste0(out_bfn,"_table.txt"),quote=FALSE,sep="\t",col.names=FALSE,row.names=TRUE)
	
	shades2=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
	hmat <- pmin(csps$merged[gene_order,], fc_max)
	hmtitle <- sprintf("Coexspressed genes in %s and %s", paste(sp1_focus,collapse=","), paste(sp2_focus,collapse=","))
	col_ord <- colnames(hmat)
	row_ord <- str_remove(rownames(hmat),"\\.\\d") #rownames(hmat)
	row_ord_1 <- annot_sp1[row_ord,1]
	row_ord_1[nchar(row_ord_1) > 30] <- paste0(substr(row_ord_1[nchar(row_ord_1) > 30],1,27),"...")
	row_ord_1 <- str_replace(row_ord_1,"\"\"","")
	row_ord_2 <- annot_sp1[row_ord,2]
	row_ord_2[nchar(row_ord_2) > 30] <- paste0(substr(row_ord_2[nchar(row_ord_2) > 30],1,27),"...")
	row_ord_2 <- str_replace(row_ord_2,"\"\"","")
	top_column_col_ha <- HeatmapAnnotation(
		which = "column",
		lab = anno_text(which = "column", just = "left", location = unit(0, "npc"), col_ord, gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	bottom_column_col_ha <- HeatmapAnnotation(
		which = "column",
		lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	left_row_col_ha <- HeatmapAnnotation(
		which = "row",
		lab = anno_text(which = "row", row_ord_2, just = "right", location = unit(1, "npc"), gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	right_row_col_ha <- HeatmapAnnotation(
		which = "row",
		lab = anno_text(which = "row", row_ord_1, gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	hm <- Heatmap(
		hmat, col = shades2, name = "expression", border = TRUE,
		rect_gp = gpar(col = "gray50", lwd = 0.2), 
		show_row_names = FALSE, show_column_names = FALSE, 
		#row_title = hmtitle, row_title_gp = gpar(fontsize = 2*label_font_size),
		cluster_rows = FALSE, cluster_columns = FALSE, 
		left_annotation = left_row_col_ha, right_annotation = right_row_col_ha, 
		bottom_annotation = bottom_column_col_ha, top_annotation = top_column_col_ha, 
		heatmap_legend_param = list(
			title = "expression\n", 
			border = TRUE,
			legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
			title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
			labels_gp = gpar(fontsize = label_font_size)
		)
	)
	
	if (!is.null(out_fn)) {
		
		png(paste0(out_bfn,"_boxplot.png"),h=500,w=2000)
		par(mar=c(15,7,2,2))
		boxplot(log2(csps$merged[sps1_selected_ids,]),las=2,pch=20,outcol=alpha("black",0.2),ylab="log2FC")
		abline(h=0)
		dev.off()
		
		height=max(length(gene_order) * height_cex,2000)
		width=ncol(csps$merged) * width_cex
		png(paste0(out_bfn,"_heatmap.png"),height=height,width=width)
		ht_opt(
			COLUMN_ANNO_PADDING=unit(5,"mm"), 
			ROW_ANNO_PADDING=unit(5,"mm"), 
			DIMNAME_PADDING=unit(5,"mm"),
			HEATMAP_LEGEND_PADDING=unit(10,"mm"),
			TITLE_PADDING=unit(10,"mm")
		)
		draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
		dev.off()
		
	}
	return(list(heatmap=hm, output_table=output_table))
}


#' Select one or more MCs in the first species and project expression 
#' of top genes across both species in comparison. 
#' 
#' @param csps cross-species comparison object as output by 
#'   `csps_create_crossspecies_object`
#' @param focus_mcs character, one or more metacell (or cell type) names (i.e. columns in `csps$merged`)
#' @param gene_annot_sp1_fn character, path to gene annotation file for the first species 
#'   in `csps`
#' @param sps1_tfs_fn character, path to the matrix of gene fold change
#'   in metacells (or cell types) for the first species in `csps` (usually metacell 
#'   footprint matrix, `mc@mc_fp`)
#' @param max_n_genes integer, number of top genes to project (default: 100)
#' @param fc_thrs numeric, fold change threshold (default: 2)
#' @param out_fn character, path to png file where expression heatmap will be saved
#'   (default: NULL)
#' @param width_cex numeric, figure height scalling
#' @param height_cex numeric, figure width scalling
#' @param res numeric (default: NA)
#' @param annotation_size numeric, height of the annotation color bar
#' @param label_font_size numeric, size of annotation labels
#' @param fc_max numeric, fc color scaling max value (default: 3)
#' 
#' @return list with following elements
#'   1) `heatmap` ComplexHeatmap object
#'   2) `expr_mat` matrix with the expression of top genes
#'   
csps_focused_coexpression=function(
	csps, focus_mcs, gene_annot_sp1_fn, sps1_tfs_fn, 
	max_n_genes=100, fc_thrs=2, 
	out_fn=NULL, width_cex=20, height_cex=20, res=NA,
	fc_max=3, annotation_size = 10, label_font_size=12
) {
	
	tfs_sp1=read.table(sps1_tfs_fn,h=TRUE,row.names=1,sep="\t",quote="",stringsAsFactors=FALSE,fill=TRUE)
	annot_sp1=read.table(gene_annot_sp1_fn,row.names=1,sep="\t",quote="",stringsAsFactors=FALSE,fill=TRUE,header=FALSE)
	
	ffall=sort(apply(csps$merged[,focus_mcs,drop=FALSE],1,median))
	#ff=tail(ffall[!grepl("\\.\\d",names(ffall))],max_n_genes)
	ff=tail(ffall,max_n_genes)
	focus_genes=names(which(ff > fc_thrs))
	
	shades2=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
	#label_col=ifelse(focus_genes %in% rownames(tfs_sp1),"red","black")
	hmtitle=sprintf("Top %s genes in %s",max_n_genes, focus_mcs)
	hmat <- pmin(csps$merged[rev(focus_genes),],fc_max)
	col_ord <- colnames(hmat)
	row_ord <- str_remove(rownames(hmat),"\\.\\d") #rownames(hmat)
	row_ord_1 <- annot_sp1[row_ord,1]
	row_ord_1[nchar(row_ord_1) > 30] <- paste0(substr(row_ord_1[nchar(row_ord_1) > 30],1,27),"...")
	row_ord_1 <- str_replace(row_ord_1,"\"\"","")
	row_ord_2 <- annot_sp1[row_ord,2]
	row_ord_2[nchar(row_ord_2) > 30] <- paste0(substr(row_ord_2[nchar(row_ord_2) > 30],1,27),"...")
	row_ord_2 <- str_replace(row_ord_2,"\"\"","")
	top_column_col_ha <- HeatmapAnnotation(
		which = "column",
		lab = anno_text(which = "column", just = "left", location = unit(0, "npc"), col_ord, gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	bottom_column_col_ha <- HeatmapAnnotation(
		which = "column",
		lab = anno_text(which = "column", col_ord, gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	left_row_col_ha <- HeatmapAnnotation(
		which = "row",
		lab = anno_text(which = "row", row_ord_2, just = "right", location = unit(1, "npc"), gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	right_row_col_ha <- HeatmapAnnotation(
		which = "row",
		lab = anno_text(which = "row", row_ord_1, gp = gpar(fontsize = label_font_size)),
		border = FALSE, show_annotation_name = FALSE, show_legend = FALSE
	)
	hm <- Heatmap(
		hmat, col = shades2, name = "expression", border = TRUE,
		rect_gp = gpar(col = "gray50", lwd = 0.2), 
		show_row_names = FALSE, show_column_names = FALSE, 
		row_title = hmtitle, row_title_gp = gpar(fontsize = 2 * label_font_size),
		cluster_rows = FALSE, cluster_columns = FALSE, 
		left_annotation = left_row_col_ha, right_annotation = right_row_col_ha, 
		bottom_annotation = bottom_column_col_ha, top_annotation = top_column_col_ha, 
		heatmap_legend_param = list(
			title = "expression\n", border = TRUE,
			legend_height = unit(6, "cm"), grid_width = unit(annotation_size,"mm"),
			title_position = "leftcenter-rot", title_gp = gpar(fontsize = label_font_size),
			labels_gp = gpar(fontsize = label_font_size)
		)
	)
	
	if (!is.null(out_fn)) {
		out_fn_prefix <- str_remove(out_fn,"\\.png")
		out_bfn=paste0(out_fn_prefix,"_",strtrim(paste(focus_mcs,collapse="_"),50),".png")
		height=max(length(focus_genes) * height_cex,2000)
		width=ncol(csps$merged) * width_cex
		
		png(out_bfn,h=height,w=width,res=NA)
		
		ht_opt(
			COLUMN_ANNO_PADDING=unit(5,"mm"), 
			ROW_ANNO_PADDING=unit(5,"mm"), 
			DIMNAME_PADDING=unit(5,"mm"),
			HEATMAP_LEGEND_PADDING=unit(10,"mm"),
			TITLE_PADDING=unit(10,"mm")
		)
		draw(hm, padding = unit(c(50, 50, 50, 50), "mm")) #bottom, left, top, right
		dev.off()
		
	}
	
	return(list(heatmap=hm, expr_mat=hmat))
}


# # # # # # # # # # # # # # # # #
#                               #
#    DOWNSTREAM FUNCTIONS       #
#    take matrix as input       #
#                               #
# # # # # # # # # # # # # # # # #

#' Identify genes with conserved expression in broad cell types of two species 
#' by calulating Jaccard index for expression in cell types
#' 
#' @param mc_fp_1,mc_fp_2 cell type gene expression matrix, with genes in rows 
#'   and cell types in columns. At least some of the column names should either 
#'   contain the pattern specified with `cts` argument, or start with common broad 
#'   cell type particle (four letters) - this will be used for overlap calculation.
#' @param cts pattern to use for matching cell typs; if NULL (default), broad cell type
#'   particle is used (i.e. the first four letters of expression matrix column names)
#' @param fc_thrs numeric, fold change threshold for genes selection
#' @param jacc_thrs numeric, Jaccard index threshold for genes selection
#' @return data.table with conserved cell types for each gene
#' 
ctConservedGenes <- function(mc_fp_1, mc_fp_2, cts=NULL, fc_thrs=1.7, jacc_thrs=0.5) {
	
	# jaccard distance for overlap in cell types
	.jacc_ct <- function(cell_types_1, cell_types_2) {
		cell_types_int <- intersect(cell_types_1,cell_types_2)
		cell_types_uni <- c(cell_types_1,cell_types_2)
		length_int <- length(cell_types_uni[cell_types_uni %in% cell_types_int])
		length_uni <- length(cell_types_uni)
		length_int / length_uni
	}
	
	# cts for all genes
	gen_ct_1 <- apply(mc_fp_1, 1, function(x) {
		id <- x > fc_thrs
		if (is.null(cts)) {
			cts <- substr(colnames(mc_fp_1)[id],1,4)
		} else {
			cts <- str_extract(colnames(mc_fp_1)[id],pattern=paste(cts,collapse="|"))
		}
		if (length(cts > 0)) {
			cts
		} else {
			NULL
		}
	})
	gen_ct_1 <- gen_ct_1[-which(lapply(gen_ct_1,is.null) == TRUE)]
	gen_ct_1 <- sapply(gen_ct_1, function(g) g[!is.na(g)], USE.NAMES=TRUE, simplify=FALSE)
	gen_ct_2 <- apply(mc_fp_2, 1, function(x) {
		id <- x > fc_thrs
		if (is.null(cts)) {
			cts <- substr(colnames(mc_fp_2)[id],1,4)
		} else {
			cts <- str_extract(colnames(mc_fp_2)[id],pattern=paste(cts,collapse="|"))
		}
		
		if (length(cts > 0)) {
			cts
		} else {
			NULL
		}
	})
	gen_ct_2 <- gen_ct_2[-which(lapply(gen_ct_2,is.null) == TRUE)]
	gen_ct_2 <- sapply(gen_ct_2, function(g) g[!is.na(g)], USE.NAMES=TRUE, simplify=FALSE)
	
	# common genes
	cg <- intersect(names(gen_ct_1),names(gen_ct_2))
	message(length(cg), " genes in intersect")
	
	# found common cts for genes
	ctj <- sapply(cg, function(g) {
		jacc <- .jacc_ct(gen_ct_1[[g]], gen_ct_2[[g]])
		if (length(jacc) > 0 & !is.na(jacc)) {
			if (jacc > jacc_thrs) {
				table(c(gen_ct_1[[g]],gen_ct_2[[g]]))
			} else {
				NULL
			}
		} else {
			NULL
		}
	}, USE.NAMES=TRUE, simplify=FALSE)
	nulls <- lapply(ctj,is.null) == TRUE
	if (any(nulls))
		ctj <- ctj[-which(nulls)]
	message(length(ctj))
	ctj_tbl <- table(sapply(ctj,length))
	message(sprintf(
		"Found %s conserved genes (Jaccard > 0.5); 
    \n%s genes with unique conserved ct;
    \n%s genes with non-unique conserved cts",
		length(ctj), ctj_tbl["1"], sum(ctj_tbl) - ctj_tbl["1"] 
	))
	dt <- rbindlist(lapply(ctj, function(x) {
		x <- as.matrix(unclass(x))
		data.table(ct=rownames(x),occurence=x[,1])
	}),idcol="gene")
	dt
}

#' Identify genes with conserved expression in selected columns of expression matrix 
#' (metacells, cell types) by expression fold change thresholding.
#' 
#' @param mc_fp cell type gene expression matrix, with genes in rows and cell types in columns
#' @param cols integer or character specifiying poisitions or names of columns to use
#' @param feature_in_thrs (default: 1.5)
#' @param feature_out_thrs (default: 1)
#' @param method,methodbg character specifying how to summarize gene expression in selected columns, 
#'   one of "absolute" (default) or "median"
#' @param abs_leakyness,abs_leakynessbg numeric, percent of selected columns or background columns 
#'   in which the expression value can be below or above the specified threshold, respectively;
#'   this is only used if `method`` or `methodbg`` is "absolute"; (default 0 and 0.05, respectively)
#' 
commonGenes <- function(
	mc_fp, cols, feature_in_thrs = 1.5, feature_out_thrs = 1,
	method="absolute", methodbg="absolute", abs_leakyness=0, abs_leakynessbg=0.05
) {
	
	# check that given cols are in the matrix and get in and out cols
	if (all(class(cols) == "integer")) {
		cols=colnames(mc_fp)[cols]
	} else if (all(class(cols) == "character")) {
		cols=intersect(cols,colnames(mc_fp))
	}
	cols_out <- setdiff(colnames(mc_fp),cols)
	
	# function to calculate leakyness
	.calc_leakyness <- function(abs_leakyness,n) {
		if (abs_leakyness < 1) {
			pmin(round(c(1 - abs_leakyness) * n), n)
		} else {
			pmin(round(n - abs_leakyness), n)
		}
	}
	# featues in and out
	if (method == "median") {
		fin=apply(mc_fp[,cols,drop=FALSE],1,median) > feature_in_thrs
		fin_inv=apply(mc_fp[,cols_out,drop=FALSE],1,median) > feature_in_thrs
	} else if (method == "absolute") {
		fin=apply(mc_fp[,cols,drop=FALSE], 1, function(x) 
			!(sum(x > feature_in_thrs) < .calc_leakyness(abs_leakyness,n=length(cols)))
		)
		fin_inv=apply(mc_fp[,cols_out,drop=FALSE], 1, function(x) 
			!(sum(x > feature_in_thrs) <  .calc_leakyness(abs_leakyness,n=length(cols_out)))
		)
	} else {
		stop("method should be either 'median' or 'absolute'")
	}
	if (methodbg == "median") {
		fout=apply(mc_fp[,cols_out,drop=FALSE],1,median) < feature_out_thrs
		fout_inv=apply(mc_fp[,cols,drop=FALSE],1,median) < feature_out_thrs
	} else if (methodbg == "absolute") {
		fout=apply(mc_fp[,cols_out,drop=FALSE], 1, function(x) 
			!(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(cols_out)))
		)
		fout_inv=apply(mc_fp[,cols,drop=FALSE], 1, function(x) 
			!(sum(x < feature_out_thrs) < .calc_leakyness(abs_leakynessbg,n=length(cols)))
		)
	} else {
		stop("methodbg should be either 'median' or 'absolute'")
	}
	
	f_in=which(fin & fout)
	f_out=which(fin_inv & fout_inv)
	
	list(genes_in=names(f_in), genes_out=names(f_out))
	
}


#' Plot violin plots of expression of selected genes in selected cell type, versus all other cell types
#' 
#' @param mc_fp expression matrix, rows are genes, columns are cell types with species prefix
#' @param ct one of colnames of mc_fp
#' @param species vector of species four-letter abbreviations
#' @param logical, whether cell type names are preceeded by species four-letter abbreviations
#' @param genes cell type specific genes 
#' @param scale.fc.max numeric, scale gene fc to this max value for plotting (default NULL)
#' @param sign.test function to calculate significance test, either `t.test` or `wilcox.test`
#' @param sign.label "pval","padj","p.adj.signif"
#' @param title logical, show cell type names as title
#' @param x.names logical, hide cell type names on x axis
#' @param x.names.replacement named character, optional replacement pattern for cell type labes on x axis
#' 
csps_plot_groupped_gene_expression <- function(
	mc_fp, ct, species=c("Spis","Nvec","Xesp","Hvul"), species_cell_type_names=TRUE, 
	genes, min.genes=3, scale.fc.max=NULL, sign=TRUE, sign.test=t.test, sign.label="p.adj.signif",
	col, title=TRUE, x.names=TRUE, x.names.angle=30, x.names.replacement=NULL,
	text.size=12, label.size=12, title.size=12, tip.length=0, step.increase=0.1
) {
	
	spreg=paste(species,collapse="|")
	ctl <- paste(ct,collapse=" ")
	message("Plotting expression for ", ctl)
	
	# transform matrix to data.table
	mc_fp_dt=melt.data.table(as.data.table(
		mc_fp,keep.rownames="gene"
	),id.vars="gene",variable.name="cell_type",value.name="fc")
	mc_fp_dt[,species:=str_extract(cell_type,spreg)]
	mc_fp_dt[,species:=factor(species,levels=unique(species))]
	if (!species_cell_type_names) 
		mc_fp_dt[,cell_type:=str_remove(cell_type,sprintf("(%s)_",spreg))]
	
	# genes to plot
	if (length(genes) > 0) {
		features_in_mc_fp=mc_fp_dt[gene %in% genes]
		features_in_mc_fp[,selected_cell_type:="other"]
		features_in_mc_fp[cell_type %in% ct, selected_cell_type:=cell_type]
		features_in_mc_fp[,selected_cell_type:=factor(selected_cell_type,levels=c(unique(ct),"other"))]
		uniqgenes <- min(unique(features_in_mc_fp[,.N,.(cell_type,species)]$N))
		if (is.infinite(uniqgenes)) uniqgenes <- 0
		if (uniqgenes > min.genes) {
			
			stattest=copy(features_in_mc_fp)
			
			# species must have at least two levels
			ss=stattest[,unique(.SD[,.(species,selected_cell_type)]),.(species)][,.N,species][N > 1]$species
			stattest=stattest[species %in% ss]
			
			# cell types to compare to "other"
			levels_ct <- levels(stattest$selected_cell_type)
			levels_compare <- levels_ct[levels_ct != "other"]
			
			dt <- stattest[,.(fc=mean(fc)),.(gene,selected_cell_type,species)]
			dt <- dt[,cell_type:=as.character(selected_cell_type)]
			st <- dt[,rstatix::pairwise_sign_test(data=.SD, formula = fc ~ cell_type, ref.group="other"),species]
			st[,species:=factor(species,levels=unique(species))]
			
			# plot
			if (!"other" %in% names(col)) {
				col <- append(col,"gray60")
				names(col)[length(col)] <- "other"
			}
			features_plot <- dt[species %in% ss]
			if (is.null(scale.fc.max)) {
				dt[, y:=max(fc) + max(fc) * 0.1,species ]
			}else {
				features_plot[,fc:=pmin(fc,scale.fc.max)]
				st[,y:=scale.fc.max + scale.fc.max * 0.1]
			}
			features_plot[,selected_cell_type_labels:=str_remove_all(selected_cell_type,sprintf("(%s)_*",spreg))]
			cell_types_labels <- features_plot$selected_cell_type_labels
			if (!is.null(x.names.replacement) & length(x.names.replacement) > 0)
				cell_types_labels <- str_replace_all(cell_types_labels,x.names.replacement)
			names(cell_types_labels) <- features_plot$selected_cell_type
			features_plot[,selected_cell_type:=droplevels(selected_cell_type)]
			selected_cell_type <- unique(features_plot$selected_cell_type)
			message("groups: ",paste(selected_cell_type,collapse=", "))
			message("colors: ",paste(col[selected_cell_type],collapse=", "))
			gp=ggplot(features_plot,aes(selected_cell_type,fc,color=selected_cell_type)) + 
				facet_grid(. ~ species,scales="free_x",space="free_x",switch="x",drop=TRUE) + #strip.position="bottom"
				geom_jitter(size=1, alpha=0.5, width=0.4) +
				geom_violin(alpha=1, lwd=0.2) +
				scale_color_manual(values=col) + scale_fill_manual(values=col) + 
				scale_x_discrete(labels=cell_types_labels) +
				geom_text(aes(y=fc * 1.2, label="")) + 
				theme(
					legend.position="none", text=element_text(size=text.size),
					axis.line.x=element_blank(), #axis.ticks.x=element_blank(), 
					axis.text.x=element_text(angle=x.names.angle,hjust=1), # vjust=0.5 if angle=90
					plot.title=element_text(size=title.size),
					strip.background=element_blank(), strip.placement="outside",
					panel.border=element_blank()
				) + 
				labs(x=NULL,y="gene fold change",title=ctl)
			if (sign) {
				gp <- gp + stat_pvalue_manual(
					st, label=sign.label, y.position="y",
					step.increase=step.increase, step.group.by="species",
					label.size = label.size / 3.88, tip.length=tip.length
				)
			}
			if (!title) {
				gp <- gp + theme(plot.title=element_blank())
			}
			if (!x.names)
				gp <- gp + theme(axis.text.x=element_blank())
		} else {
			gp <- NULL
			message(sprintf("Not enough genes (%s provided, more than %s required)", uniqgenes, min.genes))
		}
	} else {
		gp <- NULL
		message(sprintf("Not enough genes in expression data (%s provided, more than %s required)", length(genes), min.genes))
	}
	return(gp)
}


#' Plot cross-species chord diagram
#' 
#' @param mat matrix with similarities to plot on the circos; colnames and 
#'   rownames should start with four-letter species abbreviations
#' @param sp character, four-letter abbreviations of the species for which to 
#'   plot the circos, (at least some of) `colnames(mat)` and `rownames(mat)` 
#'   should start with these
#' @param revert logical, wheter to use `1/mat` for plotting links, set this to 
#'   TRUE if the values in matrix are a measure of divergence, and to FALSE if 
#'   they measure similarity
#' @param threshold numeric, between 0 and 1, a threshold to use for selecting 
#'   links to be shown on circos plot - the values in the matrix ae scaled to 
#'   [0,1] range and only those with scaled value above specified threshold are 
#'   plotted
#' @param name.suffix character, optional suffix to be appended to the plot
#'   filename which by default will include species names
#' @param outdir character, directory in which the plot will be saved
#' @param width numeric, width of the plot in inches
#' @param height numeric, height of the plot in inches
#' @param sectors.order charcter vector with names of sectors to be plotted 
#'   on the circos, should be in `c(colnames(mat)`, `rownames(mat))`
#' @param sectors.labels named character, optional labels for sectors, the names
#'   should be `sectors.order`
#' @param sectors.colors named character, optional colors for sectors, the names 
#'   should be `sectors.order`
#' @param sectors.groups named factor or character, indicates grouping of sectors, 
#'   levels of factor determine ordering on the circos, and the names should be 
#'   `sectors.order`
#' @param sectors.groups.labels named character, optional labels for groups of 
#'   sectors, the names should be `unique(sectors.groups)`
#' @param start.degree numeric
#' @param annotation.track logical
#' @param annotation.names logical
#' 
csps_plot_chord_diagram <- function(
	mat, sp, revert=FALSE, threshold=NULL, threshold.quantile=NULL,
	scale.values=TRUE, regularize.values=TRUE, reg.param=1,
	save.plot=TRUE, save.data=TRUE, name.suffix="", outdir=".", width=18, height=18, mar=c(10,5,10,5),
	sectors.order=NULL, sectors.labels=NULL, sectors.colors=NULL, grid.border=NULL,
	sectors.groups=NULL, sectors.groups.labels=NULL,
	sector.scale=TRUE, sector.width.link=FALSE, remove.empty.sectors=FALSE, self.links=FALSE, 
	sector.groups.padding=c(0.1,0,1.8,0), sector.groups.col="gray98", sector.groups.border="gray88", 
	start.degree=360, small.gap=0, big.gap=5, annotation.track=TRUE, annotation.names=TRUE,
	sector.text.cex=1, sector.groups.text.cex=2, sector.groups.text.vjust=-6,
	title.main=NULL, title.sub=NULL, title.main.cex=2, title.sub.cex=2
) {
	# functions
	col_fun = function(x,min,max) {
		f=circlize::colorRamp2(c(min,max), c("white","red"))
		f(x)
	}
	transparency_fun = function(x,n=4,scale=TRUE,pseudocount=0.0001) {
		if (scale) {
			x=(x - min(x)) / (max(x) - min(x))
		}
		x= -n * log(x + pseudocount)
		(x - min(x)) / (max(x) - min(x))
	}
	highlightSector = function(
		group.name,group,track.index=1,
		col="gray98",border="black",padding=c(0,0,0,0),
		text=group.name,text.col="black",text.cex=1,text.vjust=0.5,niceFacing=TRUE
	) {
		highlight.sector(
			sector.index=names(group[group == group.name]), track.index=1, 
			col=col, border=border, padding=padding,
			text=text, text.vjust=text.vjust, text.col=text.col, cex=text.cex,
			niceFacing=niceFacing, facing="bending.inside"
		)
	}
	# data
	if (is.null(sectors.order))
		sectors.order=colnames(mat)
	if (is.null(sectors.labels)) {
		sectors.labels=sectors.order
		names(sectors.labels)=sectors.order
	}
	if (is.null(sectors.colors)) {
		sectors.colors=rand_color(length(sectors.labels))
		names(sectors.colors)=sectors.order
	}
	dt=as.data.table(mat,keep.rownames="ct1")
	dtm=melt.data.table(dt,id.vars="ct1",variable.name="ct2",value.name="value")
	dtm[,ct1:=factor(ct1,levels=sectors.order[sectors.order %in% ct1])]
	dtm[,ct2:=factor(ct2,levels=sectors.order[sectors.order %in% ct2])]
	dtm[,sp1:=stringr::str_extract(ct1,"[A-z]{4}")]
	dtm[,sp2:=stringr::str_extract(ct2,"[A-z]{4}")]
	if (revert) {
		dtm[,value := 1 / value]
		dtm[is.infinite(value),value:=0]
	}
	#setcolorder(dtm,c("ct1","ct2","value","sp1","sp2"))
	dtm[,linkvalue:=value]
	if (regularize.values == TRUE) {
		regdt=dtm[sp1 == sp2][ct1 != ct2][,.(rg=max(value)),ct1]#[,.(rg=max(value)),.(sp1,sp2)]
		dtm[regdt,on="ct1",rg:=i.rg]
		dtm[,linkvalue:=linkvalue / (reg.param * rg)]
	}
	if (scale.values) 
		dtm[,linkvalue:=linkvalue / max(linkvalue)]
	if (self.links) {
		dtms=dtm[sp1 == sp2]
	} else {
		dtms=NULL
	}
	if (length(sp) > 4) {
		# use the first group as a reference and plot only links from it to others
		dtms=rbindlist(list(
			dtms,
			dtm[grepl(sp[1],dtm[["ct1"]])][grepl(paste(sp[2:length(sp)],collapse="|"),dtm[["ct2"]])]
		))
	} else {
		if (!is.na(sp[2])) {
			# plot links: 1-2
			dtms=rbindlist(list(
				dtms,dtm[grepl(sp[1],dtm[["ct1"]]) & grepl(sp[2],dtm[["ct2"]])]
			))
		} 
		if (!is.na(sp[3])) {
			# plot links: 1-2, 1-3, 2-3
			dtms=rbindlist(list(
				dtms,
				dtm[grepl(sp[1],dtm[["ct1"]]) & grepl(sp[3],dtm[["ct2"]])],
				dtm[grepl(sp[2],dtm[["ct1"]]) & grepl(sp[3],dtm[["ct2"]])]
			))
		}
		if (!is.na(sp[4])) {
			# plot links: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
			dtms=rbindlist(list(
				dtms,
				dtm[grepl(sp[2],dtm[["ct1"]]) & grepl(sp[4],dtm[["ct2"]])],
				dtm[grepl(sp[3],dtm[["ct1"]]) & grepl(sp[4],dtm[["ct2"]])]
			))
		}
	}
	dtcd=dtms[,.(ct1,ct2,linkvalue)][order(ct1,ct2)]
	# select links above threshold
	if (!any(!is.null(threshold),!is.null(threshold.quantile)) |
		all(is.null(threshold),is.null(threshold.quantile))) {
		stop("You must specify either threshold or threshold.quantile!")
	}
	if (!is.null(threshold.quantile)) {
		threshold=quantile(dtcd$linkvalue,threshold.quantile)
	}
	dtcds=dtcd[linkvalue > threshold]
	# add empty space at the link origin for links below threshold
	dtcds=rbindlist(list(
		dtcds,
		dtcd[!ct1 %in% dtcds$ct1][,.SD[order(linkvalue,decreasing=TRUE)][1],ct1][,.(ct1,ct2,linkvalue)]
	))
	# add empty space at the link destination for links below threshold
	if (remove.empty.sectors == FALSE)
		dtcall=rbindlist(list(
			dtcds,
			dtcd[!ct2 %in% dtcds$ct2][,.SD[order(linkvalue,decreasing=TRUE)][1],ct2][,.(ct1,ct2,linkvalue)]
		))
	dtcall=dtcall[order(ct1,ct2)]
	dtcall[,ct1:=as.character(ct1)][,ct2:=as.character(ct2)]
	# link color by link value
	# dtcall[,col:=col_fun(linkvalue,min(linkvalue),max(linkvalue))]
	# link color by originating sector, hide links below threshold
	dtcall[,col:=sectors.colors[as.character(ct1)]]
	dtcall[!(linkvalue > threshold),col:="white"]
	# remove below threshold links (i.e. make above threshold links the same width as sector)
	if (sector.width.link) {
		single_link_cts=dtcall[,.N,ct1][N == 1]$ct1
		empty_cts=dtcall[ct1 %in% single_link_cts][col == "white"]$ct1[1]
		dtcall[col == "white",ct1:=empty_cts]
	}
	# remove sectors that only have below threshold links
	if (remove.empty.sectors) {
		multiple_link_cts=dtcall[,.N,ct1][N > 1]$ct1
		dtcall=dtcall[!(ct1 %in% multiple_link_cts & col == "white")]
	}
	# set transparency
	dtcall[,transparency:=transparency_fun(linkvalue)]
	dtcall[!(linkvalue > threshold),transparency:=1]
	# set order
	dtcall[,zindex:=rank(linkvalue)]
	setcolorder(dtcall,c("ct1","ct2","linkvalue","col","transparency","zindex"))
	if (save.data)
		saveRDS(
			dtcall, 
			file.path(outdir,sprintf("circos_%s_%s.RDS",paste0(sp,collapse="_"),name.suffix))
		)
	
	# circos
	if (save.plot) {
		pdf(
			file.path(outdir,sprintf("circos_%s_%s.pdf",paste0(sp,collapse="_"),name.suffix)),
			w=width,h=height,useDingbats=TRUE
		)
		par(mar=mar)
	}
	circos.par(
		"canvas.xlim"=c(-2, 2),"canvas.ylim"=c(-2, 2),
		"track.height"=0.2, "track.margin"=c(0,0), cell.padding=c(0,1,0,1),
		"start.degree"=start.degree
	)
	annotationTrack=ifelse(annotation.track,"grid",NULL)
	chordDiagram(
		dtcall, grid.col=sectors.colors, grid.border=grid.border,
		annotationTrack=annotationTrack, small.gap=small.gap, big.gap=big.gap, 
		directional=2, #direction.type=c("arrows"), link.arr.type = "big.arrow", 
		#annotationTrackHeight = c(0.05,-0.05), 
		scale=sector.scale, col=dtcall$col, transparency=dtcall$transparency + 0.1,
		link.zindex=dtcall$zindex, group=sectors.groups, preAllocateTracks=1,
	)
	groupchord=sectors.groups[names(sectors.groups) %in% c(as.character(dtcall$ct1),as.character(dtcall$ct2))]
	for (group.name in unique(c(dtm$sp1,dtm$sp2))) {
		if (!is.null(sectors.groups.labels)) {
			sectorlabs=sectors.groups.labels[group.name]
		} else {
			sectorlabs=group.name
		}
		highlightSector(
			group.name=group.name, group=groupchord,
			padding=sector.groups.padding, col=sector.groups.col, border=sector.groups.border,
			text=sectorlabs, text.vjust=sector.groups.text.vjust, text.cex=sector.groups.text.cex
		)
	}
	circos.trackPlotRegion(track.index=1, panel.fun = function(x, y) {
		xlim = get.cell.meta.data("xlim")
		ylim = get.cell.meta.data("ylim")
		sector.name = get.cell.meta.data("sector.index")
		if (annotation.names) {
			sector.text=sectors.labels[sector.name]
		} else {
			sector.text=""
		}
		circos.text(
			mean(xlim), ylim[1] + 0.1, sector.text, 
			facing="clockwise", niceFacing=TRUE, adj=c(0,0.5), 
			col=sectors.colors[sector.name], 
			cex = sector.text.cex
		)
	}, bg.border = NA)
	circos.clear()
	if (any(!is.null(title.main),!is.null(title.sub))) {
		title(
			main=title.main, sub=title.sub,
			cex.main=title.main.cex, cex.sub=title.sub.cex
		)
	}
	if (save.plot) 
		dev.off()
	return(dtcall)
}




#' Obtain one-to-one orthology pairs based on orthology pairs file and expression correlation
#' 
#' @param sp1_fp_fn,sp2_fp_fn either a path to a gene expression matrix (rows are genes, columns are metacells/cell types/etc) or one such matrix
#' @param OG_pairs_fn either a path to a ortholog pairs table, or one such table. Requires two columns: first are genes from species 1 (`sp1_fp_fn`), second are genes from species 2 (`sp2_fp_fn`)
#' @param cor_method correlation method to identify co-correlation of genes across species (any `cor` method is valid); default is pearson
#' 
#' @return data.frame with one-to-one ortholog pairs. Column 1 and 2 are genes from species 1 and 2, and column 3 indicates whether this pair is a true one-to-one ortholog (according to input table) or whether it has been added by discarding paralogs of uncorrelated expression.
#' 
csps_ortholog_o2o_pairs_from_expression = function(
	sp1_fp_fn, sp2_fp_fn, OG_pairs_fn, cor_method="pearson"
) {

	# Read and parse expression data  
	# species 1
	extension = rev(stringr::str_split(sp1_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp1_fp_fn)) {
		sp1_fp = sp1_fp_fn
	} else if ("character" %in% class(sp1_fp_fn) & extension == ".RDS") {
		sp1_fp=readRDS(sp1_fp_fn)
	} else if ("character" %in% class(sp1_fp_fn)) {
		sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp1_fp_fn is neither a matrix object nor a path to such an object, I quit.")    
	}
	# species 2
	extension = rev(stringr::str_split(sp2_fp_fn, "\\.")[[1]])[1]
	if ("matrix" %in% class(sp2_fp_fn)) {
		sp2_fp = sp2_fp_fn
	} else if ("character" %in% class(sp2_fp_fn) & extension == ".RDS") {
		sp2_fp=readRDS(sp2_fp_fn)
	} else if ("character" %in% class(sp2_fp_fn)) {
		sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
	} else {
		stop("sp2_fp_fn is neither a matrix object nor a path to such an object, I quit.")    
	}

	# load orthologous pairs
	if ("data.frame" %in% class(OG_pairs_fn) | "matrix" %in% class(OG_pairs_fn)) {
		og_pairs = OG_pairs_fn
	} else if ("character" %in% class(sp2_fp_fn)) {
		og_pairs = data.table::fread(OG_pairs_fn, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	}	
	colnames(og_pairs) = c("sp1", "sp2")
	
	# retrieve truly one2one pairs
	og_pairs_o2o = og_pairs [ !duplicated(og_pairs$sp1) & !duplicated(og_pairs$sp2), ]
	# among these, keep pairs that are also expressed (to use as reference for expression similarity)
	og_pairs_ref = og_pairs_o2o [ 
		og_pairs_o2o$sp1 %in% rownames(sp1_fp) & og_pairs_o2o$sp2 %in% rownames(sp2_fp) ,
	]
	# list of paralogous genes from species 1 that have one (or more) orthologs in sp2; and viceversa
	list_m2o_genes_sp1 = unique(og_pairs$sp1 [ table(og_pairs$sp1) > 1 ])
	list_m2o_genes_sp2 = unique(og_pairs$sp2 [ table(og_pairs$sp2) > 1 ])
	# restrict to genes that are expressed
	list_m2o_genes_sp1 = list_m2o_genes_sp1 [ list_m2o_genes_sp1 %in% rownames(sp1_fp) ]
	list_m2o_genes_sp2 = list_m2o_genes_sp2 [ list_m2o_genes_sp2 %in% rownames(sp2_fp) ]
	# restrict to genes that have expressed orthologs AND whose orthologs are not in the one2one table
	has_expressed_og1 = unlist( lapply( list_m2o_genes_sp1, function(g) any(og_pairs[og_pairs$sp1 == g, ]$sp2 %in% rownames(sp2_fp) & !og_pairs[og_pairs$sp1 == g, ]$sp2 %in% og_pairs_o2o$sp2  ) ) )
	has_expressed_og2 = unlist( lapply( list_m2o_genes_sp2, function(g) any(og_pairs[og_pairs$sp2 == g, ]$sp1 %in% rownames(sp1_fp) & !og_pairs[og_pairs$sp2 == g, ]$sp1 %in% og_pairs_o2o$sp1  ) ) )
	list_m2o_genes_sp1 = list_m2o_genes_sp1 [ has_expressed_og1 ]
	list_m2o_genes_sp2 = list_m2o_genes_sp2 [ has_expressed_og2 ]

	# first, correlation of expression between all non-unique sp1 genes and the one2one reference set
	ge1_exp_nonuniq = t(sp1_fp[list_m2o_genes_sp1,])
	ge1_exp_refogps = t(sp1_fp[og_pairs_ref$sp1,])
	ge1_cor_m = cor(ge1_exp_nonuniq, ge1_exp_refogps, method=cor_method)

	# second, same correlation matrices for sp2
	ge2_exp_nonuniq = t(sp2_fp[list_m2o_genes_sp2,])
	ge2_exp_refogps = t(sp2_fp[ rownames(sp2_fp) %in% og_pairs_ref$sp2, ])
	ge2_cor_m = cor(ge2_exp_nonuniq, ge2_exp_refogps, method=cor_method)

	# loop through unmatched species 1 genes, trying to find best paralog among candidates
	# looping function
	find_best_paralog = function(gei, spi = "sp1", spj="sp2", spi_cor_m=ge1_cor_m, spj_cor_ref=ge2_exp_refogps, spj_exp_m=sp2_fp  ) {

		# message(gei)
		# species 1 correlation vector
		gei_cor_v = spi_cor_m[gei,]

		# correlation matrix between candidates and sp2 reference
		gej_candidates = as.data.frame(og_pairs) [ as.data.frame(og_pairs)[,spi] == gei, spj]
		gej_candidates = gej_candidates [ gej_candidates %in% rownames(spj_exp_m) ]
		gej_exp_candidates = t(spj_exp_m[ rownames(spj_exp_m) %in% gej_candidates, ])

		if (length(gej_candidates) > 1) {
			
			gej_cor_m = cor(
				gej_exp_candidates,
				spj_cor_ref,
				method=cor_method)

			# correlation matrix between candidates and sp2 vector
			gep_cor_m = cor(
				gei_cor_v,
				t(gej_cor_m),
				method=cor_method)
			
			# keep top gene	
			gep_cor_v = as.vector(gep_cor_m)
			names(gep_cor_v) = colnames(gep_cor_m)
			gep_cor_v = sort(gep_cor_v, decreasing = TRUE)
			gej_keep = names(gep_cor_v[1])

		} else {

			gej_keep = gej_candidates

		}

		return(gej_keep)
	}

	# loop around sp1
	og_pairs_ge1_o2m = data.frame(sp1 = list_m2o_genes_sp1, sp2 = NA)
	og_pairs_ge1_o2m$sp2 = unlist(lapply(
		og_pairs_ge1_o2m$sp1, 
		function(g) find_best_paralog(gei = g, spi = "sp1", spj="sp2", spi_cor_m = ge1_cor_m, spj_cor_ref = ge2_exp_refogps, spj_exp_m = sp2_fp)
	))
	# loop around sp2
	og_pairs_ge2_o2m = data.frame(sp1 = NA, sp2 = list_m2o_genes_sp2)
	og_pairs_ge2_o2m$sp1 = unlist(lapply(
		og_pairs_ge2_o2m$sp2, 
		function(g) find_best_paralog(gei = g, spi = "sp2", spj="sp1", spi_cor_m = ge2_cor_m, spj_cor_ref = ge1_exp_refogps, spj_exp_m = sp1_fp)
	))
	# concatenate
	og_pairs_extra_o2m = rbind(og_pairs_ge1_o2m, og_pairs_ge2_o2m)
	# remove pairs that are not reciprocal and unique
	og_pairs_extra_o2m_filt = og_pairs_extra_o2m [ !duplicated(og_pairs_extra_o2m$sp1) & !duplicated(og_pairs_extra_o2m$sp2) , ]

	# add to original o2o table
	og_pairs_clean = rbind(og_pairs_o2o, og_pairs_extra_o2m_filt)
	og_pairs_clean$is_true_o2o = paste(og_pairs_clean$sp1, og_pairs_clean$sp2) %in% paste(og_pairs_o2o$sp1, og_pairs_o2o$sp2)
	og_pairs_clean = og_pairs_clean [ !duplicated(og_pairs_clean$sp1) & !duplicated(og_pairs_clean$sp2) , ]

	# return pairs of orthologs (one2one if possible; and among those with paralogs, the best paralog best on expression)
	return(og_pairs_clean)

}


#' Calculate expression conservation (EC) score for pairs of orthologous genes in a cross-species iterative comparison of coexpression (ICC)
#'
#' @param mat_sp1 expression matrix of species 1 (rows are genes, columns are any conditions)
#' @param mat_sp2 expression matrix of species 2 (rows are genes, columns are any conditions)
#' @param og_pairs either a path to a ortholog pairs table, or one such table. Requires two columns: first are genes from species 1, 
#'    second are genes from species 2
#' @param niter maximum number of ICC iterations (default = 100, min is 2)
#' @param icc_thr similarity threshold to stop ICC iteration (ICC stops when difference between iteration is below this value; default = 0.05)
#' @param method ICC correlation method (default is `pearson`)
#' @param do_duplicates bool (default `FALSE`): if `TRUE`, calculate EC values for best reciprocal pairs of duplicate genes in addition 
#'    to the one-to-one ortholog pairs (which is the default)
#' @param num_o2o_ref_genes if `do_duplicates` is TRUE, select up to this number of genes for the duplicate EC scoring step. Less than 1000 is not
#'    recommended (very unstable EC values).
#' @param use_variable_o2o_genes bool; if `TRUE` (default), use only one-to-one orthologs with variable expression in both species as a reference 
#'    matrix for the duplicates' EC calculations. Variable genes will be computed from footprint matrices (default: top quartile of max footprint
#'    across cells). If input matrices are not footprint matrices, you can provide pre-computed variable genes with `vargenes_sp1` and `vargenes_sp1`
#'    instead (see below).
#' @param variable_o2o_thr if use_variable_o2o_genes is `TRUE`, use this quantile threshold to select variable genes based on the distribution of 
#'    max observed footprint across all cells (default is 0.75, i.e. the top quartile of variable genes). Increasing the size of the reference o2o
#'    matrix can result in extremely long runtimes. Decreasing it can result in poor EC score estimation for paralogs. Rough benchmarking: thr=0.75 
#'    is fairly accurate and fast (~1-2h, cor ~0.80 with whole-matrix result. A thr=0.5 is more accurate but slower (~12h, cor ~0.95).
#' @param vargenes_sp1, vargenes_sp2 (default is `NULL`): vectors of variable one-to-one orthologous genes to use as reference for the 
#'    `do_duplicates` step. Only applicable if `do_duplicates` is `TRUE`. 
#' @param nthreads_icc number of threads to use for ICC calculations (i.e. parallelisation of the correlation matrix).
#' @param nthreads_dup number of parallel processes of duplicate scoring (each of them will run with n=`nthreads_icc` threads).
#' @param cross_fc_thrs numeric, fold change threshold (default: 2, set to NULL to skip filtering by fc)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` we require in each species (default: 1, set to NULL to skip filtering by fc)
#'
#' @return a list with the following: a new csps object using ICC pairs as markers, a dataframe of unique gene pairs from species 1 and 2 with their cross-species expression conservation score, and a dataframe with paralog genes (non-unique pairs).
#'
csps_markers_icc = function(
	mat_sp1,
	mat_sp2, 
	og_pairs, 
	niter = 100, 
	icc_thr = 0.05, 
	method = "pearson", 
	do_duplicates = FALSE, 
	num_o2o_ref_genes = NULL,
	use_variable_o2o_genes = TRUE, 
	variable_o2o_thr = 0.75, 
	vargenes_sp1 = NULL, 
	vargenes_sp2 = NULL,
	cross_fc_thrs = 2, 
	do_quantile_normalisation = TRUE,
	cross_n = 1,
	nthreads_icc = 2,
	nthreads_dup = 1) {
	
	require("igraph")
	require("data.table")
	
	# register cores
	registerDoParallel(cores = nthreads_dup)

	# ensure data are matrices
	mat_sp1 = as.matrix(mat_sp1)
	mat_sp2 = as.matrix(mat_sp2)
	
	# load orthologous pairs
	if ("character" %in% class(og_pairs)) {
		og_pairs = data.table::fread(og_pairs, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	}	
	colnames(og_pairs) = c("sp1", "sp2")

	# restrict pairs to expressed genes
	og_pairs = og_pairs [ og_pairs[,1] %in% rownames(mat_sp1) & og_pairs[,2] %in% rownames(mat_sp2) , ]
	
	### Find ec values of one-to-one orthologs
	# get one to one pairs (reference)
	list_o2o_sp1 = names(which(table(og_pairs[,1]) == 1))
	list_o2o_sp2 = names(which(table(og_pairs[,2]) == 1))
	bool_o2o = og_pairs[,1] %in% list_o2o_sp1 & og_pairs[,2] %in% list_o2o_sp2
	og_pairs_o2o = og_pairs [ bool_o2o, ]

	# reference o2o matrix
	mar_sp1 = mat_sp1 [ og_pairs_o2o[,1] , ]
	mar_sp2 = mat_sp2 [ og_pairs_o2o[,2] , ]
	
	# ec values for the one to one orthologs
	message(sprintf("ICC markers | pairs of one-to-one orthologs = %i", nrow(og_pairs_o2o)))
	ecv_o2o = csps_calc_icc(mat_sp1 = mar_sp1, mat_sp2 = mar_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = TRUE, num_cores = nthreads_icc)
	# ecv_o2o = ecv_o2o_list[[1]]
	ecv_o2o [ ecv_o2o$ec_value < 0 | is.na(ecv_o2o$ec_value) , "ec_value" ] = 0
	ecv_o2o$is_o2o = 1
	
	### Find best reciprocal duplicate pairs (highest ec value)
	# loop through pairs of non-o2o homologs, add them to the original o2o matrix, and calculate ec value
	if (do_duplicates) { 
		
		# select genes for the reference o2o matrix
		if (use_variable_o2o_genes & is.null(vargenes_sp1) & is.null(vargenes_sp2)) {

			# o2o orthologs with variable expression
			max_fcs_sp1 = apply(mat_sp1, 1, max, na.rm = TRUE)
			max_fcs_sp2 = apply(mat_sp2, 1, max, na.rm = TRUE)
			vargenes_sp1 = rownames(mat_sp1) [ max_fcs_sp1 > quantile(max_fcs_sp1, variable_o2o_thr) ]
			vargenes_sp2 = rownames(mat_sp2) [ max_fcs_sp2 > quantile(max_fcs_sp2, variable_o2o_thr) ]
			# variable AND o2o orthologs
			ecv_o2o_ref = ecv_o2o [ ecv_o2o[,1] %in% vargenes_sp1 & ecv_o2o[,2] %in% vargenes_sp2 , ]
			refgenes_sp1 = ecv_o2o_ref[,1]
			refgenes_sp2 = ecv_o2o_ref[,2]
		
		} else if (!is.null(vargenes_sp1) & !is.null(vargenes_sp2)) { 
			
			ecv_o2o_ref = ecv_o2o [ ecv_o2o[,1] %in% vargenes_sp1 & ecv_o2o[,2] %in% vargenes_sp2 , ]
			refgenes_sp1 = ecv_o2o_ref[,1]
			refgenes_sp2 = ecv_o2o_ref[,2]
			
		} else if (is.null(num_o2o_ref_genes)) { 
			
			refgenes_sp1 = ecv_o2o[,1]
			refgenes_sp2 = ecv_o2o[,2]
			
		} else {
			
			ixs_ref = sample(1:nrow(ecv_o2o), min(nrow(ecv_o2o), num_o2o_ref_genes), replace = FALSE)
			refgenes_sp1 = ecv_o2o[ixs_ref,1]
			refgenes_sp2 = ecv_o2o[ixs_ref,2]
			
		}
		
		# sanity check: do we have enough genes in the reference matrix?
		message(sprintf("ICC markers | reference one-to-one orthologs for duplicate EC scoring = %i", length(refgenes_sp1)))
		if (length(refgenes_sp1) == 0) { 
			stop(sprintf("Only %i pairs of one-to-one orthologs in the reference matrix for duplicate EC calculations. I can't do this.", length(refgenes_sp1)))
		} else if (length(refgenes_sp1) < 1000 & length(refgenes_sp1) > 0) {
			warning(sprintf("Only %i pairs of one-to-one orthologs in the reference matrix for duplicate EC calculations. Too few?", length(refgenes_sp1)))
		}
		
		# get duplicates (all types)
		list_dup_sp1 = names(which(table(og_pairs[,1]) > 1))
		list_dup_sp2 = names(which(table(og_pairs[,2]) > 1))
		bool_dup = og_pairs[,1] %in% list_dup_sp1 | og_pairs[,2] %in% list_dup_sp2
		og_pairs_dup = og_pairs [ bool_dup , ]
		og_pairs_dup = og_pairs_dup[ order(og_pairs_dup[,1], og_pairs_dup[,2]) , ]
		
		# graph of duplicates
		og_pairs_dup_sps = as.matrix(og_pairs_dup)
		og_pairs_dup_sps[,1] = paste("sp1", og_pairs_dup_sps[,1])
		og_pairs_dup_sps[,2] = paste("sp2", og_pairs_dup_sps[,2])
		graph_dup = igraph::graph_from_edgelist(og_pairs_dup_sps, directed = TRUE)
		graph_dup_components = igraph::components(graph_dup)
		
		# dict of duplicates
		dict_dups = list()
		list_components = unique(graph_dup_components$membership)
		for (component in list_components) {
			dups = names(graph_dup_components$membership) [graph_dup_components$membership == component]
			dict_dups[[component]] = dups
		}

		if (length(dict_dups) > 0) {

			# log progress		
			message(sprintf("ICC markers | clusters with non-o2o homologs = %i (total = %i genes)", length(dict_dups), length(unique(unlist(dict_dups)))))
			
			# loop through duplicates in dict_dups
			# this loop is a function for easy parallelisation with plyr::adply
			loop_dups = function(i) {
				
				if (i %% (as.integer(length(dict_dups)/100)) == 0 | i == length(dict_dups)) {
					message(sprintf("ICC markers | clusters with non-o2o homologs = %i (%i/%i processed)", length(dict_dups), i, length(dict_dups)))
				}
				# get lists of genes from sp1 and sp2
				d = dict_dups[[i]]
				dg = stringr::word(d, 2)
				genes1 = dg [ stringr::word(d, 1) == "sp1" ]
				genes2 = dg [ stringr::word(d, 1) == "sp2" ]
				
				# prepare vector of ec values for this set of paralogs
				gen_dups = og_pairs [ og_pairs[,1] %in% genes1 & og_pairs[,2] %in% genes2 , ]
					
				# add test genes to reference expression matrix, sp1
				mai_sp1 = rbind( mat_sp1[ gen_dups[,1] ,], mar_sp1[refgenes_sp1,] )
				rownames(mai_sp1)[ 1:nrow(gen_dups) ] = gen_dups[,1]
				# add test genes to reference expression matrix, sp2
				mai_sp2 = rbind( mat_sp2[ gen_dups[,2] ,], mar_sp2[refgenes_sp2,] )
				rownames(mai_sp2)[ 1:nrow(gen_dups) ] = gen_dups[,2]

				# recalculate ec values using reference expression matrix and new set of duplicates
				ecv_dup_i = csps_calc_icc(mat_sp1 = mai_sp1, mat_sp2 = mai_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = FALSE, num_cores = nthreads_icc)
				ecv_dup_q = ecv_dup_i [ ecv_dup_i$sp1 %in% gen_dups[,1] & ecv_dup_i$sp2 %in% gen_dups[,2] , ]
						
				return(ecv_dup_q)
				
			}
			
			# apply loop (not parallelised internally)
			ecv_dup = plyr::adply(.data = 1:length(dict_dups), .margins = 1,  .fun = loop_dups, .parallel = TRUE, .id = "cluster", .progress = "none")
			ecv_dup = ecv_dup [ , c("sp1", "sp2", "ec_value", "cluster") ]
			
			# keep best reciprocal pairs of orthologous genes
			ecv_dup_t = data.table::as.data.table(ecv_dup)
			ecv_dup_t$ec_value [ is.na(ecv_dup_t$ec_value) ] = -1
			ecv_dup_f1 = ecv_dup_t[ ecv_dup_t[, .I[ ec_value == max(ec_value) ], by=c("sp1","cluster")]$V1 ]
			ecv_dup_f2 = ecv_dup_t[ ecv_dup_t[, .I[ ec_value == max(ec_value) ], by=c("sp2","cluster")]$V1 ]
			ecv_dup_f  = rbind(ecv_dup_f1,ecv_dup_f2)
			ecv_dup_f  = ecv_dup_f [ duplicated(ecv_dup_f) & ecv_dup_f$ec_value >= 0 , ]
			
			# break ties at random (so as to obtain a truly one-to-one table)
			ecv_dup_f = ecv_dup_f [ !duplicated(ecv_dup_f$sp1), ]
			ecv_dup_f = ecv_dup_f [ !duplicated(ecv_dup_f$sp2), ]
			
			# log progress		
			message(sprintf("ICC markers | pairs of non-o2o homologs kept = %i", nrow(ecv_dup_f)))
			
			# output contains both o2o orthologs & best reciprocal pairs of duplicated genes
			ecv_dup_f = ecv_dup_f [ , c("sp1","sp2","ec_value") ]
			ecv_dup_f$is_o2o = 0
			ecv_out = rbind(ecv_o2o, ecv_dup_f)
		
		} else{
		
			# log progress		
			message(sprintf("ICC markers | clusters with non-o2o homologs = %i (total = %i genes) | Can't score duplicates!", length(dict_dups), length(unique(unlist(dict_dups)))))
			# output will only contain one-to-one orthologs
			ecv_out = ecv_o2o
			ecv_dup = NULL
			
		}

		
	} else {
		
		# output will only contain one-to-one orthologs
		ecv_out = ecv_o2o
		ecv_dup = NULL
	
	}
	
	
	#### Output ####

	# create cross-species object
	csps = list(
		merged = cbind(mat_sp1[ ecv_out[,1], ], mat_sp2[ ecv_out[,2], ]),
		og_pairs = ecv_out[,c(1,2)],
		og_pairs_is_o2o   = ecv_out$is_o2o,
		og_pairs_ec_value = ecv_out$ec_value,
		sp1 = mat_sp1[ ecv_out[,1], ], 
		sp2 = mat_sp2[ ecv_out[,2], ],
		top_cross_sp1 = NULL, 
		top_cross_sp2 = NULL
	)
	
	# do quantile normalisation on merged matrix?
	if (do_quantile_normalisation) {
		csps$merged = quantile_normalisation(csps$merged)
	}
	
	# find covariable genes
	if (!is.null(cross_fc_thrs) & !is.null(cross_n)) {
		
		# get covariable genes
		message(sprintf("ICC markers | select genes that are variable in both species, fc %.2f in %i samples", cross_fc_thrs, cross_n))
		cross_variable_genes = csps_select_covariable_genes(
			sp1_fp = csps$sp1,
			sp2_fp = csps$sp2, 
			merged = csps$merged,
			cross_fc_thrs = cross_fc_thrs,
			cross_n = cross_n)
		
		top_cross_sp1 = cross_variable_genes$sp1
		top_cross_sp2 = cross_variable_genes$sp2
			
	} else {

		# omit getting covariable genes, get them all instead
		top_cross_sp1 = rownames(ecv_out[,1])
		top_cross_sp2 = rownames(ecv_out[,2])
		
	}
	csps$top_cross_sp1 = top_cross_sp1
	csps$top_cross_sp2 = top_cross_sp2
	
	# return data.frame with EC scores for pairs of genes
	message(sprintf("ICC markers | total num of markers = %i", nrow(ecv_out)))
	return(list(csps = csps, ec_markers = ecv_out, ec_duplicates = ecv_dup))
	
	# TODO: regress-out effect of total expression from final ec values?

}


#' Calculate expression conservation (EC) score for pairs of orthologous genes in a cross-species iterative comparison of coexpression (ICC)
#'
#' @param mat_sp1 expression matrix of species 1 (rows are genes, columns are any conditions)
#' @param mat_sp2 expression matrix of species 2 (rows are genes, columns are any conditions)
#' @param og_pairs either a path to a ortholog pairs table, or one such table. Requires two columns: first are genes from species 1, 
#'    second are genes from species 2
#' @param niter maximum number of ICC iterations (default = 100, min is 2)
#' @param icc_thr similarity threshold to stop ICC iteration (ICC stops when difference between iteration is below this value; default = 0.05)
#' @param method ICC correlation method (default is `pearson`)
#' @param do_duplicates bool (default `FALSE`): if `TRUE`, calculate EC values for best reciprocal pairs of duplicate genes in addition 
#'    to the one-to-one ortholog pairs (which is the default)
#' @param num_o2o_ref_genes if `do_duplicates` is TRUE, select up to this number of genes for the duplicate EC scoring step. Less than 1000 is not
#'    recommended (very unstable EC values).
#' @param use_variable_o2o_genes bool; if `TRUE` (default), use only one-to-one orthologs with variable expression in both species as a reference 
#'    matrix for the duplicates' EC calculations. Variable genes will be computed from footprint matrices (default: top quartile of max footprint
#'    across cells). If input matrices are not footprint matrices, you can provide pre-computed variable genes with `vargenes_sp1` and `vargenes_sp1`
#'    instead (see below).
#' @param variable_o2o_thr if use_variable_o2o_genes is `TRUE`, use this quantile threshold to select variable genes based on the distribution of 
#'    max observed footprint across all cells (default is 0.75, i.e. the top quartile of variable genes). Increasing the size of the reference o2o
#'    matrix can result in extremely long runtimes. Decreasing it can result in poor EC score estimation for paralogs. Rough benchmarking: thr=0.75 
#'    is fairly accurate and fast (~1-2h, cor ~0.80 with whole-matrix result. A thr=0.5 is more accurate but slower (~12h, cor ~0.95).
#' @param vargenes_sp1, vargenes_sp2 (default is `NULL`): vectors of variable one-to-one orthologous genes to use as reference for the 
#'    `do_duplicates` step. Only applicable if `do_duplicates` is `TRUE`. 
#' @param nthreads_icc number of threads to use for ICC calculations (i.e. parallelisation of the correlation matrix).
#' @param nthreads_dup number of parallel processes of duplicate scoring (each of them will run with n=`nthreads_icc` threads).
#' @param cross_fc_thrs numeric, fold change threshold (default: 2, set to NULL to skip filtering by fc)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` we require in each species (default: 1, set to NULL to skip filtering by fc)
#'
#' @return a list with the following: a new csps object using ICC pairs as markers, a dataframe of unique gene pairs from species 1 and 2 with their cross-species expression conservation score, and a dataframe with paralog genes (non-unique pairs).
#'
csps_markers_icc_noobj = function(
	mat_sp1,
	mat_sp2, 
	og_pairs, 
	niter = 100, 
	icc_thr = 0.05, 
	method = "pearson", 
	do_duplicates = FALSE, 
	num_o2o_ref_genes = NULL,
	use_variable_o2o_genes = TRUE, 
	variable_o2o_thr = 0.75, 
	vargenes_sp1 = NULL, 
	vargenes_sp2 = NULL,
	cross_fc_thrs = 2, 
	do_quantile_normalisation = TRUE,
	cross_n = 1,
	nthreads_icc = 2,
	nthreads_dup = 1) {
	
	require("igraph")
	require("data.table")
	
	# register cores
	registerDoParallel(cores = nthreads_dup)

	# ensure data are matrices
	mat_sp1 = as.matrix(mat_sp1)
	mat_sp2 = as.matrix(mat_sp2)
	
	# load orthologous pairs
	if ("character" %in% class(og_pairs)) {
		og_pairs = data.table::fread(og_pairs, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	}	
	colnames(og_pairs) = c("sp1", "sp2")

	# restrict pairs to expressed genes
	og_pairs = og_pairs [ og_pairs[,1] %in% rownames(mat_sp1) & og_pairs[,2] %in% rownames(mat_sp2) , ]
	
	### Find ec values of one-to-one orthologs
	# get one to one pairs (reference)
	list_o2o_sp1 = names(which(table(og_pairs[,1]) == 1))
	list_o2o_sp2 = names(which(table(og_pairs[,2]) == 1))
	bool_o2o = og_pairs[,1] %in% list_o2o_sp1 & og_pairs[,2] %in% list_o2o_sp2
	og_pairs_o2o = og_pairs [ bool_o2o, ]

	# reference o2o matrix
	mar_sp1 = mat_sp1 [ og_pairs_o2o[,1] , ]
	mar_sp2 = mat_sp2 [ og_pairs_o2o[,2] , ]
	
	# ec values for the one to one orthologs
	message(sprintf("ICC markers | pairs of one-to-one orthologs = %i", nrow(og_pairs_o2o)))
	ecv_o2o = csps_calc_icc(mat_sp1 = mar_sp1, mat_sp2 = mar_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = TRUE, num_cores = nthreads_icc)
	# ecv_o2o = ecv_o2o_list[[1]]
	ecv_o2o [ ecv_o2o$ec_value < 0 | is.na(ecv_o2o$ec_value) , "ec_value" ] = 0
	ecv_o2o$is_o2o = 1
	
	### Find best reciprocal duplicate pairs (highest ec value)
	# loop through pairs of non-o2o homologs, add them to the original o2o matrix, and calculate ec value
	if (do_duplicates) { 
		
		# select genes for the reference o2o matrix
		if (use_variable_o2o_genes & is.null(vargenes_sp1) & is.null(vargenes_sp2)) {

			# o2o orthologs with variable expression
			max_fcs_sp1 = apply(mat_sp1, 1, max, na.rm = TRUE)
			max_fcs_sp2 = apply(mat_sp2, 1, max, na.rm = TRUE)
			vargenes_sp1 = rownames(mat_sp1) [ max_fcs_sp1 > quantile(max_fcs_sp1, variable_o2o_thr) ]
			vargenes_sp2 = rownames(mat_sp2) [ max_fcs_sp2 > quantile(max_fcs_sp2, variable_o2o_thr) ]
			# variable AND o2o orthologs
			ecv_o2o_ref = ecv_o2o [ ecv_o2o[,1] %in% vargenes_sp1 & ecv_o2o[,2] %in% vargenes_sp2 , ]
			refgenes_sp1 = ecv_o2o_ref[,1]
			refgenes_sp2 = ecv_o2o_ref[,2]
		
		} else if (!is.null(vargenes_sp1) & !is.null(vargenes_sp2)) { 
			
			ecv_o2o_ref = ecv_o2o [ ecv_o2o[,1] %in% vargenes_sp1 & ecv_o2o[,2] %in% vargenes_sp2 , ]
			refgenes_sp1 = ecv_o2o_ref[,1]
			refgenes_sp2 = ecv_o2o_ref[,2]
			
		} else if (is.null(num_o2o_ref_genes)) { 
			
			refgenes_sp1 = ecv_o2o[,1]
			refgenes_sp2 = ecv_o2o[,2]
			
		} else {
			
			ixs_ref = sample(1:nrow(ecv_o2o), min(nrow(ecv_o2o), num_o2o_ref_genes), replace = FALSE)
			refgenes_sp1 = ecv_o2o[ixs_ref,1]
			refgenes_sp2 = ecv_o2o[ixs_ref,2]
			
		}
		
		# sanity check: do we have enough genes in the reference matrix?
		message(sprintf("ICC markers | reference one-to-one orthologs for duplicate EC scoring = %i", length(refgenes_sp1)))
		if (length(refgenes_sp1) == 0) { 
			stop(sprintf("Only %i pairs of one-to-one orthologs in the reference matrix for duplicate EC calculations. I can't do this.", length(refgenes_sp1)))
		} else if (length(refgenes_sp1) < 1000 & length(refgenes_sp1) > 0) {
			warning(sprintf("Only %i pairs of one-to-one orthologs in the reference matrix for duplicate EC calculations. Too few?", length(refgenes_sp1)))
		}
		
		# get duplicates (all types)
		list_dup_sp1 = names(which(table(og_pairs[,1]) > 1))
		list_dup_sp2 = names(which(table(og_pairs[,2]) > 1))
		bool_dup = og_pairs[,1] %in% list_dup_sp1 | og_pairs[,2] %in% list_dup_sp2
		og_pairs_dup = og_pairs [ bool_dup , ]
		og_pairs_dup = og_pairs_dup[ order(og_pairs_dup[,1], og_pairs_dup[,2]) , ]
		
		# graph of duplicates
		og_pairs_dup_sps = as.matrix(og_pairs_dup)
		og_pairs_dup_sps[,1] = paste("sp1", og_pairs_dup_sps[,1])
		og_pairs_dup_sps[,2] = paste("sp2", og_pairs_dup_sps[,2])
		graph_dup = igraph::graph_from_edgelist(og_pairs_dup_sps, directed = TRUE)
		graph_dup_components = igraph::components(graph_dup)
		
		# dict of duplicates
		dict_dups = list()
		list_components = unique(graph_dup_components$membership)
		for (component in list_components) {
			dups = names(graph_dup_components$membership) [graph_dup_components$membership == component]
			dict_dups[[component]] = dups
		}

		if (length(dict_dups) > 0) {

			# log progress		
			message(sprintf("ICC markers | clusters with non-o2o homologs = %i (total = %i genes)", length(dict_dups), length(unique(unlist(dict_dups)))))
			
			# loop through duplicates in dict_dups
			# this loop is a function for easy parallelisation with plyr::adply
			loop_dups = function(i) {
				
				if (i %% (as.integer(length(dict_dups)/100)) == 0 | i == length(dict_dups)) {
					message(sprintf("ICC markers | clusters with non-o2o homologs = %i (%i/%i processed)", length(dict_dups), i, length(dict_dups)))
				}
				# get lists of genes from sp1 and sp2
				d = dict_dups[[i]]
				dg = stringr::word(d, 2)
				genes1 = dg [ stringr::word(d, 1) == "sp1" ]
				genes2 = dg [ stringr::word(d, 1) == "sp2" ]
				
				# prepare vector of ec values for this set of paralogs
				gen_dups = og_pairs [ og_pairs[,1] %in% genes1 & og_pairs[,2] %in% genes2 , ]
					
				# add test genes to reference expression matrix, sp1
				mai_sp1 = rbind( mat_sp1[ gen_dups[,1] ,], mar_sp1[refgenes_sp1,] )
				rownames(mai_sp1)[ 1:nrow(gen_dups) ] = gen_dups[,1]
				# add test genes to reference expression matrix, sp2
				mai_sp2 = rbind( mat_sp2[ gen_dups[,2] ,], mar_sp2[refgenes_sp2,] )
				rownames(mai_sp2)[ 1:nrow(gen_dups) ] = gen_dups[,2]

				# recalculate ec values using reference expression matrix and new set of duplicates
				ecv_dup_i = csps_calc_icc(mat_sp1 = mai_sp1, mat_sp2 = mai_sp2, niter = niter, icc_thr = icc_thr, method = method, verbose = FALSE, num_cores = nthreads_icc)
				ecv_dup_q = ecv_dup_i [ ecv_dup_i$sp1 %in% gen_dups[,1] & ecv_dup_i$sp2 %in% gen_dups[,2] , ]
						
				return(ecv_dup_q)
				
			}
			
			# apply loop (not parallelised internally)
			ecv_dup = plyr::adply(.data = 1:length(dict_dups), .margins = 1,  .fun = loop_dups, .parallel = TRUE, .id = "cluster", .progress = "none")
			ecv_dup = ecv_dup [ , c("sp1", "sp2", "ec_value", "cluster") ]
			
			# keep best reciprocal pairs of orthologous genes
			ecv_dup_t = data.table::as.data.table(ecv_dup)
			ecv_dup_t$ec_value [ is.na(ecv_dup_t$ec_value) ] = -1
			ecv_dup_f1 = ecv_dup_t[ ecv_dup_t[, .I[ ec_value == max(ec_value) ], by=c("sp1","cluster")]$V1 ]
			ecv_dup_f2 = ecv_dup_t[ ecv_dup_t[, .I[ ec_value == max(ec_value) ], by=c("sp2","cluster")]$V1 ]
			ecv_dup_f  = rbind(ecv_dup_f1,ecv_dup_f2)
			ecv_dup_f  = ecv_dup_f [ duplicated(ecv_dup_f) & ecv_dup_f$ec_value >= 0 , ]
			
			# break ties at random (so as to obtain a truly one-to-one table)
			ecv_dup_f = ecv_dup_f [ !duplicated(ecv_dup_f$sp1), ]
			ecv_dup_f = ecv_dup_f [ !duplicated(ecv_dup_f$sp2), ]
			
			# log progress		
			message(sprintf("ICC markers | pairs of non-o2o homologs kept = %i", nrow(ecv_dup_f)))
			
			# output contains both o2o orthologs & best reciprocal pairs of duplicated genes
			ecv_dup_f = ecv_dup_f [ , c("sp1","sp2","ec_value") ]
			ecv_dup_f$is_o2o = 0
			ecv_out = rbind(ecv_o2o, ecv_dup_f)
		
		} else{
		
			# log progress		
			message(sprintf("ICC markers | clusters with non-o2o homologs = %i (total = %i genes) | Can't score duplicates!", length(dict_dups), length(unique(unlist(dict_dups)))))
			# output will only contain one-to-one orthologs
			ecv_out = ecv_o2o
			ecv_dup = NULL
			
		}

		
	} else {
		
		# output will only contain one-to-one orthologs
		ecv_out = ecv_o2o
		ecv_dup = NULL
	
	}
	
	
	#### Output ####

	# matrices to return
	mot_merged = cbind(mat_sp1[ ecv_out[,1], ], mat_sp2[ ecv_out[,2], ])
	mot_sp1 = mat_sp1[ ecv_out[,1], ]
	mot_sp2 = mat_sp2[ ecv_out[,2], ]
	
	# do quantile normalisation on merged matrix?
	if (do_quantile_normalisation) {
		mot_merged = quantile_normalisation(mot_merged)
	}
	
	# find covariable genes
	if (!is.null(cross_fc_thrs) & !is.null(cross_n)) {
		
		# get covariable genes
		message(sprintf("ICC markers | select genes that are variable in both species, fc %.2f in %i samples", cross_fc_thrs, cross_n))
		cross_variable_genes = csps_select_covariable_genes(
			sp1_fp = mot_sp1,
			sp2_fp = mot_sp2, 
			merged = mot_merged,
			cross_fc_thrs = cross_fc_thrs,
			cross_n = cross_n)
		
		top_cross_sp1 = cross_variable_genes$sp1
		top_cross_sp2 = cross_variable_genes$sp2
			
	} else {

		# omit getting covariable genes, get them all instead
		top_cross_sp1 = rownames(ecv_out[,1])
		top_cross_sp2 = rownames(ecv_out[,2])
		
	}
	
	# return data.frame with EC scores for pairs of genes
	message(sprintf("ICC markers | total num of markers = %i", nrow(ecv_out)))
	return(list(
		ec_markers = ecv_out,
		ec_duplicates = ecv_dup,
		mat_merged = mot_merged,
		mat_sp1 = mot_sp1,
		mat_sp2 = mot_sp2,
		top_cross_sp1 = top_cross_sp1,
		top_cross_sp2 = top_cross_sp2
	))
	

}


#' Calculate expression conservation scores with ICC, from symmetrical tables (same genes & order in both; rownames can be different). This is in an internal function called by `csps_markers_icc`
#'
#' @param mat_sp1 expression matrix of species 1 (rows are genes, columns are any conditions). Row order in matrix 1 and 2 must be matched (pairs of orthologs), but row names need not be.
#' @param mat_sp2 expression matrix of species 2 (rows are genes, columns are any conditions). Row order in matrix 1 and 2 must be matched (pairs of orthologs), but row names need not be.
#' @param niter maximum number of ICC iterations (default = 100, min is 2)
#' @param icc_thr similarity threshold to stop ICC iteration (ICC stops when difference between iteration is below this value; default = 0.05)
#' @param method ICC correlation method (default is `pearson`)
#' @param verbose is an adjective
#'
#' @return dataframe of unique gene pairs from species 1 and 2, and their cross-species expression conservation score
#'
csps_calc_icc = function(mat_sp1, mat_sp2, niter = 100, icc_thr = 0.05, method = "pearson", verbose = TRUE, num_cores = 2) {
	
	require("WGCNA")
	require("tgstat")
	
	# coexpression correlation matrices
	exc_sp1 = WGCNA::cor(t(mat_sp1), method = method, nThreads = num_cores)
	exc_sp2 = WGCNA::cor(t(mat_sp2), method = method, nThreads = num_cores)
	
	# iteration 0: populate with the original cross-species matrix
	# cross-species matrix (single-species matrices must be symmetrical)
	i = 1	
	if (verbose) {
		message(sprintf("ICC iteration %i" , i - 1))
	}
	exc_csp = WGCNA::cor(exc_sp1, exc_sp2, method = method, nThreads = num_cores, use = "pairwise.complete.obs")
	
	# store expression correlation matrices
	# ecv object
	icc_ecv = list()
	icc_ecv[[i]] = diag(exc_csp)
	# per-species correlations
	exp_sp1 = exc_sp1
	exp_sp2 = exc_sp2
	
	for (i in 2:niter) {
		
		# which per-gene ec values are above zero?
		ixs = which(icc_ecv[[i - 1]] > 0)
		
		# for values above zero, keep matrix values from previous iteration
		exi_sp1 = matrix(0, nrow = nrow(exp_sp1), ncol = ncol(exp_sp1))
		exi_sp1[ixs,ixs] = exp_sp1[ixs,ixs]
		exi_sp2 = matrix(0, nrow = nrow(exp_sp2), ncol = ncol(exp_sp2))
		exi_sp2[ixs,ixs] = exp_sp2[ixs,ixs]
		
		# create vector of weights: previous if positive, zero otherwise
		exi_weights = rep(0, length.out = length(icc_ecv[[i - 1]]))
		exi_weights[ixs] = icc_ecv[[i - 1]][ixs]
		
		# weighted correlation
		icc_eci = WGCNA::cor(x = exi_sp1, y = exi_sp2, weights.x = exi_weights, weights.y = exi_weights, method = method, nThreads = num_cores)
		icc_ecv[[i]] = diag(icc_eci)
		
		# calculate delta icc value between iterations
		icc_delta = sum( (icc_ecv[[i]][ixs] - icc_ecv[[i - 1]][ixs] ) ^ 2 )
		if (verbose) {
			message(sprintf("ICC iteration %i | delta = %.1e" , i - 1, icc_delta))
		}
		
		# store expression correlation matrices for next iteration
		exp_sp1 = exi_sp1
		exp_sp2 = exi_sp2

		# quit loop as soon as delta icc value drops below threshold
		if (icc_delta < icc_thr) {
			break
		}
		
	}
	
	# return dataframe with gene pairs and their ec score
	ec_values = data.frame( sp1 = rownames(exc_sp1), sp2 = rownames(exc_sp2), ec_value = icc_ecv[[i]] )
	return(ec_values)
	# return(list(ec_values, exc_csp))
	
}




#' Create cross-species comparison object. 
#' THIS FUNCTION IS DEPRECATED, DO NOT USE
#' 
#' @param sp1_fp_fn,sp2_fp_fn path to expression matrices for the analyzed species, 
#'   rows are genes and columns are cells, metacells or cell types (usually these 
#'   are metacell footprint matrices, `mc@mc_fp`); this can be either raw fold changes 
#'   UMI counts or UMI fraction matrices in either RDS or tab-separated text format.
#'   Alternatively, you can provide a matrix object (like `mc@mc_fp`) as input too.
#' @param OG_pairs_fn either a data.frame with pairs of orthologs, or a path to file to such a 
#'   file (in broccoli-style, tab-separated text format; make sure they are in the right order.
#' @param sp_names character of length 2, optional species names/abbreviations;
#'   if not specified, will use `c("sp1","sp2")`
#' @param make_sp_colnames logical, append the species abbreviations at the 
#'   begining of column names (default: TRUE)
#' @param quant_norm logical (default: TRUE)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` 
#'   we require in each species (default: 1, set to NULL to skip filtering by fc)
#' @param cross_fc_thrs numeric, fold change threshold (default: 2, set to NULL 
#'   to skip filtering by fc)
#' @param preselected_pairs either a dataframe with pre-computed gene pairs to include in the
#'   cross-species comparison matrix, or a path to one such dataframe. If NULL (default), one-to-one
#'   orthologs will be used instead (based on the `OG_pairs_fn` table).
#'
#' @return a list with the following elements:
#'   * merged, a matrix with selected orthologous genes in rows and metacells 
#'     (or cell types) from both species in columns
#'   * og_pairs, data.frame with orthologous pairs in columns
#'   * sp1, a matrix with selected orthologous genes in rows and metacells 
#'     of first species in columns
#'   * sp2, a matrix with selected orthologous genes in rows and metacells 
#'     of the second species in columns
#'   * top_cross_sp1, top orthologous genes in the first species
#'   * top_cross_sp2, top orthologous genes in the second species
#' This will be a standard object for cross-species analyses
#' 
# csps_create_crossspecies_object2 <- function(
# 	sp1_fp_fn, sp2_fp_fn, OG_pairs_fn, sp_names = c("sp1","sp2"), 
# 	make_sp_colnames = TRUE, quant_norm = TRUE, cross_n = 1, cross_fc_thrs = 2, 
# 	preselected_pairs = NULL
# ){ 
	
# 	require("data.table")
# 	require("stringr")
	
# 	### Read expression data ###
	
# 	# species 1
# 	extension = rev(stringr::str_split(sp1_fp_fn, "\\.")[[1]])[1]
# 	if ("matrix" %in% class(sp1_fp_fn)) {
# 		sp1_fp = sp1_fp_fn
# 	} else if ("character" %in% class(sp1_fp_fn) & extension == ".RDS") {
# 		sp1_fp=readRDS(sp1_fp_fn)
# 	} else if ("character" %in% class(sp1_fp_fn)) {
# 		sp1_fp=read.table(sp1_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
# 	} else {
# 		stop("sp1_fp_fn is neither a matrix object nor a path to such an object, I quit.")    
# 	}
# 	# species 2
# 	extension = rev(stringr::str_split(sp2_fp_fn, "\\.")[[1]])[1]
# 	if ("matrix" %in% class(sp2_fp_fn)) {
# 		sp2_fp = sp2_fp_fn
# 	} else if ("character" %in% class(sp2_fp_fn) & extension == ".RDS") {
# 		sp2_fp=readRDS(sp2_fp_fn)
# 	} else if ("character" %in% class(sp2_fp_fn)) {
# 		sp2_fp=read.table(sp2_fp_fn,h=TRUE,row.names=1,sep="\t",quote="",check.names=FALSE,stringsAsFactors=FALSE)
# 	} else {
# 		stop("sp2_fp_fn is neither a matrix object nor a path to such an object, I quit.")    
# 	}
	
	
# 	### Read orthology data ###
	
# 	# read orthology pairs
# 	if ("character" %in% class(OG_pairs_fn)) {
# 		og_pairs = data.table::fread(OG_pairs_fn, header=FALSE, sep="\t", stringsAsFactors=FALSE)
# 	} else {
# 		colnames(og_pairs) = c("sp1","sp2")
# 	}
	
# 	# reduce expression matrices to genes present in the orthology table	
# 	sp1 = sp_names[1]
# 	sp2 = sp_names[2]
# 	sp1_fp = sp1_fp[ intersect(rownames(sp1_fp), og_pairs[,1]), ]
# 	sp2_fp = sp2_fp[ intersect(rownames(sp2_fp), og_pairs[,2]), ]
# 	if (make_sp_colnames) {
# 		colnames(sp1_fp) = paste(sp1, colnames(sp1_fp), sep = "_")
# 		colnames(sp2_fp) = paste(sp2, colnames(sp2_fp), sep = "_")
# 	}
	
# 	# reduce orthology table to only expressed genes
# 	og_pairs = og_pairs [ og_pairs[,1] %in% rownames(sp1_fp) & og_pairs[,2] %in% rownames(sp2_fp) , ]
	
	
# 	### Select cross-species markers ###
	
# 	# if available, reduce og_pairs to pairs present in a second list	
# 	if (!is.null(preselected_pairs)) { 
		
# 		# get table with preselected pairs
# 		if ("character" %in% class(preselected_pairs)) {
# 			list_pairs_csp_t = data.table::fread(preselected_pairs, header=FALSE, sep="\t", stringsAsFactors=FALSE)
# 		} else {
# 			list_pairs_csp_t = preselected_pairs[,1:2]
# 			colnames(list_pairs_csp) = c("sp1","sp2")
# 		}
# 		# list of preselected pairs
# 		list_pairs_csp = paste(list_pairs_csp_t[,1], list_pairs_csp_t[,2], sep = " ")
# 		list_pairs_csp = list_pairs_csp [ list_pairs_csp %in% paste(og_pairs[,1], og_pairs[,2], sep = " ") ]
# 		# reduce og_pairs to preselected pairs
# 		og_pairs = og_pairs [ paste(og_pairs[,1], og_pairs[,2], sep = " ") %in% list_pairs_csp, ]
# 		message(sprintf("csps object | %s-%s | keep gene pairs among preselected pairs, n = %i", sp1, sp2, nrow(og_pairs)))
		
# 	}
	
# 	# if one2one, keep only one-to-one orthologous pairs (defined after discarding unexpressed genes)
# 	if (one2one) {
		
# 		# get one to one pairs
# 		list_o2o_sp1  = names(which(table(og_pairs[,1]) == 1))
# 		list_o2o_sp2  = names(which(table(og_pairs[,2]) == 1))
# 		bool_o2o      = og_pairs[,1] %in% list_o2o_sp1 & og_pairs[,2] %in% list_o2o_sp2
# 		og_pairs      = og_pairs [ bool_o2o, ]
# 		cross_markers = og_pairs[,1]
# 		message(sprintf("csps object | %s-%s | keep one-to-one orthologs, n = %i", sp1, sp2, nrow(og_pairs)))
		
# 	}


# 	### Merged cross-species expression matrix ###
	
# 	# reorder species-specific matrices, and merge
# 	if (class(sp1_fp) == "dgCMatrix" | class(sp2_fp) == "dgCMatrix") {
# 		merged = Matrix::Matrix( cbind( sp1_fp[og_pairs[,1],], sp2_fp[og_pairs[,2],] ), sparse=TRUE)
# 	} else {
# 		merged = data.matrix( cbind( sp1_fp[og_pairs[,1],], sp2_fp[og_pairs[,2],] ) )
# 	}
	
# 	# perform quantile normalisation
# 	if (quant_norm) {
		
# 		# adapt quant_norm parameters according to data matrix class
# 		if (any(class(merged) == "dgCMatrix")) {
# 			sparse = TRUE
# 			merged = as.matrix(merged)
# 			copymat = FALSE
# 		} else {
# 			sparse = FALSE
# 			copymat = TRUE
# 		}
		
# 		# quantile normalisation
# 		merged = tryCatch({
# 			preprocessCore::normalize.quantiles(merged, copy=copymat)
# 		}, error = function(e){
# 			warning(e)
# 			quantile_normalisation(merged)
# 		})
		
# 		# convert to sparse matrix
# 		if (sparse) {
# 			merged = Matrix::Matrix(merged, sparse=TRUE)
# 		}
		
# 	}
	
# 	# select variable genes in BOTH species
# 	if (!is.null(cross_fc_thrs) & !is.null(cross_n)) {
		
# 		# select genes that are variable in both species, defined as having fc 
# 		# above a certain threshold in at least n samples/mcs/cell types/etc.
# 		top_cross_ix = which(
# 				apply(merged[,1:ncol(sp1_fp)                 ], 1, function(x) sort(x, decreasing=TRUE)[cross_n]) > cross_fc_thrs & 
# 				apply(merged[,(ncol(sp1_fp) + 1):ncol(merged)], 1, function(x) sort(x, decreasing=TRUE)[cross_n]) > cross_fc_thrs
# 		)
# 		top_cross_sp1 = unique(og_pairs[,1][top_cross_ix])
# 		top_cross_sp2 = unique(og_pairs[,2][top_cross_ix])
# 		message(sprintf("csps object | %s-%s | retrieve variable markers, n sp1 = %i | sp2 = %i", sp1, sp2, length(top_cross_sp1),length(top_cross_sp2) ) )
		
# 	} else {
		
# 		top_cross_sp1 = unique(og_pairs[,1])
# 		top_cross_sp2 = unique(og_pairs[,2])
# 		message(sprintf("csps object | %s-%s | skip variable markers filter, n sp1 = %i | sp2 = %i", sp1, sp2, length(top_cross_sp1),length(top_cross_sp2) ) )
		
# 	}
	
# 	### Output ###
	
# 	# Return a list with both matrices
# 	out_m1 = merged[ ,1:ncol(sp1_fp) ]
# 	rownames(out_m1) = OG_pairs[[sp1]]
# 	out_m2 = merged[ ,(ncol(sp1_fp) + 1):ncol(merged) ]
# 	rownames(out_m2)=make.names(OG_pairs[[sp2]],unique=TRUE)
	
# 	ids <- unlist(lapply(top_cross,function(x) which(OG_pairs[[sp1]] == x)))
# 	top_cross_sp2=OG_pairs[ids,][[sp2]]
	
# 	# return cross-species object
# 	csps = list(
# 		merged = merged,
# 		og_pairs = OG_pairs,
# 		sp1 = out_m1, sp2 = out_m2,
# 		top_cross_sp1 = top_cross_sp1, top_cross_sp2 = top_cross_sp2
# 	)
# 	return(csps)
# }


# Select genes that are variable in two species, defined as having fc above a certain threshold in at least n samples/mcs/cell types/etc.
#' 
#' @param sp1_fp,sp2_fp,merged matrices of gene (rownames) expression over mcs/cell types/etc (colnames). Rows must be paired, i.e. correspond to comparable genes in the same order
#' @param method character, one of "min_fc" (default: keep genes with fp>`cross_fc_thrs` in at least `cross_n` cell types/metacells), "top_fc_each" (keep up to n=`genes_n` markers with fp>`cross_fc_thrs` in *each* cell type), or "top_fc_pair" (keep up to n=`genes_n` markers with fp>`cross_fc_thrs` in each pair of cross-sps cell types).
#' @param cross_fc_thrs numeric, fold change threshold (default: 2)
#' @param cross_n integer, how many metacells with a fold change of `cross_fc_thrs` we require in each species (default: 1); only relevant if `method="min_fc"`
#' @param genes_n integer, keep up to this many genes for each cell type (default: 100); only relevant if `method="top_fc_each"`
#' 
#' @return a list with two vectors of variable genes from species 1 and 2, respectively
#' 
csps_select_covariable_genes = function(
	sp1_fp, 
	sp2_fp,
	merged,
	method = "min_fc",
	cross_fc_thrs = 2,
	cross_n = 1,
	genes_n = 200) {
	
	# get submatrices from merged matrix (it's not the same as the unmerged matrices because this has been quantile normalised)
	merged_sp1 = merged[,1:ncol(sp1_fp)]
	merged_sp2 = merged[,(ncol(sp1_fp) + 1):ncol(merged) ]
	
	if (method == "min_fc") {
		
		top_cross_ix = which(
			apply(merged_sp1, 1, function(x) sort(x, decreasing=TRUE)[cross_n]) > cross_fc_thrs & 
			apply(merged_sp2, 1, function(x) sort(x, decreasing=TRUE)[cross_n]) > cross_fc_thrs
		)
		top_cross_sp1 = rownames(sp1_fp)[top_cross_ix]
		top_cross_sp2 = rownames(sp2_fp)[top_cross_ix]
	
	} else if (method == "top_fc_each") {
		
		# for each column in the merged matrix, get up to `genes_n` marker genes with a fc>thrs
		top_cross_sp1 = unique(unlist(lapply(
			1:ncol(merged), 
			function(i) {
				fpso = sort(merged[,i], decreasing = TRUE)
				gnso = names(fpso) [ fpso > cross_fc_thrs ] [ 1:genes_n ]
				return(gnso)
			}
		)))
		# reorder sp1
		top_cross_sp1 = top_cross_sp1 [ !is.na(top_cross_sp1) ]
		top_cross_sp1 = top_cross_sp1 [ order( match(top_cross_sp1, rownames(sp1_fp)) ) ]
		# get sp2
		top_cross_sp2 = rownames(sp2_fp) [ rownames(sp1_fp) %in% top_cross_sp1 ]
		
	} else if (method == "top_fc_pair") {
		
		# for each pair of columns in the merged matrix, get up to `genes_n` marker genes with a fc>thrs
		top_cross_sp1 = c()
		for (ni in 1:ncol(merged_sp1)) {
			fpso_i = sort(merged_sp1[,ni], decreasing = TRUE)
			gnso_i = names(fpso_i) [ fpso_i > cross_fc_thrs ] [ 1:genes_n ]
			for (nj in 1:ncol(merged_sp2)) {
				fpso_j = sort(merged_sp2[,nj], decreasing = TRUE)
				gnso_j = names(fpso_j) [ fpso_i > cross_fc_thrs ] [ 1:genes_n ]
				gnso_c = intersect(gnso_i, gnso_j)
				top_cross_sp1 = c(top_cross_sp1, gnso_c)
			}
		}
		# reorder sp1
		top_cross_sp1 = top_cross_sp1 [ !is.na(top_cross_sp1) ]
		top_cross_sp1 = unique(top_cross_sp1)
		top_cross_sp1 = top_cross_sp1 [ order( match(top_cross_sp1, rownames(sp1_fp)) ) ]
		# get sp2
		top_cross_sp2 = rownames(sp2_fp) [ rownames(sp1_fp) %in% top_cross_sp1 ]
		
	} else {
		
		stop("`method` has to be either `min_fc` (default) or `top_fc_each`")
		
	}
	
	# return
	return(list(sp1 = top_cross_sp1, sp2 = top_cross_sp2))

}
