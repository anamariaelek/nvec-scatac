# libraries
require("universalmotif")
require("PWMEnrich")
require("GenomicRanges")
require("IRanges")
require("rtracklayer")
require("plyr")
require('tidyr')
require("doParallel")
require("WGCNA")
require("ggraph")
require("igraph")
require("circlize")
require("ComplexHeatmap")
require("data.table")
require("stringr")
require("xml2")

# deactivate complexheatmap warnings
ht_opt$message = FALSE

##################################
## Helper and utility functions ##
##################################

# Helper function to save or return plots
#' 
#' @param output_file path to file to which the plot will be saved; if `NULL` (default), the plot is returned to stdout
#' @param height,width,res numeric, the width,  height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param EXPR expression that produces the plot (if multiple lines, enclose it in `{ ... }`)
#' 
#' @return plot (if `output_file` is `NULL`), otherwise `NULL`
#' 
plotting_function <- function(
    output_file = NULL, 
    width, height, res = NA, 
    EXP
) {
	if (!is.null(output_file)) {
		extension <- stringr::str_extract(output_file,"(png|pdf)$")
		if ( is.na(extension) ) {
			extension = "pdf"
		}
		# open graphics device
		if (extension == "png") {
			png(output_file, height = height, width = width, res=res)
		} else if (extension == "pdf" ) {
			pdf(output_file, height = height, width = width, useDingbats=TRUE)
		}
	}
	EXP
	if (!is.null(output_file)) dev.off()
}


#' Helper function to convert bytes to human-readable format
#' 
format_memory <- function(size_in_bytes) {
  units <- c("B", "KB", "MB", "GB", "TB")
  power <- min(floor(log(size_in_bytes, 1024)), length(units) - 1)
  size <- size_in_bytes / (1024^power)
  paste0(round(size, 2), " ", units[power + 1])
}

#' Helper function to load bed files as GRanges, with metadata
#' 
#' @param bed_fn path to BED file. Two possible formats: 3-col ("chr", "start", "end"), or more columns (assumes 4th-to-6th cols are "name", "score", "strand"; all other ones are metadata)
#' @param metadata_col_ix vector of indexes of metadata columns to load (5 by default, i.e. score metadata is loaded into the the GRanges object)
#' @param metadata_col_names vector of names for the metadata columns ("score" by default).
#' 
#' @return GRanges object with metadata if appropriate
#' 
mta_bed_to_granges = function(
    bed_fn,
    minimal_columns = 1:3,
    metadata_col_ix = 5,
    metadata_col_name = "score"
) {
  
  require("GenomicRanges")
  require("data.table")
  
  # load
  bed_d = data.table::fread(bed_fn, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
  
  # if file is a short bed file (only three columns), add missing standard columns (name, score, strand)
  if ( is.null(bed_d$V4) ) { 
    bed_d$V4 = NA 
  }
  if ( is.null(bed_d$V5) ) { 
    bed_d$V5 = 0 
  }
  if ( is.null(bed_d$V6) ) { 
    bed_d$V6 = "*"
  }
  
  # build complete dataframe
  bed_d = bed_d [ , c(1:6, metadata_col_ix) ]
  colnames(bed_d) = c("chr","start","end","name","score","strand", metadata_col_name)
  # add 1 to start position (bed start is 0-based)
  bed_d$start = bed_d$start + 1
  # ensure all strans are '+' '-' '*'
  bed_d$strand [ ! (bed_d$strand %in% c("+","-","*")) | is.na(bed_d$strand) ] = "*"
  
  # dataframe to GRanges
  bed_r = GenomicRanges::GRanges(
    seqnames=bed_d$chr, 
    ranges = IRanges(start = bed_d$start, end = bed_d$end), 
    strand = bed_d$strand, 
    name = bed_d$name)
  
  # add metadata columns
  mcols(bed_r)[metadata_col_name] = bed_d[,metadata_col_name]
  
  # return
  return(bed_r)
  
}


#' Helper function to output BED file from GRanges object.
#' 
#' @param gr GRanges object to export.
#' @param bed_fn path to BED file for output. If NULL, returns a `data.frame`
#' @param mcols vector of metadata columns from the GRanges object, to add to BED (in addition to default: name, score, strand).
#' 
#' @return if bed_fn is NULL, a data.frame. If bed_fn is a path to a file, BED is written there and nothing is returned.
#' 
mta_granges_to_bed = function(
    gr,
    bed_fn = NULL,
    mcols = NULL
) {
  
  # what to do if there are no names, scores or strands? (mandatory columns for a complete, well-formatted bed file)
  # name: empty
  if (is.null(mcols(gr)$name)) {
    gr_names = ""
  } else {
    gr_names = mcols(gr)$name
  }
  # score: 0
  if (is.null(mcols(gr)$score)) {
    gr_score = 0
  } else {
    gr_score = mcols(gr)$score
  }
  # strand: "*"
  if (is.null(strand(gr))) {
    gr_strand = "*"
  } else {
    gr_strand = strand(gr)
  }
  
  # create dataframe
  bed_d = data.frame(
    chr = seqnames(gr),
    start = start(gr) - 1,
    end = end(gr),
    name = gr_names,
    score = gr_score,
    strand = gr_strand
  )
  
  # add metadata if available
  for (mcol in mcols) {
    bed_d[,mcol] = mcols(gr)[[mcol]]
  }
  
  # return
  if (! is.null(bed_fn)) {
    options(scipen=10)
    write.table(bed_d, file = bed_fn, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
    options(scipen=0)
  } else {
    return(bed_d)
  }
  
}


#' Modified version of `universalmotif::read_homer` that does not recompute p-values when importing HOMER motifs
#' 
#' @param file path to HOMER motif file
#' @param skip num of lines to skip at the beginning of the homer motif file
#' @param pval_thr ignore HOMER motifs with pvalue above this threshold (default NULL, i.e. all motifs are included)
#' @param sanitise_names default TRUE; whether to sanitise motif names upon import (remove backslashes, spaces, filler text)
#' @param prefix string to add as a prefix to each imported motif (default NULL, i.e. no prefix is added)
#' 
#' @return motif object
#' 
mta_read_homer_mod = function(
	file,
	skip = 0, 
	pval_thr = NULL, 
	ignore_SeqBias = TRUE, 
	sanitise_names = TRUE, 
	prefix = NULL, 
	make_unique_names = TRUE
) {
	
	# parsing functions for modified universalmotif functions
	universalmotif_cpp <- function(motif, name = "new motif", altname = NA_character_, family = NA_character_, organism = NA_character_, alphabet = "DNA", type = NA_character_, icscore = as.numeric( c()), nsites = as.numeric( c()), pseudocount = 1.0, bkg = as.numeric( c()), bkgsites = as.numeric( c()), consensus = NA_character_, strand = "+-", pval = as.numeric( c()), qval = as.numeric( c()), eval = as.numeric( c()), extrainfo = NA_character_, isgapped = NA_integer_, gaploc = as.numeric( c()), mingap = as.numeric( c()), maxgap = as.numeric( c())) {
		.Call("_universalmotif_universalmotif_cpp", PACKAGE = "universalmotif", motif, name, altname, family, organism, alphabet, type, icscore, nsites, pseudocount, bkg, bkgsites, consensus, strand, pval, qval, eval, extrainfo, isgapped, gaploc, mingap, maxgap)
	}
	
	validObject_universalmotif <- function(motif, throw_error = TRUE) {
		.Call("_universalmotif_validObject_universalmotif", PACKAGE = "universalmotif", motif, throw_error)
	}
	parse_motifs <- function(x, y) {
		motif <- raw_lines[x:y]
		motif <- vapply(motif, function(x) as.numeric(strsplit(x, "\\t+")[[1]]), numeric(4))
		matrix(motif, ncol = 4, byrow = TRUE)
	}
	parse_meta <- function(x) {
		x <- strsplit(x, "\\t+")[[1]]
		if (grepl("/", x[2])) {
			y <- strsplit(x[2], "/")[[1]]
			y <- strsplit(y[1], "\\(")[[1]]
			x[2] <- y[1]
			family <- strsplit(y[2], "\\)")[[1]][1]
		}
		else family <- character(0)
		x2 <- strsplit(x[6], ",")[[1]]
		nsites <- strsplit(strsplit(x2[1], "T:")[[1]][2], "\\(")[[1]][1]
		bkgsites <- strsplit(strsplit(x2[2], "B:")[[1]][2], "\\(")[[1]][1]
		pval <- strsplit(x2[3], "P:")[[1]][2]
		c(name = x[2], nsites = nsites, bkgsites = bkgsites, 
		  pval = pval, threshold = x[3], family = family)
	}
	
	# read motif object
	raw_lines <- readLines(con <- file(file))
	close(con)
	if (skip > 0) 
		raw_lines <- raw_lines[-seq_len(skip)]
	raw_lines <- raw_lines[raw_lines != ""]
	headers <- which(grepl("^>", raw_lines))
	motif_starts <- headers + 1
	if (length(headers) == 1) {
		motif_stops <- length(raw_lines)
	} else {
		motif_stops <- c(headers[-1] - 1, length(raw_lines))
	}
	
	# load motifs & metadata
	motif_list <- mapply(parse_motifs, motif_starts, motif_stops, 
						 SIMPLIFY = FALSE)
	motif_meta <- lapply(raw_lines[headers], parse_meta)
	
	# construct & validate motif object
	homer2umot <- function(x, y) {
		mot <- universalmotif_cpp(
			name = gsub(" ", "_" , x[1]), 
			nsites = ifelse(is.na(as.numeric(x[2])), numeric(0), as.numeric(x[2])), 
			bkgsites = ifelse(is.na(as.numeric(x[3])), numeric(0), as.numeric(x[3])), 
			motif = t(y), 
			alphabet = "DNA", 
			type = "PPM", 
			family = x[6])
		validObject_universalmotif(mot)
		mot
	}
	
	# create motif object
	motifs <- mapply(homer2umot, motif_meta, motif_list, SIMPLIFY = FALSE)
	
	# get pvalue directly from Homer
	# Beware: this pvalues are quite different from those calculated by `universalmotif`
	motifs <- mapply(function(x, y) {
		x@pval = as.numeric(y["pval"])
		x
	}, motifs, motif_meta)
	
	# drop motifs with low pvalue
	if (!is.null(pval_thr)) {
		pvals = lapply(motif_meta, function(i) i["pval"])
		motifs = motifs [ as.numeric(pvals) <= pval_thr ]
	}
	
	# drop motifs from SeqBias
	if (ignore_SeqBias) {
		motifs = motifs [ grep("SeqBias:",plyr::laply(motifs, function(i) i@name), invert = TRUE) ]
	}
	
	# clean motif names upon imprt
	if (sanitise_names) {
		lapply(1:length(motifs), function(i) motifs[[i]]@name <<- gsub("/","", motifs[[i]]@name) )
		lapply(1:length(motifs), function(i) motifs[[i]]@name <<- gsub(" ","", motifs[[i]]@name) )
		lapply(1:length(motifs), function(i) motifs[[i]]@name <<- gsub(",BestGuess:",":", motifs[[i]]@name) )
	}
	
	# add prefix?
	if (make_unique_names & !is.null(prefix)) {
		lapply(1:length(motifs), function(i) motifs[[i]]@name <<- gsub("^", sprintf("%s_m%i_", prefix, i), motifs[[i]]@name) )
	} else if (make_unique_names &is.null(prefix)) {
		lapply(1:length(motifs), function(i) motifs[[i]]@name <<- gsub("^", sprintf("m%i_", i), motifs[[i]]@name) )
	} else if (!make_unique_names & !is.null(prefix)) {
		lapply(1:length(motifs), function(i) motifs[[i]]@name <<- gsub("^", sprintf("%s_", prefix), motifs[[i]]@name) )
	}
	
	# return
	message(sprintf("load motifs | n = %i", length(motifs)))
	if (length(motifs) == 1) 
		motifs <- motifs[[1]]
	return(motifs)
	
}


#' Read CISBP files into universalmotif format, without frills
#' 
#' @param file path to CISBP motif file
#' @param motif_id motif id, inferred from filename if not set
#' @param family family id, empty if not set
#' 
mta_read_cisbp_mod <- function (
	file, 
	motif_id = NULL, 
	family = NULL
) {
	
	# parsing functions for modified universalmotif functions
	universalmotif_cpp <- function(motif, name = "new motif", altname = NA_character_, family = NA_character_, organism = NA_character_, alphabet = "DNA", type = NA_character_, icscore = as.numeric( c()), nsites = as.numeric( c()), pseudocount = 1.0, bkg = as.numeric( c()), bkgsites = as.numeric( c()), consensus = NA_character_, strand = "+-", pval = as.numeric( c()), qval = as.numeric( c()), eval = as.numeric( c()), extrainfo = NA_character_, isgapped = NA_integer_, gaploc = as.numeric( c()), mingap = as.numeric( c()), maxgap = as.numeric( c())) {
		.Call("_universalmotif_universalmotif_cpp", PACKAGE = "universalmotif", motif, name, altname, family, organism, alphabet, type, icscore, nsites, pseudocount, bkg, bkgsites, consensus, strand, pval, qval, eval, extrainfo, isgapped, gaploc, mingap, maxgap)
	}

	# read file
	pwm = read.table(file, header = TRUE, row.names = 1)

	# do we have motif id name?
	if (is.null(motif_id)) {
		motif_id = gsub(".txt$","",basename(file))
	}
	if (is.null(family)) {
		family = ""
	}
	
	# create universalmotif object
	if ( nrow(pwm > 0) ) {
		motifs = universalmotif_cpp(
			name = motif_id, 
			nsites = NA, 
			bkgsites = NA, 
			motif = t(pwm), 
			alphabet = "DNA", 
			type = "PPM", 
			family = family)
	} else {
		warning(sprintf("empty motif %s", file))
		motifs = NULL
	}
	
	return(motifs)
}
#' Read XML files (as produced by Dimont) into universalmotif format
#' 
#' @param file path to XML file
#' @param motif_id motif id, inferred from filename if not set
#' @param base_order order of bases in the PWM
#' 
#' @return motif object
#' 
mta_read_xml <- function(file, motif_id = NULL, base_order = c("A", "C", "G", "T")) {

	# Parsing functions for modified universalmotif functions
	universalmotif_cpp <- function(motif, name = "new motif", altname = NA_character_, family = NA_character_, organism = NA_character_, alphabet = "DNA", type = NA_character_, icscore = as.numeric( c()), nsites = as.numeric( c()), pseudocount = 1.0, bkg = as.numeric( c()), bkgsites = as.numeric( c()), consensus = NA_character_, strand = "+-", pval = as.numeric( c()), qval = as.numeric( c()), eval = as.numeric( c()), extrainfo = NA_character_, isgapped = NA_integer_, gaploc = as.numeric( c()), mingap = as.numeric( c()), maxgap = as.numeric( c())) {
		.Call("_universalmotif_universalmotif_cpp", PACKAGE = "universalmotif", motif, name, altname, family, organism, alphabet, type, icscore, nsites, pseudocount, bkg, bkgsites, consensus, strand, pval, qval, eval, extrainfo, isgapped, gaploc, mingap, maxgap)
	}
	
	validObject_universalmotif <- function(motif, throw_error = TRUE) {
		.Call("_universalmotif_validObject_universalmotif", PACKAGE = "universalmotif", motif, throw_error)
	}

	# Read the XML file
    xml_data <- xml2::read_xml(file)

    # Extract the PWM positions (loop through all <pos> tags)
    positions <- xml2::xml_find_all(xml_data, ".//pos")

    # Initialize an empty matrix to store the PWM values
    num_positions <- as.integer(xml_text(xml_find_first(xml_data, ".//length")))
    num_bases <- length(base_order)
    
    # Create a matrix to store the PWM values
    pwm_matrix <- matrix(NA, nrow = num_bases, ncol = num_positions)

    # Fill the PWM matrix
    for (i in 1:num_positions) {
      # Extract the values for each position (base probabilities)
      base_values <- xml_find_all(xml_data, paste0(".//pos[@val='", i-1, "']/pos"))
      
      for (j in 1:length(base_values)) {
        base_value <- as.numeric(xml_text(base_values[j]))
        pwm_matrix[j, i] <- base_value
      }
    }

    # Set row names for bases
    rownames(pwm_matrix) <- base_order

    # Set column names for positions
    colnames(pwm_matrix) <- seq_len(num_positions)

	# Get the motif name
	if (is.null(motif_id)) {
		motif_id <- gsub(".xml$", "", basename(file))
	}
	
	# construct & validate motif object
	mot <- universalmotif_cpp(
		name = motif_id, 
		nsites =  numeric(0), 
		bkgsites =  numeric(0), 
		motif = pwm_matrix, 
		alphabet = "DNA", 
		type = "PPM"
	)
	validObject_universalmotif(mot)
	
    return(mot)

}

#' Parse Homer results
#'
#' @param fn path to homerResults or knownResults folder
#' @param eta numeric, constant to add for calculating log2 fold change
#' 
#' @return data.table with motifs results
#' 
mta_parse_homer <- function(res_dir, eta = 1) {

  # parse homerResults folder to get names of best hits for denovo motifs
  if (grepl("homerResults$", res_dir)) {
	ms <- list.files(res_dir, pattern = "motif\\d+.motif", full.names = TRUE)
  } else if (grepl("knownResults$", res_dir)) {
	ms <- list.files(res_dir, pattern = "known\\d+.motif", full.names = TRUE)
  }
  dt <- rbindlist(sapply(ms, function(x) {
    m <- universalmotif::read_homer(x)
    data.table(
      name = m@name,
      consensus = m@consensus,
      logodds_threshold = m@extrainfo[["logodds.threshold"]],
      pval = m@pval,
      fg_count = m@nsites,
      bg_count = m@bkgsites
    )[]
  }, simplify = FALSE, USE.NAMES = TRUE), idcol = "filepath")
  # dt[, name := str_remove(name, ".+(?=\\,BestGuess),BestGuess:")]
  dt[, fc := (fg_count + eta) / (bg_count + eta)]
  dt[, .(
    name, consensus, fc, pval,
    fg_count, bg_count, filepath
  )]
}


#' Trim granges object by genome coordinates
#' 
#' @param gr_object either a path to bed file with ranges, or a preloaded 
#' GRanges object to trim
#' @param index_object either a path to a genome index file (`.fai`) with chr 
#' names and lengths, or a preloaded data.frame with these two columns
#' 
mta_trim_granges = function(gr_object, index_object) {
  
  # read genome index
  if ("character" %in% class(index_object)) {
    gix = read.table(index_object, stringsAsFactors = FALSE)
    gix = gix[,1:2]
    colnames(gix) = c("chr","length")
  } else if ("data.frame" %in% class(index_object)) {
    gix = index_object
    gix = gix[,1:2]
    colnames(gix) = c("chr","length")
  } else {
    stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`)", index_object))
  }
  
  # read peaks
  if ("character" %in% class(gr_object)) {
    # load
    pka_r = mta_bed_to_granges(gr_object)
    names(pka_r) = pka_r$name
  } else if ("GRanges" %in% class(gr_object)) {
    pka_r = gr_object
    names(pka_r) = pka_r$name
  } else {
    stop(sprintf("`gr_object` %s has to be either a GRanges object, or a path to a compatible BED file", gr_object))
  }
  
  # add chr lengths to genome annotation and peak annotation
  seqlevels(pka_r)  <- as.character(gix$chr)
  seqlengths(pka_r) <- as.integer(gix$length)
  
  # trim
  pka_r_trimmed = GenomicRanges::trim(pka_r)
  
  # return 
  return(pka_r_trimmed)
  
}


#' Liftover annotation bed file from concatenated to original coordinates
#' 
#' @param faif character, path to chromosome size (fasta.fai) of the oiginal assembly, 
#' which was used to construct concatenated genome, it should have been sorted alphanumerically by scaffold
#' @param bedf character, path to input annotation bed or bed-like file, 
#' it should have at least the first three columns: chr, start, end
#' @param bedfout character, output filename, if NULL will use the original filename with '-liftover' appended to it before extension
#'
mta_liftover_concat_bed_files <- function(faif, bedf, bedfout=NULL) {
  
  # load fasata sizes file
  fai <- fread(faif)
  setnames(fai, colnames(fai)[1:2], c("chr_fai","size"))
  
  # load bed file
  bed <- fread(bedf)
  j <- pmin(7,ncol(bed))
  setnames(bed, colnames(bed)[1:j], c("chr","start","end","name","width","strand","score")[1:j])
  setkey(bed,start,end)
  
  # output bed file
  if (is.null(bedfout))
    bedfout <- str_replace(bedf,".bed","-liftover.bed")
  
  # calculate coordinates of original scaffolds in concatenated assembly
  tot <- cumsum(fai$size)
  fai[,start_fai:=data.table::shift(size,n=1L)]
  fai[is.na(start_fai),start_fai:=1]
  fai[,start_fai:=cumsum(start_fai)]
  fai[,end_fai:=cumsum(size)]
  setkey(fai,start_fai,end_fai)
  
  # overlap = get original scaffold for each bed interval
  fovl <- foverlaps(bed, fai, type="within")
  fovl <- fovl[!is.na(chr_fai)]
  
  # calculate coordinates of original intervals in originla assembly
  fovl[,start_new:=start-start_fai]
  fovl[,end_new:=end-start_fai]
  fovl[,width_new:=end_new-start_new]
  
  # save
  cols <- c("chr_fai","start_new","end_new","name","width","strand","score")[1:j]
  ordcols <- c(cols, setdiff(colnames(bed), c("chr","start","end",cols)))
  bed_new <- fovl[,..ordcols]
  fwrite(bed_new, bedfout, sep="\t", col.names = FALSE)
  
  return(bed_new)
}


#####################
## Peaks functions ##
#####################

#' Match ATAC peaks to genes based on coordinates and evidence-based promoter regions (optionally)
#' 
#' @param gff_object either a path to a GFF file with gene annotations, or a preloaded `GRanges` object
#' @param peak_object either a path to a BED file with ATAC peaks, or a preloaded `GRanges` object
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param feature_to_match type of feature in the GFF that will be matched to the ATAC peaks (default: "gene"; "transcript" can be used for GTFs, etc.)
#' @param exclude_genes vector of genes to exclude (e.g. nonexpressed genes); defaults to NULL
#' @param max_tss_dist max distance away from a gene's TSS to check for ATAC peaks that can be assigned to it (default: 20000 bp)
#' @param min_overlap min number of overlapping bases between a gene/gene extended region and an ATAC peak so as to match them(default: 0 bp, i.e. any overlap is sufficient)
#' @param promoter_upstream,promoter_downstream distance upstream/downstream of a TSS region that's used to define promoters (inferred from `gff_object`)
#' @param promoter_object optionally, an additional set of promoters can be provided (will be merged with GFF-defined promoters). Either a path to a BED file, or a preloaded `GRanges` object
#' 
#' @return data.frame with many-to-many matches between genes and peaks
#' 
mta_match_peaks_to_genes = function(
	gff_object, peak_object, index_object, 
	list_genes = NULL, 
	feature_to_match = "gene", 
	feature_field = "gene_id", 
	exclude_genes = NULL,     # genes to exclude from gff file (can be nonexpressed genes, for example)
	max_tss_dist = 20000,    # peaks beyond this distance should never be assigned to a particular gene
	min_overlap = 0, # min overlap between a peak and a gene/gene extended region to link them
	max_gap = 1,
	promoter_upstream = 200, # these regions are used to define promoters around TSS: peaks within this region get automatically assigned to the relevant gene (and NEVER to any other gene)
	promoter_downstream = 50,
	promoter_object = NULL # a bed file with promoter regions, so that (i) we can directly assign peaks to promoters (and assign them to a gene as if these were TSS, see above), and (ii) we can enforce a rule that prevents peaks to be assigned to a gene if there's a promoter between them
) {
	
	require("rtracklayer")
	require("GenomicRanges")
	require("IRanges")
	
	# temporarily suppress warnings
	options(warn=-1)
	on.exit(options(warn = 0))
	
	### Read data ###
	
	message(sprintf("%s | Loading data...", Sys.time()))
	
	# read gff 
	if ("character" %in% class(gff_object)) {
		gff_r = rtracklayer::readGFFAsGRanges(gff_object)
	} else if ("GRanges" %in% class(gff_object)) {
		gff_r = gff_object
	} else {
		stop(sprintf("`gff_object` %s has to be either a GRanges object from a GFF, or a path to a GFF", gff_object))
	}
	
	# read genome index
	if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`)", index_object))
	}
	
	# read peaks
	if ("character" %in% class(peak_object)) {
		# load
		pka_r = mta_bed_to_granges(peak_object)
		names(pka_r) = pka_r$name
	} else if ("GRanges" %in% class(peak_object)) {
		pka_r = peak_object
		names(pka_r) = pka_r$name
	} else {
		stop(sprintf("`peak_object` %s has to be either a GRanges object, or a path to a compatible BED file", peak_object))
	}
	
	# if available, load promoter regions
	if ("character" %in% class(promoter_object)) {
		# load
		pro_r = mta_bed_to_granges(promoter_object)
		names(pro_r) = pro_r$name
	} else if ("GRanges" %in% class(promoter_object)) {
		pro_r = promoter_object
		names(pro_r) = pro_r$name
	} else if (!is.null(promoter_object)) {
		stop(sprintf("`promoter_object` %s has to be either a GRanges object, or a path to a compatible BED file", promoter_object))
	}    
	
	
	### Create range objects ###
	# get genes (or any other feature to match)
	ann_r = gff_r[gff_r$type == feature_to_match,]
	names(ann_r) = as.data.frame(ann_r)[,feature_field]
	ann_r <- ann_r[!is.na(names(ann_r))]
	ann_r$name = names(ann_r)	

	# add chr lengths to genome annotation and peak annotation
	seqlevels(ann_r)  <- as.character(gix$chr)
	seqlengths(ann_r) <- as.integer(gix$length)
	seqlevels(pka_r)  <- as.character(gix$chr)
	seqlengths(pka_r) <- as.integer(gix$length)
	if (!is.null(promoter_object)) {
		seqlevels(pro_r)  <- as.character(gix$chr)
		seqlengths(pro_r) <- as.integer(gix$length)
	}
	
	# exclude genes if necessary
	if (!is.null(exclude_genes)) {
		ann_r = ann_r [ !names(ann_r) %in% exclude_genes ]
	}
	
	# log
	message(sprintf("%s | Considering %i features...", Sys.time(), length(ann_r)))
	
	# get tss (strand-aware)
	ann_tss_r = GenomicRanges::resize(ann_r, width = 1)
	
	# get promoter regions
	ann_pro_r = GenomicRanges::promoters(ann_tss_r, upstream = promoter_upstream, downstream = promoter_downstream) 
	ann_pro_r = GenomicRanges::trim(ann_pro_r)
	# ann_pro_r = GenomicRanges::GRanges(seqnames = seqnames(ann_pro_r), ranges = ranges(ann_pro_r))
	ann_pro_r$matched_gene = names(ann_pro_r)
	
	# if we have coordinates of promoters:
	# - get promoter regions overlapping with a TSS, and extend gene feature-inferred promoter to the full extend of the peak (better promoter)
	# - also keep promoter regions NOT overlapping with TSS, which will be used to define boundaries around genes. Peaks beyond these boundaries 
	#   won't be assigned to these genes.
	if (!is.null(promoter_object)) {
		
		# which predefined promoters overlap with GFF-inferred promoter regions?
		ovs_pro = IRanges::findOverlapPairs(ann_pro_r, pro_r,  select = "all", type = "any")
		
		# expand GFF-inferred promoters for each gene to include predefined peaks (punion = parallel union)
		ovs_pro_r = GenomicRanges::punion(ovs_pro@first, ovs_pro@second)
		ovs_pro_r = GenomicRanges::trim(ovs_pro_r)
		ovs_pro_r = unique(ovs_pro_r)
		ovs_pro_r$matched_gene = names(ovs_pro_r)
		# add these extended promoters to original ranges
		ann_pro_r [ names(ovs_pro_r) ] = ovs_pro_r
		message(sprintf("%s | %s promoters assigned to genes", Sys.time(), length(ovs_pro_r)))
		
		# get promoters not assigned to any gene
		ovs_nopro = IRanges::findOverlapPairs(ann_pro_r, pro_r,  select = "all", type = "any")
		ovs_nopro_r = pro_r [ ! names(pro_r) %in% names(ovs_nopro@second) ]
		if (length(ovs_nopro_r)>1) {
		  ovs_nopro_r$matched_gene = "unknown"
		}
		# add unassigned promoters to original ranges
		ann_pro_r = c(ann_pro_r, ovs_nopro_r)
		ann_pro_r = GenomicRanges::trim(ann_pro_r)
		message(sprintf("%s | %s promoters not assigned to genes", Sys.time(), length(ovs_nopro_r)))
		
	}
	
	# FIRST: find ATAC peaks that are overlapping promoter regions. 
	# These regions get assigned to these genes, straight away.
	pro_r = ann_pro_r [ mcols(ann_pro_r)$matched_gene %in% names(ann_r)]
	ovs_pka_pro = IRanges::findOverlapPairs(pka_r, pro_r, select="all", minoverlap=min_overlap, maxgap=max_gap)
	ovs_pka_pro_t = data.table(
		gene = names(ovs_pka_pro@second),
		peak = names(ovs_pka_pro@first),
		dist = distance(ovs_pka_pro@first, ann_tss_r[ ovs_pka_pro@second$name ], ignore.strand=FALSE)
	)
	ovs_pka_pro_t = unique(ovs_pka_pro_t)
	message(sprintf("%s | %s peaks overlapping promoters", Sys.time(), length(unique(ovs_pka_pro_t$peak))))
	
	# for each peak assigned to two genes, keep closest TSS unique assignment 
	mult_peaks <- ovs_pka_pro_t[,.N,peak][N>1]$peak
	ovs_pka_pro_c <- ovs_pka_pro_t[ peak %in% mult_peaks ]
	setorder(ovs_pka_pro_c, peak, -dist)
	ovs_pka_pro_c[, min_dist:=min(dist),peak ]
	ovs_pka_pro_c <- ovs_pka_pro_c[ dist==min_dist ]
	
	# combine
	ovs_pka_pro_c[,c("dist","min_dist"):=NULL]
	ovs_pka_pro_t[,dist:=NULL]
	ovs_pka_pro_t <- rbindlist(list(ovs_pka_pro_t[!peak %in% mult_peaks], ovs_pka_pro_c))
	class(ovs_pka_pro_t) <- "data.frame" # because others don't use data.table :(
	ovs_pka_pro_t_list = paste(ovs_pka_pro_t$gene, ovs_pka_pro_t$peak)
	message(sprintf("%s | %i peak-gene assignments in promoters", Sys.time(), nrow(ovs_pka_pro_t)))
	
	# SECOND: find ATAC peaks that are located within each gene's body.
	ovs_pka_gen = IRanges::findOverlapPairs(pka_r, ann_r, select="all", minoverlap=min_overlap, ignore.strand=TRUE)
	ovs_pka_gen_t = data.frame(
	  gene = names(ovs_pka_gen@second),
	  peak = names(ovs_pka_gen@first)
	)
	ovs_pka_gen_t <- unique(ovs_pka_gen_t)
	
	
	# THIRD: find ATAC peaks that are located within each gene's extended region.
	# These extended regions are defined by a wider area around the gene (e.g. +/- 20kb), 
	# but are constrained by the presence of promoters of other genes. These other 
	# promoters can be defined from the GFF, or from a precomputed track of promoter 
	# coordinates (e.g. from CHIP-seq H3K4me3 data).
	
	message(sprintf("%s | Getting extended regions genes, this may take a while", Sys.time()))
	
	# get extended regions of interest (TSS +/- distance)
	ann_ext_r = GenomicRanges::resize(ann_tss_r, width = width(ann_tss_r) + (max_tss_dist * 2), fix = "center")
	ann_ext_r = GenomicRanges::trim(ann_ext_r)

	# ------------------------------------------------------------------------------------------------------
	# This version is faster but follow and preceede don't consider overlapping TSS 
	# for limiting extended regions
	#
	# # limit the regions by presence of nearest non-overlaping TSS
	# fl <- follow(ann_r, ann_tss_r, ignore.strand=TRUE)
	# fl_i <- !is.na(fl)
	# start(ann_ext_r[fl_i]) <- pmax(start(ann_ext_r[fl_i]), start(ann_tss_r[fl[fl_i]]))
	# pr <- precede(ann_r, ann_tss_r, ignore.strand=TRUE)
	# pr_i <- !is.na(pr)
	# end(ann_ext_r[pr_i]) <- pmin(end(ann_ext_r[pr_i]), end(ann_tss_r[pr[pr_i]]))
	# # get peaks in this region
	# ovs_pka_ext = IRanges::findOverlaps(pka_r, ann_ext_r, select="all", minoverlap=min_overlap, maxgap=max_gap)
	# ovs_pka_ext_t = data.table(
	#   gene = ann_ext_r[subjectHits(ovs_pka_ext)]$name,
	#   peak = pka_r[queryHits(ovs_pka_ext)]$name
	# )
	# ovs_pka_ext_t = unique(ovs_pka_ext_t)
	
	# ------------------------------------------------------------------------------------------------------
	# This version takes into account overlapping TSS for limiting extended regions
	# but it is super slow because of looping over genes
	# 
	# for each gene, restrict extended region by presence of nearest upstream and downstream prooters
	ovs_pka_ext_l <- lapply(1:length(ann_ext_r$name), function(i) {
	  if (i%%1000==0) message(sprintf("%s | Matching gene %i out of %i",Sys.time(),i,length(ann_ext_r$name)))
	  gene <- ann_ext_r$name[i]
	  ovl_genes <- ""
	  # exclude genes which overlap TSS of the gene of interest
	  ann_r_g <- ann_r[ann_r$name!=gene]
	  ann_tss_r_g <- ann_tss_r[ann_tss_r$name==gene]
	  ovl <- findOverlaps(ann_tss_r_g, ann_r_g, ignore.strand=TRUE)
	  if (length(ovl)>0) {
	    ovl_genes <- ann_r_g[subjectHits(ovl)]$name
	  }
	  # exclude genes whose TSS overlaps the gene of interest
	  ann_r_g <- ann_r[ann_r$name==gene]
	  ann_tss_r_g <- ann_tss_r[ann_tss_r$name!=gene]
	  ovl <- findOverlaps(ann_tss_r_g, ann_r_g, ignore.strand=TRUE)
	  if (length(ovl)>0) {
	    ovl_genes <- c(ovl_genes, ann_tss_r_g[queryHits(ovl)]$name)
	  }
	  ann_tss_r_o <- ann_tss_r[!ann_tss_r$name %in% c(gene,ovl_genes)]
	  ann_pro_r_g <- narrow(promoters(ann_r_g, upstream = 0, downstream = 50), start = 25, width = 1)
	  # find closest downstream and upstream TSS
	  prec_i <- unique(precede(ann_pro_r_g, ann_tss_r_o, ignore.strand=TRUE))
	  if (!is.na(prec_i)) {
	    prec_g <- ann_tss_r_o[prec_i]
	    x <- ifelse(strand(prec_g)=="+", end(prec_g), start(prec_g))
	  } else {
	    x <- ifelse(strand(ann_r_g)=="+", start(ann_ext_r[i]), end(ann_ext_r[i]))
	  }
	  foll_i <- unique(follow(ann_pro_r_g, ann_tss_r_o, ignore.strand=TRUE))
	  if (!is.na(foll_i)) {
	    foll_g <- ann_tss_r_o[foll_i]
	    y <- ifelse(strand(foll_g)=="+", start(foll_g), end(foll_g))
	  } else {
	    y <- ifelse(strand(ann_r_g)=="+", end(ann_ext_r[i]), start(ann_ext_r[i]))
	  }
	  ext_r <- GRanges(seqnames(ann_pro_r_g), IRanges(start=min(x,y), end = max(x,y)))
	  # get peaks in this region
	  ovs_pka_ext = IRanges::findOverlaps(pka_r, ext_r, select="all", minoverlap=min_overlap, maxgap=max_gap)
	  ovs_pka_ext_t = data.table(
	    gene = gene,
	    peak = pka_r[queryHits(ovs_pka_ext)]$name
	  )
	  ovs_pka_ext_t = unique(ovs_pka_ext_t)
	  ovs_pka_ext_t
	})
	ovs_pka_ext_t <- rbindlist(ovs_pka_ext_l)
	ovs_pka_ext_t <- ovs_pka_ext_t[!is.na(peak)]
	
	# ------------------------------------------------------------------------------------------------------
	
	# explicitly exclude promoters from the extended region data!
	ovs_pka_ext_t <- ovs_pka_ext_t[! ovs_pka_ext_t$peak %in% ovs_pka_pro_t$peak, ]
	message(sprintf("%s | %s peaks overlapping extended regions",Sys.time(),length(unique(ovs_pka_ext_t$peak))))
			
	# filtering peaks further than max distance
	ovs_pka_ext_t$dist_to_tss = GenomicRanges::distance( ann_tss_r [ ovs_pka_ext_t$gene ] , pka_r [ ovs_pka_ext_t$peak ], ignore.strand = TRUE) 
	ovs_pka_ext_t = ovs_pka_ext_t [ ovs_pka_ext_t$dist_to_tss <= max_tss_dist , ]
	ovs_pka_ext_t$dist_to_tss <- NULL

	# merge gene-peak assignments from promoters and extended regions
	ovs_pka_g = rbind(ovs_pka_pro_t, ovs_pka_gen_t)
	ovs_pka_t = rbind(ovs_pka_g, ovs_pka_ext_t)
	ovs_pka_t = unique(ovs_pka_t)
	ovs_pka_t$dist_to_tss = GenomicRanges::distance( ann_tss_r [ ovs_pka_t$gene ] , pka_r [ ovs_pka_t$peak ], ignore.strand = TRUE) 

	# ordering
	ovs_pka_t = ovs_pka_t [ order(ovs_pka_t$peak), ]
	ovs_pka_t = ovs_pka_t [ order(ovs_pka_t$gene), ]
	ovs_pka_t = ovs_pka_t [ !duplicated(ovs_pka_t), ]
	ovs_pka_t = ovs_pka_t[ !is.na(ovs_pka_t$peak), ]
	
	# Peak annotations
	# get coordinates
	ovs_pka_t$chr   = as.vector(seqnames(pka_r [ ovs_pka_t$peak ]))
	ovs_pka_t$start = start(pka_r [ ovs_pka_t$peak ])
	ovs_pka_t$end   = end(pka_r [ ovs_pka_t$peak ])
	ovs_pka_t$strand   = as.character(strand(pka_r [ ovs_pka_t$peak ]))
	# are peaks located in intergenic regions?
	ovs_pka_t$is_promoter   = paste(ovs_pka_t$gene, ovs_pka_t$peak) %in% ovs_pka_pro_t_list
	ovs_pka_t$is_intergenic = ! overlapsAny(pka_r [ ovs_pka_t$peak ], ann_r)
	message(sprintf("%s | %i peak-gene assignments outside promoters", Sys.time(), sum(!ovs_pka_t$is_promoter) ))

	# output
	# return(list(matches = ovs_pka_t, broken_ranges = ovs_extpro_r,broken_dist = dist_ext_to_pro, final_ranges = ann_ext_fil_r))
	return(ovs_pka_t)
	
}

#' Match ATAC peaks to genes based on coordinates and evidence-based promoter regions (optionally)
#' (this is the original function written by Xavi)
#' 
#' @param gff_object either a path to a GFF file with gene annotations, or a preloaded `GRanges` object
#' @param peak_object either a path to a BED file with ATAC peaks, or a preloaded `GRanges` object
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param feature_to_match type of feature in the GFF that will be matched to the ATAC peaks (default: "gene"; "transcript" can be used for GTFs, etc.)
#' @param exclude_genes vector of genes to exclude (e.g. nonexpressed genes); defaults to NULL
#' @param max_tss_dist max distance away from a gene's TSS to check for ATAC peaks that can be assigned to it (default: 20000 bp)
#' @param min_overlap min number of overlapping bases between a gene/gene extended region and an ATAC peak so as to match them(default: 0 bp, i.e. any overlap is sufficient)
#' @param promoter_upstream,promoter_downstream distance upstream/downstream of a TSS region that's used to define promoters (inferred from `gff_object`)
#' @param promoter_object optionally, an additional set of promoters can be provided (will be merged with GFF-defined promoters). Either a path to a BED file, or a preloaded `GRanges` object
#' 
#' @return data.frame with many-to-many matches between genes and peaks
#' 
mta_match_peaks_to_genes_old = function(
	gff_object, peak_object, index_object, 
	list_genes = NULL, 
	feature_to_match = "gene", 
	feature_field = "gene_id", 
	exclude_genes = NULL,     # genes to exclude from gff file (can be nonexpressed genes, for example)
	max_tss_dist = 20000,    # peaks beyond this distance should never be assigned to a particular gene
	min_overlap = 0, # min overlap between a peak and a gene/gene extended region to link them
	promoter_upstream = 200, # these regions are used to define promoters around TSS: peaks within this region get automatically assigned to the relevant gene (and NEVER to any other gene)
	promoter_downstream = 50,
	promoter_object = NULL # a bed file with promoter regions, so that (i) we can directly assign peaks to promoters (and assign them to a gene as if these were TSS, see above), and (ii) we can enforce a rule that prevents peaks to be assigned to a gene if there's a promoter between them
) {
	
	require("rtracklayer")
	require("GenomicRanges")
	require("IRanges")
	
	# temporarily suppress warnings
	options(warn=-1)
	on.exit(options(warn = 0))
	
	### Read data ###
	
	message("match peaks & genes | Load data...")
	
	# read gff 
	if ("character" %in% class(gff_object)) {
		gff_r = rtracklayer::readGFFAsGRanges(gff_object)
	} else if ("GRanges" %in% class(gff_object)) {
		gff_r = gff_object
	} else {
		stop(sprintf("`gff_object` %s has to be either a GRanges object from a GFF, or a path to a GFF", gff_object))
	}
	
	# read genome index
	if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`)", index_object))
	}
	
	# read peaks
	if ("character" %in% class(peak_object)) {
		# load
		pka_r = mta_bed_to_granges(peak_object)
		names(pka_r) = pka_r$name
	} else if ("GRanges" %in% class(peak_object)) {
		pka_r = peak_object
		names(pka_r) = pka_r$name
	} else {
		stop(sprintf("`peak_object` %s has to be either a GRanges object, or a path to a compatible BED file", peak_object))
	}
	
	# if available, load promoter regions
	if ("character" %in% class(promoter_object)) {
		# load
		pro_r = mta_bed_to_granges(promoter_object)
		names(pro_r) = pro_r$name
	} else if ("GRanges" %in% class(promoter_object)) {
		pro_r = promoter_object
		names(pro_r) = pro_r$name
	} else if (!is.null(promoter_object)) {
		stop(sprintf("`promoter_object` %s has to be either a GRanges object, or a path to a compatible BED file", promoter_object))
	}    
	
	
	### Create range objects ###
	# get genes (or any other feature to match)
	ann_r = gff_r[gff_r$type == feature_to_match,]
	names(ann_r) = as.data.frame(ann_r)[,feature_field]
	
	# add chr lengths to genome annotation and peak annotation
	seqlevels(ann_r)  <- as.character(gix$chr)
	seqlengths(ann_r) <- as.integer(gix$length)
	seqlevels(pka_r)  <- as.character(gix$chr)
	seqlengths(pka_r) <- as.integer(gix$length)
	if (!is.null(promoter_object)) {
		seqlevels(pro_r)  <- as.character(gix$chr)
		seqlengths(pro_r) <- as.integer(gix$length)
	}
	
	# exclude genes if necessary
	if (!is.null(exclude_genes)) {
		ann_r = ann_r [ !names(ann_r) %in% exclude_genes ]
	}
	
	# log
	message(sprintf("match peaks & genes | Considering %i features...", length(ann_r)))
	
	# get tss (strand-aware)
	ann_tss_r = GenomicRanges::resize(ann_r, width = 1)
	
	# get promoter regions
	ann_pro_r = GenomicRanges::promoters(ann_tss_r, upstream = promoter_upstream, downstream = promoter_downstream) 
	ann_pro_r = GenomicRanges::trim(ann_pro_r)
	# ann_pro_r = GenomicRanges::GRanges(seqnames = seqnames(ann_pro_r), ranges = ranges(ann_pro_r))
	ann_pro_r$matched_gene = names(ann_pro_r)
	
	# if we have coordinates of promoters:
	# - get promoter regions overlapping with a TSS, and extend gene feature-inferred promoter to the full extend of the peak (better promoter)
	# - also keep promoter regions NOT overlapping with TSS, which will be used to define boundaries around genes. Peaks beyond these boundaries 
	#   won't be assigned to these genes.
	if (!is.null(promoter_object)) {
		
		# which predefined promoters overlap with GFF-inferred promoter regions?
		ovs_pro = IRanges::findOverlapPairs(ann_pro_r, pro_r,  select = "all", type = "any")
		
		# expand GFF-inferred promoters for each gene to include predefined peaks (punion = parallel union)
		ovs_pro_r = GenomicRanges::punion(ovs_pro@first, ovs_pro@second)
		ovs_pro_r = GenomicRanges::trim(ovs_pro_r)
		ovs_pro_r = unique(ovs_pro_r)
		ovs_pro_r$matched_gene = names(ovs_pro_r)
		# add these extended promoters to original ranges
		ann_pro_r [ names(ovs_pro_r) ] = ovs_pro_r
		message(sprintf("match peaks & genes | %s promoters assigned to genes", length(ovs_pro_r)))
		
		# get promoters not assigned to any gene
		ovs_nopro = IRanges::findOverlapPairs(ann_pro_r, pro_r,  select = "all", type = "any")
		ovs_nopro_r = pro_r [ ! names(pro_r) %in% names(ovs_nopro@second) ]
		ovs_nopro_r$matched_gene = "unknown"
		# add unassigned promoters to original ranges
		ann_pro_r = c(ann_pro_r, ovs_nopro_r)
		ann_pro_r = GenomicRanges::trim(ann_pro_r)
		message(sprintf("match peaks & genes | %s promoters not assigned to genes", length(ovs_nopro_r)))
		
	}
	
	# FIRST: find ATAC peaks that are overlapping promoter regions. 
	# These regions get assigned to these genes, straight away.
	ovs_pka_pro = IRanges::findOverlapPairs(pka_r, ann_pro_r [ mcols(ann_pro_r)$matched_gene %in% names(ann_r) ], select="all", minoverlap=min_overlap)
	ovs_pka_pro_t = data.frame(
		gene = names(ovs_pka_pro@second),
		peak = names(ovs_pka_pro@first)
	)
	ovs_pka_pro_t_list = paste(ovs_pka_pro_t$gene, ovs_pka_pro_t$peak)
	message(sprintf("match peaks & genes | %i peak-gene assignments in promoters", nrow(ovs_pka_pro_t)))
	
	
	# SECOND: find ATAC peaks that are located within each gene's extended region.
	# These extended regions are defined by a wider area around the gene (e.g. +/- 20kb), 
	# but are constrained by the presence of promoters of other genes. These other 
	# promoters can be defined from the GFF, or from a precomputed track of promoter 
	# coordinates (e.g. from CHIP-seq H3K4me3 data).
	
	# get extended regions of interest (TSS +/- distance)
	ann_ext_r = GenomicRanges::resize(ann_tss_r, width = width(ann_tss_r) + (max_tss_dist * 2), fix = "center")
	ann_ext_r = GenomicRanges::trim(ann_ext_r)
	
	# subtract promoters from extended regions    
	sub_extpro_r = GenomicRanges::setdiff(ann_ext_r, ann_pro_r, ignore.strand=TRUE)
	# remap fragmented extended regions to their genes of origin
	ovs_extpro = GenomicRanges::findOverlaps(sub_extpro_r, ann_ext_r, select = "all", minoverlap = 1L)
	ovs_extpro_r = sub_extpro_r [ ovs_extpro@from ]
	names(ovs_extpro_r) = names(ann_tss_r [ ovs_extpro@to ])
	strand(ovs_extpro_r) = strand(ann_tss_r [ ovs_extpro@to ])
	# find which extended region fragment is closest the promoter of each gene
	dist_ext_to_pro = data.frame(
		gene = names(ovs_extpro_r),
		dist_to_pro = GenomicRanges::distance(ovs_extpro_r, ann_pro_r [ names(ovs_extpro_r) ]),
		r = ovs_extpro_r
	)
	dist_ext_to_pro_order = as.numeric(as.factor(dist_ext_to_pro$gene))
	dist_ext_to_pro = dist_ext_to_pro [ order(dist_ext_to_pro_order, dist_ext_to_pro$dist_to_pro, decreasing = FALSE), ]
	rownames(dist_ext_to_pro) = 1:nrow(dist_ext_to_pro)
	# keep extended areas adjacent to the promoter of each gene
	# dist_ext_to_pro_f = dist_ext_to_pro [ dist_ext_to_pro$dist_to_pro == 0, ]
	# get extended fragment closest to promoter
	bool_closest = ! duplicated(dist_ext_to_pro$gene)
	dist_ext_to_pro_f1 = dist_ext_to_pro [ bool_closest, ]
	# get extended fragment 2nd-closest to promoter
	dist_ext_to_pro_f2 = dist_ext_to_pro [ ! bool_closest, ]
	dist_ext_to_pro_f2 = dist_ext_to_pro_f2 [ ! duplicated(dist_ext_to_pro_f2$gene), ]
	# merge & reorder
	dist_ext_to_pro_f = rbind(dist_ext_to_pro_f1, dist_ext_to_pro_f2)
	dist_ext_to_pro_f_order = as.numeric(as.factor(dist_ext_to_pro_f$gene))
	dist_ext_to_pro_f = dist_ext_to_pro_f [ order(dist_ext_to_pro_f_order, dist_ext_to_pro_f$dist_to_pro, decreasing = FALSE), ]
	
	# rebuild extended regions
	ann_ext_fil_r = GenomicRanges::GRanges(
		seqnames=dist_ext_to_pro_f$r.seqnames, 
		ranges = IRanges(start = dist_ext_to_pro_f$r.start, end = dist_ext_to_pro_f$r.end), 
		strand = dist_ext_to_pro_f$r.strand)
	names(ann_ext_fil_r) = dist_ext_to_pro_f$gene
	
	# finally, find overlap between extended regions and atac peaks
	ovs_pka_ext = IRanges::findOverlapPairs(pka_r, ann_ext_fil_r, select="all", minoverlap=min_overlap)
	ovs_pka_ext_t = data.frame(
		gene = names(ovs_pka_ext@second),
		peak = names(ovs_pka_ext@first)
	)
	
	
	# THIRD: merge gene-peak assignments from promoters and extended regions
	ovs_pka_t = rbind(ovs_pka_pro_t, ovs_pka_ext_t) 
	ovs_pka_t$dist_to_tss = GenomicRanges::distance( ann_tss_r [ ovs_pka_t$gene ] , pka_r [ ovs_pka_t$peak ] ) 
	ovs_pka_t = ovs_pka_t [ ovs_pka_t$dist_to_tss <= max_tss_dist , ]
	ovs_pka_t = ovs_pka_t [ order(ovs_pka_t$gene), ]
	ovs_pka_t = ovs_pka_t [ !duplicated(ovs_pka_t) , ]
	
	# Peak annotations
	# get coordinates
	ovs_pka_t$chr   = as.vector(seqnames(pka_r [ ovs_pka_t$peak ]))
	ovs_pka_t$start = start(pka_r [ ovs_pka_t$peak ])
	ovs_pka_t$end   = end(pka_r [ ovs_pka_t$peak ])
	# are peaks located in intergenic regions?
	ovs_pka_t$is_promoter   = paste(ovs_pka_t$gene, ovs_pka_t$peak) %in% ovs_pka_pro_t_list
	ovs_pka_t$is_intergenic = ! pka_r [ ovs_pka_t$peak ] %over% ann_r
	message(sprintf("match peaks & genes | %i peak-gene assignments outside promoters", sum(!ovs_pka_t$is_promoter) ))
	
	# output
	# return(list(matches = ovs_pka_t, broken_ranges = ovs_extpro_r,broken_dist = dist_ext_to_pro, final_ranges = ann_ext_fil_r))
	return(ovs_pka_t)
	
}

#' Refine peak to gene assignment based on gene-peak correlation and peak-peak co-accessibility
#' 
#' @param genes_peaks data.frame, output from `mta_match_peaks_to_genes()`
#' @param genes_peaks_corr data.frame with three columns, peak, gene and their corelation
#' @param peaks_peaks data.frame with three columns, peak1, peak2 and their co-accessibility
#' @param coaccess_thr numeric, threshold for peak-peak coaccasibility (should be between 0 and 1)
#' @param delta_corr_thr numeric, absolute decrease of ranked peak-to-gene correlation values
#'   after which the peak-to-gene assignments are removed
#' @param verbose
#' 
mta_refine_peaks_to_genes_by_coaccessibility <- function(
  genes_peaks, genes_peaks_corr, peaks_peaks,
  coaccess_thr = 0.5, delta_corr_thr = 0.2,
  verbose = TRUE
) {

  # peak to gene assignments
  setDT(genes_peaks)
	peak_asign_orig <- unique(genes_peaks[, .(gene, peak)])
	setnames(peak_asign_orig, c("gene2", "peak2"))
    
	# peak to gene correlation
  setDT(genes_peaks_corr) 
	setnames(genes_peaks_corr, c("peak", "gene", "corr"))
    
	# combine
	genes_comb <- merge.data.table(
	  genes_peaks, genes_peaks_corr,
    by = c("peak", "gene"),
    all.x = TRUE, sort = FALSE
  )
	genes_comb[is.na(corr), corr := 0]

  # peak to peak correlation
	setDT(peaks_peaks)
	setnames(peaks_peaks, c("peak1", "peak2", "coaccess"))

  # for each peak assigned to more than one gene,
  # check best correlated genes it is assigned to
  # and coaccessible peaks and genes they are assigned to

  mult_peaks <- genes_comb[,.N, peak][N > 1]$peak
	message(sprintf(
		"There are %s%% (%s / %s) peaks assigned to more than one gene.",
		format(length(mult_peaks) / length(unique(genes_comb$peak)) * 100, digits = 4),
		length(mult_peaks), length(unique(genes_comb$peak))
	))
	genes_comb[, coaccess_corr := numeric()]
  # add coaccessibility to multiply-assigned peaks
	message("Getting co-accessibile peaks for multiply-assigned peaks, this might take a while...")
  for (i in seq_along(mult_peaks)) {
    if ((verbose & i %% 1000 == 0) | i == length(mult_peaks))
	    message(sprintf("... %s/%s", i, length(mult_peaks)))
    p <- mult_peaks[i]
    # genes that a peak is initially assigned to
    gs <- genes_comb[peak == p]$gene
    # coaccessible peaks
    cp <- peaks_peaks[peak1 == p][coaccess > coaccess_thr]$peak2
    # correlation
    gp <- genes_comb[peak %in% cp & gene %in% gs]
    if (nrow(gp) > 0) {
      gc <- gp[, .(coaccess_corr = mean(corr, na.rm = TRUE)), gene][order(-coaccess_corr)]
      gv <- gc[, setNames(coaccess_corr, gene)]
      genes_comb[peak == p, coaccess_corr := gv[gene]]
    }
  }

  # calculate correlations
	# genes_comb[is.na(coaccess_corr), coaccess_corr := corr]
	genes_comb[, mean_corr := mean(
      c(corr, coaccess_corr), na.rm = TRUE
  ), by = seq_len(nrow(genes_peaks))]
	genes_comb[, peak := factor(peak, levels = unique(genes_comb$peak))]
  setorder(genes_comb, peak, -mean_corr)
  genes_comb[, gene_rank := 1:nrow(.SD), by = peak]
  genes_comb[, corr_diff := .SD$mean_corr[1] - mean_corr, by = peak]
  genes_comb[, keep := TRUE][corr_diff > delta_corr_thr, keep := FALSE]
	rmvd_peaks <- unique(genes_comb[keep == FALSE]$peak)
	message(sprintf(
		"Removed assignments for %s%% (%s / %s) peaks initially assigned to multiple genes.",
		format(length(rmvd_peaks) / length(mult_peaks) * 100, digits = 4),
		length(rmvd_peaks), length(mult_peaks)
	))

	# check if initial coaccessibility groups are geting split
	# by gene assignment?

  # get coaccessibility groups and assignments
  setnames(peaks_peaks, "peak1", "peak")
  peaks_peaks_genes <- merge.data.table(
    peaks_peaks, unique(genes_comb[, .(peak, gene, keep)]),
    by = "peak", all.x = TRUE, sort = FALSE, allow.cartesian = TRUE
  )
	peaks_peaks_coacces <- merge.data.table(
		peaks_peaks_genes, peak_asign_orig, by = "peak2",
		all.x = TRUE, sort = FALSE, allow.cartesian = TRUE
	)
	setcolorder(peaks_peaks_coacces, c("peak", "peak2", "coaccess", "gene", "gene2"))
	
	# select only peaks initially assigned to same genes
	coacces_dt <- peaks_peaks_coacces[gene == gene2][, gene2 := NULL]

	# select only co-accessible pairs of peaks
	coacces_dt <- coacces_dt[coaccess > coaccess_thr]
	
	# split co-accessible groups within gene if needed
	message("Getting co-accessibile groups for genes, this might take a while...")
	test_peaks <- unique(coacces_dt$peak)
	for (i in seq_along(test_peaks)) {
		if ((verbose & i %% 1000 == 0) | i == length(test_peaks))
			message(sprintf("... %s/%s", i, length(test_peaks)))
		p1 <- test_peaks[i]
		gs <- coacces_dt[peak == p1]$gene
		for (g in gs) {
			p2 <- coacces_dt[peak == p1 & gene == g]$peak2
			pd <- coacces_dt[peak %in% c(p1, p2) & gene == g]
			pg <- paste(c(sort(unique(c(pd$peak, pd$peak2))), g), collapse = "_")
			coacces_dt[peak == p1 & gene == g, coaccess_group := pg]
		}
	}

	# either all of the co-accessible pairs should be kept, or none
	split_coaccess_group <- unique(coacces_dt[, .(keep, coaccess_group)])[, .N, .(coaccess_group)][N > 1]$coaccess_group

	# which ones to save
	coacces_st <- unique(coacces_dt[coaccess_group %in% split_coaccess_group][keep == FALSE][, .(peak, gene)])
	savd_peaks <- unique(coacces_st$peak)
	message(sprintf(
		"Of those, %s%% (%s / %s) are co-accessible with other peaks assigned to the same gene, we save them.",
		format(length(savd_peaks) / length(rmvd_peaks) * 100, digits = 4),
		length(savd_peaks), length(rmvd_peaks)
	))
	
	genes_coac <- merge.data.table(
	  genes_comb, coacces_st[, saved_by_coaccess := TRUE],
		by = c("peak", "gene"), all = TRUE, sort = FALSE
	)
	genes_coac[is.na(saved_by_coaccess), saved_by_coaccess := FALSE]
	genes_coac[saved_by_coaccess == TRUE, keep := TRUE]
	mult_peaks_new <- genes_coac[keep == TRUE][,.N, peak][N > 1]$peak
	message(sprintf(
		"Now there are %s%% (%s / %s) peaks assigned to more than one gene.",
		format(length(mult_peaks_new) / length(unique(genes_coac$peak)) * 100, digits = 4),
		length(mult_peaks_new), length(unique(genes_coac$peak))
	))

	genes_coac
}


#' Classify promoters as constitutive (CP), specific (SP) or alternative (AP)
#'  
#' @param scatac_peaks GRanges or data.frame with scATAC peaks assigned to genes and cell types
#'   (i.e. it has to have "cell_type" and "gene" columns)
#' @param scrna_peaks GRanges or data.frame with 5' scRNA peaks assigned to cell types
#'   (i.e. it has to have a "cell_type" column)
#' @param chip_peaks GRanges or data.frame with ChIP H3K4me3 peaks assigned to genes
#'   (i.e. it has to have a "gene" column)
#' @param genes GRanges or data.frame with genes annotation
#' @param max_dist_to_tss integer, peaks beyond this distance will never be assigned as promoters
#' @param extend_scrna integer, how much to extend scRNA peaks for overlap with scATC peaks
#' @param extend_chip integer, how much to extend ChIP peaks for overlap with scATC peaks
#' @param class_cols named vector of colors for plotting
#' @param vrerbose logical
mta_class_promoters <- function(
  scatac_peaks,
  scrna_peaks = NULL,
  chip_peaks = NULL,
  pairs_dt = NULL,
  genes,
  max_dist_to_tss = 3000,
  extend_scrna = 250,
  extend_chip = 250,
  class_cols = c(
    "scATAC peaks" = "#7570b3",
    "5' scRNA" = "#e66101",
    "H3K4me3" = "#e6ab02",
    "AP" = "#d01c8b",
    "SP" = "#37b8aa",
    "CP" = "#3935d4"
  ),
  verbose = FALSE
) {
  
  require(ggraph)
  require(tidygraph)
  require(ArchR)

  # LOAD DATA
  
  ## convert to GRanges
  if ("data.frame" %in% class(scatac_peaks)) {
    scatac_gr <- makeGRangesFromDataFrame(scatac_peaks, keep.extra.columns = TRUE)
  } else if ("GRanges" %in% class(scatac_peaks)) {
    scatac_gr <- scatac_peaks
  }
  if ("data.frame" %in% class(scrna_peaks)) {
    scrna_gr <- makeGRangesFromDataFrame(scrna_peaks, keep.extra.columns = TRUE)
  } else if ("GRanges" %in% class(scrna_peaks)) {
    scrna_gr <- scrna_peaks
  } else if (is.null(scrna_peaks)) {
    scrna_gr <- GenomicRanges::GRanges()
  }
  if ("data.frame" %in% class(chip_peaks)) {
    chip_gr <- makeGRangesFromDataFrame(chip_peaks, keep.extra.columns = TRUE)
  } else if ("GRanges" %in% class(chip_peaks)) {
    chip_gr <- chip_peaks
  } else if (is.null(chip_peaks)) {
    chip_gr <- GenomicRanges::GRanges()
  }
  if ("data.frame" %in% class(genes)) {
    genes_gr <- makeGRangesFromDataFrame(genes, keep.extra.column = TRUE)
  } else if ("GRanges" %in% class(genes)) {
    genes_gr <- genes
  }
  ## chromosomes
  seqlvl <- seqlevels(scatac_gr)
  seqlevels(scrna_gr, pruning.mode="coarse") <- seqlvl
  seqlevels(chip_gr, pruning.mode="coarse") <- seqlvl
  
  ## add gene strand to scATAC peaks
  genes_str <- as.data.frame(genes_gr)[, c("gene","strand")]
  peaks_str <- genes_str[match(scatac_gr$gene,genes_str$gene),]
  peaks_str[is.na(peaks_str$strand),"strand"] <- "*"
  strand(scatac_gr) <- peaks_str$strand
  
  ## overlap scATAC with scRNA peaks
  scrna_gr <- ArchR::extendGR(scrna_gr, upstream = extend_scrna, downstream = 0)
  ovl <- findOverlaps(query=scatac_gr, subject=scrna_gr, minoverlap=10, ignore.strand=FALSE)
  
  ### scATAC peaks overlapping scRNA peaks
  ovl_gr <- scatac_gr[queryHits(ovl)]
  meta_df <- mcols(scrna_gr[subjectHits(ovl)])
  mcols(ovl_gr)[,"scRNA_peak"] <- meta_df$peak
  mcols(ovl_gr)[,"scRNA_strand"] <- strand(scrna_gr[subjectHits(ovl)])
  mcols(ovl_gr)[,"scRNA_cell_type"] <- meta_df$cell_type
  
  ### scATAC peaks not overlapping scRNA peaks
  nonovl <- setdiff(seq_along(scatac_gr), queryHits(ovl))
  nonovl_gr <- scatac_gr[nonovl]
  
  ### all scATAC peaks with scRNA info
  all_gr <- sort(c(nonovl_gr,ovl_gr))
  seqlevels(all_gr) <- seqlvl
  all_gr <- sort(all_gr)
  
  ### keep only entries where scRNA peak matches scATAC cell type
  if (!is.null(scrna_peaks)) {
    if (any(is.na(all_gr$cell_type))) all_gr[is.na(all_gr$cell_type)]$cell_type <- ""
    if (any(is.na(all_gr$scRNA_cell_type)))  all_gr[is.na(all_gr$scRNA_cell_type)]$scRNA_cell_type <- ""
    all_gr$id <- paste(all_gr$gene, all_gr$peak, all_gr$cell_type)
    ovl_gr <- all_gr[all_gr$cell_type == all_gr$scRNA_cell_type | is.na(all_gr$scRNA_peak)]
    novl_gr <- all_gr[all_gr$cell_type != all_gr$scRNA_cell_type]
    novl_gr <- novl_gr[!novl_gr$id %in% ovl_gr$id]
    novl_gr$scRNA_peak <- ""
    novl_gr$scRNA_cell_type <- ""
    peaks_gr <- c(ovl_gr, novl_gr)
  } else {
    peaks_gr <- all_gr
    peaks_gr$scRNA_peak <- ""
  }
  
  ## add ChIP peaks overlap
  gs <- unique(peaks_gr$gene)
  gs <- gs[gs %in% peaks_gr$gene]
  gs <- gs[gs!=""]
  if (!is.null(chip_peaks)) {
    message(Sys.time()," | ","Overlap scATAC and ChIP peaks per gene")
    chip_gr <- extendGR(chip_gr, upstream = extend_chip, downstream = extend_chip)
    if (is.null(pairs_dt)) {
      pairs_dt <- rbindlist(lapply(1:length(gs), function(i) {
        if (verbose & i%%1000==0) message(Sys.time()," | ",i,"/",length(gs))
        g <- gs[i]
        g_gr <- peaks_gr[peaks_gr$gene==g]
        c_gr <- chip_gr[chip_gr$gene==g]
        ovl <- findOverlaps(query=g_gr, subject=c_gr, minoverlap=1, ignore.strand=TRUE)
        if (length(ovl)>0) {
          ovl_peaks <- unique(g_gr[queryHits(ovl)]$peak)
          data.table(peak=ovl_peaks)[,H3K4me3_peak:=c_gr$peak][,gene:=g][]
        }
      }))
    } else {
      setnames(pairs_dt, c("peak", "H3K4me3_peak", "gene"))
    }
    all_dt <- merge.data.table(
      as.data.table(mcols(peaks_gr)),
      pairs_dt,
      by=c("gene", "peak"), all.x=TRUE
    )
  } else {
    all_dt <- as.data.table(mcols(peaks_gr))[, H3K4me3_peak := ""]
  }

  ## prepare input data
  dt <- all_dt[gene != ""][cell_type != ""]
  dt[is.na(H3K4me3_peak), H3K4me3_peak := ""]
  dt[is.na(scRNA_peak), scRNA_peak := ""]
  dt[, strand := as.vector(strand(genes_gr[match(dt$gene, genes_gr$gene)]))]
  dt[,c("scRNA_cell_type","scRNA_strand","id"):=NULL]
  genes_gr <- genes_gr[match(dt$gene, genes$gene)]

  ## peaks that are downstream from gene body should not be considered potential promoters
  promoters_gr <- promoters(genes[match(dt$gene, genes$gene)], upstream = 0, downstream = 1)
  minus_strand_genes <- as.vector(strand(genes_gr) == "-")
  tss_coord <- start(genes_gr)
  tss_coord[minus_strand_genes] <- end(genes_gr)[minus_strand_genes]
  dt[, tss := tss_coord]
  tes_coord <- end(genes_gr)
  tes_coord[minus_strand_genes] <- start(genes_gr)[minus_strand_genes]
  dt[, tes:=tes_coord]
  
  peaks_gr <- scatac_gr[match(dt$peak,scatac_gr$peak)]
  #minus_strand_peaks <- as.vector(strand(peaks_gr)=="-")
  pst <- start(peaks_gr)
  #pst[minus_strand_peaks] <- end(peaks_gr)[minus_strand_peaks]
  pen <- end(peaks_gr)
  #pen[minus_strand_peaks] <- start(peaks_gr)[minus_strand_peaks]
  dt[,peak_start:=pst]
  dt[,peak_end:=pen]
  
  dt[strand=="+", is_downstream_from_gene:=peak_start>tes]
  dt[strand=="-", is_downstream_from_gene:=peak_start<tes]
  
  #dt[strand=="+", is_downstream_from_tss:=peak_start<tss]
  #dt[strand=="-", is_downstream_from_tss:=peak_end>tss]
  
  ## distance peak-TSS
  dt$dist_to_tss <- GenomicRanges::distance(
    promoters_gr, 
    peaks_gr, 
    ignore.strand=TRUE
  ) 
  
  ## distance peak-H3K4me3
  ids <- which(dt$gene %in% chip_gr$gene)
  chip_mid <- chip_gr[match(dt[ids]$gene,chip_gr$gene)]
  hw <- width(chip_mid)/2
  start(chip_mid) <- start(chip_mid)+hw
  end(chip_mid) <- end(chip_mid)-hw
  
  peaks_mid <- scatac_gr[match(dt[ids]$peak,scatac_gr$peak)]
  hw <- width(peaks_mid)/2
  start(peaks_mid) <- start(peaks_mid)+hw
  end(peaks_mid) <- end(peaks_mid)-hw
  
  dt[ids, peak_mid:=start(peaks_mid)]
  dt[ids, chip_mid:=start(chip_mid)]
  dt[,strand:=as.character(strand(genes[match(dt$gene,genes$gene)]))]
  ids_dist <- distance(peaks_mid, chip_mid, ignore.strand=TRUE)
  dt[ids, dist_to_k4me3:=ids_dist]
  dt[strand=="+" & peak_mid<chip_mid, dist_to_k4me3:=-1*dist_to_k4me3]
  dt[strand=="-" & peak_mid>chip_mid, dist_to_k4me3:=-1*dist_to_k4me3]
  dt[,c("peak_mid","chip_mid"):=NULL]
  
  ## decision tree
  n_peaks <- length(unique(dt$peak))
  n_peaks_all <- length(unique(all_dt$peak))
  n_genes <- length(unique(dt$gene))
  n_cell_types <- length(unique(dt$cell_type))
  message(sprintf(
    "%s | %i peaks (%.2f%%) assigned to %i genes and %i cell types",
    Sys.time(),n_peaks,n_peaks/n_peaks_all*100,n_genes,n_cell_types
  ))
  message(Sys.time()," | ","Classify peaks assigned to genes")
  
  gs <- unique(dt$gene)
  outl <- lapply(1:length(gs), function(i) {
    
    if (verbose & i%%1000==0) message(Sys.time()," | ",i,"/",length(gs))
    g <- gs[i]

	# peaks table
    dtg <- unique(dt[gene==g])
    if ("peak_adult" %in% colnames(dtg)) dtg[,peak_adult:=NULL]
	if ("peak_gastrula" %in% colnames(dtg)) dtg[,peak_gastrula:=NULL]
	dtg <- unique(dtg)
	
    # table to save tree path
    dtt <- data.table(gene=g)
    
    # don't consider peaks downstream from gene and further than max distance
    dtg_d <- dtg[is_downstream_from_gene==TRUE | dist_to_tss>max_dist_to_tss][,promoter:="NO"]
    dtg <- dtg[!peak %in% dtg_d$peak]
    
    # count peaks
    dtg[,n_peaks:=length(unique(dtg$peak))]
    dtg[,n_peaks_in_cell_type:=length(unique(.SD$peak)), cell_type]
    dtg[,n_cell_type:=length(unique(.SD$cell_type)), peak]
    
    if (nrow(dtg)>1) {
     
      # 1) one peak in each cell type?
      if (max(dtg$n_peaks_in_cell_type)==1) {
        dtt[,n1_one_peak:=TRUE]
        
        # 2a) same peak in all available cell types?
        if (all(dtg$n_peaks==1)) {
          dtt[,n2a_same_peak_across_cell_types:=TRUE]
          
          # 3a) in all cell types?
          if (all(dtg$n_cell_type==n_cell_types)) {
            dtt[,n3a_all_cell_types:=TRUE]
            dtg[,promoter:="CP"]
          } else {
            dtt[,n3a_all_cell_types:=FALSE]
            dtg[,promoter:="SP"]
          }
            
        } else {
          dtt[,n2a_same_peak_across_cell_types:=FALSE]
          
          # 3b) any peak with 5' overlap?
          if (any(dtg$scRNA_peak!="")) {
            dtt[,n3b_any_peak_overlaps_rna:=TRUE]
            
            # 4a) different peaks with 5' overlap?
            if (length(unique(dtg[scRNA_peak!=""]$peak))>1) {
              dtt[,n4a_different_peaks_overlap_rna:=TRUE]
              dtg[,promoter:="AP+"]
            } else {
              dtt[,n4a_different_peaks_overlap_rna:=FALSE]
              dtg[,promoter:="AP"]
            }
            
          } else {
            dtt[,n3b_any_peak_overlaps_rna:=FALSE]
            dtg[,promoter:="AP"]
          }
          
        }
        
      } else {
        dtt[,n1_one_peak:=FALSE]
        
        # 2b) any peak with H3K4me3 overlap?
        if (any(dtg$H3K4me3_peak!="")) {
          
          dtt[,n2b_any_peak_overlaps_h3k4me3:=TRUE]
          
          # 3c) many peaks with H3K4me3 overlap?
          if (length(unique(dtg[H3K4me3_peak!=""]$peak))>1) {
            dtt[,n3c_different_peaks_overlap_h3k4me3:=TRUE]
            
            # 4b) any peak with 5' scRNA overlap?
            if (nrow(dtg[H3K4me3_peak!="" & scRNA_peak!=""])>0) {
              dtt[,n4b_any_peak_overlaps_rna:=TRUE]
              
              # 5a) many peaks with 5' scRNA overlap in all available cell types?
              if (length(unique(dtg[H3K4me3_peak!="" & scRNA_peak!=""]$peak))>1) {
                dtt[,n5a_many_peaks_overlap_rna:=TRUE]
                promoter_peaks <- unique(dtg[H3K4me3_peak!="" & scRNA_peak!=""][dist_to_tss<max_dist_to_tss]$peak)
                
                # 6a) any alternative peaks?
                alternative_peaks <- unlist(sapply(promoter_peaks, function(p) {
                  ct <- dtg[peak == p, cell_type]
                  if (!any(setdiff(promoter_peaks, p) %in% dtg[cell_type %in% ct]$peak)) 
                    p
                }, USE.NAMES = TRUE, simplify = FALSE))
                if (length(alternative_peaks)>1) {
                  dtt[,n6a_alternative:=TRUE]
                  dtg[peak %in% alternative_peaks, promoter:="AP+"]
                } else {
                  dtt[,n6a_alternative:=FALSE]
                  # select the one peak which is in more cell types
                  promoter_peak <- dtg[peak %in% promoter_peaks,][order(-n_cell_type)]$peak[1]
                  
                  # 7a) in all cell types?
                  if (all(dtg[peak==promoter_peak]$n_cell_type==n_cell_types)) {
                    dtt[,n7a_all_cell_types:=TRUE]
                    dtg[peak==promoter_peak, promoter:="CP"]
                  } else {
                    dtt[,n7a_all_cell_types:=FALSE]
                    dtg[peak==promoter_peak, promoter:="SP"]
                  }
                }
                
              } else {
                dtt[,n5a_many_peaks_overlap_rna:=FALSE]
                
                # 6b) in all cell types?
                promoter_peak <- unique(dtg[H3K4me3_peak!="" & scRNA_peak!=""]$peak)
                if (all(dtg[peak==promoter_peak]$n_cell_type==n_cell_types)) {
                  dtt[,n6b_all_cell_types:=TRUE]
                  dtg[peak==promoter_peak, promoter:="CP"]
                } else {
                  dtt[,n6b_all_cell_types:=FALSE]
                  dtg[peak==promoter_peak, promoter:="SP"]
                }
                
              }
              
            } else {
              dtt[,n4b_any_peak_overlaps_rna:=FALSE]
              
              # 5b) any peak across all available cell types?
              n_ct <- length(unique(dtg$cell_type))
              if (any(dtg[H3K4me3_peak!=""]$n_cell_type==n_ct)) {
                dtt[,n5b_in_all_available_cell_types:=TRUE]
                promoter_peaks <- unique(dtg[H3K4me3_peak!=""][n_cell_type==n_ct]$peak)
                
                # select closest upstream of H3K4me3
                promoter_peak <- dtg[peak %in% promoter_peaks][order(-dist_to_k4me3)]$peak[1]
              
                # 6c) in all cell types?
                if (all(dtg[peak==promoter_peak]$n_cell_type==n_cell_types)) {
                  dtt[,n6c_all_cell_types:=TRUE]
                  dtg[peak==promoter_peak, promoter:="CP"]
                } else {
                  dtt[,n6c_all_cell_types:=FALSE]
                  dtg[peak==promoter_peak, promoter:="SP"]
                }
                
              } else {
                dtt[,n5b_in_all_available_cell_types:=FALSE]
                promoter_peaks <- dtg[H3K4me3_peak!=""][dist_to_tss<max_dist_to_tss]$peak
                
                # 6e) any alternative peaks?
                alternative_peaks <- unlist(sapply(promoter_peaks, function(p) {
                  ct <- dtg[peak == p, cell_type]
                  if (!any(setdiff(promoter_peaks, p) %in% dtg[cell_type %in% ct]$peak)) 
                    p
                }, USE.NAMES = TRUE, simplify = FALSE))
                
                if (length(alternative_peaks)>1) {
                  dtt[,n6e_alternative:=TRUE]
                  dtg[peak %in% alternative_peaks, promoter:="AP"]
                } else {
                  dtt[,n6e_alternative:=FALSE]
                  # select closest upstream of H3K4me3
                  promoter_peak <- dtg[peak %in% promoter_peaks][dist_to_k4me3<0][order(-dist_to_k4me3)]$peak[1]
                  dtg[peak == promoter_peak, promoter:="SP"]
                }
                
              }
  
            }
            
          } else {
            dtt[,n3c_different_peaks_overlap_h3k4me3:=FALSE]
            
            # 4c) in all cell types?
            if (nrow(dtg[H3K4me3_peak!="" & n_cell_type==n_cell_types])>0) {
              dtt[,n4c_all_cell_types:=TRUE]
              dtg[H3K4me3_peak!="",promoter:="CP"]
            } else {
              dtt[,n4c_all_cell_types:=FALSE]
              dtg[H3K4me3_peak!="",promoter:="SP"]
            }
            
          }
          
        } else {
          dtt[,n2b_any_peak_overlaps_h3k4me3:=FALSE]
          
          # 3d) any peak with 5' scRNA overlap?
          if (nrow(dtg[scRNA_peak!=""])>0) {
            dtt[,n3d_any_peak_overlaps_rna:=TRUE]
            promoter_peaks <- unique(dtg[scRNA_peak!=""]$peak)
            
            # 4d) many peaks with 5' scRNA overlap?
            if (length(unique(dtg[scRNA_peak!=""]$peak))>1) {
              dtt[,n4d_many_peaks_overlap_rna:=TRUE]
              promoter_peaks <- unique(dtg[scRNA_peak!=""]$peak)

              # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
              # 
              # # 5c) any alternative peaks?
              # alternative_peaks <- unlist(sapply(promoter_peaks, function(p) {
              #   ct <- dtg[peak == p, cell_type]
              #   if (!any(setdiff(promoter_peaks, p) %in% dtg[cell_type %in% ct]$peak))
              #     p
              # }, USE.NAMES = TRUE, simplify = FALSE))
              # if (length(alternative_peaks)>1) {
              #   dtt[,n5c_alternative:=TRUE]
              #   dtg[peak %in% alternative_peaks, promoter:="AP"]
              # } else {
              #   dtt[,n5c_alternative:=FALSE]
              #   # select the one peak which is in more cell types
              #   promoter_peak <- dtg[peak %in% promoter_peaks][order(-n_cell_type)]$peak[1]
              # 
              #   # 6d) in all cell types?
              #   if (all(dtg[peak==promoter_peak]$n_cell_type==n_cell_types)) {
              #     dtt[,n6d_all_cell_types:=TRUE]
              #     dtg[peak==promoter_peak, promoter:="CP"]
              #   } else {
              #     dtt[,n6d_all_cell_types:=FALSE]
              #     dtg[peak==promoter_peak, promoter:="SP"]
              #   }
              # }
              # 
              # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
              
              # 5c) any peak in all considered cell types?
              ncts <- length(unique(dtg$cell_type))
              if (nrow(dtg[peak %in% promoter_peaks][n_cell_type == ncts])>0) {
                dtt[,n5c_alternative:=FALSE]
                promoter_peak <- dtg[peak %in% promoter_peaks][n_cell_type == ncts]$peak[1]

                # 6d) in all cell types?
                if (all(dtg[peak==promoter_peak]$n_cell_type==n_cell_types)) {
                  dtt[,n6d_all_cell_types:=TRUE]
                  dtg[peak==promoter_peak, promoter:="CP"]
                } else {
                  dtt[,n6d_all_cell_types:=FALSE]
                  dtg[peak==promoter_peak, promoter:="SP"]
                }
              } else {
                # if a cell type has more than one potential promoter
                # select the more frequent one (n_cell_type)
                # alternatively, the one closest to TSS (dist_to_tss)
                apdt <- dtg[peak %in% promoter_peaks][scRNA_peak!=""][order(-n_cell_type)][
                  , .SD[1], cell_type][
                    , id := paste(peak, cell_type)]
                promoter_peak <- unique(apdt$peak)
                if (length(promoter_peak) == 1) {
                  dtt[,n5c_alternative:=FALSE]
                  dtg[peak==promoter_peak, promoter:="SP"]
                } else {
                  dtt[,n5c_alternative:=TRUE]
                  promoter_ids <- apdt$id
                  for (apid in promoter_ids) {
                    ap <- strsplit(apid, " ")[[1]]
                    appk <- ap[1]
                    apct <- ap[2]
                    dtg[peak %in% appk & cell_type == apct, promoter:="AP"]
                  }
                  # cell types that are not supported by 5' RNA but AP peak(s) open
                  promoterless_cts <- setdiff(
                    unique(dtg[is.na(promoter)]$cell_type),
                    unique(dtg[!is.na(promoter)]$cell_type)
                  )
                  if (length(promoterless_cts) > 0) {
                    dtg[cell_type %in% promoterless_cts & peak %in% promoter_peaks, promoter:="AP"]
                    promoterless_cts <- setdiff(
                      unique(dtg[is.na(promoter)]$cell_type),
                      unique(dtg[!is.na(promoter)]$cell_type)
                    )
                    # any cell types still left without a promoter?
                    if (length(promoterless_cts) > 0) {
                      apdt <- dtg[cell_type %in% promoterless_cts][order(-n_cell_type, dist_to_tss)][
                        , .SD[1], cell_type][
                          , id := paste(peak, cell_type)]
                      promoter_ids <- apdt$id
                      for (apid in promoter_ids) {
                        ap <- strsplit(apid, " ")[[1]]
                        appk <- ap[1]
                        apct <- ap[2]
                        dtg[peak %in% appk & cell_type == apct, promoter:="AP"]
                      }
                    }
                  }
                }
                dtt[,n6d_all_cell_types:=FALSE]
              }

              # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #   
            } else {
              dtt[,n4d_many_peaks_overlap_rna:=FALSE]
              promoter_peak <- unique(dtg[scRNA_peak!=""]$peak)
              
              # 5d) in all cell types?
              if (all(dtg[peak==promoter_peak]$n_cell_type==n_cell_types)) {
                dtt[,n5d_all_cell_types:=TRUE]
                dtg[peak==promoter_peak, promoter:="CP"]
              } else {
                dtt[,n5d_all_cell_types:=FALSE]
                dtg[peak==promoter_peak, promoter:="SP"]
              }
              
            }
            
              
          } else {
            dtt[,n3d_any_peak_overlaps_rna:=FALSE]
            
            # 4e) any peak across all available cell types?
            n_ct <- length(unique(dtg$cell_type))
            if (any(dtg$n_cell_type==n_ct)) {
              promoter_peaks <- dtg[n_cell_type==n_ct]$peak
              dtt[,n4e_in_all_available_cell_types:=TRUE]
              # if multiple, select closest to TSS
              promoter_peak <- dtg[peak %in% promoter_peaks][order(dist_to_tss)]$peak[1]
              
              # 5e) in all cell types?
              if (all(dtg[peak==promoter_peak]$n_cell_type==n_cell_types)) {
                dtt[,n5e_all_cell_types:=TRUE]
                dtg[peak==promoter_peak, promoter:="CP"]
              } else {
                dtt[,n5e_all_cell_types:=FALSE]
                dtg[peak==promoter_peak, promoter:="SP"]
              }
              
            } else {
              dtt[,n4e_in_all_available_cell_types:=FALSE]
              
              # select the peak which is closest to TSS
              promoter_peaks <- unique(dtg[order(dist_to_tss)][,.SD[1],cell_type][dist_to_tss<max_dist_to_tss]$peak)
  
              # 5f) any alternative peaks?
              alternative_peaks <- unlist(sapply(promoter_peaks, function(p) {
                ct <- dtg[peak == p, cell_type]
                if (!any(setdiff(promoter_peaks, p) %in% dtg[cell_type %in% ct]$peak)) 
                  p
              }, USE.NAMES = TRUE, simplify = FALSE))
              
              if (length(alternative_peaks)>1) {
                dtt[,n5f_alternative:=TRUE]
                dtg[peak %in% alternative_peaks, promoter:="AP"]
              } else {
                dtt[,n5f_alternative:=FALSE]
                # select the peak which is in most cell types
                promoter_peak <- dtg[peak %in% promoter_peaks][order(-n_cell_type)][,.SD[1]]$peak
                dtg[peak == promoter_peak, promoter:="SP"]
              }
            }
            
          }
          
        }
        
      }
    } else {
      
      dtg[,promoter:="SP"]
      
    }
    
    dtg <- rbindlist(list(dtg,dtg_d), fill=TRUE, use.names=TRUE)
    
    # cell types without promoters?
    dtg[is.na(promoter),promoter:="NO"]
    dtg[,promoter:=factor(promoter,levels = c("NO","AP","AP+","SP","CP"))]
    proms <- dtg[order(-promoter), .SD[1], cell_type][,.(cell_type,promoter)]
    if (any(proms$promoter == "NO")) {
      orphan_cell_types <- unique(proms[promoter=="NO"]$cell_type)
      orphan_proms_dt <- dtg[cell_type %in% orphan_cell_types]
      # select closest to TSS
      orphan_proms <- orphan_proms_dt[!(is_downstream_from_gene==TRUE | dist_to_tss>max_dist_to_tss)][order(dist_to_tss),.SD[1],cell_type]
      if (nrow(orphan_proms)>0) {
        for  (i in 1:nrow(orphan_proms)) {
          dtg[cell_type == orphan_proms[i]$cell_type & peak == orphan_proms[i]$peak, promoter:="SP"]
        }
        dtg[promoter!="NO",promoter:="SP"]
        # remove cell types that only have peak(s) too far from TSS
        orphans_rm <- orphan_proms_dt[!cell_type %in% orphan_proms$cell_type]
        if (nrow(orphans_rm)>0) {
          dtg <- dtg[!cell_type %in% orphans_rm$cell_type]
        }
      }
    }
    
    list(dtg, dtt)
    
  })
  
  # parse annotation results
  gdt_l <- lapply(outl, function(i) i[[1]])
  gdt <- rbindlist(gdt_l, use.names=TRUE, fill=TRUE)
  gdt[is.na(promoter), promoter:="NO"]
  
  # parse counting results
  dtt_l <- lapply(outl, function(i) i[[2]])
  dtt <- rbindlist(dtt_l, use.names=TRUE, fill=TRUE)

  # summarize tree nodes
  ndt <- melt.data.table(dtt, id.vars="gene", variable.name="node")[,.N,.(node,value)][,value:=factor(value,levels=c(TRUE,FALSE))][order(node,value)]
  ndt <- ndt[!is.na(value)]
  ndt[,level:=str_extract(node,"(?<=n)\\d")]
  ndt[,perc:=N/sum(ndt[1:2]$N)]
  ndt[,node_N:=sum(.SD$N),node]
  ndt[,node_perc:=sum(.SD$perc),node]
  ndt[,node_short:=str_extract(node,"n\\d[a-z]*")]
  setnames(ndt,"node","node_long")
  setnames(ndt,"node_short","node")
  setcolorder(ndt, c("node","node_long"))
  ndt[
  ,desc:=str_remove(node_long,"n\\d[a-z]*_")][
    ,desc:=str_replace_all(desc,"_"," ")][
      ,desc:=paste0(str_to_sentence(desc),"?")][
        ,desc:=str_replace_all(desc,"h3k4me3","H3K4me3")][
          ,desc:=str_replace_all(desc," rna"," 5' scRNA")]
  
  # tree elements
  ## edges
  edges_dt <- do.call("rbind", list(
    c("n1","n2a", "YES"),
    c("n1","n2b", "NO"),
    c("n2a","n3a", "YES"),
    c("n2a","n3b", "NO"),
    c("n2b","n3c", "YES"),
    c("n2b","n3d", "NO"),
    c("n3b","n4a", "YES"),
    c("n3d","n4d", "YES"),
    c("n3d","n4e", "NO"),
    c("n3c","n4b", "YES"),
    c("n3c","n4c", "NO"),
    c("n4b","n5a", "YES"),
    c("n4b","n5b", "NO"),
    c("n4d","n5c", "YES"),
    c("n4d","n5d", "NO"),
    c("n4e","n5e", "YES"),
    c("n4e","n5f", "NO"),
    c("n5a","n6a", "YES"),
    c("n5a","n6b", "NO"),
    c("n5b","n6c", "YES"),
    c("n5b","n6e", "NO"),
    c("n5c","n6d", "NO"),
    c("n6a","n7a", "NO")
  ))
  edges_dt <- as.data.table(as.data.frame(edges_dt))
  setnames(edges_dt, c("node","child","decision"))
  tips_dt <- copy(ndt)[,child:=paste0(node,"_",value)][,.(node,child)]
  tips_dt <- tips_dt[
    node %in% c("n3a","n4a","n4c","n5d","n5f","n5e","n6b","n6c","n6d","n6e","n7a") |
      (node == "n3b" & child == "n3b_FALSE") |
      (node == "n5c" & child == "n5c_TRUE") |
      (node == "n6a" & child == "n6a_TRUE") 
  ]
  tips_dt[,decision:=ifelse(grepl("TRUE",child), "YES", "NO")]
  edges_dt <- rbindlist(list(edges_dt, tips_dt))
  ## nodes
  nodes_dt <- data.table(node = unique(c(edges_dt$node, edges_dt$child)))
  nodes_inner_dt <- unique(ndt[,.(node,node_long,node_N,node_perc,desc)])
  nodes_outer_dt <- merge.data.table(
    tips_dt[,.(child)], 
    copy(ndt)[,child:=paste(node,value,sep="_")][,node_long_child:=paste(node_long,value,sep="_")], 
    by="child", all.x=TRUE
  )
  proms <- c(
    "n3a_TRUE" = "CP",
    "n3a_FALSE" = "SP",
    "n3b_FALSE" = "AP",
    "n4a_TRUE" = "AP+",
    "n4a_FALSE" = "AP",
    "n4c_TRUE" = "CP",
    "n4c_FALSE" = "SP",
    "n5c_TRUE" = "AP",
    "n5d_TRUE" = "CP",
    "n5d_FALSE" = "SP",
    "n5e_TRUE" = "CP",
    "n5e_FALSE" = "SP",
    "n5f_TRUE" = "AP",
    "n5f_FALSE" = "SP",
    "n6a_TRUE" = "AP+",
    "n6b_TRUE" = "CP",
    "n6b_FALSE" = "SP",
    "n6c_TRUE" = "CP",
    "n6c_FALSE" = "SP",
    "n6d_TRUE" = "CP",
    "n6d_FALSE" = "SP",
    "n6e_TRUE" = "AP",
    "n6e_FALSE" = "SP",
    "n7a_TRUE" = "CP",
    "n7a_FALSE" = "SP"
  )
  nodes_outer_dt[,node_long_child:=proms[child]]
  nodes_outer_dt <- nodes_outer_dt[,.(child,node_long_child,N,perc,node_long_child)]
  all_nodes_dt <- rbindlist(list(nodes_inner_dt, nodes_outer_dt), use.names=FALSE)
  nodes_dt <- merge.data.table(nodes_dt, all_nodes_dt, by="node", all.x=TRUE)
  nodes_dt[is.na(node_N),node_N:=0][is.na(node_perc),node_perc:=0]
  nodes_dt[,class:="scATAC peaks"][
    grepl("scRNA",desc),class:="5' scRNA"][
      grepl("H3K4me3",desc),class:="H3K4me3"][
        grepl("TRUE|FALSE",node),class:=str_extract(desc,"[A-Z]+")]
  
  ## graph
  graph_dt <- tbl_graph(nodes = nodes_dt, edges = edges_dt, directed = TRUE)
  gg <- ggraph(graph_dt, layout = 'dendrogram', circular = FALSE) + 
    geom_edge_elbow(aes(color=decision)) +
    geom_node_point(aes(size = node_N, fill = class), color = "black", shape = 21) + 
    scale_size_continuous(range = c(1, 30), guide="none") +
    scale_fill_manual(values = class_cols) +
    scale_edge_colour_manual(values = c("YES"="#31a354", "NO"="#d7191c")) +
    geom_node_text(aes(label = sprintf("%.0f%%\n(%i)", node_perc*100, node_N)), repel = FALSE) +
    geom_node_text(aes(label = str_wrap(desc, width = 25)), nudge_y = -0.35, repel = FALSE) +
    guides(fill = guide_legend(override.aes = list(size=5))) +
    theme_void()
    
  # genes per nodes
  setnames(dtt, c("gene", str_extract(colnames(dtt)[-1],"n\\d[a-z]*")))
  
  # output
  res_list_filtered <- list(
    annotaiton = gdt,
    tree_graph = gg,
    tbl_graph = graph_dt,
    nodes = ndt,
    nodes_genes = dtt
  )
  return(res_list_filtered)
}

#' Calculate gene scores as weighted sum of peaks
#' 
#' @param genes_peaks_table data.frame, output from `mta_match_peaks_to_genes()`
#' @param gff_object either a path to a GFF file with gene annotations, or a preloaded `GRanges` object
#' 	it should contain gene name in "symbol" column
#' @param peak_object either a path to a BED file with ATAC peaks, or a preloaded `GRanges` object
#' @param peaks_mat matrix, rows are peaks, columns are cells
#' @param weight_peaks logical, whether to weight peaks by distance and variability (if available)
#' @param cells_groups data.frame with at least two columns, the first one cells and the second groups for calculating peak variability across groups
mta_gene_scores <- function(
    genes_peaks_table, 
    gff_object, 
    peak_object, 
    peaks_mat,
	weight_peaks = TRUE,
    cells_groups = NULL
) {
  
  require(ineq)
  
  # read gff 
  if ("character" %in% class(gff_object)) {
    gff_r = rtracklayer::readGFFAsGRanges(gff_object)
  } else if ("GRanges" %in% class(gff_object)) {
    gff_r = gff_object
  } else {
    stop(sprintf("`gff_object` %s has to be either a GRanges object from a GFF, or a path to a GFF", gff_object))
  }
  
  # read peaks
  if ("character" %in% class(peak_object)) {
    # load
    pka_r = mta_bed_to_granges(peak_object)
    names(pka_r) = pka_r$name
  } else if ("GRanges" %in% class(peak_object)) {
    pka_r = peak_object
    names(pka_r) = pka_r$name
  } else {
    stop(sprintf("`peak_object` %s has to be either a GRanges object, or a path to a compatible BED file", peak_object))
  }
  
  # add gene and peak info to table used to calculate weights
  message("Matching peaks")
  assigned_peaks <- pka_r[match(genes_peaks_table$peak,pka_r$name)]
  assigned_genes <- gff_r[match(genes_peaks_table$gene,gff_r$symbol)]
  assigned_tss <- resize(assigned_genes, width = 1)
  DT <- as.data.table(genes_peaks_table)
  setnames(DT, "seqnames", "chr")
  DT[,gene_chr:=as.character(seqnames(assigned_genes))]
  DT[,gene_start:=start(assigned_genes)]
  DT[,gene_end:=end(assigned_genes)]
  DT[,gene_length:=width(assigned_genes)]
  DT[,gene_strand:=as.character(strand(assigned_genes))]
  DT[,gene_tss:=start(assigned_tss)]
  DT[,peak_chr:=as.character(seqnames(assigned_peaks))]
  DT[,peak_start:=start(assigned_peaks)]
  DT[,peak_end:=end(assigned_peaks)]
  DT[,peak_width:=width(assigned_peaks)]
  DT[,peak_midpoint:=peak_start + peak_width/2]
  DT[,is_gene_body:=countOverlaps(
	GRanges(chr, IRanges(start, end)),
	GRanges(gene_chr, IRanges(gene_start, gene_end))
  ) > 0]
  
  # peaks distance weights: exponentially decaying function
  if (weight_peaks==TRUE) {
	DT[,w_dist:=exp(1)^(-1*dist_to_tss/5000)+exp(1)^-1]
  } else {
	DT[,w_dist:=1]
  }
  
  # inside gene body set distance weight to 1
  DT[is_gene_body==TRUE,w_dist:=1]

  # peaks variability weight: e ^ z-standardized Gini scores
  if (!is.null(cells_groups)) {
	setDT(cells_groups)
	sum_peak_list <- tapply(cells_groups[[1]], cells_groups[[2]], function(x) {
			xj <- match(intersect(x,colnames(peaks_mat)),colnames(peaks_mat))
			Matrix::rowSums(peaks_mat[,xj,drop=FALSE])
	}, simplify = FALSE)
	sum_peak <- do.call(cbind,sum_peak_list)
	gini <- apply(sum_peak,1,function(x) ineq(x,type="Gini"))
	gini_scaled <- (gini-mean(gini,na.rm=TRUE))/sd(gini,na.rm=TRUE)
	DT[,peak_gini:=gini[peak]]
	DT[,peak_gini_scaled:=gini_scaled[peak]]
	DT[is.na(peak_gini_scaled),peak_gini_scaled:=0]
	DT[,w_var:=exp(1)^(peak_gini_scaled)]
	# final weights
	DT[,weight:=w_dist*w_var]
  } else {
	# final weights
	DT[,weight:=w_dist]
  }

  # annotations by genomic location
  DT[,rel_dist_to_tss:=dist_to_tss]
  DT[gene_chr==peak_chr & (gene_strand=="+" & peak_midpoint<gene_start) | (gene_strand=="-" & peak_midpoint>gene_start), rel_dist_to_tss:=rel_dist_to_tss*-1]
  DT[is_intergenic==TRUE,peak_location_annotation:="intergenic"]
  DT[is_promoter==TRUE,peak_location_annotation:="promoter"]
  DT[is_gene_body==TRUE,peak_location_annotation:="gene_body"]
  DT[is.na(peak_location_annotation), peak_location_annotation:="other"]
  DT[,peak_location_annotation:=factor(peak_location_annotation,levels=c("promoter","gene_body","intergenic","other"))]
  
  # weighted sum
  setorder(DT,gene,start)
  genes <- unique(DT$gene)
  gene_score_list <- lapply(1:length(genes), function(i) {
    if (i%%100==0) message(i,"/",length(genes))
    g <- genes[i]
    dt <- DT[gene==g]
    if (!all(dt$peak %in% rownames(peaks_mat))) {
      warning("Some peaks are not found in peak matrix!")
      dt <- dt[peak %in% rownames(peaks_mat)]
    }
    x <- match(dt$peak, rownames(peaks_mat))
    Matrix::colSums(t(t(peaks_mat[x, , drop=FALSE]) * dt$weight))
  })
  names(gene_score_list) <- genes
  gene_scores_mat <- do.call(rbind,gene_score_list)
  
  return(list(
    genes_peaks_table = DT,
    genes_scores_matrix = gene_scores_mat
  ))
}

#' Find marker peaks for each cell type based on gene expression and save BED files
#' 
#' @param fp footprint table (columns are metacells, cell types, etc; rows are genes)
#' @param gene_peak_table data.frame matching genes (first column) to peaks (second column). Can be produced with `mta_match_peaks_to_genes`.
#' @param peaks_object either a path to a BED file with ATAC peaks, or a preloaded `GRanges` object.
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param fixed_peak_size fixed width to apply to the peaks in `peaks_object` (default is 250 bp). If NULL, this is omitted and ranges are used as given.
#' @param threshold threshold to apply to `fp` table, to find genes in the foreground (default is 1.5)
#' @param threshold_bg threshold to apply to `fp` table, to find genes in the background (default is NULL, which means that all peaks not in the foreground are used as background)
#' @param n_top_genes max number of genes to keep for each cell type (default is 500)
#' @param cell_types vector of cell types (columns in `fp`) to find markers from (if NULL, all columns are used)
#' @param save_bed default TRUE; whether to save bed files with foreground and background peaks
#' @param bed_prefix string: prefix of bed files (a suffix specifying cell type names and fg/bg is always added)
#' 
#' @return a list where each entry is a cell type, containing a list of two `GRanges` objects in each case: the foreground and background peaks.
#' 
mta_marker_peaks_per_ct = function(
	fp,
	gene_peak_table,
	peaks_object,
	index_object,
	fixed_peak_size = 250,
	threshold = 1.5,
	threshold_bg = NULL,
	n_top_genes = 500,
	cell_types = NULL,
	save_bed = TRUE,
	bed_prefix = "marker_peaks"
) {
	
	require("rtracklayer")
	require("GenomicRanges")
	
	# list of cell types
	if (is.null(cell_types)) {
		cell_types = colnames(fp)
	}
	
	# table matching genes to peaks
	gene_peak_table = gene_peak_table[,c(1,2)]
	colnames(gene_peak_table) = c("gene", "peak")
	
	# read genome index
	if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`)", index_object))
	}
	
	# load peaks bed file
	if ("character" %in% class(peaks_object)) {
		pka_r = mta_bed_to_granges(peaks_object)
		names(pka_r) = pka_r$name
	} else if ("GRanges" %in% class(peaks_object)) {
		pka_r = peaks_object
		names(pka_r) = pka_r$name
	} else {
		stop(sprintf("`peaks_object` %s has to be either a GRanges object, or a path to a compatible BED file", peaks_object))
	}
	# resize bed regions
	if (!is.null(fixed_peak_size)) {
		pka_r = GenomicRanges::resize(pka_r, width = fixed_peak_size, fix = "center")
		pka_r = GenomicRanges::trim(pka_r)
	}
	
	# loop through cell types, finding marker genes and matching them to peaks
	pka_markers_list = list()
	for (cti in cell_types) {
		
		# get gene markers
		fg_markers = data.frame(gene = rownames(fp) [ fp[,cti] >= threshold ], fp = fp[,cti] [ fp[,cti] >= threshold ])
		fg_markers = fg_markers[ order(fg_markers$fp, decreasing = TRUE), ]
		n_initial_fg_markers = nrow(fg_markers)
		
		if (n_initial_fg_markers > 0) {
		
			# if we have specified a bg threshold, use this to define the background
			# otherwise, background will be all peaks not in foreground (which will include peaks from non-expressed genes)
			if (!is.null(threshold_bg)) {
				bg_markers = data.frame(gene = rownames(fp) [ fp[,cti] < threshold_bg ], fp = fp[,cti] [ fp[,cti] < threshold_bg ])
			}
			
			# top N markers only?
			if (!is.null(n_top_genes)) {
				fg_markers = head(fg_markers, n_top_genes)
			}
			
			# identify fg
			pka_markers_fg = merge(fg_markers, gene_peak_table, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
			pka_markers_bed_fg = pka_r [ pka_markers_fg$peak ]

			# identify bg
			if (!is.null(threshold_bg)) {
				pka_markers_bg = merge(bg_markers, gene_peak_table, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
				pka_markers_bed_bg = pka_r [ pka_markers_bg$peak ]
			} else {
				pka_markers_bed_bg = pka_r [ ! names(pka_r) %in% pka_markers_fg$peak ]
			}
			
			# keep gene and footprint info for fg
			mcols(pka_markers_bed_fg)$gene = pka_markers_fg$gene
			mcols(pka_markers_bed_fg)$fp   = pka_markers_fg$fp
			
			if (save_bed) {
				
				# save peaks as bed (foreground and background)
				# bed files contain equal-sized regions centered around peak
				dir.create(dirname(bed_prefix), showWarnings = FALSE)
				mta_granges_to_bed(gr = pka_markers_bed_fg, bed_fn = sprintf("%s.%s.fg.bed", bed_prefix, cti), mcols = c("gene","fp"))
				mta_granges_to_bed(gr = pka_markers_bed_bg, bed_fn = sprintf("%s.%s.bg.bed", bed_prefix, cti))
				
			}
			
			# log
			message(sprintf("marker peaks | %s | %i marker genes have peaks (n=%i), out of %i genes in this cell type (pp=%.1f)", cti, length(unique(pka_markers_bed_fg$gene)), length(pka_markers_bed_fg), n_initial_fg_markers, 100 * length(unique(pka_markers_bed_fg$gene)) / n_initial_fg_markers  ))
			
			# save granges
			pka_markers_list[[cti]] = list()
			pka_markers_list[[cti]]["fg"] = pka_markers_bed_fg
			pka_markers_list[[cti]]["bg"] = pka_markers_bed_bg
		
		} else { 
			
			message(sprintf("marker peaks | %s | 0 marker genes in this cell type, skip", cti, n_initial_fg_markers))
			
		}
		
	}
	
	# return list of dataframes with foreground peaks
	return(pka_markers_list)
	
}


#' Find marker peaks for each cell type based on gene expression at the metacell level and save BED files
#' 
#' @param fp footprint table (columns are metacells, cells, etc; rows are genes)
#' @param cell_type_vector named vector of cell types for each cell, metacell, etc in the `fp` table (same length as columns in `fp`, names are colnames in `fp`)
#' @param gene_peak_table data.frame matching genes (first column) to peaks (second column). Can be produced with `mta_match_peaks_to_genes`.
#' @param peaks_object either a path to a BED file with ATAC peaks, or a preloaded `GRanges` object.
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param fraction_of_mcs a gene is considered part of the foreground for a given cell type if it's part of the foreground of at least this fraction of metacells
#' @param bg_method a gene is considered part of the background for a given cell type if it's never found in the foreground of any of its metacells ("never_in_fg" mode), or if it's not in the aggregated set of fg peaks ("not_in_ct_fg")
#' @param fixed_peak_size fixed width to apply to the peaks in `peaks_object` (default is 250 bp). If NULL, this is omitted and ranges are used as given.
#' @param threshold threshold to apply to `fp` table, to find genes in the foreground (default is 1.5)
#' @param n_top_genes max number of genes to keep for each cell type (default is 500)
#' @param cell_types vector of cell types (values in `cell_type_vector`) to find markers from (if NULL, all columns are used)
#' @param save_bed default TRUE; whether to save bed files with foreground and background peaks
#' @param bed_prefix string: prefix of bed files (a suffix specifying cell type names and fg/bg is always added)
#' 
#' @return a list where each entry is a cell type, containing a list of two `GRanges` objects in each case: the foreground and background peaks.
#' 
mta_marker_peaks_per_ct_using_mcs = function(
	fp,
	cell_type_vector,
	gene_peak_table,
	peaks_object,
	index_object,
	fraction_of_mcs = 0.5,
	bg_method = "never_in_fg",
	fixed_peak_size = 250,
	threshold = 1.5,
	n_top_genes = 500,
	save_bed = TRUE,
	cell_types = NULL,
	bed_prefix = "marker_peaks"
) {
	
	require("rtracklayer")
	require("GenomicRanges")
	
	# list of cell types
	if (is.null(cell_types)) {
		cell_types = unique(cell_type_vector)
	}
	
	# table matching genes to peaks
	gene_peak_table = gene_peak_table[,c(1,2)]
	colnames(gene_peak_table) = c("gene", "peak")
	
	# read genome index
	if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`)", index_object))
	}
	
	# load peaks bed file
	if ("character" %in% class(peaks_object)) {
		pka_r = mta_bed_to_granges(peaks_object)
		names(pka_r) = pka_r$name
	} else if ("GRanges" %in% class(peaks_object)) {
		pka_r = peaks_object
		names(pka_r) = pka_r$name
	} else {
		stop(sprintf("`peaks_object` %s has to be either a GRanges object, or a path to a compatible BED file", peaks_object))
	}
	# resize bed regions
	if (!is.null(fixed_peak_size)) {
		pka_r = GenomicRanges::resize(pka_r, width = fixed_peak_size, fix = "center")
		pka_r = GenomicRanges::trim(pka_r)
	}
	
	# loop through cell types, finding marker genes and matching them to peaks
	pka_markers_list = list()
	for (cti in cell_types) {
		
		mcs_in_cti = names(cell_type_vector) [ cell_type_vector == cti ]
		
		fg_markers_v = c()
		for (mci in mcs_in_cti) {
			
			# get gene markers
			fg_markers_i = data.frame(gene = rownames(fp) [ fp[,mci] >= threshold ], fp = fp[,mci] [ fp[,mci] >= threshold ])
			fg_markers_i = fg_markers_i[ order(fg_markers_i$fp, decreasing = TRUE), ]
			
			# top N markers only?
			if (!is.null(n_top_genes)) {
				fg_markers_i = head(fg_markers_i, n_top_genes)
			}
			
			# keep
			fg_markers_v = c(fg_markers_v, fg_markers_i$gene)
			
		}
		
		# which genes appear in most metacells?
		fg_markers_t = table(fg_markers_v)
		fg_markers_g = names(which(fg_markers_t >= floor(length(mcs_in_cti) * fraction_of_mcs) ))
		n_initial_fg_markers = length(fg_markers_t)
		
		# identify fg
		fg_markers = data.frame(gene = fg_markers_g)
		pka_markers_fg = merge(fg_markers, gene_peak_table, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
		pka_markers_bed_fg = pka_r [ pka_markers_fg$peak ]

		# identify bg
		if (bg_method == "never_in_fg") {
			# which genes are never part of the foreground? this is the background
			bg_markers_g = rownames(fp) [ ! rownames(fp) %in% fg_markers_v ]
		} else if (bg_method == "not_in_ct_fg") {
			# which genes are NOT part of the foreground? this is the background
			bg_markers_g = rownames(fp) [ ! rownames(fp) %in% fg_markers_g ]
		}
		bg_markers = data.frame(gene = bg_markers_g)
		pka_markers_bg = merge(bg_markers, gene_peak_table, by.x = "gene", by.y = "gene", all.x = FALSE, all.y = FALSE)
		pka_markers_bed_bg = pka_r [ pka_markers_bg$peak ]
		
		# keep gene info for fg
		mcols(pka_markers_bed_fg)$gene = pka_markers_fg$gene
		
		if (save_bed) {
			
			# save peaks as bed (foreground and background)
			# bed files contain equal-sized regions centered around peak
			dir.create(dirname(bed_prefix), showWarnings = FALSE)
			mta_granges_to_bed(gr = pka_markers_bed_fg, bed_fn = sprintf("%s.%s.fg.bed", bed_prefix, cti), mcols = c("gene"))
			mta_granges_to_bed(gr = pka_markers_bed_bg, bed_fn = sprintf("%s.%s.bg.bed", bed_prefix, cti))
			
		}
		
		# log
		message(sprintf("marker peaks | %s | %i marker genes have peaks (n=%i), out of %i genes in this cell type (pp=%.1f)", cti, length(unique(pka_markers_bed_fg$gene)), length(pka_markers_bed_fg), n_initial_fg_markers, 100 * length(unique(pka_markers_bed_fg$gene)) / n_initial_fg_markers  ))
		
		# save granges
		pka_markers_list[[cti]] = list()
		pka_markers_list[[cti]]["fg"] = pka_markers_bed_fg
		pka_markers_list[[cti]]["bg"] = pka_markers_bed_bg
		
	}
	
	# return list of dataframes with foreground peaks
	return(pka_markers_list)
	
}


#############################
## Motif merging functions ##
#############################

#' Plot logos of merged motifs from `universalmotif`
#' 
#' @param motifs_merged motif object after merging (e.g. output of `universalmotif::merge_similar`)
#' @param motifs_original motif object before merging (e.g. input of `universalmotif::merge_similar`)
#' 
#' @return outputs to plotting device
mta_plot_merged_motifs = function(motifs_merged, motifs_original) {
	
	require("stringr")
	require("universalmotif")
	
	# get names of motifs
	merged_names = unlist(lapply(1:length(motifs_merged), function(i) motifs_merged[[i]]@name ))
	merged_names_list = stringr::str_split(merged_names, "/")
	
	# keep merged
	ix_merged = which(lengths(merged_names_list) > 1)
	merged_names_list = merged_names_list [ ix_merged ]
	merged_names = merged_names [ ix_merged ]
	motifs_merged = motifs_merged [ ix_merged ]
	
	# original names
	original_names = unlist(lapply(1:length(motifs_original), function(i) motifs_original[[i]]@name ))
	
	# loop to find original indexes of the motifs that were merged to merge
	for (m in 1:length(merged_names_list)) {
		# which motifs to plot (merged 1st, original later)
		ms = merged_names_list[[m]]
		ixs_o = which(original_names %in% ms)
		# plot
		if (length(ixs_o) > 0 ) {
			print(universalmotif::view_motifs(c(motifs_merged[m], motifs_original[ixs_o])))
		}
	}
	
}


#' Merge motifs by genome-wide score similarity
#'
#' @param motifs list of motif matrices (`universalmotif` format). Will be converted to `PWM` format if not already in this format.
#' @param genome_object either a path to a genome file (`.fasta`), or a preloaded `DNAStringSet` object (from `Biostrings::readDNAStringSet`)
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param correlation_matrix if you have a pre-computed correlation matrix, you can specify it here and avoid costly computations. Only relevant if `clus_method="hclust`. Rownames and colnames are motif names in `motifs`.
#' @param score_matrix if you have a pre-computed score matrix, you can specify it here and avoid costly computations. Rownames are motifs, and columns are their scores over random genomic bins.
#' @param bin_width width of the genomic bins used to calculate genome-wide max alignment scores for each motif (default is 250bp)
#' @param given_gr GRanges object defining specific regions to score. Default is NULL, i.e. the whole genome is used. 
#' @param subsample_fraction fraction of the genome to include in the calculations of genome-wide alignment score distributions (default = 0.10, 10%). If `given_gr` are provided, this is applied to the size of the ranges instead.
#' @param method correlation method ("pearson" or "spearman", default is "pearson").
#' @param clus_method clustering method, either "hclust" (see `universalmotif::merge_similar`) or "dbscan"
#' @param threshold clustering threshold applied to `hclust`, applied if `clus_method="hclust"`
#' @param eps_dbscan epsilon value for `dbscan`, applied if `clus_method="dbscan"`
#' @param nthreads num of threads to use for motif scoring, for parallelisation over motifs.
#'
#' @return a list with four elements: a data.frame with motifs and their cluster membership, a matrix of motif-motif score correlation, a matrix of the scores found per bin for each motif, and the list of bins (`GRanges` object).
#'
mta_merge_motifs_by_gw_score_correlation = function(
	motifs, 
	genome_object,
	index_object,
	correlation_matrix = NULL,
	score_matrix = NULL,
	bin_width = 250,
	given_gr = NULL,
	subsample_fraction = 0.10,
	threshold = 0.95,
	method = "pearson",
	clus_method = "hclust",
	eps_dbscan = 1,
	nthreads = 2
) {
	
	require("universalmotif")
	require("GenomicRanges")
	require("Biostrings")
	require("PWMEnrich")
	require("doParallel")
	require("plyr")
	require("WGCNA")
	
	# register cores
	registerDoParallel(cores = nthreads)
	
	#### Load data ####
	# read genome index
	if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`)", index_object))
	}
	
	# read genome fasta
	if ("character" %in% class(genome_object)) {
		gen = Biostrings::readDNAStringSet(genome_object, format = "fasta")
	} else if ("DNAStringSet" %in% class(genome_object)) {
		gen = genome_object
	} else {
		stop(sprintf("`genome_object` %s has to be either a DNAStringSet or a path to a fasta file", DNAStringSet))
	}
	
	
	#### Genome-wide scores ####
	# get genome sequences from subsampled bins
	if (is.null(given_gr)) {
		
		# if no ranges are provided, use whole genome
		gix_seqlengths = gix$length
		names(gix_seqlengths) = gix$chr
		sub_bin = GenomicRanges::tileGenome(gix_seqlengths, tilewidth = bin_width * 10)
		sub_bin = sample(sub_bin, size = floor(length(sub_bin) * subsample_fraction), replace = FALSE)
		sub_gen = BSgenome::getSeq(gen, sub_bin)
		sub_gen = unlist(sub_gen)
		names(sub_gen) = 1:length(sub_gen)
		message(sprintf("merge by motif-motif correlation | %i subsampled bins (f=%.2fpc of genome, %.1fMbp)", length(sub_gen), subsample_fraction * 100, sum(width(sub_gen)) / 1e6 ))
		
	} else {
		
		# if ranges are provided, use them instead
		sub_bin = sample(given_gr, size = floor(length(given_gr) * subsample_fraction), replace = FALSE)
		sub_gen = BSgenome::getSeq(gen, sub_bin)
		names(sub_gen) = 1:length(sub_gen)
		message(sprintf("merge by motif-motif correlation | %i given ranges (f=%.2fpc of genome, %.1fbp)", length(sub_gen), subsample_fraction * 100, sum(width(sub_gen)) / 1e6 ))
		
	}
	
	# find distribution of scores over equally sized windows
	# first	create genome bins 
	sub_seqlengths = width(sub_gen)
	names(sub_seqlengths) = 1:length(sub_gen)
	sub_bin_r = GenomicRanges::tileGenome(sub_seqlengths, tilewidth = bin_width)
	sub_bin_r = sub_bin_r [ width(sub_bin_r) == bin_width ]
	
	# motif list
	names_motifs = sapply(1:length(motifs), function(i) motifs@listData[[i]]@name)
	
	# for each motif, find score distribution in genome-wide peaks
	# loop over motifs
	registerDoParallel(cores = nthreads)
	fun_loop_scoring = function(i) {
		
		# log
		if (i == 1 | i %% 100 == 0 | i == length(motifs)) {
			message(sprintf("motif scoring | motif %s/%s", i, length(motifs)))
		}
		
		# motif score distribution
		sub_gen_r = monaLisa::findMotifHits(
			query = motifs[[i]],
			subject = sub_gen,
			min.score = 0, 
			BPPARAM = BiocParallel::MulticoreParam(nthreads)
		)
		
		# overlap between bins and motif alignment sites
		sub_bin_ovs = GenomicRanges::findOverlaps(sub_bin_r, sub_gen_r)
		
		# find max score per bin
		sub_bin_scores_d = data.frame(from_bin = sub_bin_ovs@from, score = mcols(sub_gen_r [ sub_bin_ovs@to ])$score, motif = mcols(sub_gen_r [ sub_bin_ovs@to ])$pwmname)
		sub_bin_scores_d$from_bin = factor(sub_bin_scores_d$from_bin, levels = 1:length(sub_bin_r))

		# store
		if (nrow(sub_bin_scores_d) > 0) {
			sub_bin_scores_v = aggregate(score ~ from_bin, data = sub_bin_scores_d, FUN = max, drop = FALSE)[,2]
			names(sub_bin_scores_v) = 1:length(sub_bin_r)
			sub_bin_scores_v [ is.na(sub_bin_scores_v) ] = 0
		} else {
			sub_bin_scores_v = rep(0, length(sub_bin_r))
			names(sub_bin_scores_v) = 1:length(sub_bin_r)
		}
		
		return(list(scores = sub_bin_scores_v, motif = names_motifs[i]))
	}
	
	# do we have a precomputed motif score matrix?
	if (is.null(score_matrix)) {
		
		# execute parallelised loop to get vectors of scores in bins
		bin_scores_l = plyr::alply(.data = 1:length(motifs), .margins = 1,  .fun = fun_loop_scoring, .parallel = TRUE, .inform = TRUE)
		# as matrix
		bin_scores_m = t(sapply(1:length(motifs), function(i) bin_scores_l[[i]]$scores))
		rownames(bin_scores_m) = sapply(1:length(motifs), function(i) bin_scores_l[[i]]$motif)

	} else {
		
		message(sprintf("merge by motif-motif correlation | load precomputed motif-bin score matrix..."))
		bin_scores_m = score_matrix
		sub_bin_r = NULL
	}
	
	# clusters dataframe
	if (!is.null(clus_method)) {
		
		if (clus_method == "hclust") {
			
			message(sprintf("merge by motif-motif correlation | %s clustering, thr = %.2f", clus_method, threshold))
			# do we have a precomputed correlation matrix?
			if (is.null(correlation_matrix)) {			
				bin_scores_cor = WGCNA::cor(t(bin_scores_m), method = method)
			} else {
				message(sprintf("merge by motif-motif correlation | load precomputed motif-motif correlation matrix..."))
				bin_scores_cor = correlation_matrix
			}
			bin_scores_hcl = hclust(as.dist(1 - bin_scores_cor), method = "ward.D2")
			bin_scores_clu = cutree(bin_scores_hcl, h = 1 - threshold)
			mot_merge_d = data.frame(
				motif = rownames(bin_scores_m),
				cluster = bin_scores_clu [ rownames(bin_scores_m) ]
			)
			
		} else if (clus_method == "dbscan") {
		
			message(sprintf("merge by motif-motif correlation | %s clustering, eps = %.2f", clus_method, eps_dbscan))
			mot_clu = dbscan::dbscan(bin_scores_m, eps = eps_dbscan, minPts = 2)
			mot_merge_d = data.frame(
				motif = rownames(bin_scores_m),
				cluster = mot_clu$cluster
			)
			max_num_clu = max(mot_clu$cluster) + 1
			fin_num_clu = max_num_clu + length(mot_merge_d$cluster [ mot_clu$cluster == 0 ]) - 1
			mot_merge_d$cluster [ mot_clu$cluster == 0 ] = max_num_clu:fin_num_clu
			bin_scores_cor = NULL
		
		} else {
		
			stop(sprintf("Clustering method %s is invalid", clus_method))
		
		}
		
		# reorder by cluster size
		mot_merge_d = transform(mot_merge_d, freq = ave(seq(nrow(mot_merge_d)), cluster, FUN = length))
		mot_merge_d = mot_merge_d [ order(-mot_merge_d$freq, mot_merge_d$cluster) , ]
		mot_merge_d$freq = NULL
		message(sprintf("merge motifs by similarity | motif clusters n = %i (%i motifs)", length(unique(mot_merge_d$cluster)), nrow(mot_merge_d)))
	
	} else if (is.null(clus_method)) { 
		
		message(sprintf("merge by motif-motif correlation | omit..."))
		bin_scores_cor = NULL
		mot_merge_d = NULL
		
	}
		
	# return
	return(list(clusters = mot_merge_d, score_correlation = bin_scores_cor, scored_matrix = bin_scores_m, scored_ranges = sub_bin_r))
	
}


#' Merge motifs by position matrix similarity
#'
#' @param motifs list of motif matrices (`universalmotif` format).
#' @param matrix_class class of motifs in list (`universalmotif` classes, default is `ICM`).
#' @param dist_method distance calculation method, from `universalmotif::merge_similar`.
#' @param clus_method clustering method, either "hclust" (see `universalmotif::merge_similar`) or "dbscan"
#' @param threshold,min_overlap parameters for `universalmotif::merge_similar`.
#' @param do_heatmap whether to return a heatmap of motif-motif similarity
#' @param col_vec vector of colors for heatmap
#'
#' @return a list with three elements: a data.frame with motifs and their cluster membership, the merged motifs object from `universalmotif`, and the matrix of pairwise motif similarities.
#'
mta_merge_motifs_by_similarity = function(
	motifs,
	matrix_class = "ICM",
	threshold = 0.95,
	dist_method = "WPCC",
	clus_method = "hclust",
	min.overlap = 6,
	min.mean.ic = 0.25,
	min.position.ic = 0,
	normalise.scores = TRUE,
	eps_dbscan = 1,
	do_heatmap = FALSE,
	col_vec = c("white","#d6e72e","#6fb600","#003f4d")
) {
	
	require("universalmotif")
	require("dbscan")
	
	# convert motifs to adequate format
	message(sprintf("merge motifs by similarity | merge %s motifs, n = %i", matrix_class, length(motifs)))
	mot_mod  = universalmotif::convert_type(motifs, matrix_class, pseudocount = 1)
	
	# merge
	message(sprintf("merge motifs by similarity | %s similarity matrix", dist_method))
	mot_mod_m = mta_merge_similar_umot_mod(
		mot_mod,
		threshold = threshold, 
		method = dist_method, 
		min.overlap = min.overlap, 
		min.mean.ic = min.mean.ic, 
		min.position.ic = min.position.ic, 
		normalise.scores = normalise.scores
	)
	
	# clusters dataframe
	if (clus_method == "hclust") {
		
		message(sprintf("merge motifs by similarity | %s clustering, thr = %.2f", clus_method, threshold))
		# get vector of motifs
		mot_mod_m_names = plyr::laply(mot_mod_m$motifs, function(i) i@name)
		names_original_motifs = strsplit(mot_mod_m_names, "/")
		# get vector of merged motifs (merging represents motif clusters)
		names_merged_motifs = unlist(plyr::llply(1:length(names_original_motifs), function(i) { rep(mot_mod_m_names[i], length(names_original_motifs[[i]])) }))
		mot_merge_d = data.frame(
			motif = unlist(names_original_motifs),
			cluster = as.numeric(as.factor(names_merged_motifs))
		)
		
	} else if (clus_method == "dbscan") {
		
		message(sprintf("merge motifs by similarity | %s clustering, eps = %.2f", clus_method, eps_dbscan))
		mot_clu = dbscan::dbscan(mot_mod_m$matrix, eps = eps_dbscan, minPts = 2)
		mot_merge_d = data.frame(
			motif = rownames(mot_mod_m$matrix),
			cluster = mot_clu$cluster
		)
		max_num_clu = max(mot_clu$cluster) + 1
		fin_num_clu = max_num_clu + length(mot_merge_d$cluster [ mot_clu$cluster == 0 ]) - 1
		mot_merge_d$cluster [ mot_clu$cluster == 0 ] = max_num_clu:fin_num_clu
		
	} else {
		stop(sprintf("Clustering method %s is invalid", clus_method))
	}
	
	# reorder by cluster size
	mot_merge_d = transform(mot_merge_d, freq = ave(seq(nrow(mot_merge_d)), cluster, FUN = length))
	mot_merge_d = mot_merge_d [ order(-mot_merge_d$freq, mot_merge_d$cluster) , ]
	mot_merge_d$freq = NULL
	# log
	message(sprintf("merge motifs by similarity | motif clusters n = %i (%i motifs)", length(unique(mot_merge_d$cluster)), nrow(mot_merge_d)))
	
	# reorder similarity matrix
	mot_mod_m$matrix = mot_mod_m$matrix [ mot_merge_d[,1] , mot_merge_d[,1] ]
	
	# plot similarity heatmap?
	if (do_heatmap) {
		
		require("ComplexHeatmap")
		ht_opt$message = FALSE
		
		# plot motif-motif similarities
		col_fun = circlize::colorRamp2(seq(0.5, 1, length.out = length(col_vec)), col_vec)
		# vector of clusters (row/col annotation)
		cluster_vector = as.character(mot_merge_d$cluster)
		# categorical colors for clusters (get recycled)
		catcol_lis = c("magenta4","firebrick1","orange","khaki1","springgreen2","darkgreen","deepskyblue","cadetblue1","mediumblue","darkviolet","violet")
		catcol_vec = rep(catcol_lis, length.out = length(unique(cluster_vector)))
		names(catcol_vec) = unique(cluster_vector)
		cluster_colors = catcol_vec [ cluster_vector ]
		names(cluster_colors) = cluster_colors
		# quantitative colorscale for heatmap
		catcol_fun = circlize::colorRamp2(1:length(catcol_vec), catcol_vec)
		names(cluster_colors) = cluster_vector
		# plot
		hmmat = mot_mod_m$matrix
		hm = ComplexHeatmap::Heatmap(
			hmmat, 
			name = dist_method,
			use_raster = TRUE,
			cluster_rows = FALSE,
			cluster_columns = FALSE,
			row_title = sprintf("%i motifs in %i clusters", nrow(hmmat), length(unique(cluster_vector))),
			column_title = sprintf("%i motifs in %i clusters", nrow(hmmat), length(unique(cluster_vector))),
			show_row_names = FALSE,
			show_column_names = FALSE,
			top_annotation = ComplexHeatmap::HeatmapAnnotation(clusters = cluster_vector, col = list(clusters = cluster_colors), which = "column", show_legend = FALSE),
			left_annotation = ComplexHeatmap::HeatmapAnnotation(clusters = cluster_vector, col = list(clusters = cluster_colors), which = "row", show_legend = FALSE),
			col = col_fun)
		
		
	} else {
		
		hm = NULL
		
	}
	
	# return
	return(list(clusters = mot_merge_d, merged_motifs = mot_mod_m$motifs, matrix = mot_mod_m$matrix, matrix_hm = hm))
	
}


# Modified version of the `universalmotif::merge_similar` function that will return the similarity matrix
mta_merge_similar_umot_mod = function(motifs,
									  threshold = 0.95, threshold.type = "score.abs", method = "PCC",
									  use.type = "PPM", min.overlap = 6, min.mean.ic = min.mean.ic, min.position.ic = min.position.ic,
									  tryRC = TRUE, relative_entropy = FALSE, normalise.scores = FALSE,
									  score.strat.compare = "a.mean",
									  score.strat.merge = "sum", nthreads = 1) {
	
	# internal functions
	collapse_cpp <- function(x) {
		.Call("_universalmotif_collapse_cpp", PACKAGE = "universalmotif", x)
	}
	internal_convert <- function(motifs, class = NULL) {
		if (is.null(class)) {
			CLASS <- class(motifs)
			CLASS_PKG <- attributes(CLASS)$package
			CLASS_IN <- collapse_cpp(c(CLASS_PKG, "-", CLASS))
			CLASS_IN
		} else {
			if (length(class) == 1 && class[1] != "universalmotif-universalmotif") {
				tryCatch(motifs <- universalmotif::convert_motifs(motifs, class),
						 error = function(e) message("motifs converted to class 'universalmotif'"))
				
			} else if (length(class) > 1)
				message("motifs converted to class 'universalmotif'")
			motifs
		}
	}
	sort_by_ic <- function(x) {
		ICs <- vapply(x, function(x) x@icscore, numeric(1))
		x[order(ICs, decreasing = TRUE)]
	}
	
	# prepare input
	if (is.list(motifs)) {
		CLASS_IN <- vapply(motifs, internal_convert, character(1))
	} else { 
		CLASS_IN <- internal_convert(motifs)
	}
	motifs <- convert_motifs(motifs)
	if (!is.list(motifs)) motifs <- list(motifs)
	
	if (length(motifs) == 1) {
		return(internal_convert(motifs[[1]], unique(CLASS_IN)))
	}
	
	# For now, just merge based on score -- add p-value merging later.
	
	if (score.strat.compare == "sum")
		stop("`score.strat.compare` cannot be \"sum\".", call. = FALSE)
	if (method %in% c("ALLR", "ALLR_LL"))
		stop(wmsg("`method` cannot be \"ALLR\", \"ALLR_LL\" as the resulting similarity",
				  " matrix cannot be converted to a distance matrix."), call. = FALSE)
	
	comp.mat <- compare_motifs(
		motifs, method = method,
		use.type = use.type, min.overlap = min.overlap, min.mean.ic = min.mean.ic,
		tryRC = tryRC, relative_entropy = relative_entropy,
		normalise.scores = normalise.scores, min.position.ic = min.position.ic,
		score.strat = score.strat.compare, nthreads = nthreads)
	comp.mat.o = comp.mat
	
	if (anyNA(comp.mat)) {
		stop(wmsg("Clustering is not possible when NA values are present in the ",
				  "comparison matrix. Lower `min.mean.ic` and/or `min.position.ic` until ",
				  "no more NA values are generated, or omit the problematic motifs."),
			 call. = FALSE)
	}
	
	if (method %in% c("PCC", "SW", "WPCC")) {
		threshold <- comp.mat[1] - threshold
		comp.mat <- comp.mat[1] - comp.mat
	}
	
	comp.clust <- hclust(as.dist(comp.mat))
	comp.groups <- cutree(comp.clust, h = threshold)
	
	final.mots <- unname(tapply(motifs, comp.groups,
								function(x) merge_motifs(sort_by_ic(x),
														 method = method, use.type = use.type, min.overlap = min.overlap,
														 min.mean.ic = min.mean.ic, tryRC = tryRC,
														 relative_entropy = relative_entropy, normalise.scores = normalise.scores,
														 min.position.ic = min.position.ic, score.strat = score.strat.merge)))
	
	internal_convert(final.mots, unique(CLASS_IN))
	
	# return
	return(list(motifs = final.mots, matrix = comp.mat.o))
	
}

#' Find most representative motif in each cluster
#'
#' @param motifs list of motif matrices prior to merging (`universalmotif` format).
#' @param clusters_df data.frame with two columns: first are motif names (present in `motifs` object), second are cluster memberships
#' @param criterion one of: "max_ic", ""max_ic_min_fg", "max_nsites", "min_pval" (default: "max_ic_min_fg", i.e. the representative motif will be the one with highest IC, selected from those with a minimum number of foreground counts (`nsites`). If there are no motifs with that number of `nsites`, the highest IC is selected.).
#' @param min_nsites numeric (default 10), min number of nsites (foreground sites) required to select a motif as representative in the "max_ic_min_fg" criterion (otherwise, ignored)
#' @return a list with two elements: a data.frame with motifs with representative motifs per cluster, and the filtered list of motifs
#'
mta_merge_motifs_find_best = function(
	motifs,
	clusters_df,
	criterion = "max_nsites",
	min_nsites = 10) {
	
	require("universalmotif")
	require("plyr")
	
	mot_name = plyr::laply(motifs, function(i) i@name)
	
	# apply selection criterion
	if (criterion == "max_nsites") {
		
		# get number of observed sites per motifs
		mof_d = data.frame(
			motif = plyr::laply(motifs, function(i) i@name),
			nsites = plyr::laply(motifs, function(i) i@nsites)
		)
		# merge with merged clusters df
		mof_m = merge(mof_d, clusters_df, by.x = "motif", by.y = "motif", all.x = TRUE, all.y = FALSE)
		# keep top member of each cluster
		mof_m = mof_m [ order(mof_m$cluster, mof_m$nsites, decreasing = TRUE ), ]
		mof_f = mof_m [ !duplicated(mof_m$cluster) , ]
		mof_f = mof_f [ order(mof_f$cluster) , ]
		# list of motifs to keep
		mot_keep = mof_f$motif
		
	} else if (criterion == "max_ic") {
		
		# get number IC scores per motif
		mof_d = data.frame(
			motif = plyr::laply(motifs, function(i) i@name),
			icscore = plyr::laply(motifs, function(i) i@icscore)
		)
		# merge with merged clusters df
		mof_m = merge(mof_d, clusters_df, by.x = "motif", by.y = "motif", all.x = TRUE, all.y = FALSE)
		# keep top member of each cluster
		mof_m = mof_m [ order(mof_m$cluster, mof_m$icscore, decreasing = TRUE ), ]
		mof_f = mof_m [ !duplicated(mof_m$cluster) , ]
		mof_f = mof_f [ order(mof_f$cluster) , ]
		# list of motifs to keep
		mot_keep = mof_f$motif
		
	} else if (criterion == "max_ic_min_fg") {
		
		# get number IC scores per motif
		mof_d = data.frame(
			motif = plyr::laply(motifs, function(i) i@name),
			nsites = plyr::laply(motifs, function(i) i@nsites),
			high_nsite = as.numeric(plyr::laply(motifs, function(i) i@nsites) > min_nsites),
			icscore = plyr::laply(motifs, function(i) i@icscore)
		)
		# merge with merged clusters df
		mof_m = merge(mof_d, clusters_df, by.x = "motif", by.y = "motif", all.x = TRUE, all.y = FALSE)
		# keep top member of each cluster
		mof_m = mof_m [ order(mof_m$cluster, mof_m$high_nsite, mof_m$icscore, decreasing = TRUE ), ]
		mof_f = mof_m [ !duplicated(mof_m$cluster) , ]
		mof_f = mof_f [ order(mof_f$cluster) , ]
		# list of motifs to keep
		mot_keep = mof_f$motif
		
	} else if (criterion == "min_pval") {
		
		# get number IC scores per motif
		mof_d = data.frame(
			motif = plyr::laply(motifs, function(i) i@name),
			pval = plyr::laply(motifs, function(i) i@pval)
		)
		# merge with merged clusters df
		mof_m = merge(mof_d, clusters_df, by.x = "motif", by.y = "motif", all.x = TRUE, all.y = FALSE)
		# keep top member of each cluster
		mof_m = mof_m [ order(mof_m$cluster, mof_m$pval, decreasing = FALSE ), ]
		mof_f = mof_m [ !duplicated(mof_m$cluster) , ]
		mof_f = mof_f [ order(mof_f$cluster) , ]
		# list of motifs to keep
		mot_keep = mof_f$motif
		
	} else {
		
		stop(sprintf("selection criterion %s is invalid", criterion))
		
	}
	
	# subset original motifs
	mot_f = motifs [ mot_name %in% mot_keep ]
	mot_keep_o = mot_name [ mot_name %in% mot_keep ]
	message(sprintf("best motif selection | %i representative motifs selected by %s", length(mot_keep), criterion))
	
	# output
	return(list(keep_motifs = mot_keep_o, motifs = mot_f))
	
}


#' Motif archetypes
#' 
#' @param motifs named list of motif matrices (`universalmotif` format)
#' @param sim_mat matrix, motif-motif similarity, rows and column names should be names of `motifs`
#' @param type character, one of `c('PCM', 'PPM', 'PWM')`
#' @param clusters named character, motif cluster memberships
#' @param recluster logical, whether to recluster the motifs similarity matrix
#' @param min_cluster_similarity numeric, minimum similarity between all motifs in archetyping cluster
#' @param hclust_method character, only applicable if `recluster=TRUE`
#' @param dist_method character, only applicable if `recluster=TRUE`
#' @param cuts, numeric, where to cut the hierarchical clustering tree to generate clusters, only applicable if `recluster=TRUE`
#' @param block_filter logical, whether to block-filter motifs before generating archetypes
#'     see `?mta_filter_by_ic_block` for more information
#' @param bkg numeric, background nucleotides frequencies
#' @param pseudocount numeric, pseudocount for convering PPM to ICM
#' @param len_threshold numeric, exclude motifs shorter than this value from archetyping  
#'     (default: NULL, no length filtering is done before archetyping)
#' @param IC_threshold numeric, positions with IC lower than this value in consensus motif 
#'     will be trimmed off the end os consensus motif
#' @param occupancy_threshold numeric, if occupancy_threshold < 0, a percentage of motifs 
#'     and if occupancy_threshold > 0, an number of motifs, that any position needs to be 
#'     present in, to be retained for consensus calculation
#' 
#' 
mta_merge_archetype <- function(
	motifs,
	sim_mat,
	type = "PPM",
	clusters = NULL,
	recluster = FALSE, 
	hclust_method = "complete",
	dist_method = "euclidean",
	cuts = 1:5,
	min_cluster_similarity = 0.5,
	block_filter = FALSE,
	bkg = rep(0.25,4), 
	pseudocount = 0.0001, 
	IC_threshold = 0.5, 
	len_threshold = NULL,
	occupancy_threshold = 1,
	verbose = FALSE,
	...
) {
	
	# clustering 
	ord <- rownames(sim_mat)
	if (recluster==TRUE) {
		if (verbose==TRUE) message("Reclustering motifs")
		hc <- hclust(dist(sim_mat, method = dist_method), method = hclust_method)
		ord <- hc$labels[hc$order]
		cuts_scores <- sapply(cuts, function(h) {
			clusters <- cutree(hc, k=h)
			cl_scores <- sapply(unique(clusters),function(x) {
				ms <- names(clusters[clusters==x])
				within_cl <- median(sim_mat[ms,ms], na.rm=TRUE)
				between_cl <- median(unlist(sim_mat[!(rownames(sim_mat)%in%ms),ms]), unlist(sim_mat[ms,!(colnames(sim_mat)%in%ms)]), na.rm=TRUE)
				if (is.na(between_cl)) between_cl <- 1
				within_cl/between_cl
			})
			mean(cl_scores, na.rm=TRUE)
		})
		k <- cuts[which.max(cuts_scores)]
		clusters <- cutree(hc, k=k)
	}
	
	# helper function to split cluster
	split_subcluster <- function(sim_mat_sub, min_cluster_similarity, dist_method, hclust_method, verbose) {
		min_sim_mat_sub <- min(sim_mat_sub)
		k <- 2
		while (min_sim_mat_sub < min_cluster_similarity) {
			if (verbose==TRUE) 
				message("Splitting cluster into ", k, " because minimum similarity (",min(sim_mat_sub),") is below threshold (",min_cluster_similarity,").")
			hc <- hclust(dist(sim_mat_sub, method = dist_method), method = hclust_method)
			ord <- hc$labels[hc$order]
			spl <- cutree(hc, k=k)
			newcl <- paste0(".",spl)
			names(newcl) <- names(spl)
			min_sim_mat_sub <- min(unlist(lapply(unique(newcl), function(nc) min(sim_mat_sub[newcl==nc,newcl==nc]))))
			k <- k+1
		}
		newcl
	}
	
	
	# archetyping function
	aha <- function(x) {
		
		ppms <- motifs[x]
		
		# block filtering
		if (block_filter) {
			ids <- suppressMessages(
				mta_filter_by_ic_block(motifs = ppms, ic_thr=0.5, len_uniblock=4, len_multiblock=3)
			)
			ppms <- ppms[ids]
			if (verbose==TRUE) message(length(ppms), " motif(s) after block filtering")
		}
		
		# length filtering
		if (!is.null(len_threshold)) {
			mlens <- lapply(ppms, function(x) ncol(x@motif))
			ppms <- ppms[mlens>len_threshold]
			if (verbose==TRUE) message(length(ppms), " motif(s) after length filtering (", len_threshold, ")")
		} else {
			len_threshold=3
		}
		
		# archetyping
		if (length(ppms) == 0) {
			
			list()
			
		} else if (length(ppms) == 1) {
			
			ppms <- ppms[[1]]
			
			ppms <- mta_trim_by_ic(ppms, IC_threshold)
			
			consensus_ic <- tryCatch(
				universalmotif::convert_type(ppms, "ICM", pseudocount = pseudocount, relative_entropy = TRUE), 
				error=function(e) universalmotif::convert_type(ppms, "ICM", relative_entropy = TRUE)
			)
			
			list(
				ppms = ppms,
				ppms_aligned = ppms,
				ppm_consensus = ppms,
				ic_consensus = consensus_ic
			)
			
		} else {
			
			if (verbose==TRUE) message("Aligning ", length(ppms), " motifs")
			arr_list <- view_motifs(
				motifs = ppms,
				use.type = type, 
				method = "PCC", 
				normalise.scores = FALSE,
				return.raw = TRUE
			)
			
			# consensus
			arr <- array(
				data = unlist(arr_list), 
				dim = c(4, ncol(arr_list[[1]]),length(arr_list)), 
				dimnames = list(c("A","C","G","T"),NULL,names(arr_list))
			)
			out_mat <- rowSums(arr,dims=2)
			consensus_ppm <- apply(out_mat,2,function(x) x/sum(x))
			consensus_ppm[is.na(consensus_ppm)] <- 0
			
			# occupancy 
			occ_list <- lapply(arr_list, function(x) colSums(x)>0)
			occ_arr <- array(
				data = unlist(occ_list), 
				dim = c(1, length(occ_list[[1]]),length(occ_list)), 
				dimnames = list(NULL,NULL,names(arr_list))
			)
			occup_mat <- rowSums(occ_arr, dims=2)
			if (occupancy_threshold>1 | occupancy_threshold==1) {
				filt_occup <- occup_mat > pmin(occupancy_threshold,length(occ_list)-1)
			} else if (occupancy_threshold<1) {
				filt_occup <- occup_mat > occupancy_threshold*length(occ_list)
			}
			
			# IC
			consensus_ppm_motif <- universalmotif::create_motif(input=consensus_ppm, alphabet="DNA", type="PPM", bkg=bkg)
			filt_ic <- mta_trim_by_ic(consensus_ppm_motif, IC_threshold=IC_threshold, return_positions=TRUE)
			
			# filter columns of consensus matrix by IC and occupancy
			keep_cols <- filt_ic & filt_occup
			
			if (sum(keep_cols) > len_threshold) {
				
				col_start <- which(keep_cols==TRUE)[1]
				col_end <- rev(which(keep_cols==TRUE))[1]
				consensus_ppm <- suppressMessages(consensus_ppm[,col_start:col_end])
				consensus_ppm <- suppressMessages(universalmotif::create_motif(
					input=consensus_ppm, alphabet="DNA", type=type, bkg=bkg
				))
				consensus_ppm@name <- paste(sapply(ppms, function(m) m@name), collapse="__")
				
				# convert final consensus to IC matrix
				consensus_ic <- tryCatch(
					universalmotif::convert_type(consensus_ppm, "ICM", pseudocount = pseudocount, relative_entropy = TRUE), 
					error=function(e) universalmotif::convert_type(consensus_ppm, "ICM", relative_entropy = TRUE)
				)
				
				list(
					ppms = ppms,
					ppms_aligned = arr_list,
					ppm_consensus = consensus_ppm,
					ic_consensus = consensus_ic
				)
				
			} else {
				
				if (verbose) message("No motifs above threshold!")
				NULL
				
			}
			
		}
		
	}

	# catch unnamed motif list
	if (is.null(names(motifs))) {
	  names(motifs) <- sapply(motifs, function(x) x@name)
	}
	
	# archetype motif per cluster
	arch_list <- lapply(unique(clusters), function(i) {
		
		# get motifs similaritys
		x <- names(clusters)[clusters==i]
		if (verbose==TRUE) message("\n", length(x), " motif(s) in cluster ", i)
		sim_mat_sub <- sim_mat[x,x]
		
		# split if necessary
		if (min(sim_mat_sub) < min_cluster_similarity) {
			if (verbose==TRUE) 
				message("Splitting cluster ",i," because minimum similarity (",min(sim_mat_sub),") is below threshold (",min_cluster_similarity,").")
			splits <- split_subcluster(sim_mat_sub, min_cluster_similarity, dist_method, hclust_method, verbose)
			subclusters <- structure(paste0(i,splits), names=names(splits))
			lapply(unique(subclusters), function(ii) {
				xx <- names(subclusters)[subclusters==ii]
				if (verbose==TRUE) message("\n", length(xx), " motif(s) in subcluster ", ii)
				aha(xx)
			})
		} else {
			list(aha(x))
		}
		
	})
	
	# prepare output
	arch_list = unlist(arch_list, recursive = FALSE)
	message(sprintf(
		"list of %i archetypes including %i motifs (%i archetype(s) fail filters)", 
		length(arch_list) - sum(sapply(arch_list, is.null)), 
		sum(lengths(arch_list)), 
		sum(sapply(arch_list, is.null))
	))
	
	# ignore empty archetypes
	arch_list = arch_list [ ! sapply(arch_list, is.null) ]
	
	# return
	return(arch_list)
	
}

#' Generate dictionary for arcetyping results
#' 
#' Maps archetypes to original motifs and generates short informative names for archetypes. 
#' If annotation files are provided, these annotations are included in the output mappings.
#' 
#' @param arch list, output of `mta_merge_archetype()`
#' @param TF_annotation_file character, path to tsv file with TF annotations 
#' Gene ids should be in the first column, followed by additional annotations in other columns.
#' @param TF_motifs_file character, path to tsv file with TF genes to motif annotations. 
#' Genes should be in the first and motifs in the second column
#' @param TF_family_annotation_file character, path to tsv file with two columns, the first containing TF family and the second TF class 
#' (multiple rows per TF family), or a comma-separated list of classes belonging to this family (one row per TF family).
#' This file is used to grep-match TFs to families.
#' @param CisBP_family_annotation_file character, path to tsv file with TF information from CisBP, 
#' it should contain at least Motif_ID, TF_Name, and Family_Name columns. This is used to annotate motifs that are assigned to genes that do not have 
#' TF annotations (orthogroup).
#' @param CisBP_TF_family_mapping_file character, path to tsv file containing mapping between 
#' CisBP family ids (first column) and TF family ids (second column)
#' @param clean_names, logical, whether to construct shorter and cleaner archetype names
#' @param short_name_target_char numeric, target maximum characters for constructing short archetype name
#' 
#' @return data.table containing motif archetypes and original motifs mapped to provided annotations
#' 
mta_archetype_dictionary <- function(
	arch,
	TF_annotation_file = NULL,
	TF_motifs_file = NULL,
	TF_family_annotation_file = "data/gene_families_searchinfo.tsv",
	CisBP_family_annotation_file = "data/CisBP_2021_08_11_10_18_am_TF_Information.txt",
	CisBP_TF_family_mapping_file = "data/CisBP_Tf_mapping.tsv",
	clean_names = TRUE,
	short_name_target_char = 35
) {
	
	# TF annotations
	if (!is.null(TF_annotation_file)) {
		tf_annot <- fread(TF_annotation_file)
		tf_annot <- unique(tf_annot)
		setkey(tf_annot, gene)
	}
	
	# TF family annotations
	if (!is.null(TF_family_annotation_file)) {
		tf_fam_annot <- fread(TF_family_annotation_file)[,1:2]
		setnames(tf_fam_annot, c("tf_family","tf_class"))
		if (any(grepl(",",tf_fam_annot$tf_class))) {
			tf_fam_annot <- as.data.table(tidyr::separate_rows(tf_fam_annot, "tf_class", sep=","))[tf_class!=""]
		}
		tf_fam_vector <- structure(tf_fam_annot$tf_family, names=tf_fam_annot$tf_class)
	} else {
		tf_fam_vector <- NULL
	}
	
	# TF gene to motif annotations
	if (!is.null(TF_motifs_file)) {
		tf_m_dt <- fread(TF_motifs_file, select = 1:2)
		setnames(tf_m_dt, c("gene", "motif"))
		setkey(tf_m_dt, gene)
		tf_m_dt <- as.data.table(
			tidyr::separate_rows(tf_m_dt, motif, sep = ",", convert = TRUE)
		)
		tf_m_dt <- tf_m_dt[!is.na(gene)]
		tf_m_dt <- tf_m_dt[gene != ""]
		tf_m_dt <- unique(tf_m_dt)
	}
	
	if (all(!is.null(TF_annotation_file), !is.null(TF_motifs_file))) {
		# merge motif and TF annotations
		mdt <- merge.data.table(tf_m_dt, tf_annot, by="gene", all.x = TRUE, allow.cartesian = TRUE)
		mdt <- unique(mdt)
		
	} else if (all(!is.null(TF_annotation_file), is.null(TF_motifs_file))) {
		# use only TF annotationa
		mdt <- tf_annot
	} else if (all(is.null(TF_annotation_file), !is.null(TF_motifs_file))) {
		# use only CisBP
		mdt <- tf_m_dt
	} else {
		mdt <- NULL
	}
	
	# get archetype motifs names
	arch_list <- lapply(arch, function(x) x$ppm_consensus)
	arch_nums <- paste0("ARCH",seq_along(arch_list))
	arch_names <- make.unique(sapply(arch_list, function(x) x@name), sep="__")
	arch_lens <- sapply(arch_list, function(x) ncol(x@motif))
	arch_numots <- sapply(arch, function(x) length(x[[1]]))
	dict <- data.table(
	  archetype = arch_names, 
	  archetype_num = arch_nums,
	  archetype_length = arch_lens, 
	  archetype_num_motifs = arch_numots
	)
	
	# separate entry for original motifs
	dict[, motif:=archetype]
	dict <- tidyr::separate_rows(dict, motif, sep="__", convert = TRUE)
	setDT(dict)
	
	# add gene annotations
	if (!is.null(mdt)) {
		dict <- merge.data.table(dict, mdt, by="motif", all.x=TRUE, sort = FALSE)
		for (j in colnames(mdt))
			data.table::set(dict,which(is.na(dict[[j]])),j,"")
	}
	
	# grep TF family in motif name - this is for motifs that are not from CisBP and not assigned to a gene
	if (!is.null(tf_fam_vector)) {
		all_motifs <- unique(dict$motif)
		if (length(all_motifs)>100)
			message(sprintf("Mapping families to %i motifs, this might take a few minutes",length(all_motifs)))
		all_motifs_fam <- rbindlist(lapply(all_motifs, function(mot) {
			unique(rbindlist(
				lapply(names(tf_fam_vector), function(cls) {
					mot_fam_grep <- grep(cls, mot, ignore.case = TRUE, value = TRUE)
					mot_fam_dt <- data.table(motif=mot_fam_grep)[,tf_class:=cls][,tf_family:=tf_fam_vector[tf_class]]
					unique(mot_fam_dt)
				})
			))
		}))
		dict <- merge.data.table(dict, all_motifs_fam, by="motif", all.x = TRUE, sort=FALSE)
		for (j in c("tf_family","tf_class"))
			data.table::set(dict,which(is.na(dict[[j]])),j,"")
	}
	
	# extract TF family annotations directly from orthogroup
	message("Getting families from OG annotation for TFs")
	if ("og" %in% colnames(dict)) {
		dict[og!="",tf_family:=str_remove(og,"\\.(?<=\\.).+")]                                                                           
	}
	
	# TF family annotations from CisBP
	message("Getting families from CisBP")
	if (!is.null(CisBP_family_annotation_file)) {
		cisbp_ann <- fread(CisBP_family_annotation_file, select=c("Motif_ID","TF_Name","Family_Name"))
		cisbp_ann <- tidyr::separate_rows(cisbp_ann, Motif_ID, sep=",", convert = TRUE)
		setDT(cisbp_ann)
		setnames(cisbp_ann,c("Motif_ID","TF_Name","Family_Name"),c("motif","gene","cisbp_family"))
		# if mappings to our TF familes are provided, use them
		if (!is.null(CisBP_TF_family_mapping_file)) {
			cisbp_tf_map <- fread(CisBP_TF_family_mapping_file)
			setnames(cisbp_tf_map, c("cisbp_family","tf_family","tf_class"))
			cisbp_ann <- merge.data.table(cisbp_ann,cisbp_tf_map,by="cisbp_family")
			cisbp_ann <- cisbp_ann[,.(motif,gene,tf_family,tf_class)]
		} else {
			# otherwise, use CisBP annotations
			setnames(cisbp_ann,"cisbp_family","tf_family")
		}
		tf_family_vector <- structure(cisbp_ann$tf_family, names=cisbp_ann$motif)
		# add TF family to cisBP motifs that are missing it
		dict[grepl("^M0\\d{4}_[12].00",motif) & tf_family=="", 
			 tf_family := tf_family_vector[motif]]
	} else {
		dict[, tf_name:=og]
	}
	setcolorder(dict,intersect(c(
	  "archetype", "archetype_num", "archetype_name", "archetype_name_short",
	  "archetype_length", "archetype_num_motifs",
	  "motif", "source", 
	  "gene", "og", "pfam", "gene_name", "gene_list",
	  "tf_family", "tf_class"
	), colnames(dict)))
	return(dict)
}

#' Create list of archetypes, and a table with mappings to the original motifs
#'
#' Return a flat list of archetypes in `universalmotif` format
#'
#' @param arch list, output of `mta_merge_archetype()` (i.e. a nested list of archetypes and their member ppms)
#' @param use_simple_names whether to return list with simple archetype names
#' @param prefix a prefix for the simple archetype names (default = "ARCH")
#'
#' @return list with two values: a table with motif-to-archetype mappings, and a list of archetypes in `universalmotif` format.
#' 
mta_archetype_list <- function(arch, use_simple_names = TRUE, prefix = "ARCH") {
	
	# get archetype motifs names
	arch_list <- lapply(arch, function(x) x$ppm_consensus)
	arch_nums <- sapply(arch, function(x) length(x[[1]]))
	arch_names <- sapply(arch_list, function(x) x@name)
	arch_lens <- sapply(arch_list, function(x) ncol(x@motif))
	dict <- data.table(archetype = arch_names, archetype_length = arch_lens, archetype_num_motifs = arch_nums)
	
	# separate entry for original motifs
	dict[, motif:=archetype]
	dict <- tidyr::separate_rows(dict, motif, sep="__", convert = TRUE)
	setDT(dict)
	
	dict[, archetype_num:=.GRP,.(archetype)]
	
	# use simple names in table?
	if (use_simple_names) {
		dict[, archetype:=paste(prefix, archetype_num, sep = "_")]
	}
	
	# reorder
	dict = dict [ , c("motif","archetype","archetype_length","archetype_num_motifs") ]
	
	# reassign names to input list?
	if (use_simple_names) {
		for (x in 1:length(arch_list)) {
			arch_list[[x]]@name <- paste(prefix, x, sep = "_")
		}
	} else {
		arch_list = arch_list
	}
	
	# return
	return(list(dict = dict, archetypes = arch_list))
	
}


#' Plot archetype motifs
#' 
#' Plots aligned archetype motif logos
#'
#' @param arch list, output of `mta_merge_archetype()`
#' @param dict optional, data.table with mappings of archetypes to motifs (output of `mta_archetype_dictionary()`)
#' @param dict_archetype_name optional, name of the column in `dict` that contains archetype names (default: "archetype_name_short")
#' Original archetype name will be matched to 'archetype' column in `dict`
#' @param dict_motif_name optional, name of the column in `dict` that contains motif names (default: "tf_cisbp_name"). 
#' Original names will be matche to 'motif' column in `dict`.
#' @param type see `use.type` in `view_motifs`
#' @param skip_singletons logical, skip archetypes that are one single motif (default: FALSE)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default) the plot is returned to stdout
#' @param height,width,res numeric, the width, height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' 
#' @return NULL
#' 
mta_plot_archetype <- function(
	arch,
	dict = NULL,
	dict_archetype_name = "archetype_name_short",
	dict_motif_name = "tf_cisbp_name",
	type = "ICM",
	skip_singletons = FALSE,
	output_file = NULL,
	height = 7, width = 4, res = NA
) {
	arch_num <- sapply(arch, function(x) length(x$ppms))
	if (skip_singletons==TRUE) arch <- arch[arch_num>1]
	plot_list <- lapply(seq_along(arch), function(i) {
	  archetype_result = arch[[i]]
		motifs <- archetype_result$ppms
		if ("universalmotif" %in% class(motifs)) motifs <- list(motifs)
		consensus_motif <- archetype_result$ppm_consensus
		if (!is.null(dict)) {
		  # match motifs names
		  motifs <- sapply(motifs, function(x) {
		    mot_id <- match(x@name, dict$motif)
		    mot_nm <- dict[[I(dict_motif_name)]][mot_id]
		    if (is.na(mot_nm) | mot_nm == "") mot_nm <- x@name
		    x@name <- mot_nm
		    x
		  }, simplify = FALSE, USE.NAMES = FALSE)
		  # match archetype name
		  arch_id <- match(consensus_motif@name, dict$archetype)
			arch_name <- dict[[I(dict_archetype_name)]][arch_id]
			if (is.na(arch_name) | arch_name == "") arch_name <- consensus_motif@name
			consensus_motif@name <- arch_name
		}
		consensus_motif@alphabet <- "DNA"
		consensus_motif@pseudocount <- 1
		motifs[["Consensus"]] <- consensus_motif
		tryCatch({
		  view_motifs(
		    motifs = rev(motifs),
			use.type = type,
		    dedup.names=TRUE,
		    method = "PCC", 
		    normalise.scores = FALSE,
		    relative_entropy = TRUE
		  )
		}, error = function(e) message(sprintf("Failed to plot ARCH%s\n%s",i,e)))
	})
	plotting_function(output_file, width = width, height = height , res = res, print(plot_list))
}

#' Plot archetyping results heatmap
#' 
#' Plots a heatmap of pairwise motif similarities, with motifs grouped by archetypes, and annotated by their source (Homer known, denovo and CisBP).
#'
#' @param sim_mat matrix, motif-motif similarity, rows and column names should be names of `motifs`
#' @param arch list, output of `mta_merge_archetype()`
#' @param dict optional, data.table with mappings of archetypes to motifs (output of `mta_archetype_dictionary()`)
#' @param output_file path to file to which the plot will be saved; if `NULL` (default) the plot is returned to stdout
#' @param height,width,res numeric, the width, height and resolution of plot to be saved (in pixels if png, in inches if pdf)
#' @param return_mat logical, whether to return a matrix of pairwise motif similarities grouped by archetypes and ordered for plotting
#' 
#' @return NULL; if `return_raw=TRUE` returns a matrix of pairwise motif similarities grouped by archetypes and ordered for plotting
#' 
mta_plot_archetype_heatmap <- function(
	sim_mat,
	arch,
	dict = NULL,
	cisbp_dta = NULL, 
	output_file = NULL,
	height = 14, width = 15, res = NA,
	return_mat = FALSE
) {
	
	# ordering
	message("Ordering archetypes for plotting")
	ord <- unlist(lapply(arch, function(x) {
		if (length(x$ppms) == 1) {
			x$ppms@name
		} else {
			unlist(lapply(x$ppms, function(m) m@name),  use.names=FALSE)
		}
	}), use.names=FALSE)
	if (!is.null(dict)) {
	  names(arch) <- unique(dict$archetype_name)
	} else {
	  names(arch) <- paste0("ARCH",seq_along(arch))
	}
	clusts <- unlist(lapply(names(arch), function(i) {
	  ar <- arch[[i]]
	  rep(i, length(ar$ppms))
	}), use.names=FALSE)
	cls <- unique(clusts)
	clmat <- matrix(1, nrow = length(cls), ncol = length(cls), dimnames = list(cls, cls))
	for (cl in cls) {
		cl_m <- ord[clusts==cl]
		other_cl <- setdiff(cls, cl)
		for (clo in other_cl) {
			clo_m <- ord[clusts==clo]
			clmat[cl,clo] <- mean(sim_mat[cl_m,clo_m])
		}
	}
	cl_hc <- hclust(dist(clmat, method = dist_method), method = hclust_method)
	cl_ord <- cl_hc$labels[cl_hc$order]
	cls_ord <- unlist(lapply(cl_ord, function(cl) rep(cl, length(ord[clusts==cl]))), use.names=FALSE)
	m_ord <- unlist(lapply(cl_ord, function(cl) ord[clusts==cl]), use.names=FALSE)
	
	# infer motif annotation
	message("Inferring motif source annotation")
	sdict <- dict[,.(motif,source)][grep("CisBP",source), source := "CisBP"]
	sdict <- sdict[match(m_ord,motif)]
	group <- sdict$source
	names(group) <- m_ord
	group_counts <- table(group)[c("Homer_denovo","Homer_known","CisBP")]
	group_names <- sprintf("%s (%s)",names(group_counts),group_counts)
	# group_colors <- structure(RColorBrewer::brewer.pal(3,"Set1")[1:3],names=group_names)
	group_colors <- structure(c("#e41a1c", "#ff7f00", "#377eb8"),names=group_names)
	group_repl <- structure(group_names,names=names(group_counts))
	group <- str_replace_all(group, group_repl)
	
	# infer motif family
	wdict <- dict[match(cl_ord,archetype_name)]
	if (!is.null(cisbp_dta)) {
		wdict <- merge.data.table(
			wdict, cisbp_dta[, .(Motif_ID, Family_Name)],
			by.x="motif", by.y="Motif_ID", all.x=TRUE, sort=FALSE
		)
		wdict[tf_family=="" & !is.na(Family_Name), tf_family := Family_Name]
		wdict[, Family_Name := NULL]
	}
	wdict_sum <- wdict[, .N, .(archetype_name, tf_family)]
	wdict_sum[, archetype_name := factor(archetype_name, levels=cl_ord)]
	wdict_sum <- wdict_sum[order(archetype_name, -N)][, .SD[tf_family!=""][1], by=archetype_name]
	wdict[, tf_family := NULL]
	wdict <- merge.data.table(wdict, wdict_sum, by="archetype_name", all.x=TRUE, sort=FALSE)
	wdict <- wdict[match(m_ord, motif), ]
	wdict[is.na(tf_family), tf_family := "Unknown"]
	family <- wdict$tf_family
	
	family_counts <- table(family)
	show_family_counts = FALSE
	family_names <- names(family_counts)
	if (show_family_counts) {
		family_names <- sprintf("%s (%s)", names(family_counts), family_counts)
		family_repl <- structure(family_names,names=paste0("^", names(family_counts), "$"))
		family <- str_replace_all(family, family_repl)
	}
	family_colors <- colorRampPalette(brewer.pal(brewer.pal.info["Set3", "maxcolors"], "Set3"))(length(family_counts))
	names(family_colors) <- family_names
	family_colors[grep("Unknown", names(family_colors))] <- "white"


	# ordered and clustered similarity matrix
	mat1 <- sim_mat[m_ord,m_ord]
	mat1 <- pmax(mat1,0)
	mat2 <- rbind(cls_ord)
	nc <- ncol(mat2)
	change_clust <- unname(which(sapply(2:nc, function(i) mat2[,i] != mat2[,i - 1])))
	lead_change_clust <- dplyr::lead(change_clust)
	lead_change_clust[length(lead_change_clust)] <- dim(mat1)[1]
	lab_cord <- dplyr::lead(change_clust) + round((lead_change_clust-change_clust)/2)
	lab_cord <- c(ceiling(change_clust[1]/2), lab_cord)
	col_vals <- c(0,0.5,1) 
	col <- colorRamp2(col_vals, c("#67a9cf", "#f7f7f7", "#ef8a62"))
	
	# motif source annotation
	motif_dt <- data.table(
		motif_name=rownames(mat1),
		group=group,
		family=family
	)
	annotation_dt <- as.data.frame(motif_dt[match(m_ord,motif_name), .(group,family)])
	ha_fun <- function(which, show_legend=TRUE) {
		HeatmapAnnotation(
			which = which, 
			show_legend = show_legend,
			df = annotation_dt,
			col = list(
				"group" = group_colors,
				"family" = family_colors
			)
		)
	}
	
	# archetype names annotation
	anno_lab <- substr(cl_ord,1,80)
	ha_cluster <- HeatmapAnnotation(
		which="row", 
		annotation_width = unit(100, "mm"),
		cluster = anno_mark(
			at = lab_cord,
			labels = anno_lab,
			labels_gp = gpar(fontsize=4), 
			lines_gp = gpar(lwd=0.1)
		)
	)
	
	# archetype family names annotation
	nms <- names(table(mat2[1,]))[table(mat2[1,]) > 10]
	lbs <- wdict[match(nms, archetype_name), ]$tf_family
	ids <- sapply(nms, function(nm) {
		range <- which(mat2[1, ] == nm)
		round(min(range) + ((max(range) - min(range)) / 2))
	})
	ha_cluster <- HeatmapAnnotation(
		which="row", 
		annotation_width = unit(100, "mm"),
		cluster = anno_mark(
			at = ids,
			labels = lbs,
			labels_gp = gpar(fontsize=4), 
			lines_gp = gpar(lwd=0.1)
		)
	)
	
	# heatmap
	message("Building hetmap object")
	hm <- Heatmap(
		mat1, 
		name = "similarity",
		col = col,
		cluster_rows = FALSE, 
		cluster_columns = FALSE,
		border = TRUE,
		show_row_names = FALSE, 
		show_column_names = FALSE,
		row_names_gp = gpar(fontsize = 3),
		column_names_gp = gpar(fontsize = 3),
		left_annotation = ha_fun("row"),
		top_annotation = ha_fun("column", show_legend = FALSE),
		right_annotation = ha_cluster, 
		column_title = sprintf("%i motifs in %i archetypes", ncol(mat1), length(arch))
	)
	
	# Plotting
	message("Plotting")
	plotting_function(output_file, width = width, height = height , res = res, {
		draw(hm, padding = unit(c(5, 2, 5, 2), "cm"))
		decorate_heatmap_body("similarity", {
			for (x in seq_along(change_clust)) {
				i <- change_clust[x]
				j <- change_clust[x+1]
				if (is.na(j)) j <- dim(mat2)[2]
				h <- change_clust[x-1]
				if (length(h)==0) h <- 1
				grid.lines(
					x = i/nc,
					y = c(1-(j/nc) , 1-(h/nc)),
					gp = gpar(lty = 1, lwd = 0.5)
				)
				grid.lines(
					x = c(j/nc,h/nc),
					y = 1-(i/nc),
					gp = gpar(lty = 1, lwd = 0.5)
				)
			}
		})
	})
	
	if (return_mat==TRUE) mat1
}


#' Map acessible genes to motif of the corresponding TF class
#' 
#' Takes as input gene expression and accessibility data, and maps TFs with no previously assigned motifs to 
#' the motif from the corresponding TF class which is most accessible in the selected cluster.
#' 
#' @param mcfp matrix, metacell footprint
#' @param mc integer or character, column or column name in `mcfp` to select
#' @param fc_thrs numeric, expression threshold to use for selecting TFs in `mcfp`
#' @param accessibility matrix, single-cell (columns) accessibility scores for genes (rows)
#' @param accessibility_metadata data.frame with at least the following columns: cell and `cluster_column`
#' @param cluster_column character column name in `accessibility_metadata` with cluster assignment
#' @param cluster character, ids of the clusters (in `cluster_column`) corresponding to selected metacell `mc`, this is used for selecting cells to calculate motif scores
#' @param chromvar matrix, single-cell (columns) chromvar z scores for motifs (rows)
#' @param dict data.frame, assignment of motifs to genes, with annotations (output of `mta_archetype_dictionary()`)
#' It should contain at leas the following columns: cell and `motif_column`
#' @param motif_column character, column name in `dict` with motif names (they should be the same as rownames of `chromvar`)
#' @param assign_nonorphan_motifs logical, whether, after assigning orphan motifs, assign best non-orphan motif in a givenTF class
#' for genes that are still without a motif
#' @param TF_annotation_file character, path to tsv file with TF annotations 
#' Gene ids should be in the first column, followed by additional annotations in other columns: ortogroups and pfam domains
#' @param CisBP_family_annotation_file character, path to tsv file with TF information from CisBP, 
#' it should contain at least Motif_ID and Family columns. This is used to annotate motifs that are assigned to genes that do not have 
#' TF annotations (orthogroup).
#' @param CisBP_TF_family_mapping_file character, path to tsv file containing mapping between 
#' CisBP family ids (first column) and TF family ids (second column)
#' 
#' @return list
#'  * `dict` updated dictionary with newly assigned genes
#'  * `assignments` data.frame with new assignments of genes to orphan motifs
#'  * `assignments_non_orphan` data.frame with new assignments of genes to non-orphan motifs
#'  
mta_assign_motifs_chromvar <- function(
	mcfp, mc, fc_thrs = 1, avg_mc_fun = mean,
	accessibility, accessibility_metadata, cluster_column, cluster, avg_cluster_fun = mean,
	chromvar, 
	dict, motif_column = "archetype_name", assign_nonorphan_motifs = TRUE,
	TF_annotation_file,
	TF_family_annotation_file = "data/gene_families_searchinfo.tsv",
	CisBP_family_annotation_file = "data/CisBP_2021_08_11_10_18_am_TF_Information.txt",
	CisBP_TF_family_mapping_file = "data/CisBP_Tf_mapping.tsv"
) {
	
	# tf annotation data
	tf_annotation <- fread(TF_annotation_file)
	setnames(tf_annotation, c("gene","og","pfam"))
	tf_fam_dt <- fread(TF_family_annotation_file, select=1:2)
	setnames(tf_fam_dt, c("family","class"))
	# add tf family
	tf_annotation[,tf_family:=str_remove(og,"(?<=\\.).+")]                                                                           
	tf_annotation[,tf_family:=str_remove(tf_family,"\\.$")]
	# add tf class
	tf_fam_dt <- tidyr::separate_rows(tf_fam_dt, class, sep=",", convert = TRUE)
	setDT(tf_fam_dt)
	tf_classes <- lapply(tf_annotation$gene, function(g) {
		x = tf_annotation[gene==g]$og
		fam = unique(tf_annotation[og == x]$tf_family)
		cs = tf_fam_dt[family==fam]$class
		ch = sapply(cs, function(p) grepl(p, x, ignore.case = TRUE))
		paste(names(ch)[ch==TRUE], collapse=",")
	})
	names(tf_classes) <- tf_annotation$gene
	tf_annotation[, tf_class:=tf_classes[gene]]
	
	# select TF with FC above threshold in selected metacell
	tf_genes <- intersect(rownames(mcfp),tf_annotation[[1]])
	tffp <- mcfp[tf_genes, mc]
	if (length(mc)>1) {
		tffp <- apply(tffp, 1, avg_mc_fun)
		names(tffp) <- tf_genes
	}
	tfs <- names(tffp[tffp>fc_thrs])
	message(sprintf(
		"prepare genes | %s TFs with expression value over %s in %s", 
		length(tfs), fc_thrs, paste(mc, collapse = ",")
	))
	motifs_to_genes <- copy(dict); setDT(motifs_to_genes); 
	motif_column_order <- unique(motifs_to_genes[[motif_column]])
	tfs_with_motif <- sum(tfs %in% motifs_to_genes$gene)
	tfs_wo_motif <- sum(!tfs %in% motifs_to_genes$gene)
	message(sprintf("prepare genes | TFs with assigned motifs: %s (%.0f%%)", tfs_with_motif, tfs_with_motif/(tfs_with_motif+tfs_wo_motif)*100))
	message(sprintf("prepare genes | TFs without assigned motifs: %s (%.0f%%)", tfs_wo_motif, tfs_wo_motif/(tfs_with_motif+tfs_wo_motif)*100))
	
	# select cluster cells corresponding to selected expression metacell
	mdt <- copy(accessibility_metadata); setDT(mdt)
	setnames(mdt,cluster_column,"selected_cluster_column")
	selected_cluster_cells <- mdt[selected_cluster_column %in% cluster]$cell
	other_cells <- setdiff(mdt$cell,selected_cluster_cells)
	message(sprintf(
		"accessibility | %s cells in cluster %s; %s cells in background", 
		length(selected_cluster_cells), paste(cluster, collapse=","), length(other_cells)
	))
	
	# calculate motif fold change in selected cluster
	message(sprintf(
		"accessibility | motif accessibility in cluster %s", 
		paste(cluster, collapse=",")
	))
	chromvar_fg <- apply(chromvar[,selected_cluster_cells], 1, avg_cluster_fun)
	chromvar_bg <- apply(chromvar[,other_cells], 1, avg_cluster_fun)
	chromvar_fc <- chromvar_fg / chromvar_bg
	chromvar_dt <- data.table(motif=names(chromvar_fc), chromvar_fc, chromvar_fg, chromvar_bg)
	setnames(chromvar_dt,"motif",motif_column)
	motif_dt <- merge.data.table(motifs_to_genes,chromvar_dt,by=I(motif_column),sort=FALSE)
	
	# calculate gene accessibility fold change in selected cluster
	message(sprintf(
		"accessibility | gene accessibility in cluster %s", 
		paste(cluster, collapse=",")
	))
	gene_score_fg <- apply(accessibility[,selected_cluster_cells], 1, avg_cluster_fun)
	gene_score_bg <- apply(accessibility[,other_cells], 1, avg_cluster_fun)
	gene_score_fc <- gene_score_fg / gene_score_bg
	gene_dt <- data.table(gene=names(gene_score_fc), gene_score_fc, gene_score_fg, gene_score_bg)
	gene_dt[,gene:=str_replace_all(gene,"-","_")]
	gene_dt <- gene_dt[gene %in% tf_annotation[[1]]]
	gene_dt <- merge.data.table(tf_annotation, gene_dt, by="gene",sort=FALSE)
	
	# genes not assigned to motifs
	acc_gene <- gene_dt[gene %in% tfs]$gene
	miss_gene <- acc_gene[!acc_gene %in% motifs_to_genes$gene]
	missing_gene <- gene_dt[gene %in% miss_gene]
	message(sprintf("accessibility | accessible TFs without assigned motif: %s", length(miss_gene)))
	
	# motifs not assigned to gene
	motif_dt[, num_genes:=nrow(.SD[gene!=""]), get(motif_column)]
	missing_mot <- motif_dt[num_genes==0]
	message(sprintf("prepare motifs | motifs without a gene: %s", length(unique(missing_mot[[motif_column]]))))
	
	# helper function to assign motifs
	.assignmot <- function(tf_fam, mis_g, mis_m) {
		# select tfs and motifs by family
		mg <- mis_g[tf_family==tf_fam][gene_score_fg>0][order(-gene_score_fc)]
		mm <- mis_m[tf_family==tf_fam][chromvar_fg>0][chromvar_fc>0][order(-chromvar_fc)]
		if (nrow(mm)>0) {
			# split multiple class entries and match by class
			mg <- tidyr::separate_rows(mg, tf_class, sep=",", convert = TRUE); setDT(mg)
			mm <- tidyr::separate_rows(mm, tf_class, sep=",", convert = TRUE); setDT(mm)
			# if there are direct class matches, pair those
			ovl_cl <- intersect(unique(mg$tf_class), unique(mm$tf_class))
			if (length(ovl_cl)>0) {
				ovl_cl_dt <- rbindlist(lapply(ovl_cl, function(cli) {
					# genes in this class
					genes <- mg[tf_class==cli]$gene
					# best chromvar scoring motif
					motif <- mm[tf_class==cli][order(-chromvar_fc)][1][[motif_column]]
					data.table(gene=genes)[,I(motif_column):=motif]
				}))
				ovl_cl_genes <- ovl_cl_dt$gene
			} else {
				ovl_cl_dt <- data.table()
				ovl_cl_genes <- character()
			}
			# for genes without direct class overlap, assign best chromvar scoring motif for this family
			genes <- unique(mg[!gene %in% ovl_cl_genes]$gene)
			motif <- mm[order(-chromvar_fc)][1][[motif_column]]
			noovl_cl_dt <- data.table(gene=genes)[,I(motif_column):=motif]
			# results
			rbindlist(list(ovl_cl_dt, noovl_cl_dt))[,tf_family:=tf_fam]
		}
	}
	
	# assign orphan motifs by tf family to genes with no motif
	tf_fams <- unique(missing_gene$tf_family)
	tf_fams_assign <- rbindlist(lapply(tf_fams, FUN = .assignmot, mis_g = missing_gene, mis_m = missing_mot))
	adt <- copy(tf_fams_assign)[,tf_family:=NULL]
	message(sprintf("done | assigned orphan motifs to %s genes", nrow(tf_fams_assign)))
	
	# for genes with no hits in orphan motifs, get best motif per family, regardless if it is orphan
	if (assign_nonorphan_motifs==TRUE) {
		missing_gene_2 <- missing_gene[!gene %in% tf_fams_assign$gene]
		tf_fams_2 <- unique(missing_gene_2$tf_family)
		tf_fams_assign_2 <- rbindlist(lapply(tf_fams_2, FUN = .assignmot, mis_g = missing_gene_2, mis_m = motif_dt))
		adt_2 <- copy(tf_fams_assign_2)[,tf_family:=NULL]
		adt <- rbindlist(list(adt,adt_2))
		message(sprintf("done | assigned non-orphan motifs to %s genes", nrow(tf_fams_assign_2)))
	} else {
		tf_fams_assign_2 <- NULL
	}
	
	# update dictionary
	# the function should return one-to-one gene to motifs assignments
	motifs_to_genes_updated <- mta_update_dict(
		assignments = adt, 
		dict = motifs_to_genes,
		motif_columnm = motif_columnm,
		tf_annotation = tf_annotation
	)
	
	# results list 
	res_list <- list(
		dict = motifs_to_genes_updated,
		assignments = tf_fams_assign,
		assignments_non_orphan = tf_fams_assign_2
	)
	
	# sanity checks
	if (any(res_list$assignments$gene %in% dict$gene))
		warning("Some of the newly assigned genes are already in the input dicionary, you should inspect this.")
	if (any(res_list$assignments_non_orphan$gene %in% dict$gene))
		warning("Some of the genes assigned to non-orphan motifs are already in the input dicionary, you should inspect this.")
	if (any(res_list$assignments_non_orphan$gene %in% res_list$assignments$gene))
		warning("Some of the newly assigned genes are also assigned to non-orphan motifs, you should inspect this.")
	
	return(res_list)
	
}

#' Helper function to update dictonary with new gene-motif assignments
#' @param assignment data.frame, should contain at least gene and motif_columnm
#' @param dict data.frame, dictionary to update
#' @param motif_column character, column name in both `assignment` and `dict` with motif names
#' @param tf_annotation data.frame with TF annotations to update in the dictionary
mta_update_dict <- function(assignments, dict, motif_columnm, tf_annotation) {
	setDT(dict)
	add_to_dict <- lapply(unique(assignments$gene), function(ng) {
		# gene annotation
		nv = structure(tf_annotation[gene==ng][1], names=colnames(tf_annotation))
		# new motif for gene
		nm = unique(assignments[gene==ng][[motif_column]])
		nl <- rbindlist(lapply(nm, function(nc) {
			# if there are already genes with this motif, keep them
			og = dict[get(motif_column)==nc,]$gene
			# if no genes with this motif, edit in place
			if (all(og=="")) {
				dict[get(motif_column)==nc, names(nv):=nv]
				dict_add <- NULL
				# otherwise keep new entries for dictionary
			} else {
				dict_add <- dict[get(motif_column)==nc,][,names(nv):=nv]
			}
			dict_add
		}))
	})
	dict_add <- rbindlist(add_to_dict)
	if (nrow(dict_add) > 0) {
		motifs_to_genes_updated <- rbindlist(list(dict,dict_add))
		motifs_to_genes_updated[,I(motif_column):=factor(
			motifs_to_genes_updated[[motif_column]],
			levels=motif_column_order
		)]
		setorderv(motifs_to_genes_updated, cols=motif_column)
	} else {
		motifs_to_genes_updated <- dict
	}
	return(motifs_to_genes_updated)
}

###############################
## Motif filtering functions ##
###############################

#' Apply IC block filtering to a list of motifs
#'
#' A block of at least 4 consecutive bases with IC >= 0.5 (ungapped motif)
#' or at least two blocks of at least 3 consecutive bases with IC >= 0.5
#' (gapped motif). This definitition is from Huber and Bulyk (2006)
#' 
#' @params motifs list of motif matrices (`universalmotif` format)
#' @params ic_thr threshold of total information content per motif position, which will be used to identify motif stretches with sufficient IC (default=0.5)
#' @params len_uniblock minimum length of a single stretch necessary to define a valid block (default>=4)
#' @params len_multiblock minimum length of a pair of stretches necessary to define a valid block (default>=3)
#'
#' @return index of motifs to keep
#'
mta_filter_by_ic_block = function(motifs, ic_thr=0.5, len_uniblock=4, len_multiblock=3) {
	
	require("universalmotif")
	
	# convert motifs to IC matrix
	icm = tryCatch(
		universalmotif::convert_type(motifs, "ICM"), 
		error=function(e) universalmotif::convert_type(motifs, "ICM", relative_entropy = TRUE)
	)
	
	# calculate sum of IC values per motif
	icm_sum = lapply(1:length(motifs), function(i) colSums(icm[[i]]) )
	
	# find which motifs have stretches of IC values above threshold, using
	# 2 definitions: a single long block > thr, or two shorter blocks > thr
	keep = lapply(1:length(motifs), function(i) {
		
		# which positions in motif are above ic threshold?
		i_above = as.numeric(icm_sum[[i]] >= ic_thr)
		# how many consecutive bases are above or below the ic thr?
		i_above_runs = rle(i_above)
		# does this motif contain at least one viable block of len 4, or two viable blocks of len 3?
		has_uniblock = any(i_above_runs$values == 1 & i_above_runs$lengths >= len_uniblock)
		has_multiblock = sum(i_above_runs$values == 1 & i_above_runs$lengths >= len_multiblock) > 1
		# return
		return(list(has_uniblock = has_uniblock, has_multiblock=has_multiblock) )
		
	} )
	
	# which motifs survive one or more of these two criteria?
	keep_m = as.data.frame(matrix(unlist(keep), byrow = TRUE, ncol = 2, nrow = length(motifs)))
	ixs_keep = which(apply(keep_m, 1, any))
	
	# return
	return(ixs_keep)
	
}

#' Apply IC filtering to positions in a motif
#'
#' Filter columns of IC matrix for which the total IC is lover than specified threshold.
#' 
#' @param motif or a list of motifs (`universalmotif` format)
#' @param IC_threshold numeric, IC threshold for filtering motif positions
#' @param return_positions logical, if the function should only return which positions are passing IC filtering in each motif, instead of trimming motifs (default: `FALSE`)
#'
#' @return trimmed motif; if `return_positions=TRUE` boolean vector indicating which motifs in the original list pass the threshold
#' 
mta_trim_by_ic <- function(motifs, IC_threshold, return_positions=FALSE) {
	
	require("universalmotif")
	
	if ("universalmotif" %in% class(motifs)) motifs <- list(motifs)
	
	# convert motifs to IC matrix
	consensus_ic <- tryCatch(
		universalmotif::convert_type(motifs, "ICM", pseudocount = pseudocount, relative_entropy = TRUE), 
		error=function(e) universalmotif::convert_type(motifs, "ICM", relative_entropy = TRUE)
	)
	if ("universalmotif" %in% class(consensus_ic)) consensus_ic <- list(consensus_ic)
	
	# filter positions per IC
	motifs_trimmed <- lapply(1:length(motifs), function(i) {
		
		ic_mat <-consensus_ic[[i]]@motif
		filt_ic <- apply(ic_mat,2,sum) > IC_threshold
		
		if (return_positions==TRUE) {
			filt_ic
		} else {
			# trim motif
			col_start <- which(filt_ic==TRUE)[1]
			col_end <- which(filt_ic==TRUE)[length(which(filt_ic==TRUE))]
			motif <- motifs[[i]]
			motif_trimmed_mat <- motif@motif[,col_start:col_end]
			motif_trimmed <- suppressWarnings(universalmotif::create_motif(
				input=motif_trimmed_mat, name=motif@name, pseudocount = motif@pseudocount, type=motif@type, bkg=motif@bkg
			))
			motif_trimmed@alphabet <- motif@alphabet
			motif_trimmed
		}
		
	})
	
	if (length(motifs_trimmed)==1) motifs_trimmed <- motifs_trimmed[[1]]
	
	return(motifs_trimmed)
}


#' Apply nsites filtering to a list of motifs
#' 
#' @params motifs list of motif matrices (`universalmotif` format)
#' @params min_fg_sites min number of foreground sites, default is 10
#' @params min_bg_sites min number of background sites, default is 0 (no filter)
#'
#' @return index of motifs to keep
#'
mta_filter_by_nsites = function(motifs, min_fg_sites = 10, min_bg_sites = 0) {
	
	require("universalmotif")
	
	# which blocks pass foreground and background filters
	ixs_pass_fg = which(unlist(lapply(1:length(motifs), function(n) motifs[[n]]@nsites)   >= min_fg_sites))
	ixs_pass_bg = which(unlist(lapply(1:length(motifs), function(n) motifs[[n]]@bkgsites) >= min_bg_sites))
	
	# return motifs that pass both filters
	ixs_keep = intersect(ixs_pass_fg, ixs_pass_bg)
	
	# return
	return(ixs_keep)
	
}


##########################
## Motif scan functions ##
##########################

#' Calculate genome-wide motif alignment scores and report quantile thresholds from genome-wide distribution
#'
#' @param motifs list of motif matrices (`universalmotif` format). Will be converted to `PWM` format if not already in this format.
#' @param genome_object either a path to a genome file (`.fasta`), or a preloaded `DNAStringSet` object (from `Biostrings::readDNAStringSet`)
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param given_gr GRanges object defining specific regions to analyse. If set, these regions will be used to perform the scanning step. Default is NULL, i.e. the whole genome is used (which can be very slow! Try at least to focus on promoters.). Motif scoring is not constrained by this parameter.
#' @param bin_width width of the genomic bins used to calculate genome-wide max alignment scores for each motif (default is 250bp)
#' @param subsample_fraction fraction of the genome to include in the calculations of genome-wide alignment score distributions (default = 0.10, 10%).
#' @param score_quantiles vector of quantiles to report in the score quantile data.frame (default is `c(0, 0.25, 0.5, 0.75, 0.95, 0.98, 0.99)`)
#' @param score_quantile_thr (default is 0.98): quantile threshold used so as to filter genome-wide motif alignments and report high-confidence hits
#' @param max_chrom_len cut chromosomes to this length (default 1e6 bp) and concatenate fragments (will reduce memory usage). If NULL, scan full chromosomes.
#' @param do_gw_scan whether to do a genome-wide scan of the input motifs. Defaults to TRUE. If FALSE, only subsampled genome-wide score distributions and thresholds are provided.
#' @param nthreads_mot,nthreads_chr num of threads to use for motif scoring & genome-wide scanning, for parallelisation over motifs and chromosomes respectively. Total number of motif/chromosome threads (n*n) should not exceed the total number of available threads, obviously.
#'
#' @return a list with three elements: a filtered `GRanges` object containing the coordinates of high-confidence hits (above threshold), a matrix with per-motif score quantiles, and a matrix with empirical score distributions.
#'
mta_gw_motif_score = function(
	motifs, 
	genome_object,
	index_object = NULL,
	bin_width = 250,
	subsample_fraction = 0.10,
	score_quantiles = c(0, 0.25, 0.5, 0.75, 0.95, 0.98, 0.99),
	score_quantile_thr = 0.98,
	max_chrom_len = 1e6,
	do_gw_scan = TRUE,
	given_gr = NULL,
	nthreads_mot = 2,
	nthreads_chr = 2) {
	
	require("universalmotif")
	require("GenomicRanges")
	require("Biostrings")
	require("PWMEnrich")
	require("doParallel")
	require("plyr")
	
	# register cores
	registerDoParallel(cores = nthreads_mot)
	
	#### Load data ####
	# read genome index
	if (is.null(index_object)) {
		if (!"BSgenome" %in% class(genome_object)) stop("index_object can only be NULL if `genome_object` is BSgenome")
		gix <- data.frame(
			chr = seqnames(genome_object), 
			length = seqlengths(genome_object)
		)
	} else if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`), or NULL if `genome_object` is BSgenome", index_object))
	}
	
	# read genome fasta
	if ("character" %in% class(genome_object)) {
		gen = Biostrings::readDNAStringSet(genome_object, format = "fasta")
	} else if ("DNAStringSet" %in% class(genome_object)) {
		gen = genome_object
	} else if ("BSgenome" %in% class(genome_object)) {
		gen = DNAStringSet(sapply(
			seqnames(genome_object), 
			function(chr) genome_object[[chr]], 
			USE.NAMES = TRUE, 
			simplify = FALSE
		))
	} else {
		stop(sprintf("`genome_object` %s has to be either a DNAStringSet, BSgenome or a path to a fasta file", DNAStringSet))
	}
	
	# if given_gr is given, ignore `max_chrom_len`
	if (!is.null(given_gr)) {
		max_chrom_len = max(width(given_gr))
		warning(sprintf("motif scoring | given_gr is set, ignoring max_chrom_len %i", max_chrom_len))
	}
	
	
	#### Genome-wide scores ####
	# get genome sequences from subsampled bins
	gix_seqlengths = gix$length
	names(gix_seqlengths) = gix$chr
	sub_bin = GenomicRanges::tileGenome(gix_seqlengths, tilewidth = bin_width * 10)
	sub_bin = sample(sub_bin, size = floor(length(sub_bin) * subsample_fraction), replace = FALSE)
	sub_gen = BSgenome::getSeq(gen, sub_bin)
	sub_gen = unlist(sub_gen)
	names(sub_gen) = 1:length(sub_gen)
	message(sprintf("motif scoring | genome subsampled at f=%.2fpc (%ibp)", subsample_fraction * 100, sum(width(sub_gen))))
	
	# find distribution of scores over equally sized windows
	# first	create genome bins 
	sub_seqlengths = width(sub_gen)
	names(sub_seqlengths) = 1:length(sub_gen)
	sub_bin_r = GenomicRanges::tileGenome(sub_seqlengths, tilewidth = bin_width)
	message(sprintf("motif scoring | equally sized genome bins n = %i", length(sub_bin)))
	
	# scan motifs in subsampled genome
	# convert PPM motifs to PWM, in PWMEnrich format
	mot_pwm  = universalmotif::convert_type(motifs, "PWM", pseudocount = 1)
	mot_pwme = universalmotif::convert_motifs(mot_pwm, "PWMEnrich-PWM")
	
	# for each motif, calculate max score in each bin and store vector
	if (length(mot_pwm) == 1) {
		mot_pwm_names = mot_pwm@name
		mot_pwme = list(mot_pwme)
	} else {
		mot_pwm_names = sapply(1:length(mot_pwm), function(i) mot_pwm[[i]]@name)
	}
	# define parallelised loop function
	message(sprintf("motif scoring | genome-wide scores of motifs n = %i", length(mot_pwm_names)))
	fun_loop_scoring = function(i) {
		
		# PWMEnrich scoring function
		message("motif scoring | ", i)
		sub_score_i_r = mta_pwme_scan(pwme = mot_pwme[[i]], gen = sub_gen, pwme_name = mot_pwm_names[i], max_chrom_len = max_chrom_len, nthreads = nthreads_chr)
		
		# overlap between bins and motif alignment sites
		sub_bin_ovs = GenomicRanges::findOverlaps(sub_bin_r, sub_score_i_r)
		
		# find max score per bin
		sub_bin_scores_d = data.frame(from_bin = sub_bin_ovs@from, score = mcols(sub_score_i_r [ sub_bin_ovs@to ])$score)
		sub_bin_scores_d_max = aggregate(score ~ from_bin, data = sub_bin_scores_d, FUN = max)
		sub_bin_scores = rep(-Inf, length(sub_bin_r))
		sub_bin_scores [ sub_bin_scores_d_max$from_bin ] = sub_bin_scores_d_max$score
		
		# save
		return(sub_bin_scores)
		
	}
	# calculate genome-wide vectors of motif scores
	gw_bin_scores_l = plyr::alply(.data = 1:length(mot_pwm_names), .margins = 1,  .fun = fun_loop_scoring, .parallel = TRUE, .inform = TRUE)
	# to matrix
	gw_bin_scores_m = matrix(unlist(gw_bin_scores_l), byrow = FALSE, ncol = length(gw_bin_scores_l))
	colnames(gw_bin_scores_m) = mot_pwm_names
	
	# get quantile distribution and quantile threshold
	gw_score_quantile_l = lapply(gw_bin_scores_l, function(i) quantile(i, score_quantiles))
	gw_score_quantile_t = lapply(gw_bin_scores_l, function(i) quantile(i, score_quantile_thr))
	gw_score_quantile_m = matrix(unlist(gw_score_quantile_l), byrow = TRUE, ncol = length(score_quantiles))
	colnames(gw_score_quantile_m) = paste("q",score_quantiles, sep ="")
	rownames(gw_score_quantile_m) = mot_pwm_names
	
	
	#### Scan motifs ####
	
	if (do_gw_scan) {
		
		gw_scan = mta_gw_scan_to_granges(
			mot_pwm_list = mot_pwme, 
			genome_object = gen,
			mot_pwm_names = mot_pwm_names, 
			thresholds = unlist(gw_score_quantile_t), 
			max_chrom_len = max_chrom_len, 
			nthreads_chr = nthreads_chr, 
			nthreads_mot = nthreads_mot, 
			given_gr = given_gr)
		
	} else {
		
		message(sprintf("motif scanning | omit"))
		gw_scan = NULL
		
	}
	
	# return thresholds & filtered hits
	return(list(
		gw_scan = gw_scan,
		score_distribution = gw_bin_scores_m,
		score_quantiles = gw_score_quantile_m
	))
	
}

#' Convert from `universalmotif` list to `monaLisa` list
#'
#' @param motifs list of motif matrices (`universalmotif` format).
#'
#' @return a list of motifs in `monaLisa` format.
#'
mta_convert_umot_to_monalisa = function(motifs) {
	mot_l = suppressMessages(universalmotif::convert_motifs(motifs, class = "TFBSTools-PWMatrix"))
	mot_l = do.call(TFBSTools::PWMatrixList, mot_l)
	for (motif_ix in 1:length(mot_l@listData)) {
  		mot_l@listData[[motif_ix]]@ID = mot_l@listData[[motif_ix]]@name
	}
	return(mot_l)
}


#' Calculate genome-wide motif alignment scores and report quantile thresholds from genome-wide distribution, using monaLisa (faster than universalmotif!)
#'
#' @param motifs list of motif matrices in monaLisa format
#' @param genome_object either a path to a genome file (`.fasta`), or a preloaded `DNAStringSet` object (from `Biostrings::readDNAStringSet`)
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param given_gr GRanges object defining specific regions to analyse. If set, these regions will be used to perform the scanning step. Default is NULL, i.e. the whole genome is used (which can be very slow! Try at least to focus on promoters.). Motif scoring is not constrained by this parameter.
#' @param bin_width width of the genomic bins used to calculate genome-wide max alignment scores for each motif (default is 250bp)
#' @param subsample_fraction fraction of the genome to include in the calculations of genome-wide alignment score distributions (default = 0.10, 10%).
#' @param score_quantiles vector of quantiles to report in the score quantile data.frame (default is `c(0, 0.25, 0.5, 0.75, 0.95, 0.98, 0.99, 0.995, 0.999, 1.0)`)
#' @param score_quantile_thr (default is 0.95): quantile threshold used so as to filter genome-wide motif alignments and report high-confidence hits
#' @param do_gw_scan whether to do a genome-wide scan of the input motifs. Defaults to TRUE. If FALSE, only subsampled genome-wide score distributions and thresholds are provided.
#' @param nthreads_mot,nthreads_chr num of threads to use for motif scoring & genome-wide scanning, for parallelisation over motifs and chromosomes respectively. Total number of motif/chromosome threads (n*n) should not exceed the total number of available threads, obviously.
#'
#' @return a list with three elements: a filtered `GRanges` object containing the coordinates of high-confidence hits (above threshold), a matrix with per-motif score quantiles, and a matrix with empirical score distributions.
#'
mta_gw_motif_score_monalisa = function(
	motifs, 
	genome_object,
	index_object = NULL,
	bin_width = 250,
	subsample_fraction = 0.10,
	score_quantiles = c(0, 0.25, 0.5, 0.75, 0.95, 0.98, 0.99, 0.995, 0.999, 1.0),
	score_quantile_thr = 0.95,
	do_gw_scan = TRUE,
	given_gr = NULL,
	nthreads = 2) {
	
	require("monaLisa")
	require("GenomicRanges")
	require("Biostrings")
	require("doParallel")
	require("plyr")
	
	#### Load data ####
	
	# read genome index
	message(sprintf("motif scoring | load data..."))
	if (is.null(index_object)) {
		if (!"BSgenome" %in% class(genome_object)) stop("index_object can only be NULL if `genome_object` is BSgenome")
		gix <- data.frame(
			chr = seqnames(genome_object), 
			length = seqlengths(genome_object)
		)
	} else if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`), or NULL if `genome_object` is BSgenome", index_object))
	}
	
	# read genome fasta
	if ("character" %in% class(genome_object)) {
		gen = Biostrings::readDNAStringSet(genome_object, format = "fasta")
	} else if ("DNAStringSet" %in% class(genome_object)) {
		gen = genome_object
	} else if ("BSgenome" %in% class(genome_object)) {
		gen = DNAStringSet(sapply(
			seqnames(genome_object), 
			function(chr) genome_object[[chr]], 
			USE.NAMES = TRUE, 
			simplify = FALSE
		))
	} else {
		stop(sprintf("`genome_object` %s has to be either a DNAStringSet, BSgenome or a path to a fasta file", DNAStringSet))
	}
	
	# if given_gr is given, ignore `max_chrom_len`
	if (!is.null(given_gr)) {
		max_chrom_len = max(width(given_gr))
		warning(sprintf("motif scoring | given_gr is set, ignoring max_chrom_len %i", max_chrom_len))
	}
	
	# motif list
	names_motifs = sapply(1:length(motifs), function(i) motifs@listData[[i]]@name)
	
	
	#### Genome-wide scores ####
	
	# get genome sequences from subsampled bins
	gix_seqlengths = gix$length
	names(gix_seqlengths) = gix$chr
	sub_bin = GenomicRanges::tileGenome(gix_seqlengths, tilewidth = bin_width )
	sub_bin = sample(sub_bin, size = floor(length(sub_bin) * subsample_fraction), replace = FALSE)
	sub_gen = Biostrings::getSeq(gen, sub_bin)
	sub_gen = unlist(sub_gen)
	names(sub_gen) = 1:length(sub_gen)
	message(sprintf("motif scoring | genome subsampled at f=%.2fpc (%.1fMbp)", subsample_fraction * 100, sum(width(sub_gen)) / 1e6 ))
	
	# find distribution of scores over equally sized windows
	# first	create genome bins 
	sub_seqlengths = width(sub_gen)
	names(sub_seqlengths) = 1:length(sub_gen)
	sub_bin_r = GenomicRanges::tileGenome(sub_seqlengths, tilewidth = bin_width)
	message(sprintf("motif scoring | equally sized genome bins n = %i", length(sub_bin)))
	
	# for each motif, find score distribution in genome-wide peaks
	# loop over motifs
	registerDoParallel(cores = nthreads)
	fun_loop_scoring = function(i) {
		
		# log
		if (i == 1 | i %% 100 == 0 | i == length(motifs)) {
			message(sprintf("motif scoring | motif %s/%s", i, length(motifs)))
		}
		
		# motif score distribution
		sub_gen_r = monaLisa::findMotifHits(
			query = motifs[[i]],
			subject = sub_gen,
			min.score = 0, 
			BPPARAM = BiocParallel::MulticoreParam(nthreads)
		)
		
		# overlap between bins and motif alignment sites
		sub_bin_ovs = GenomicRanges::findOverlaps(sub_bin_r, sub_gen_r)
		
		# find max score per bin
		sub_bin_scores_d = data.frame(from_bin = sub_bin_ovs@from, score = mcols(sub_gen_r [ sub_bin_ovs@to ])$score, motif = mcols(sub_gen_r [ sub_bin_ovs@to ])$pwmname)

		# store
		if (nrow(sub_bin_scores_d) > 0) {
			sub_bin_scores_d_max = aggregate(score ~ from_bin, data = sub_bin_scores_d, FUN = max)
			gw_bin_scores_q = quantile(sub_bin_scores_d_max [ , "score" ], score_quantiles)
			names(gw_bin_scores_q) = as.character(score_quantiles)
		} else {
			gw_bin_scores_q = rep(0, length(score_quantiles))
			names(gw_bin_scores_q) = as.character(score_quantiles)
		}
		
		#### Scan motifs ####
		
		if (do_gw_scan) {
			
			if (is.null(given_gr)) {
				
				# scan whole genome
				gw_scan = monaLisa::findMotifHits(
					query = motifs[[i]],
					subject = gen,
					min.score = gw_bin_scores_q [ as.character(score_quantile_thr) ])
				
				# rename some columns
				mcols(gw_scan)$motif       = mcols(gw_scan)$pwmname	
				mcols(gw_scan)$motif_score = mcols(gw_scan)$score
				mcols(gw_scan)$pwmname = NULL
				mcols(gw_scan)$pwmid = NULL
				mcols(gw_scan)$matchedSeq = NULL
				mcols(gw_scan)$score = NULL
					
			} else {
				
				# get sequences to scan
				given_seq = Biostrings::getSeq(gen, given_gr)
				names(given_seq) = paste("s", 1:length(given_seq), sep = "")
				given_scan = monaLisa::findMotifHits(
					query = motifs[[i]],
					subject = given_seq,
					min.score = gw_bin_scores_q [ as.character(score_quantile_thr) ])
					
				# convert given coordinates to genomic coordinates
				given_scan_d = data.frame(sequence = given_scan@seqnames,  seq_start = start(given_scan), seq_stop = end(given_scan), motif = mcols(given_scan)$pwmname, motif_score = mcols(given_scan)$score)
				given_geno_d = data.frame(chr = given_gr@seqnames        , start = start(given_gr),       stop = end(given_gr),       sequence = paste("s", 1:length(given_gr), sep = ""))
				given_scan_d = merge(given_scan_d, given_geno_d, by = "sequence")
				
				# ranges of genomic coordinates 
				gw_scan = GenomicRanges::GRanges(
					seqnames = given_scan_d$chr,
					ranges = IRanges::IRanges(start = given_scan_d$start + given_scan_d$seq_start, end = given_scan_d$start + given_scan_d$seq_stop),
					motif = given_scan_d$motif,
					motif_score = given_scan_d$motif_score
				)

			}
			
		} else {
			
			gw_scan = NULL
			
		}
		
		return(list(gw_scan = gw_scan, score_quantiles = gw_bin_scores_q, motif = names_motifs[i]))
	}
	
	# execute parallelised loop
	gw_scan_l = plyr::alply(.data = 1:length(motifs), .margins = 1,  .fun = fun_loop_scoring, .parallel = TRUE, .inform = TRUE)
	
	# concatenate granges
	message(sprintf("motif scoring | done, conncatenating"))
	gw_scan_l_r = sapply(1:length(motifs), function(i) gw_scan_l[[i]]$gw_scan )
	gw_scan_l_r = unlist(as(gw_scan_l_r, "GRangesList"))
	
	# concatenate score distributions
	gw_scan_l_q = t(sapply(1:length(motifs), function(i) gw_scan_l[[i]]$score_quantiles ))
	rownames(gw_scan_l_q) = names_motifs
	
	# return thresholds & filtered hits
	return(list(
		gw_scan = gw_scan_l_r,
		score_quantiles = gw_scan_l_q
	))
	
}

#' Calculate genome-wide motif alignment scores with an empirical FDR
#'
#' @param motifs list of motif matrices in monaLisa format
#' @param genome_object either a path to a genome file (`.fasta`), or a preloaded `DNAStringSet` object (from `Biostrings::readDNAStringSet`)
#' @param index_object either a path to a genome index file (`.fai`) with chr names and lengths, or a preloaded data.frame with these two columns
#' @param fg_r GRanges object defining specific regions to analyse.
#' @param bg_r GRanges object defining specific regions to use as background. If NULL, background is sampled from the rest of the genome (I recommend you do this).
#' @param bg_not_in_fg if TRUE and bg_r is NULL (i.e. use all gneome) select as background regions not in the foreground (rather than all regions). Default is TRUE.
#' @param bin_width width of the genomic bins used to calculate genome-wide max alignment scores for each motif (default is 250bp)
#' @param gcc_intervals how many GC content bins to use? (breaks are calculated to create bins of approximately equivalent sizes).
#' @param gcc_breaks instead of using N gcc intervals, use these user-provided breaks instead
#' @param min_score_reporting calculate p-values for motif alignments with >= X% of the maximum possible score (default is 0.6, i.e. 60%)
#' @param subsample_fraction fraction of the genome to use as background (default is no subsampling, i.e. subsample_fraction=1, 100%; I recommend not lowering this value too much because, after GC binning, too few bins may remain in the background for meaningful pvalue estimations).
#' @param nthreads num of threads to use
#' @param continuous_output if set to a folder, the function will save one bed file per motif as soon as they're scanned, so that you don't have to wait. Default is NULL (disabled).
#' @param verbose default is FALSE.
#'
#' @return a filtered `GRanges` object containing the coordinates of high-confidence motif hits (above threshold).
#'
mta_gw_motif_score_monalisa_fdr = function(
	motifs, 
	genome_object,
	fg_r,
	bg_r = NULL,
	index_object = NULL,
	bin_width = 250,
	bg_not_in_fg = TRUE,
	gcc_intervals = 20,
	gcc_breaks = NULL,
	min_score_reporting = 0.6,
	subsample_fraction = 1,
	nthreads = 2,
	continuous_output = NULL,
	verbose = FALSE) {
	
	require("monaLisa")
	require("GenomicRanges")
	require("Biostrings")
	require("doParallel")
	require("plyr")
	
	#### Load data ####
	
	# read genome index
	message(sprintf("motif scanning | load data..."))
	if (is.null(index_object)) {
		if (!"BSgenome" %in% class(genome_object)) stop("index_object can only be NULL if `genome_object` is BSgenome")
		gix <- data.frame(
			chr = seqnames(genome_object), 
			length = seqlengths(genome_object)
		)
	} else if ("character" %in% class(index_object)) {
		gix = read.table(index_object, stringsAsFactors = FALSE)
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else if ("data.frame" %in% class(index_object)) {
		gix = index_object
		gix = gix[,1:2]
		colnames(gix) = c("chr","length")
	} else {
		stop(sprintf("`index_object` %s has to be either a data.frame or a path to a genome index (`.fai`), or NULL if `genome_object` is BSgenome", index_object))
	}
	
	# read genome fasta
	if ("character" %in% class(genome_object)) {
		gen = Biostrings::readDNAStringSet(genome_object, format = "fasta")
	} else if ("DNAStringSet" %in% class(genome_object)) {
		gen = genome_object
	} else if ("BSgenome" %in% class(genome_object)) {
		gen = DNAStringSet(sapply(
			seqnames(genome_object), 
			function(chr) genome_object[[chr]], 
			USE.NAMES = TRUE, 
			simplify = FALSE
		))
	} else {
		stop(sprintf("`genome_object` %s has to be either a DNAStringSet, BSgenome or a path to a fasta file", DNAStringSet))
	}
	
	# motif list
	names_motifs = sapply(1:length(motifs), function(i) motifs@listData[[i]]@name)
	
	# get genome sequences from subsampled regions
	gix_seqlengths = gix$length
	names(gix_seqlengths) = gix$chr

	
	## Foreground ##

	# genome bins for foreground (either given, or create from genome)
	fg_seq = Biostrings::getSeq(gen, fg_r)
	if (!is.null(mcols(fg_r)$name)) {
		names(fg_seq) = mcols(fg_r)$name
	} else {
		names(fg_seq) = 1:length(fg_seq)
	}
	message(sprintf("motif scanning | get foreground regions, n=%i", length(fg_r)))	
	
	# GC content in foreground
	fg_seq_gcc = as.numeric(Biostrings::letterFrequency(fg_seq, letters = c("CG")) / width(fg_seq))
	if (is.null(gcc_breaks)) {
		gcc_breaks = quantile(fg_seq_gcc, 1 / gcc_intervals *  c(0:gcc_intervals) , include.lowest = TRUE) # these are used also for bg
	} else {
		gcc_intervals = length(gcc_breaks) - 1
	}
	fg_seq_gcc_category = cut(fg_seq_gcc, breaks = gcc_breaks)
	names(fg_seq_gcc_category) = names(fg_seq)
	message(sprintf("motif scanning | bin foreground regions by GCC in %i intervals, n=%i", gcc_intervals, length(fg_r)))	


	## Background ##

	# genome bins for background (either given, or create from genome)
	if (!is.null(bg_r)) {
		bg_seq = Biostrings::getSeq(gen, bg_r)
		names(bg_seq) = 1:length(bg_seq)
	} else if (bg_not_in_fg) {
		bg_r = unlist(GenomicRanges::tileGenome(gix_seqlengths, tilewidth = bin_width))
		ovs_bgfg = GenomicRanges::findOverlaps(fg_r, bg_r, type = "any", ignore.strand = TRUE)
		bg_r = bg_r [ !1:length(bg_r) %in% ovs_bgfg@to ]
		bg_seq = Biostrings::getSeq(gen, bg_r)
		names(bg_seq) = 1:length(bg_seq)
	} else {
		bg_r = unlist(GenomicRanges::tileGenome(gix_seqlengths, tilewidth = bin_width))
		bg_seq = Biostrings::getSeq(gen, bg_r)
		names(bg_seq) = 1:length(bg_seq)
	}
	message(sprintf("motif scanning | get background regions, n=%i", length(bg_r)))	
	
	# subsample background?
	if (subsample_fraction < 1) {
		
		ixs_subsample = sample(1:length(bg_r), size = subsample_fraction * length(bg_r), replace = TRUE)
		bg_r = bg_r [ ixs_subsample ]
		bg_seq = bg_seq [ ixs_subsample ]
		names(bg_seq) = 1:length(bg_seq)
		message(sprintf("motif scanning | subsample background regions to %.1fpc, n=%i", subsample_fraction*100, length(bg_r)))
		
	}
	
	# GC content in background
	bg_seq_gcc = as.numeric(Biostrings::letterFrequency(bg_seq, letters = c("CG")) / width(bg_seq))
	bg_seq_gcc_category = cut(bg_seq_gcc, breaks = gcc_breaks)
	names(bg_seq_gcc_category) = names(bg_seq)
	message(sprintf("motif scanning | bin background regions by GCC in %i intervals, n=%i", gcc_intervals, length(bg_r)))
	
	# warning if categories are too small
	if (min(table(bg_seq_gcc_category)) < 1e3) {
		message(sprintf("motif scanning | smallest GC bin in background has n=%i elements, too small for pval calculations?", min(table(bg_seq_gcc_category))) )
	}

	## Scan function ##

	# for each motif, find score distribution in genome-wide peaks
	# loop over motifs
	registerDoParallel(cores = nthreads)
	fun_loop_scoring = function(i) {
		
		# log
		if (i == 1 | i %% 10 == 0 | i == length(motifs)) {
			message(sprintf("motif scanning | motif %s (%s/%s) | %s", names_motifs[i], i, length(motifs), Sys.time()))
		} else if (verbose) {
			message(sprintf("motif scanning | motif %s (%s/%s) | %s", names_motifs[i], i, length(motifs), Sys.time()))
		}
		
		# score motif in background
		bg_r_hits = monaLisa::findMotifHits(
			query = motifs[[i]],
			subject = bg_seq,
			min.score = 0, 
			BPPARAM = BiocParallel::MulticoreParam(nthreads)
		)
		# record GC content of each bin
		mcols(bg_r_hits)$gcc_category = bg_seq_gcc_category [ as.character(seqnames(bg_r_hits)) ]
		bg_r_hits = bg_r_hits [ !is.na(mcols(bg_r_hits)$gcc_category) ]
		
		# obtain distributions of scores per GC content bin
		bg_r_hits_d = data.frame(
			peak = as.character(seqnames(bg_r_hits)),
			score = mcols(bg_r_hits)$score,
			gcc_category = mcols(bg_r_hits)$gcc_category)
		bg_r_hits_d = bg_r_hits_d [ order(bg_r_hits_d$peak, -bg_r_hits_d$score), ]
		bg_r_hits_d = bg_r_hits_d [ !duplicated(bg_r_hits_d$peak), ]
		bg_r_hits_l = list()
		for (gci in unique(bg_r_hits_d$gcc_category)) {
			bg_r_hits_l[[gci]] = bg_r_hits_d [ bg_r_hits_d$gcc_category == gci, "score" ]
		}
		
		# scan motif in foreground
		max_possible_score = sum(apply(motifs[[i]]@profileMatrix, 2, max))
		fg_r_hits = monaLisa::findMotifHits(
			query = motifs[[i]],
			subject = fg_seq,
			min.score = max_possible_score * min_score_reporting,
			BPPARAM = BiocParallel::MulticoreParam(nthreads)
		)
		
		# are there any hits?
		if (length(fg_r_hits) > 0) {
			# record GC content of each bin
			mcols(fg_r_hits)$gcc_category = fg_seq_gcc_category [ as.character(seqnames(fg_r_hits)) ]
			fg_r_hits = fg_r_hits [ !is.na(mcols(fg_r_hits)$gcc_category) ]
			
			# top score per fg region
			fg_r_hits_d = data.frame(
				peak = as.character(seqnames(fg_r_hits)),
				start = start(fg_r_hits),
				end = end(fg_r_hits),
				strand = strand(fg_r_hits),
				motif = mcols(fg_r_hits)$pwmname,
				score = mcols(fg_r_hits)$score,
				matchedSeq = mcols(fg_r_hits)$matchedSeq,
				gcc_category = mcols(fg_r_hits)$gcc_category)
			fg_r_hits_d = fg_r_hits_d [ order(fg_r_hits_d$peak, -fg_r_hits_d$score), ]
			fg_r_hits_d = fg_r_hits_d [ !duplicated(fg_r_hits_d$peak), ]
			rownames(fg_r_hits_d) = 1:nrow(fg_r_hits_d)

			# get pval based on bg
			fg_r_hits_d$pval = NA
			for (gci in unique(fg_r_hits_d$gcc_category)) {
				fg_r_hits_l_gci = fg_r_hits_d [ fg_r_hits_d$gcc_category == gci, "score" ]
				bg_r_hits_l_gci = bg_r_hits_l[[gci]]
				fg_pvi = sapply(fg_r_hits_l_gci, function(g) {
					max( 1/length(bg_r_hits_l_gci) , 1 - rank(c(g, bg_r_hits_l_gci), ties.method = "average")[1] / (length(bg_r_hits_l_gci) + 1) )
				} )
				fg_r_hits_d[fg_r_hits_d$gcc_category == gci, "pval"] = fg_pvi
			}
			fg_r_hits_p = fg_r_hits_d$pval
			names(fg_r_hits_p) = fg_r_hits_d$peak
			
			# convert region coordinates to genomic coordinates
			given_scan_d = data.frame(sequence = fg_r_hits_d$peak, seq_start = fg_r_hits_d$start, strand = fg_r_hits_d$strand, seq_stop = fg_r_hits_d$end, motif = fg_r_hits_d$motif, motif_score = fg_r_hits_d$score, motif_matched_seq = as.character(fg_r_hits_d$matchedSeq))
			given_geno_d = data.frame(chr = fg_r@seqnames, start = start(fg_r), stop = end(fg_r), sequence = names(fg_seq))
			given_scan_d = merge(given_scan_d, given_geno_d, by = "sequence")
			
			# ranges of genomic coordinates 
			gw_scan = GenomicRanges::GRanges(
				seqnames = given_scan_d$chr,
				ranges = IRanges::IRanges(start = given_scan_d$start + given_scan_d$seq_start, end = given_scan_d$start + given_scan_d$seq_stop),
				strand = given_scan_d$strand,
				motif = given_scan_d$motif,
				motif_score = given_scan_d$motif_score,
				motif_pval = fg_r_hits_p [ given_scan_d$sequence ],
				motif_matched_seq = given_scan_d$motif_matched_seq,
				region = given_scan_d$sequence
			)
			
			if (!is.null(continuous_output)) {
				dir.create(continuous_output, showWarnings = FALSE)
				mta_granges_to_bed(gw_scan, bed_fn = sprintf("%s/motif.%s.bed", continuous_output, names_motifs[i]), mcols = c("motif","motif_score","motif_pval","motif_matched_seq","region"))
			}

			return(list(motif = names_motifs[i], gw_scan = gw_scan))
			
		} else {
			
			# if nothing to report...
			message(sprintf("motif scanning | motif %s (%s/%s) | shows zero hits at score >= %.1f!", names_motifs[i], i, length(motifs), max_possible_score * min_score_reporting))
			if (!is.null(continuous_output)) {
				dir.create(continuous_output, showWarnings = FALSE)
				file.create(sprintf("%s/motif.%s.bed", continuous_output, names_motifs[i]))
			}

			# return empty ranges
			gw_scan = GenomicRanges::GRanges(seqnames = NULL, ranges = NULL, strand = NULL, motif = NULL, motif_score = NULL, motif_pval = NULL, motif_matched_seq = NULL, region = NULL)
			return(list(motif = names_motifs[i], gw_scan = gw_scan))
			
		}
		
	}
	

	## Run scan function in loop ##

	# execute parallelised loop
	gw_scan_l = plyr::alply(.data = 1:length(motifs), .margins = 1,  .fun = fun_loop_scoring, .parallel = TRUE, .inform = TRUE)

	## Output ##

	# concatenate granges
	gw_scan_l_r = sapply(1:length(motifs), function(i) gw_scan_l[[i]]$gw_scan )
	gw_scan_l_r = unlist(as(gw_scan_l_r, "GRangesList"))
	
	
	# return thresholds & filtered hits
	return(gw_scan = gw_scan_l_r)
	
}


#' Genome-wide scan of motif presence
#' 
#' @param mot_pwm_list list of motifs of `universalmotif` or `PWM` class (from `PWMEnrich`).
#' @param genome_object either a path to a genome file (`.fasta`), or a preloaded `DNAStringSet` object (from `Biostrings::readDNAStringSet`)
#' @param mot_pwm_names motif names to populate output `GRanges` object (default = NULL i.e. derived from `mot_pwm_list` if possible).
#' @param thresholds vector of reporting thresholds for each motif. If NULL, I will apply a fractional threshold to the IC score of each motif (which requires `mot_pwm_list` to be in "universalmotif" format).
#' @param fractional_threshold if thresholds=NULL, I will apply this fractional threshold (default = 0.8) to the IC score of each motif (which requires `mot_pwm_list` to be in "universalmotif" format). If set to NULL, all alignments will be recorded (can be huge!).
#' @param nhits integer, number of top motif hits per bin to report, if NULL (default), all hits above threshold will be reported.
#' @param max_chrom_len cut chromosomes to this length (default 1e6 bp) and concatenate fragments (will reduce memory usage). If NULL, scan full chromosomes.
#' @param given_gr `GRanges` object defining specific regions to analyse. If set, these regions will be used to perform the genome-wide scan. Default is NULL, i.e. the whole genome is used (which can be very slow! Try at least to focus on promoters.)
#' @param nthreads_mot,nthreads_chr num of threads to use for motif scoring & genome-wide scanning, for parallelisation over motifs and chromosomes respectively. Total number of motif/chromosome threads (n*n) should not exceed the total number of available threads, obviously.
#' 
#' @return `GRanges` object with motif presence coordinates (only those above threshold)
#' 
mta_gw_scan_to_granges = function(
	mot_pwm_list, 
	genome_object, 
	mot_pwm_names = NULL, 
	thresholds = NULL, 
	fractional_threshold = 0.8, 
	nhits = NULL,
	max_chrom_len = 1e6, 
	given_gr = NULL, 
	nthreads_chr = 1, 
	nthreads_mot = 1) {
	
	# search for bona fide hits in the genome
	message(sprintf("motif scanning | num motifs n = %i", length(mot_pwm_list)))
	message(sprintf("motif scanning | num chromosomes n = %i", length(genome)))
	
	# register cores
	registerDoParallel(cores = nthreads_mot)
	
	# convert PPM motifs to PWM, in PWMEnrich format
	if (! "PWM" %in% class(mot_pwm_list[[1]])) {
		mot_pwm  = universalmotif::convert_type(mot_pwm_list, "PWM", pseudocount = 1)
		mot_pwme = universalmotif::convert_motifs(mot_pwm, "PWMEnrich-PWM")
	} else {
		mot_pwme = mot_pwm_list
	}
	
	# get motif names if not already provided
	if (is.null(mot_pwm_names) & "universalmotif" %in% class(mot_pwm_list[[1]])) {
		mot_pwm_names = plyr::laply(mot_pwm_list, function(i) i@name)
	} else if (is.null(mot_pwm_names) & "PWM" %in% class(mot_pwm_list[[1]])) {
		mot_pwm_names = plyr::laply(mot_pwm_list, function(i) i$name)
	}
	
	# get thresholds if not already provided
	if (!is.null(thresholds)) {
		message(sprintf("motif scanning | using motif-specific thresholds"))
	} else if (is.null(thresholds) & "universalmotif" %in% class(mot_pwm_list[[1]]) & !is.null(fractional_threshold)) {
		message(sprintf("motif scanning | no thresholds provided, using %.2f * IC score of each motif", fractional_threshold))
		thresholds = plyr::laply(mot_pwm_list, function(i) i@icscore * fractional_threshold)
	} else if (is.null(fractional_threshold)) {
		if (is.null(nhits)) {
			message(sprintf("motif scanning | no thresholds provided, reporting all motif alignments"))
		} else if (nhits>0) {
			message(sprintf("motif scanning | no thresholds provided, reporting top %.0f motif alignments per bin", nhits))
		}
		thresholds = rep(-Inf, length(mot_pwm_list))
	} else {
		stop(sprintf("If you don't provide reporting thresholds, you need to provide input motifs in `universalmotif` format so that I can calculate an ad-hoc threshold based on its IC score and a `fractional_threshold` value %.2f", fractional_threshold))
	}
	
	# read genome fasta
	if ("character" %in% class(genome_object)) {
		gen = Biostrings::readDNAStringSet(genome_object, format = "fasta")
	} else if ("DNAStringSet" %in% class(genome_object)) {
		gen = genome_object
	} else {
		stop(sprintf("`genome_object` %s has to be either a DNAStringSet or a path to a fasta file", DNAStringSet))
	}
	
	# per-motif loop
	fun_loop_scanning = function(i) {
		
		# get threshold for this motif
		mot_thr = thresholds[i]
		mot_thr [ is.na(mot_thr) ] = -Inf
		
		# PWMEnrich scoring function
		message(sprintf("motif scanning | genome-wide scan %s (thr = %.2f)", mot_pwm_names[i], mot_thr))
		if (is.null(given_gr)) {
			gw_scan_i_r = mta_pwme_scan(pwme = mot_pwme[[i]], gen = gen, pwme_name = mot_pwm_names[i], max_chrom_len = max_chrom_len, nthreads = nthreads_chr)
		} else {
			gw_scan_i_r = mta_pwme_scan(pwme = mot_pwme[[i]], gen = gen, pwme_name = mot_pwm_names[i], max_chrom_len = max_chrom_len, nthreads = nthreads_chr, given_gr = given_gr)
		}
		
		# filter by reporting score
		mcols(gw_scan_i_r)$score [ is.na(mcols(gw_scan_i_r)$score) ] = -Inf
		gw_scan_i_r = gw_scan_i_r [ mcols(gw_scan_i_r)$score >= thresholds[[i]] , ]
		if (nhits>0) {
		  gw_scan_i_dt <- setDT(as.data.frame(gw_scan_i_r))
		  setorder(gw_scan_i_dt,-score)
		  gw_scan_i_r_dt <- gw_scan_i_dt[,.SD[1:pmin(nhits,.N)],bin_name]
		  setorder(gw_scan_i_r_dt,seqnames,start,end)
		  gw_scan_i_r <- makeGRangesFromDataFrame(gw_scan_i_r_dt, keep.extra.columns = TRUE)
		}
		
		# return
		return(gw_scan_i_r)
		
	}
	
	# list of hits per motif (as GRanges)
	gw_scan_l = plyr::alply(.data = 1:length(mot_pwm_names), .margins = 1,  .fun = fun_loop_scanning, .parallel = TRUE)
	gw_scan = unlist(as(gw_scan_l, "GRangesList"))
	
	# return
	return(gw_scan)
	
}


#' Convenience function to scan a multi-chromosome genome with PWMEnrich and obtain a GRanges object
#'
#' @param pwme a PWM motif in PWMEnrich format
#' @param gen a genome object in `DNAStringSet` format, contanining chromosome names.
#' @param pwme_name name of the motif (default is NULL, i.e. try to retrieve it from `pwme` object)
#' @param max_chrom_len cut chromosomes to this length (default 1e6 bp) and concatenate fragments (will reduce memory usage). If NULL, scan full chromosomes.
#' @param given_gr GRanges object defining specific regions to score. Default is NULL, i.e. the whole genome is used. This option is incompatible with `max_chrom_len=<int>` (i.e. it must be NULL)
#' @param strand_fun one of "max" or "mean", the two possible functions used by `PWMEnrich::scanWithPWM()` to summarise values over the two strands (default is "max")
#' @param nthreads num of threads to use for parallel processing of chromosomes/chromosome fragments.
#'
#' @return a `GRanges` object containing per-base scores.
#'
mta_pwme_scan = function(pwme, gen, pwme_name = NULL, max_chrom_len = 1e6, given_gr = NULL, strand_fun = "max", nthreads = 1) {
	
	require("PWMEnrich")
	require("plyr")
	require("doParallel")
	
	# register cores
	registerDoParallel(cores = nthreads)
	
	if (!"PWM" %in% class(pwme)) {
		stop(sprintf("Motif object `pwme` should be class \"PWM\" (PWMEnrich library), but it's class %s", class(pwme)))
	}
	if (is.null(pwme_name)) {
		pwme_name = pwme$name
	}
	
	# define chromosome regions to examine
	gen_seqlen = width(gen)
	names(gen_seqlen) = names(gen)
	if (!is.null(max_chrom_len) & is.null(given_gr)) {
		# bin genome
		gw_bin = GenomicRanges::tileGenome(gen_seqlen, tilewidth = max_chrom_len, cut.last.tile.in.chrom = TRUE)
		# add small bins covering inter-bin regions (to avoid intra-chr NA stretches)
		gw_bin_mid = gw_bin
		start(gw_bin_mid) = end(gw_bin_mid) - 50
		end(gw_bin_mid) =   end(gw_bin_mid) + 50
		gw_bin_mid = GenomicRanges::trim(gw_bin_mid)
		gw_bin_mid = gw_bin_mid [ width(gw_bin_mid) > 51, ]
		gw_bin = c(gw_bin, gw_bin_mid)
		gw_bin = sort(gw_bin)
	} else if (!is.null(given_gr)) {
		# use given intervals as bins
		gw_bin = given_gr
	} else {
		# skip genome binning
		gw_bin = GenomicRanges::tileGenome(gen_seqlen, tilewidth = max(gen_seqlen), cut.last.tile.in.chrom = TRUE)
	}
	
	# loop around chromosomes (or fragments of chromosomes)
	fun_loop_chr = function(i) {
		if (i %% 10000 == 0) message("bin", i)
		# create granges
		c = as.character(seqnames(gw_bin)[i])
		gw_scan_c_r = GenomicRanges::GRanges(
			seqnames = c, 
			ranges = IRanges(start = start(gw_bin)[i]:end(gw_bin)[i], end = start(gw_bin)[i]:end(gw_bin)[i]), 
			strand = "*",
			name = pwme_name)
		# if bin is named, keep bin name as metadata col
		if (!is.null(gw_bin$name[i])) {
			mcols(gw_scan_c_r)$bin_name = gw_bin$name[i]
		}
		# calculate scores in this chromosome
		if ( width(gw_bin[i]) >= ncol(pwme$pwm) ) {
			# if chr is of sufficient length...
			gw_scan_c = PWMEnrich::scanWithPWM(
				pwme, 
				dna = gen[[c]][ start(gw_bin)[i]:end(gw_bin)[i] ], 
				both.strands = FALSE, strand.fun = strand_fun)
			mcols(gw_scan_c_r)$score = gw_scan_c
		} else {
			# else, NA
			mcols(gw_scan_c_r)$score = NA
		}
		return(gw_scan_c_r)
	}
	
	message("Looping over ", length(gw_bin), " bins")
	gw_scan_l = plyr::alply(.data = 1:length(gw_bin), .margins = 1,  .fun = fun_loop_chr, .parallel = TRUE, .progress = "none")
	gw_scan_l = as(gw_scan_l [ ! plyr::laply(gw_scan_l, is.null) ], "GRangesList")
	
	# create single granges object
	gw_scan_r = unlist(gw_scan_l)
	gw_scan_r = sort(gw_scan_r)
	
	# drop duplicated positions that have been recorded as NA during parallelisation
	if (!is.null(max_chrom_len)) {
		pos_duplicated = start(gw_scan_r) [ which(duplicated(gw_scan_r)) ]
		ix_to_keep = which( !( is.na(mcols(gw_scan_r)$score) & start(gw_scan_r) %in% pos_duplicated ) )
		gw_scan_r = gw_scan_r [ ix_to_keep ]
		gw_scan_r = unique(gw_scan_r)
	}
	
	# return
	return(gw_scan_r)
	
}


#' Motif enrichment in a set of foreground peaks
#' 
#' @param sites_object either a path to a BED file with motif alignment coordinates, or a preloaded `GRanges` object; It ahould have `name` metadata column
#' @param fg_object,bg_object either a path to a BED file with foreground/background regions respectively (e.g. ATAC peaks, or gene promoter regions), or preloaded `GRanges` objects
#' @param peaks_object either a path to a BED file with peaks/promoter regions, or a preloaded `GRanges` object. Its `name` metadata column must match lists provided with `fg_list` and `bg_list`.
#' @param fg_list,bg_list vectors of names of peaks/regions in the `peaks_object` that make up the foreground and background respectively. If `bg_list` is unset (i.e. NULL), bakcground is defined as anything not in the `fg_list`.
#' @param thresholds_vector named vector of thresholds for each motif (names are motif names; all )
#' @param label, a string identifying the cell type enrichment you're working with (useful if you're running this function for many cell types in a loop).
#' @param nthreads num threads to parallelise enrichment calculations
#' @param pval_adjust pvalue adjustment method; you can use either "empirical_fdr" or any method from `p.adjust` (defaults to "BH"; use NULL to omit).
#' 
#' @return data.frame with enrichment statistics for each motif in this set of background/foreground peaks
#' 
mta_motif_enrichment_test = function(
	sites_object,
	fg_object = NULL,
	bg_object = NULL,
	peaks_object = NULL,
	fg_list = NULL,
	bg_list = NULL,
	thresholds_vector = NULL,
	label = NULL,
	nthreads = 2,
	pval_adjust = "BH") {
	
	# register cores
	registerDoParallel(cores = nthreads)
	
	#### Load data ####
	
	# read motif coordinates
	if ("character" %in% class(sites_object)) {
		sit_r = mta_bed_to_granges(sites_object)
		names(sit_r) = sit_r$name
	} else if ("GRanges" %in% class(sites_object)) {
		sit_r = sites_object
		names(sit_r) = sit_r$name
	} else {
		stop(sprintf("`sites_object` %s has to be either a GRanges object, or a path to a compatible BED file", sites_object))
	}
	message(sprintf("motif enrichment | load motif BED, n = %i sites", length(sit_r)))
	
	# filter motif sites if necessary
	if (!is.null(thresholds_vector)) {
		if (any(!sit_r$name %in% names(thresholds_vector))) {
			print("Missing motif thresholds:")
			print( head(sit_r$name[which(!sit_r$name %in% names(thresholds_vector))]) )
			stop("Some motifs in the `sites_object` are not present in the `thresholds_vector`! Can't proceed with motif filtering.")
		}
		sit_r = sit_r [ sit_r$score >= thresholds_vector [ sit_r$name ]  ]
		message(sprintf("motif enrichment | filter motif BED by score, n = %i sites kept", length(sit_r)))
	}
	
	# Read foreground and background. Two possibilities:
	# 1) provide bed/granges object with fg and bg regions
	# 2) provide a single bed/granges object with regions of interest, and two vectors of fg and bg regions from this object
	if (!is.null(fg_object) & !is.null(bg_object)) {
		
		# read regions: foreground
		if ("character" %in% class(fg_object)) {
			fg_r = mta_bed_to_granges(fg_object)
		} else if ("GRanges" %in% class(fg_object)) {
			fg_r = fg_object
		} else {
			stop(sprintf("`fg_object` %s has to be either a GRanges object, or a path to a compatible BED file", fg_object))
		}
		if (!"name" %in% colnames(mcols(fg_r)) | length(unique(fg_r$name)) < length(fg_r)) 
			fg_r$name <- paste0("fg_r_",1:length(fg_r))
		names(fg_r) = fg_r$name
		message(sprintf("motif enrichment | load fg BED, n = %i regions", length(fg_r)))
		
		# read regions: background
		if ("character" %in% class(bg_object)) {
			bg_r = mta_bed_to_granges(bg_object)
		} else if ("GRanges" %in% class(bg_object)) {
			bg_r = bg_object
		} else {
			stop(sprintf("`bg_object` %s has to be either a GRanges object, or a path to a compatible BED file", bg_object))
		}
		if (!"name" %in% colnames(mcols(bg_r)) | length(unique(bg_r$name)) < length(bg_r)) 
			bg_r$name <- paste0("bg_r_",1:length(bg_r))
		names(bg_r) = bg_r$name
		message(sprintf("motif enrichment | load bg BED, n = %i regions", length(bg_r)))
		
	} else if (!is.null(peaks_object) & !is.null(fg_list)) {
		
		# read regions: all peaks
		if ("character" %in% class(peaks_object)) {
			pks_r = mta_bed_to_granges(peaks_object)
		} else if ("GRanges" %in% class(peaks_object)) {
			pks_r = peaks_object
		} else {
			stop(sprintf("`peaks_object` %s has to be either a GRanges object, or a path to a compatible BED file", peaks_object))
		}
		if (!"name" %in% colnames(mcols(pks_r)) | length(unique(pks_r$name)) < length(pks_r)) 
			pks_r$name <- paste0("pks_r_",1:length(pks_r))
		names(pks_r) = pks_r$name
		message(sprintf("motif enrichment | load peaks BED, n = %i regions", length(pks_r)))
		
		# create fg and bg ranges from given lists of names
		fg_r = pks_r [ names(pks_r) %in% fg_list ]
		if (is.null(bg_list)) {
			bg_r = pks_r [ ! names(pks_r) %in% fg_list ]
		} else {
			bg_r = pks_r [ names(pks_r) %in% bg_list ]
		}
		
	} else {
		
		stop("I can't build foreground/background ranges! I need to have either `fg_object`+`bg_object` OR `peaks_object`+`fg_list`+`bg_list`.")
		
	}
	
	
	#### Intersect ####
	
	# intersect foreground/background with motif coordinates
	mot_fg_ovs = GenomicRanges::findOverlaps(sit_r, fg_r)
	mot_bg_ovs = GenomicRanges::findOverlaps(sit_r, bg_r)
	message(sprintf("motif enrichment | sites in fg regions n = %i", length(mot_fg_ovs)))
	message(sprintf("motif enrichment | sites in bg regions n = %i", length(mot_bg_ovs)))
	
	# list of motifs
	list_mot = sort(unique(names(sit_r)))
	# list_pks = unique(c( names(fg_r), names(bg_r) ))
	
	
	#### Motif counts & FC in fg/bg ####
	
	# matrix of foregrounds
	mat_fg = data.frame(
		mot = names(sit_r) [ mot_fg_ovs@from ],
		pks = names(fg_r) [ mot_fg_ovs@to ]
	)
	mat_fg$mot = factor(mat_fg$mot, levels = list_mot)
	mat_fg = unique(mat_fg)
	mat_fg_m = as.matrix(table(mat_fg$mot, mat_fg$pks, useNA = "ifany"))
	
	# matrix of backgrounds
	mat_bg = data.frame(
		mot = names(sit_r) [ mot_bg_ovs@from ],
		pks = names(bg_r) [ mot_bg_ovs@to ]
	)
	mat_bg$mot = factor(mat_bg$mot, levels = list_mot)
	mat_bg = unique(mat_bg)
	mat_bg_m = as.matrix(table(mat_bg$mot, mat_bg$pks, useNA = "ifany"))
	
	# merge matrices
	ixs_fg = 1:ncol(mat_fg_m)
	ixs_bg = (ncol(mat_fg_m) + 1) : ( ncol(mat_fg_m) + ncol(mat_bg_m) )
	mat = cbind(mat_fg_m, mat_bg_m)
	
	#### Enrichment ####
	
	# frequencies of each motif in fg/bg
	mot_fg_f = rowSums(mat[,ixs_fg,drop=FALSE]) / length(fg_r)
	mot_bg_f = rowSums(mat[,ixs_bg,drop=FALSE]) / length(bg_r)
	
	# per-motif enrichment table
	enri_mot = data.frame(
	  motif = rownames(mat),
	  fc = mot_fg_f / mot_bg_f,
	  f_fg = mot_fg_f,
	  f_bg = mot_bg_f,
	  n_fg = rowSums(mat[,ixs_fg,drop=FALSE]),
	  n_bg = rowSums(mat[,ixs_bg,drop=FALSE])
	)
	
	# hypergeometric test
	message(sprintf("motif enrichment | hypergeometric test..."))
	phyper_fun = function(i, ixs_fg, ixs_bg, fg_size, bg_size) { 
	  fg_counts = sum(i[ixs_fg] > 0)
	  bg_counts = sum(i[ixs_bg] > 0)
	  # nou meu
	  phyper(fg_counts - 1, m = fg_counts + bg_counts, n = fg_size + bg_size - fg_counts - bg_counts, k = fg_size, lower.tail = FALSE)
	}
	enri_mot$pval = apply(mat, 1, FUN = function(i) 
	  phyper_fun(i, ixs_fg=ixs_fg, ixs_bg=ixs_bg, fg_size=length(fg_r), bg_size=length(bg_r))
	)
	
	# adjust pvalues?
	if (!is.null(pval_adjust)) {
		
		message(sprintf("motif enrichment | adjust pval %s", pval_adjust))
		
		if (pval_adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
			
			enri_mot$padj = p.adjust(enri_mot$pval, method = pval_adjust)
			
		} else if (pval_adjust == "empirical_fdr") {
			
			# number of permutation iterations
			n_iter = 100
			# permute pvalues
			fun_permute = function(i) {
				map = mat[ , sample(ncol(mat), replace = FALSE) ]
				pvp = apply(map, 1, function(i) phyper_fun(i, ixs_fg=ixs_fg, ixs_bg=ixs_bg) )
			}
			p_perm = unlist(plyr::alply(1:n_iter, 1, .fun = fun_permute, .parallel = TRUE))
			
			# sort original pvalues
			p_test = enri_mot$pval
			names(p_test) = rownames(enri_mot)
			p_test = sort(p_test)
			n_test = length(p_test)
			# rank concatenated list of pvalues
			# the rank of each p-value (minus its rank in the real p-values) represents how many randomised cases were better
			p_rank = rank( c(p_test, p_perm) ) [ 1:n_test ]    
			# the FDR is the fraction of false discoveries (n_fd) from the real pvalues
			n_fd   = p_rank - 1:n_test
			p_fdr  = n_fd / ( n_iter * n_test )
			# get adjusted pvalues
			p_tadj = p_test
			p_tadj [ order(p_tadj) ] = p_fdr
			# reorder adjusted pvals and add them to original table
			enri_mot$padj = p_tadj [ rownames(enri_mot) ]
			
		} else {
			
			warning(sprintf("unknown pval_adjust method %s, return NA"), pval_adjust)
			enri_mot$padj = NA
			
		}
	}
	
	# do we have a label for this run?
	if (!is.null(label)) {
		enri_mot$label = label
	} else if (is.null(label) & !is.null(fg_object) & "character" %in% class(fg_object)) {
		enri_mot$label = fg_object
	} else { 
		enri_mot$label = NA
	}
	
	# return
	return(enri_mot)
	
}

#' Max motif scores per peak
#' 
#' @param peak_object either a path to a BED file with peaks, or a preloaded `GRanges` object. `name` field is peak name.
#' @param peak_gene_mapping named vector mapping peak ids (names) to genes (values). If provided, a matrix of max scores per gene is provided (alongside the usual max scores per peak).
#' @param motif_scan_object either a path to a BED file with motif genome-wide scan, or a preloaded `GRanges` object. `name` field is motif name. `score` field is motif alignment score at that position.
#' @param motif_scan_name_ix,motif_scan_score_ix index of columns in the BED file that will be used to find motif name and motif score, respectively (by default, 7th and 8th). Integer.
#'
#' @return matrix with max scores of each motif per peak
#'
mta_per_peak_score_profiles = function(
    peak_object,
    motif_scan_object,
    peak_gene_mapping = NULL,
    motif_scan_name_ix = 7,
    motif_scan_score_ix = 8
) {
  
  # read peaks
  if ("character" %in% class(peak_object)) {
    pka_r = mta_bed_to_granges(peak_object)
    names(pka_r) = pka_r$name
  } else if ("GRanges" %in% class(peak_object)) {
    pka_r = peak_object
    names(pka_r) = pka_r$name
  } else {
    stop(sprintf("`peak_object` %s has to be either a GRanges object, or a path to a compatible BED file", peak_object))
  }
  
  # read motif scan object
  if ("character" %in% class(motif_scan_object)) {
    bed_r = mta_bed_to_granges(motif_scan_object, metadata_col_ix = c(motif_scan_name_ix, motif_scan_score_ix), metadata_col_name = c("motif_name","motif_score"))
    mcols(bed_r)$score = mcols(bed_r)$motif_score
    mcols(bed_r)$name  = mcols(bed_r)$motif_name
    names(pka_r) = pka_r$name
  } else if ("GRanges" %in% class(motif_scan_object)) {
    bed_r = motif_scan_object
    if (is.null(bed_r$score) | is.null(bed_r$name)) {
      stop("GRanges object lacks motif `score` or `name` metadata fields")
    }
  } else {
    stop(sprintf("`motif_scan_object` %s has to be either a GRanges object, or a path to a compatible BED file", peak_object))
  }
  
  # overlap between peaks and motif alignment sites
  message("max score per peak | find overlaps...")
  ovs_pka_bed = GenomicRanges::findOverlaps(pka_r, bed_r)
  
  # find max score per peak
  message("max score per peak | sort peak by motif score...")
  ovs_pka_bed_d = data.frame(
    peak = mcols(pka_r [ ovs_pka_bed@from ])$name,
    motif = mcols(bed_r [ ovs_pka_bed@to ])$name,
    score = mcols(bed_r [ ovs_pka_bed@to ])$score)
  ovs_pka_bed_d = ovs_pka_bed_d [ order(ovs_pka_bed_d$peak, ovs_pka_bed_d$motif, -ovs_pka_bed_d$score), ]
  ovs_pka_bed_d_max = ovs_pka_bed_d [ !duplicated(paste(ovs_pka_bed_d$peak, ovs_pka_bed_d$motif)), ]
  ovs_pka_bed_d_max$motif = factor(ovs_pka_bed_d_max$motif, levels = sort(unique(ovs_pka_bed_d_max$motif)))
  
  # as matrix
  message("max score per peak | build peak matrix...")
  mat_pm_s = reshape2::dcast(ovs_pka_bed_d_max, motif ~ peak, value.var = "score", drop = TRUE)
  rownames(mat_pm_s) = mat_pm_s$motif
  mat_pm_s$motif = NULL
  # mat_pm_s [ is.na(mat_pm_s) ] = 0
  
  # do it for genes as well?
  if (!is.null(peak_gene_mapping)) {
    
    # add gene		
    message("max score per peak | build gene matrix...")
    ovs_pka_bed_d_max$gene = peak_gene_mapping [ ovs_pka_bed_d_max$peak ]
    # as matrix
    mat_gm_s = reshape2::dcast(ovs_pka_bed_d_max [ !is.na(ovs_pka_bed_d_max$gene), ], motif ~ gene, value.var = "score", drop = TRUE, fun.aggregate = max)
    rownames(mat_gm_s) = mat_gm_s$motif
    mat_gm_s$motif = NULL
    # mat_gm_s [ is.na(mat_gm_s) ] = 0
    
  } else { 
    
    mat_gm_s = NULL
    
  }
  
  return(list(peak_score_matrix = mat_pm_s, gene_score_matrix = mat_gm_s))
  
}


#' Determine binding of motifs in peaks (aka TOBIAS BINDetect)
#' 
#' The function takes precalculated motifs scores and footprint bigwig files
#' from TOBIAS ScoreBigwig as input to estimate bound TF binding sites.
#' 
#' @param mts_fn data.frame or path to a file with with motif hits, 
#'   it should contain at least the following columns: 
#'   sqnames, start, end (coordinates of the motif binding site), 
#'   and motif (motif name)
#' @param ftp_fn path to bigwig file(s) with footprint scores
#' @param pks_smpl numeric, fraction of peaks from each sample to use 
#'   for background distribution calculation
#' @param pval p value threshold for footprint scores
#' @param seed integer, seed for reproducibility
#' @param normalize_motif_scores logical, quantile normalize motif scores 
#'   (this does not affect binding detection)
#' 
#' @return list with the following elements:
#'   motifs_df, motifs scores data frame, with added footprint binding data
#'   background_dist, list of per-sample background distributions
#'   
mta_binding_detect <- function(
  mts_fn, ftp_fn,
  pks_smpl = 1e6, pval = 1e-4, seed = 1950,
  normalize_motif_scores = FALSE
) {
  
  set.seed(seed)

  # motif scores
  if ("character" %in% class(mts_fn)) {
    mts <- fread(mts_fn)
  } else if ("data.frame"  %in% class(mts_fn)) {
    mts <- as.data.table(mts_fn)
  }
  mts_gr <- makeGRangesFromDataFrame(mts, keep.extra.columns = TRUE)
  mts_gr <- unique(mts_gr)
  # take middle position of each motif 
  # https://github.com/loosolab/TOBIAS/issues/97#issuecomment-951647898
  mts_gr <- resize(mts_gr, width = 1, fix = "center")

  # footprint scores
  if (length(ftp_fn) == 1) {
    if (!grepl("(bigwig|bw)$", ftp_fn)) {
      ftp_fn <- list.files(ftp_fn, pattern = "(bw|bigwig|bigWig)$", full.names = TRUE)
    }
  }
  if (is.null(names(ftp_fn))) {
	names(ftp_fn) <- str_remove(basename(ftp_fn), ".(bw)|(bigwig)")
  }
  ftp_list <- sapply(names(ftp_fn), function(smp) {
	rtracklayer::import(ftp_fn[smp])
  }, USE.NAMES = TRUE, simplify = FALSE)
  
  # overlap footprints and motif scores
  message("motif scores | getting footprints for motif hits...")
  bnd_dt <- rbindlist(lapply(names(ftp_fn), function(smp) {
	message("motif scores | ", smp)
	ftp_gr <- ftp_list[[smp]]
    ovl <- findOverlaps(query = ftp_gr, subject = mts_gr)
    ovl_score <- ftp_gr[queryHits(ovl)]$score
    ovl_mots <- mts_gr[subjectHits(ovl)]
    ovl_mots$sample <- smp
    ovl_mots$ftp_score <- ovl_score
    ovl_dt <- as.data.table(ovl_mots)

	# top score per motif
	# setorder(ovl_dt, -ftp_score)
	# ovl_dt <- ovl_dt[, .SD[1], .(peak, sample, motif)][]

	# quantile normalize motif scores distribution
	if (normalize_motif_scores == TRUE) {
		bnd_dtc <- dcast.data.table(ovl_dt, peak ~ motif, value.var = "score")
		bnd_mt <- data.matrix(bnd_dtc[,-1])
		bnd_norm_mt <- preprocessCore::normalize.quantiles(bnd_mt)
		colnames(bnd_norm_mt) <- colnames(bnd_mt)
		rownames(bnd_norm_mt) <- bnd_dtc[[1]]
		bnd_norm_dt <- melt.data.table(
			as.data.table(bnd_norm_mt, keep.rownames = "peak"),
			id.vars = "peak", variable.name = "motif", value.name = "norm_score"
		)
		ovl_dt <- merge.data.table(
			ovl_dt, bnd_norm_dt,
			by = c("peak", "motif"),
			all.x = TRUE,
			sort = FALSE
		)
		setcolorder(ovl_dt, c("seqnames", "start", "end", "width", "strand", "score", "norm_score"))
	}
	ovl_dt
	
  }))
  bnd_dt[, seqnames:= factor(seqnames, levels = levels(seqnames(mts_gr)))]
  setorder(bnd_dt, seqnames, start, end)

  # sample peaks for background dist
  message("footprints | sampling footprint scores for background distribution...")
  ftp_dt <- rbindlist(lapply(names(ftp_list), function(smp) {
	ftp_gr <- ftp_list[[smp]]
	as.data.table(ftp_gr)[, sample := smp]
  }))
  setnames(ftp_dt, "score", "ftp_score")
  bgr_dt <- ftp_dt[, .SD[sample(
    seq_len(.N), size = pks_smpl, replace = ifelse(pks_smpl < .N, FALSE, TRUE)
  )], sample][, .(ftp_score, sample)]

  # quantile normalize background distribution
  bgr_dt[, peak := paste0("bgp", 1:.N), sample]
  bgr_dtc <- dcast.data.table(bgr_dt, peak ~ sample, value.var = "ftp_score")[,-1]
  bgr_mt <- data.matrix(bgr_dtc)
  bgr_norm_mt <- preprocessCore::normalize.quantiles(bgr_mt)
  colnames(bgr_norm_mt) <- colnames(bgr_mt)
  
  # determine footprint threshold at given pvalue
  message("footprints | determining binding thresholds at pvalue ", pval)
  bkg_dist <- fitdistrplus::fitdist(as.vector(bgr_norm_mt), "norm")
  bnd_thr <- qnorm(
	1 - pval,
	mean = bkg_dist$estimate[1],
	sd = bkg_dist$estimate[2]
  ) 

  # label motif binding sites as bound or unbound
  bnd_dt[, ftp_threshold := bnd_thr]
  bnd_dt[, bound := FALSE]
  bnd_dt[ftp_score > ftp_threshold, bound := TRUE]

  out_list <- list(
    motifs_df = bnd_dt,
	background_dt = bgr_dt,
    background_dist = bkg_dist
  )
  return(out_list)
}


#' Differential binding analysis (aka TOBIAS BINDetect)
#' 
#' @param motifs_df data.frame with motif binding scores,
#'   from the output of `mta_binding_detect()`
#' @param background_dt data.frame with background binding scores,
#'   from the output of `mta_binding_detect()`
#' @param comparisons character, samples to compare separated by `_vs_` (e.g. `treatment_vs_control`),
#'   you can indicate more than one comaprison
#' @param pseudo_score numeric, small constant to add to binding footprint scores for log2FC calculation,
#'   if NA, minimum value of footprint score is used
#' @param seed
#' 
mta_binding_compara <- function(
	motifs_df, background_dt, comparisons, pseudo_score = NA, seed = 1950
) {
	set.seed(seed)
	sapply(comparisons, function(compara) {
		message("Starting comparison ", compara)
		lvs <- strsplit(compara, "_vs_")[[1]]
		lv1 <- lvs[1]
		lv2 <- lvs[2]
		# select samples to compare
		cdt <- motifs_df[sample %in% lvs]
		bdt <- background_dt[sample %in% lvs]
		# transform data
		if (is.na(pseudo_score)) {
			pseudo_score <- min(c(
				cdt[!is.infinite(ftp_score)]$ftp_score,
				bdt[!is.infinite(ftp_score)]$ftp_score
			), na.rm = TRUE)
		}
		if (all(c("peak", "gene") %in% colnames(cdt))) {
			fo <- seqnames + start + end + motif + peak + gene + bound ~ sample
		} else if ("peak" %in% colnames(cdt)) {
			fo <- seqnames + start + end + motif + peak + bound ~ sample
		} else if ("gene" %in% colnames(cdt)) {
			fo <- seqnames + start + end + motif + gene + bound ~ sample
		} else {
			fo <- seqnames + start + end + motif + bound ~ sample
		}
		cdt <- dcast.data.table(cdt, fo, value.var = "ftp_score")
		cdt[is.na(get(lv1)), I(lv1) := 0]
		cdt[is.na(get(lv2)), I(lv2) := 0]
		bdt <- dcast.data.table(bdt, peak ~ sample, value.var = "ftp_score")
		bdt[is.na(get(lv1)), I(lv1) := 0]
		bdt[is.na(get(lv2)), I(lv2) := 0]
		# calculate fold change per TFBS
		cdt[, log2FC := log2(
			(cdt[[lv1]] + pseudo_score) /
			(cdt[[lv2]] + pseudo_score)
		)]
		bdt[, log2FC := log2(
			(bdt[[lv1]] + pseudo_score) /
			(bdt[[lv2]] + pseudo_score)
		)]
		# calculate differential binding score per TF
		cdt[, mean := mean(log2FC), motif]
		cdt[, sd := sd(log2FC), motif]
		bmn <- bdt[, .(mean = mean(log2FC))]
		bsd <- bdt[, .(sd = sd(log2FC))]
		cdt[, diff_bind_score := (abs(mean) - abs(bmn)) / ((sd + bsd) / 2), motif]
		cdt[, diff_bind_pval := ks.test(
			.SD$log2FC, sample(bdt$log2FC, 100),
			alternative = "two.sided"
		)["p.value"], motif]
		cdt[, diff_bind_qval := p.adjust(pval, method = "fdr", n = nrow(.SD)), motif]
		cdt
	}, simplify = FALSE, USE.NAMES = TRUE)	
}

#' Function to plot GenomicRanges, primarily for inspecting `mta_reduce_motif_hits()` results
#' 
#' @param gr GRanges object
#' @param label column name in the GRanges object to use as label
#' @param color column name in the GRanges object to use as color
#' @param xlims numeric vector with x-axis limits
#' 
mta_plot_granges <- function(gr, label = NULL, color = NULL, xlims = c(NA, NA)) {
  
  gr_dt <- as.data.frame(gr)
  setDT(gr_dt)
  setorder(gr_dt, seqnames, start, end)
  gr_dt[, id := as.factor(.I)]
  
  gr_gp <- ggplot(gr_dt, aes(x = start, xend = end, y = id, yend = id))
  
  if (!is.null(color)) {
    gr_gp <- gr_gp + 
      geom_segment(aes(color = get(color)), linewidth = 2) +
      scale_color_viridis_c(name = color)
  } else {
    gr_gp <- gr_gp + 
      geom_segment(color = "black", linewidth = 2)
  }
  
  gr_gp <- gr_gp +
    geom_vline(
      xintercept = unique(gr_dt$start), 
      color = "grey", linewidth = 0.05, linetype = "dotted"
    ) +
    geom_vline(
      xintercept = unique(gr_dt$end), 
      color = "grey", linewidth = 0.05, linetype = "dotted"
    ) +
    scale_x_continuous(
      breaks = unique(sort(c(gr_dt$start, gr_dt$end))),
      limits = xlims
    ) +
    theme_minimal() +
    labs(x = "genomic position", y = "") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(fill = NA)
    )
  
  if (!is.null(label)) {
    gr_gp <- gr_gp +
      geom_text_repel(
        aes(x = (start + end) / 2, y = id, label = get(label)), 
        vjust = -1
      )
  }
  
  gr_gp
}

#' Function to select non-overlapping hits.
#' 
#' For each peak in the input GRanges object, this function selects the groups of overlapping peaks, 
#' and for each of them then selects one top scored motif hit.
#' 
#' @param mta_gr GRanges object with motif hits, it should have "peak", "motif" and `order_col` metadata columns
#' @param reciprocal_overlap numeric, minimum reciprocal overlap between peaks to merge them
#' @param order_col character, column name in the GRanges object metadata to use for selecting top motif hit per overlapping group
#' @param order_decrease logical, if TRUE, order the hits in decreasing order (i.e. max value of `order_col` is selected)
#' @param seqnames character vector, if specified, restrict analysis only to given chromosome
#' @param verbose logical, print function memory usage at various steps
#' 
#' @return list with the following elements:
#'  select_hits, top scored motif hits per reduced group of overlapping peaks
#'  reduce_hits, reduced groups of overlapping peaks
#'  inputs_hits, input hits object
#'  
mta_reduce_motif_hits <- function(
    mta_gr, 
    reciprocal_overlap = 0.5,
    order_col = "motif_score",
    order_decrease = TRUE,
    seqnames = NULL, 
    verbose = FALSE
) {
  
  require(doParallel)
  require(foreach)
  require(pryr)
  require(pbapply)
  
  if (!is.null(seqnames)) {
    mta_gr_all <- mta_gr[seqnames(mta_gr) %in% seqnames,]
  } else {
    mta_gr_all <- mta_gr
  }
  strand(mta_gr_all) <- "*"
  
  # get unique motifs and peaks
  mts <- unique(mta_gr$motif)
  pks <- unique(mta_gr$peak)
  message(
    "Working with hits from ", 
    length(mts), 
    " motifs in ", 
    length(pks), 
    " peaks ",
    ifelse(is.null(seqnames), "", paste("on chromosome(s)", paste(seqnames, sep = ", ")))
  )
  
  # log memory usage
  if (verbose) message("Initial memory usage: ", format_memory(mem_used()))
  
  # extract clusters from hits object
  extractClustersFromSelfHits <- function(hits)
  {
    stopifnot(is(hits, "Hits"))
    stopifnot(queryLength(hits) == subjectLength(hits))
    hits <- GenomicRanges::union(hits, t(hits))
    qh <- queryHits(hits)
    sh <- subjectHits(hits)
    cid <- seq_len(queryLength(hits))  # cluster ids
    while (TRUE) {
      h <- Hits(qh, cid[sh],
                queryLength(hits), subjectLength(hits))
      cid2 <- pmin(cid, selectHits(h, "first"))
      if (identical(cid2, cid))
        break
      cid <- cid2
    }
    unname(splitAsList(seq_len(queryLength(hits)), cid))
  }
  
  # merge ranges that are "connected" (directly or indirectly)
  # via a hit (or several hits)
  mergeConnectedRanges <- function(x, hits)
  {
    stopifnot(is(x, "GenomicRanges"))
    stopifnot(is(hits, "Hits"))
    stopifnot(queryLength(hits) == subjectLength(hits))
    stopifnot(queryLength(hits) == length(x))
    clusters <- extractClustersFromSelfHits(hits)
    ans <- range(extractList(x, clusters))
    if (any(elementNROWS(ans) != 1L))
      stop(wmsg("some connected ranges are not on the same ",
                "chromosome and strand, and thus cannot be ",
                "merged"))
    ans <- unlist(ans)
    mcols(ans)$revmap <- clusters
    ans
  }
  
  # keep only hits that achieve % overlap
  hits <- findOverlaps(mta_gr_all)
  x <- mta_gr_all[queryHits(hits)]
  y <- mta_gr_all[subjectHits(hits)]
  relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
  hits <- hits[relative_overlap >= reciprocal_overlap]
  
  # merge the ranges that are connected via one or more hits in 'hits'.
  mta_gr_red <- mergeConnectedRanges(mta_gr_all, hits)
  
  # select top scored motif for each reduced group of overlapping peaks
  pboptions(type = "timer")
  mta_gr_top <- pblapply(mta_gr_red$revmap, function (is) {
    mta_gr_is <- mta_gr_all[is]
    mta_gr_is[order(mcols(mta_gr_is)[,order_col], decreasing = order_decrease), ][1]
  })
  mta_gr_top <- do.call("c", mta_gr_top)
  
  # final memory usage
  if (verbose) message("Final memory usage: ", format_memory(mem_used()))
  
  # outputs
  return(list(
    select_hits = mta_gr_top, 
    reduce_hits = mta_gr_red,
    inputs_hits = mta_gr_all
  ))
  
}