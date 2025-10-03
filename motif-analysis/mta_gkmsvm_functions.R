#' Modfied genNullSeqs to generate null sequences matched by GC content, given a set of positive regions
#' 
#' Generates null sequences (negative set) with matching GC content as the input bed file for positive set regions. Modified from `gkmSVM::genNullSeqs`. Works without `BSgenome`.
#' 
#' @param inputBedFN positive set regions
#' @param genome_object either a path to a genome file (`.fasta`), or a preloaded `DNAStringSet` object (from `Biostrings::readDNAStringSet`)
#' @param outputBedFN utput file name for the null sequences genomic regions. Default='negSet.bed'
#' @param outputPosFastaFN output file name for the positive set sequences. Default='posSet.fa'
#' @param outputNegFastaFN output file name for the negative set sequences. Default='negSet.fa'
#' @param xfold controls the desired number of sequences in the negative set. Default=1 (same number as in positive set)
#' @param GC_match_tol tolerance for difference in GC content. Default=0.02
#' @param length_match_tol tolerance for difference in relative sequence length. Default=0.02
#' @param batchsize number of candidate random sequences tested in each trial. Default=5000
#' @param nMaxTrials maximum number of trials. Default=20.
#' 
gkm_mod_genNullSeqs = function (
	inputBedFN, genome_object, select_background_from_inputBedFN = NULL, outputBedFN = "negSet.bed", 
    outputPosFastaFN = "posSet.fa", outputNegFastaFN = "negSet.fa", 
    xfold = 1, GC_match_tol = 0.02, 
    length_match_tol = 0.02, batchsize = 5000, nMaxTrials = 20) {
    
	requireNamespace("GenomicRanges", quietly = TRUE)
	requireNamespace("Bsgenome", quietly = TRUE)      # this is not used, but loading it allows using the Biostrings::getSeq function with DNAstringset objects
	requireNamespace("rtracklayer", quietly = TRUE)
	requireNamespace("BiocGenerics", quietly = TRUE)
	requireNamespace("Biostrings", quietly = TRUE)
	requireNamespace("IRanges", quietly = TRUE)
	requireNamespace("S4Vectors", quietly = TRUE)

	# load genome object
	if ("character" %in% class(genome_object)) {
		genome = Biostrings::readDNAStringSet(genome_object, format = "fasta")
	} else if ("DNAStringSet" %in% class(genome_object)) {
		genome = genome_object
	} else {
		stop(sprintf("`genome_object` %s has to be either a DNAStringSet or a path to a fasta file", DNAStringSet))
	}

	seqnams = names(genome)
	chrlens = width(genome)
	names(chrlens) = seqnams
	chrpos = cumsum(as.numeric(chrlens))
	pmax = max(chrpos)
	chrpos = c(chrpos, 1e+12)
	chrpos0 = c(0, chrpos)
	ichrA = as.character(names(chrlens))
	
	# functions
	getichrpos = function(ipos) {
		j = order(ipos)
		ipos = sort(ipos)
		ci = 1
		res = rep(NA, length(ipos))
		for (i in 1:length(ipos)) {
			while (ipos[i] > chrpos[ci]) {
				ci = ci + 1
			}
			res[j[i]] = ci
		}
		return(res)
	}
	generateRandomGenSeqs = function(seqlens) {
		rpos = sample(pmax, length(seqlens), replace = TRUE)
		ichr1 = getichrpos(rpos)
		ichr2 = getichrpos(rpos + seqlens)
		jj = which(ichr1 != ichr2)
		while (length(jj) > 0) {
			rpos[jj] = sample(pmax, length(jj), replace = TRUE)
			ichr1 = getichrpos(rpos)
			ichr2 = getichrpos(rpos + seqlens)
			jj = which(ichr1 != ichr2)
		}
		chr = ichrA[ichr1]
		start = rpos - chrpos0[ichr1]
		names <- chr
		ranges <- IRanges::IRanges(start = start, width = seqlens)
		strand <- BiocGenerics::strand(sample(c("+", "-"), length(names), replace = TRUE))
		gr <- GenomicRanges::GRanges(seqnames = names, ranges = ranges, strand = strand)
	}
	
	# import foreground bed
	inBed = mta_bed_to_granges(inputBedFN)
	# restrict bed regions to genome limits
	seqlengths(inBed) <- as.integer( chrlens [ levels(seqnames(inBed)) ] )
	inBed = GenomicRanges::trim(inBed)
	inBed = unique(inBed)
	message(sprintf("gkmsvm | n=%i unique regions in foreground", length(inBed)))
	# bed as dataframe
	inbed = GenomicRanges::as.data.frame(inBed)
	
	jj = which(is.na(match(as.character(inbed$seqnames), as.character(seqnams))))
	if (length(jj) > 0) {
		cat(paste("ERROR: Chromosome name not recognized for", length(jj), "sequences.\n"))
		cat(unique(as.character(inbed$seqnames[jj])))
		return(NULL)
	}
	jj = which(inbed$end > width(genome)[as.character(inbed$seqnames)])
	if (length(jj) > 0) {
		cat("ERROR: Region outside chromosome. (Check the genome version) \n")
		print(inbed[jj, ])
		return(NULL)
	}
	
	# functions to match GC content and find appropriate background
	gcContent <- function(seqs) {
		alf <- Biostrings::alphabetFrequency(seqs, as.prob = TRUE)
		gc = rowSums(alf[, c("G", "C"), drop = FALSE])
	}
	matchSeqs = function(gc1, gc2, len1, len2, 
		gc_th = 0.02, len_th = 0.02, rpt_th = 0.02) {
		len_th = len_th * len1
		i1 = order(gc1)
		i2 = order(gc2)
		gc1 = gc1[i1]
		gc2 = gc2[i2]
		len1 = len1[i1]
		len2 = len2[i2]
		gc2 = c(gc2, 1e+10)
		len_th = len_th[i1]
		m2 = 1
		N = length(i1)
		N2 = length(i2)
		mtc1 = rep(NA, N)
		mtc2 = rep(0, length(i2))
		for (i in 1:N) {
			gc1i = gc1[i]
			len1i = len1[i]
			len_thi = len_th[i]
			while (gc1i - gc2[m2] > gc_th) {
				m2 = m2 + 1
			}
			if (m2 <= N2) {
				m2b = m2
				while (gc2[m2b] - gc1i <= gc_th) {
				if ((mtc2[m2b] == 0) & (abs(len1i - len2[m2b]) <= len_thi) ) {
					mtc2[m2b] = i
					mtc1[i] = m2b
					if (m2b == m2) {
					m2 = m2 + 1
					}
					break
				}
				m2b = m2b + 1
				}
			}
			else {
				break
			}
		}
		mtc1 = i2[mtc1]
		res = rep(NA, N)
		res[i1] = mtc1
		return(res)
	}

	# prepare foreground sequences
	inSeqs = Biostrings::getSeq(genome, inBed)
	seqlens = inbed$width
	
	# foreground gc content
	inGC = gcContent(inSeqs)
	
	# mask regions that we don't want to consider for the background? (optional)
	if (!is.null(select_background_from_inputBedFN)) {
		bg_to_use_r = mta_bed_to_granges(select_background_from_inputBedFN)
		seqlengths(bg_to_use_r) <- as.integer( chrlens [ levels(seqnames(bg_to_use_r)) ] )
		bg_to_use_r = GenomicRanges::trim(bg_to_use_r)
		bg_to_use_r = unique(bg_to_use_r)
		message(sprintf("gkmsvm | n=%i unique regions to inspect for background", length(bg_to_use_r)))
	}
	
	# prepare background output object
	nout = round(nrow(inbed) * xfold)
	outbed = matrix(ncol = ncol(inbed), nrow = nout)
	outSeq = rep(inSeqs, length = nout)
	colnames(outbed) = colnames(inbed)
	unmatched = 1:length(outSeq)
	desGC = rep(inGC, length = nout)
	desLens = rep(seqlens, length = nout)
	
	# find appropriate background in N iterations
	for (iter in 1:nMaxTrials) {
		if (length(unmatched) > 0) {
			# message(sprintf("gkmsvm | find GC-matching background, iter %i/%i", iter, nMaxTrials))
			
			if (!is.null(select_background_from_inputBedFN)) {
				if (iter == 1) { bg_to_use_r_keep = bg_to_use_r }
				sample_ixs = sample(1:length(bg_to_use_r_keep), min(batchsize, length(bg_to_use_r)))
				# bg_to_use_r_keep = bg_to_use_r_keep [ -sample_ixs ]
				rndBed = bg_to_use_r_keep [ sample_ixs ]
				rndBed$score = NULL
				rndBed$name  = NULL
			} else {
				rndBed = generateRandomGenSeqs(seqlens = rep(desLens[unmatched], length.out = batchsize))
			}
			
			rndbed = GenomicRanges::as.data.frame(rndBed)
			rndSeqs = Biostrings::getSeq(genome, rndBed)
			rndGC = gcContent(rndSeqs)
			rndRpt = rndBed
			desRpt = rndBed
			mtc = matchSeqs(desGC[unmatched], rndGC, desLens[unmatched], BiocGenerics::width(rndBed), gc_th = GC_match_tol, len_th = length_match_tol, rpt_th = repeat_match_tol)
			jj = which(!is.na(mtc))
			if (length(jj) > 0) {
				outbed[unmatched[jj], 1:5] = as.matrix(rndbed[mtc[jj], ])
				outSeq[unmatched[jj], ] = rndSeqs[mtc[jj], ]
				unmatched = unmatched[-jj]
			}
		}
	}
	
	# prepare output
	if (length(unmatched) > 0) {
		outbed = outbed[-unmatched, ]
		outbed = unique(outbed)
		outSeq = outSeq[-unmatched, ]
	}
	message(sprintf("gkmsvm | n=%i GC-matched background regions found", length(outSeq)))
	
	# write output
	outbed = gsub(" ", "", outbed)
	utils::write.table(as.matrix(outbed[, 1:3]), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, file = outputBedFN)
	if (requireNamespace("seqinr", quietly = TRUE)) {
		outseqnams = paste(outbed[, 1], outbed[, 2], outbed[, 3], "neg", 1:nrow(outbed), sep = "_")
		seqinr::write.fasta(sequences = sapply(as.character(outSeq), strsplit, ""), names = outseqnams, file.out = outputNegFastaFN)
		inseqnams = paste(as.character(inbed[, 1]), inbed[, 2], inbed[, 3], "pos", 1:nrow(inbed), sep = "_")
		seqinr::write.fasta(sequences = sapply(as.character(inSeqs), strsplit, ""), names = inseqnams, file.out = outputPosFastaFN)
	}
	
}
