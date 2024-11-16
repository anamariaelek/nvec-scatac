setwd("/users/asebe/aelek/proj/scATAC_nvec_v2/Nematostella_scATAC")

library(data.table)

messagetimed <- function(...) {
  message(sprintf("%s %s", Sys.time(), sprintf(...)))
}

count_co_occurrences <- function(dt, min_ovl = NULL) {
  
  # Convert to data.table
  setDT(dt)
  
  # If we don't want to check overlap between pairs of motifs
  # we can use matrix operations to count co-occurrences 
  # (this is much faster!)
  if (is.null(min_ovl)) {

    # Count the number of motif hits in each peak
    mta_hits_dc <- dcast.data.table(dt, peak ~ motif, fun.aggregate = length)
    mta_hits <- as.matrix(mta_hits_dc[, -1])
    rownames(mta_hits) <- mta_hits_dc[[1]]
    
    # Convert to a binary indicator matrix (1 if hit count > 0, otherwise 0)
    mta_hits_ind <- (mta_hits > 0) * 1
    
    # Count co-occurrences of all pairs of motifs using matrix multiplication
    count_matrix <- crossprod(mta_hits_ind)
    
    # Set diagonal to the number of motifs with at least 2 hits
    diag(count_matrix) <- colSums(mta_hits >= 2)
    
    # Return count matrix with motif names
    rownames(count_matrix) <- colnames(count_matrix) <- colnames(mta_hits_ind)
  
  # If we want to check for overlap between every pair of motifs
  # we have to loop through each peak and count pairs of motifs
  } else {

    # Function to check if motifs don't overlap
    motifs_nonoverlap <- function(end1, start2) {
      (end1 - start2) < min_ovl
    }

    # Initialize co-occurrence matrix
    motifs <- unique(dt$motif)
    co_occurrence_matrix <- matrix(
      0,
      nrow = length(motifs), ncol = length(motifs), 
      dimnames = list(motifs, motifs)
    )
    
    # Group data by peak
    grouped_df <- split(dt, by = "peak")

    # Save results to data table
    ds <- data.table()
    
    for (peak in names(grouped_df)) {
      
      message(sprintf(
        "%s Processing %s (%d/%d)", 
        Sys.time(), peak, which(names(grouped_df) == peak), length(names(grouped_df))
      ))
      
      group_df <- grouped_df[[peak]]
      unique_motifs <- unique(group_df$motif)
      ms <- CJ(m1 = unique_motifs, m2 = unique_motifs)[m1 <= m2] 
      
      for (pair in seq_len(nrow(ms))) {

        m1 <- ms$m1[pair]
        m2 <- ms$m2[pair]
        
        # different motifs:
        if (m1 != m2) {

          # filter instances and use vectorized overlap check
          motif1 <- group_df[motif == m1]
          motif2 <- group_df[motif == m2]

          if (nrow(motif1) > 0 & nrow(motif2) > 0) {
            # sort all motif pair instances so that motif1 is always the one with lower start
            m <- CJ(i = seq_len(nrow(motif1)), j = seq_len(nrow(motif2)))
            s <- data.table(
              m1 = motif1$motif_occurrence, m2 = motif2$motif_occurrence,
              m1_start = motif1$start[m$i],
              m2_start = motif2$start[m$j],
              m1_end = motif1$end[m$i],
              m2_end = motif2$end[m$j]
            )
            s[m1_start <= m2_start, ':='(start1 = m1_start, start2 = m2_start, end1 = m1_end, end2 = m2_end)]
            s[m2_start < m1_start, ':='(start1 = m2_start, start2 = m1_start, end1 = m2_end, end2 = m1_end)]
            s[, novl := motifs_nonoverlap(end1, start2)]
            s <- s[, .(novl = all(novl)), .(m1, m2)]
            novl <- any(s$novl)
          } else {
            novl <- FALSE
          }

        # same motifs
        } else {

          motif <- group_df[motif == m1]

          if (nrow(motif) > 1) {
            # sort all motif pair instances so that motif1 is always the one with lower start
            m <- CJ(i = seq_len(nrow(motif)), j = seq_len(nrow(motif)))
            m <- m[i != j]
            s <- data.table(
              m1 = motif$motif_occurrence, m2 = motif$motif_occurrence,
              m1_start = motif$start[m$i],
              m2_start = motif$start[m$j],
              m1_end = motif$end[m$i],
              m2_end = motif$end[m$j]
            )
            s[m1_start <= m2_start, ':='(start1 = m1_start, start2 = m2_start, end1 = m1_end, end2 = m2_end)]
            s[m2_start < m1_start, ':='(start1 = m2_start, start2 = m1_start, end1 = m2_end, end2 = m1_end)]
            s[, novl := motifs_nonoverlap(end1, start2)]
            s <- s[, .(novl = all(novl)), .(m1, m2)]
            novl <- any(s$novl)
          } else {
            novl <- FALSE
          }

        }

        # update matrix
        ms[pair, nonoverlapping := novl]
        if (novl == TRUE) {
          co_occurrence_matrix[m1, m2] <- co_occurrence_matrix[m1, m2] + 1
          co_occurrence_matrix[m2, m1] <- co_occurrence_matrix[m2, m1] + 1
        }

      }
      
      ms[, peak := peak]
      ds <- rbindlist(list(ds, ms), use.names = TRUE, fill = TRUE)
    }
    
  }

  return(list(co_occurrence_matrix, ds))
}


# load data
messagetimed("Loading data")
syn_dir <- "results/Syntax"
arc_id <- "PPM-PCC-0.8-IC0.5-5bp"
q <- 0.95
mta_novl_dt <- readRDS(
  file.path(syn_dir, sprintf("motif-hits-novl-%s-mona-q-%s.rds", arc_id, q))
)

# first arg passed to this script is a bin
bin <- commandArgs(trailingOnly = TRUE)[1]
messagetimed(sprintf("Processing bin %s", bin))

# where to save results
out_fn <- file.path(syn_dir, "bins", sprintf("co-occurrences-%s.rds", bin))
messagetimed(sprintf("Results will be saved to %s", out_fn))

# filter data to bin
bin_pks <- readLines(file.path(syn_dir, "bins", sprintf("peaks-%s.txt", bin)))
mta_novl_dt <- mta_novl_dt[peak %in% bin_pks]

# calculation
mta_novl_counts <- count_co_occurrences(mta_novl_dt, min_ovl = 3)

# save results
messagetimed(sprintf("Saving results to %s", out_fn))
saveRDS(mta_novl_counts, out_fn)

# log off
messagetimed("Done")