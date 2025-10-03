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
  for(i in 1:nrow(minM)) {
    for(j in 1:ncol(minM)) {
      minM[i,j] <- min(modlen1[i], modlen2[j])
    }
  }
  
  outM <- intM/minM
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
  for(i in 1:ncol(M1)) {
    unionM[i,] <- colSums((M1[,i] + M2) > 0)
  }
  # jaccard
  outM <- intersectM/unionM
  #outM[is.na(outM)] <- 0
  rownames(outM) <- colnames(M1)
  colnames(outM) <- colnames(M2)
  outM
}

#' Similarity index-like metric
#' 
si <- function(M1,M2) {
  Gt <- overlap(M1,M2)
  Ga <- colSums(M1)
  Gb <- colSums(M2)
  si <- 1 - sqrt((1-Gt/Gb)*(1-Gt/Ga))
  si[is.na(si)] <- 0
  return(si)
}

#' Jensen-Shannon Divergence
#library(philentropy)
#library(gtools)
#mc_umifrac <- apply(mc_umifrac, 2, function(x) x/1000)
#jsd <- JSD(as.matrix(t(mc_umifrac[var_genes,])),unit="log2")
#colnames(jsd) <- rownames(jsd) <- colnames(mc_umifrac)
#jsd_dist <- as.matrix(1-sqrt(jsd))

#' KLD between columns (metacells, cell types) of matrix (footprint, genes in rows)
calcKLD=function(mc_fp) {
  nc=ncol(mc_fp)
  kld_mat=matrix(NA,nrow=nc,ncol=nc)
  kld_list=vector("list",nc)
  for (i in 1:nc) {
    message(sprintf("%s of %s",i,nc))
    for (j in 1:nc) {
      if (!i>j) {
        kld=LaplacesDemon::KLD(px=mc_fp[,i],py=mc_fp[,j])
        kld_mat[i,j]=kld$sum.KLD.px.py
        kld_mat[j,i]=kld$sum.KLD.py.px
        kld_list[[i]][[j]]=kld$KLD.px.py
        kld_list[[j]][[i]]=kld$KLD.py.px
      }
    }
  }
  colnames(kld_mat)=colnames(mc_fp);rownames(kld_mat)=colnames(mc_fp)
  list(mat=kld_mat,probs=kld_list)
}


