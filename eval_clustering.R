### Heper functions ###
#########################################################################################
# Function to match cluster labels with manual gating (reference standard) population 
# labels and calculate precision, recall, and F1 score
#
# Matching criterion: Hungarian algorithm
#
# Use this function for data sets with multiple populations of interest
#
# Lukas Weber, August 2016
#########################################################################################
#library(clue)
# arguments:
# - clus_algorithm: cluster labels from algorithm
# - clus_truth: true cluster labels
# (for both arguments: length = number of cells; names = cluster labels (integers))
helper_match_evaluate_multiple <- function(clus_algorithm, clus_truth) {
  
  # number of detected clusters
  n_clus <- length(table(clus_algorithm))
  
  # remove unassigned cells (NA's in clus_truth)
  unassigned <- is.na(clus_truth)
  clus_algorithm <- clus_algorithm[!unassigned]
  clus_truth <- clus_truth[!unassigned]
  if (length(clus_algorithm) != length(clus_truth)) warning("vector lengths are not equal")
  
  tbl_algorithm <- table(clus_algorithm)
  tbl_truth <- table(clus_truth)
  
  # detected clusters in rows, true populations in columns
  pr_mat <- re_mat <- F1_mat <- matrix(NA, nrow = length(tbl_algorithm), ncol = length(tbl_truth))
  
  for (i in 1:length(tbl_algorithm)) {
    for (j in 1:length(tbl_truth)) {
      i_int <- as.integer(names(tbl_algorithm))[i]  # cluster number from algorithm
      j_int <- as.integer(names(tbl_truth))[j]  # cluster number from true labels
      
      true_positives <- sum(clus_algorithm == i_int & clus_truth == j_int, na.rm = TRUE)
      detected <- sum(clus_algorithm == i_int, na.rm = TRUE)
      truth <- sum(clus_truth == j_int, na.rm = TRUE)
      
      # calculate precision, recall, and F1 score
      precision_ij <- true_positives / detected
      recall_ij <- true_positives / truth
      F1_ij <- 2 * (precision_ij * recall_ij) / (precision_ij + recall_ij)
      
      if (F1_ij == "NaN") F1_ij <- 0
      
      pr_mat[i, j] <- precision_ij
      re_mat[i, j] <- recall_ij
      F1_mat[i, j] <- F1_ij
    }
  }
  
  # put back cluster labels (note some row names may be missing due to removal of unassigned cells)
  rownames(pr_mat) <- rownames(re_mat) <- rownames(F1_mat) <- names(tbl_algorithm)
  colnames(pr_mat) <- colnames(re_mat) <- colnames(F1_mat) <- names(tbl_truth)
  
  # match labels using Hungarian algorithm applied to matrix of F1 scores (Hungarian
  # algorithm calculates an optimal one-to-one assignment)
  
  # use transpose matrix (Hungarian algorithm assumes n_rows <= n_cols)
  F1_mat_trans <- t(F1_mat)
  
  if (nrow(F1_mat_trans) <= ncol(F1_mat_trans)) {
    # if fewer (or equal no.) true populations than detected clusters, can match all true populations
    labels_matched <- clue::solve_LSAP(F1_mat_trans, maximum = TRUE)
    # use row and column names since some labels may have been removed due to unassigned cells
    labels_matched <- as.numeric(colnames(F1_mat_trans)[as.numeric(labels_matched)])
    names(labels_matched) <- rownames(F1_mat_trans)
    
  } else {
    # if fewer detected clusters than true populations, use transpose matrix and assign
    # NAs for true populations without any matching clusters
    labels_matched_flipped <- clue::solve_LSAP(F1_mat, maximum = TRUE)
    # use row and column names since some labels may have been removed due to unassigned cells
    labels_matched_flipped <- as.numeric(rownames(F1_mat_trans)[as.numeric(labels_matched_flipped)])
    names(labels_matched_flipped) <- rownames(F1_mat)
    
    labels_matched <- rep(NA, ncol(F1_mat))
    names(labels_matched) <- rownames(F1_mat_trans)
    labels_matched[as.character(labels_matched_flipped)] <- as.numeric(names(labels_matched_flipped))
  }
  
  # precision, recall, F1 score, and number of cells for each matched cluster
  pr <- re <- F1 <- n_cells_matched <- rep(NA, ncol(F1_mat))
  names(pr) <- names(re) <- names(F1) <- names(n_cells_matched) <- names(labels_matched)
  
  for (i in 1:ncol(F1_mat)) {
    # set to 0 if no matching cluster (too few detected clusters); use character names 
    # for row and column indices in case subsampling completely removes some clusters
    pr[i] <- ifelse(is.na(labels_matched[i]), 0, pr_mat[as.character(labels_matched[i]), names(labels_matched)[i]])
    re[i] <- ifelse(is.na(labels_matched[i]), 0, re_mat[as.character(labels_matched[i]), names(labels_matched)[i]])
    F1[i] <- ifelse(is.na(labels_matched[i]), 0, F1_mat[as.character(labels_matched[i]), names(labels_matched)[i]])
    
    n_cells_matched[i] <- sum(clus_algorithm == labels_matched[i], na.rm = TRUE)
  }
  
  # means across populations
  mean_pr <- mean(pr)
  mean_re <- mean(re)
  mean_F1 <- mean(F1)
  median_F1 <- median(F1)
  
  return(list(n_clus = n_clus, 
              pr = pr, 
              re = re, 
              F1 = F1, 
              labels_matched = labels_matched, 
              n_cells_matched = n_cells_matched, 
              mean_pr = mean_pr, 
              mean_re = mean_re, 
              mean_F1 = mean_F1,
              median_F1 = median_F1,
              Fmeasure = Fmeasure(clus_algorithm, clus_truth)))
}


#########################################################################################
# Function to match cluster labels with manual gating (reference standard) population 
# labels and calculate precision, recall, and F1 score
#
# Matching criterion: maximum F1 score
#
# Use this function for data sets with a single (e.g. rare) population of interest
#
# Lukas Weber, August 2016
#########################################################################################
# arguments:
# - clus_algorithm: cluster labels from algorithm
# - clus_truth: true cluster labels (1 = rare cluster of interest, 0 = all others)
# (for both arguments: length = number of cells; names = cluster labels (integers))
helper_match_evaluate_single <- function(clus_algorithm, clus_truth) {
  
  # number of detected clusters
  n_clus <- length(table(clus_algorithm))
  
  tbl_algorithm <- table(clus_algorithm)
  tbl_truth <- table(clus_truth)
  
  pr_mat <- re_mat <- F1_mat <- matrix(NA, nrow = length(tbl_algorithm), ncol = 1)
  
  for (i in 1:length(tbl_algorithm)) {
    i_int <- as.integer(names(tbl_algorithm))[i]  # cluster number from algorithm
    
    j_int <- 1  # true cluster number of the rare population of interest
    
    true_positives <- sum(clus_algorithm == i_int & clus_truth == j_int, na.rm = TRUE)
    detected <- sum(clus_algorithm == i_int, na.rm = TRUE)
    truth <- sum(clus_truth == j_int, na.rm = TRUE)
    
    # calculate precision, recall, and F1 score
    precision_ij <- true_positives / detected
    recall_ij <- true_positives / truth
    F1_ij <- 2 * (precision_ij * recall_ij) / (precision_ij + recall_ij)
    
    if (F1_ij == "NaN") F1_ij <- 0
    
    pr_mat[i, j_int] <- precision_ij
    re_mat[i, j_int] <- recall_ij
    F1_mat[i, j_int] <- F1_ij
  }
  
  # put back cluster labels (note some row names may be missing due to removal of unassigned cells)
  rownames(pr_mat) <- rownames(re_mat) <- rownames(F1_mat) <- names(tbl_algorithm)
  colnames(pr_mat) <- colnames(re_mat) <- colnames(F1_mat) <- "1"  # one column only
  
  # match label (single cluster only) using highest F1 score
  # use row names since some labels may have been removed due to unassigned cells
  labels_matched <- as.numeric(rownames(F1_mat)[apply(F1_mat, 2, which.max)])
  names(labels_matched) <- "1"  # one column only
  
  # precision, recall, F1 score, and number of cells for single matched cluster
  # use character names for row and column indices in case subsampling completely removes some clusters
  pr <- pr_mat[as.character(labels_matched), "1"]
  re <- re_mat[as.character(labels_matched), "1"]
  F1 <- F1_mat[as.character(labels_matched), "1"]
  
  n_cells_matched <- sum(clus_algorithm == labels_matched, na.rm = TRUE)
  
  return(list(n_clus = n_clus, 
              pr = pr, 
              re = re, 
              F1 = F1, 
              labels_matched = labels_matched, 
              n_cells_matched = n_cells_matched))
}

### R version reference FlowCAP f-measure sum
Fmeasure<-function(pred, ref){
  
  K = unique(pred)
  K = sort(K)
  C = unique(ref)
  C = sort(C)
  m = length(K)
  n = length(C)
  
  M<-matrix(NA , n, m)
  Pr<-matrix(NA, n, m)
  Re<-matrix(NA, n, m)
  Fmat<-matrix(NA, n, m)
  
  C_card<-vector("double", n)
  K_card<-vector("double", m)
  
  for(i in seq(n)){
    C_card[i] = sum(ref == C[i])
    for(j in seq(m)){
      K_card[j] = sum(pred == K[j])
      M[i,j] = sum((ref==C[i]) & (pred==K[j]))
      Pr[i,j] = M[i,j]/K_card[j]
      Re[i,j] = M[i,j]/C_card[i]
      if((Pr[i,j] + Re[i,j]) == 0.0){
        Fmat[i,j] = 0
      }else{
        Fmat[i,j] = 2.0*Pr[i,j]*Re[i,j]/(Pr[i,j] + Re[i,j])
      }
    }
  }
  
  C_card_sum = sum(C_card)
  Ffinal<-vector("double", n)
  Fsum<-vector("double", n)
  
  for(i in seq(n)){
    Ffinal[i] = max(Fmat[i,])
    Fsum[i] = Ffinal[i]*C_card[i]/C_card_sum
  }
  Ftotal = sum(Fsum)
  
  return(Ftotal)
}

eval_clustering<-function(dataset = "Nilsson_rare", method = "umap_dbs", all_data = T, repeats = 1L, ret_model = F, use_cuml = F, seed = NULL) {
  dataset<-match.arg(dataset, c("Levine_32dim","Levine_13dim","Samusik_01","Samusik_all","Nilsson_rare","Mosmann_rare", "MO_LST_all", "all"), several.ok = T)
  if ("all" %in% dataset) dataset<-c("Levine_32dim","Levine_13dim","Samusik_01","Samusik_all","Nilsson_rare","Mosmann_rare", "MO_LST_all")
  method<-match.arg(method, c("umap_dbs","flowSOM", "all"))
  if (method == "all") method<-c("umap_dbs","flowSOM")
  if (!is.numeric(repeats) || length(repeats) != 1 || repeats != round(repeats) || repeats < 0) stop("repeats must be integer > 0")
  if (is.null(seed)) fast_sgd = T else fast_sgd = F
  metres<-list()
  if ("umap_dbs" %in% method) {
    result<-list()
    for (d in dataset) {
      repres<-list()
      for (i in 1:repeats) {
        ff<-ffs[[d]]
        label<-flowCore::exprs(ff)[,marker_label[[d]]]
        if (all_data){
          data<-flowCore::exprs(ff)[,marker_cols[[d]]]
        } else data<-flowCore::exprs(ff)[!is.nan(label),marker_cols[[d]]]
        t1<-Sys.time()
        mymodel<-sumapc(data = data,
                        maxevts = umapdbs_options[[d]]$maxevts,
                        maxlvl = umapdbs_options[[d]]$maxlvl,
                        minpts = umapdbs_options[[d]]$minpts,
                        clust_options = umapdbs_options[[d]]$clust_options,
                        multi_thread = T,
                        fast_sgd = fast_sgd,
                        myqueue = NULL,
                        ret_model = T,
                        verbose = F,
                        use_cuml = use_cuml,
                        seed = seed)
        t2<-Sys.time()
        if (all_data) {
          clus_algorithm<-mymodel$cluster[!is.nan(label)]
        } else {
          clus_algorithm<-mymodel$cluster
        }
        clus_truth<-label[!is.nan(label)]
        names(clus_algorithm)<-clus_algorithm
        names(clus_truth)<-clus_truth
        if (is_rare[[d]]) {
          res<-helper_match_evaluate_single(clus_algorithm, clus_truth)
        } else res<-helper_match_evaluate_multiple(clus_algorithm, clus_truth)
        if (ret_model) model<-mymodel else model<-"none"
        repres[[paste0("repeat",i)]]<-list(res = res, model = model, runtime = t2-t1, clus_algorithm = clus_algorithm, clus_truth = clus_truth)
      }
      result[[d]]<-repres
    }
    metres[["umap_dbs"]]<-result
  }
  if ("flowSOM" %in% method) {
    result<-list()
    for (d in dataset) {
      repres<-list()
      for (i in 1:repeats) {
        ff<-ffs[[d]]
        label<-flowCore::exprs(ff)[,marker_label[[d]]]
        if (all_data){
          fft<-ff[,marker_cols[[d]]]
        } else fft<-ff[!is.nan(label),marker_cols[[d]]]
        gridsize<-10
        if (d == "Mosmann_rare") gridsize = 20
        t1<-Sys.time()
        flowSOM.res<-FlowSOM::ReadInput(fft, compensate=FALSE, transform=FALSE, scale=FALSE)
        flowSOM.res<-FlowSOM::BuildSOM(flowSOM.res, xdim=gridsize, ydim=gridsize, silent = TRUE)
        flowSOM.res<-FlowSOM::BuildMST(flowSOM.res, silent = TRUE, tSNE = FALSE)
        metacl<-FlowSOM::metaClustering_consensus(flowSOM.res$map$codes, k=40, seed = seed)
        clusters<-metacl[flowSOM.res$map$mapping[,1]] 
        t2<-Sys.time()
        if (all_data) {
          clus_algorithm<-clusters[!is.nan(label)]
        } else {
          clus_algorithm<-clusters
        }
        clus_truth<-label[!is.nan(label)]
        names(clus_algorithm)<-clus_algorithm
        names(clus_truth)<-clus_truth
        if (is_rare[[d]]) {
          res<-helper_match_evaluate_single(clus_algorithm, clus_truth)
        } else res<-helper_match_evaluate_multiple(clus_algorithm, clus_truth)
        if (ret_model) model<-flowSOM.res else model<-"none"
        repres[[paste0("repeat",i)]]<-list(res = res, model = model, runtime = t2-t1, clus_algorithm = clus_algorithm, clus_truth = clus_truth)
      }
      result[[d]]<-repres
    }
    metres[["flowSOM"]]<-result
  }
  return(metres)
}

tabulate_results<-function(res) {
  num_met<-length(res)
  mets<-names(res)
  num_ds<-length(res[[1]])
  ds<-names(res[[1]])
  restab<-matrix(NA_real_, nrow = num_ds, ncol = 4*num_met)
  cnames<-c()
  for (m in mets) {
    cnames<-c(cnames, paste(m, c("meanF1", "medianF1", "Fmeasure", "Runtime"), sep = "_"))
  }
  rownames(restab)<-ds
  colnames(restab)<-cnames
  for (d in ds){
    for (m in mets) {
      if (is_rare[[d]]) {
        restab[d, paste(m, "meanF1", sep="_")]<-res[[m]][[d]]$repeat1$res$F1
      } else {
        restab[d, paste(m, "meanF1", sep="_")]<-res[[m]][[d]]$repeat1$res$mean_F1
        restab[d, paste(m, "medianF1", sep="_")]<-res[[m]][[d]]$repeat1$res$median_F1
        restab[d, paste(m, "Fmeasure", sep="_")]<-res[[m]][[d]]$repeat1$res$Fmeasure
      }
      restab[d, paste(m, "Runtime", sep="_")]<-as.numeric(res[[m]][[d]]$repeat1$runtime, units = "secs")
    }
  }
  return(restab)
}
### End Helper funcitons


setwd("../../aval_clustering")

library(future)
future::plan(list(
  future::tweak(future::multisession, workers = 2),
  future::tweak(future::multisession, workers = future::availableCores() %/% 2)
))
datasets <- c(
  "Levine_32dim",
  "Levine_13dim",
  "Samusik_01",
  "Samusik_all",
  "Nilsson_rare",
  "Mosmann_rare",
  "MO_LST_all"
)
marker_cols <- list(
  Levine_32dim = 5:36, 
  Levine_13dim = 1:13, 
  Samusik_01   = 9:47, 
  Samusik_all  = 9:47, 
  Nilsson_rare = c(5:7, 9:18), 
  Mosmann_rare = c(7:9, 11:21),
  MO_LST_all = 13:23
)
marker_label<-list(
  Levine_32dim = "label", 
  Levine_13dim = "label", 
  Samusik_01   = "label", 
  Samusik_all  = "label", 
  Nilsson_rare = "label", 
  Mosmann_rare = "label",
  MO_LST_all = "cluster1_root"
)
files <- list(
  Levine_32dim = "Levine_32dim.fcs", 
  Levine_13dim = "Levine_13dim.fcs", 
  Samusik_01   = "Samusik_01.fcs", 
  Samusik_all  = "Samusik_all.fcs", 
  Nilsson_rare = "Nilsson_rare.fcs", 
  Mosmann_rare = "Mosmann_rare.fcs",
  MO_LST_all = "MO_LST_transformed_gated - all.fcs"
)
ffs<-lapply(files, function(x) flowCore::read.FCS(x, transformation = F, truncate_max_range = F))
is_rare <- list(
  Levine_32dim = F, 
  Levine_13dim = F, 
  Samusik_01   = F, 
  Samusik_all  = F, 
  Nilsson_rare = T, 
  Mosmann_rare = T,
  MO_LST_all = F
)
umapdbs_options <- list(
  Levine_32dim = list(
    minpts=100,
    maxevts=100000,
    maxlvl=3,
    clust_options=list(method = "sdbscan", mineps = 1, mindens = 0.05, bw = .05, nbins = 4, mvpratio = 0.5)
  ), 
  Levine_13dim = list(
    minpts=100,
    maxevts=10000,
    maxlvl=3,
    clust_options=list(method = "sdbscan", mineps = 1, mindens = 0.05, bw = .05, nbins = 4, mvpratio = 0.6)
  ), 
  Samusik_01   = list(
    minpts=100,
    maxevts=10000,
    maxlvl=3,
    clust_options=list(method = "sdbscan", mineps = 1, mindens = 0.05, bw = .05, nbins = 4, mvpratio = 0.5)
  ), 
  Samusik_all  = list(
    minpts=100,
    maxevts=10000,
    maxlvl=3,
    clust_options=list(method = "sdbscan", mineps = 1, mindens = 0.05, bw = .05, nbins = 5, mvpratio = 0.4)
  ), 
  Nilsson_rare = list(
    minpts=50,
    maxevts=10000,
    maxlvl=3,
    clust_options=list(method = "sdbscan", mineps = 1, mindens = 0.05, bw = .05, nbins = 4, mvpratio = 0.5)
  ), 
  Mosmann_rare = list(
    minpts=100,
    maxevts=100000,
    maxlvl=3,
    clust_options=list(method = "sdbscan", mineps = 1, mindens = 0.05, bw = .05, nbins = 5, mvpratio = 0)
  ),
  MO_LST_all = list(
    minpts=100,
    maxevts=100000,
    maxlvl=3,
    clust_options=list(method = "sdbscan", mineps = 1, mindens = 0.05, bw = .1, nbins = 15, mvpratio = 0)
  )
)



# evalclust<-eval_clustering(dataset = "Levine_13dim", method = "umap_dbs", all_data = T, repeats = 1, seed = 1000)
# evalclust2<-eval_clustering(dataset = "Levine_13dim", method = "flowSOM", all_data = T, repeats = 1, seed = 1000)
# evalclust3<-eval_clustering(dataset = "Mosmann_rare", method = "umap_dbs", all_data = F)
# evalclust4<-eval_clustering(dataset = "Mosmann_rare", method = "flowSOM", all_data = F)
# evalclust5<-eval_clustering(dataset = "Mosmann_rare", method = "all", all_data = T)
# 
# evalclust<-eval_clustering_cuml(dataset = "Levine_32dim", method = "umap_dbs", all_data = T, repeats = 5)
# evalclust2<-eval_clustering_cuml(dataset = "Levine_32dim", method = "flowSOM", all_data = T, repeats = 1, seed = 1000)

evalclust_final<-eval_clustering(dataset = "all", method = "all", all_data = T, repeats = 1, seed = 1000)
# Run cuml version on linux environment (or Windows 11 + WSL2) with RAPIDS/cuml installed and working (requires a NVidia GPU with compute capability 6.0+)
# Try cuml<-reticulate::import("cuml") to know if the environment is ok.
evalclust_final2<-eval_clustering(dataset = "all", method = "umap_dbs", all_data = T, repeats = 1, use_cuml = T, seed = 1000)
names(evalclust_final2)<-"umap_dbs_cuml"
evalclust_all<-(c(evalclust_final, evalclust_final2))
res_all<-as.data.frame(tabulate_results(evalclust_all))
write.csv2(res_all[,c(1,5,9,2,6,10,3,7,11,4,8,12)], "final_results.csv")

res_fastsgd<-as.data.frame(tabulate_results(evalclust_fastsgd))
write.csv2(res_fastsgd[,c(1,5,2,6,3,7,4,8)], "fastsgd_results.csv")
