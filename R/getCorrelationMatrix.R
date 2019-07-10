#' Produce CpG correlation matrix
#'
#' This function takes input as dataframe with columns of probe data and rows are samples
#' and produce a correlation matrix
#' @param m m is a data frame with columns are probe data and rows are samples
#' @keywords mvalues betavalues
#' @export
#' @examples
#' m2beta(m)
#' @return returns a dataframe

getCorrelationMatrix <- function(chr,data_cases_controls,r_low=0.1,r_high=0.3,dist_cutoff=100000000,one_prob_per_CoRSIV=F){

  if (!requireNamespace("arrangements", quietly = TRUE)) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  library(arrangements)


  probe_loc <- probe_id #read.csv("./CoRSIV_ESS_SIV_CG_sites_clusters_hg38.csv",header = T,stringsAsFactors = F)
  if(one_prob_per_CoRSIV){
  probe_loc = probe_loc[!duplicated(probe_loc$CG),]
  probe_loc = probe_loc[!duplicated(probe_loc$cluster_id),]
  }

  probs <- as.vector(as.character(probe_loc$CG))
  selected_probes <- intersect(colnames(data_cases_controls),probs)

  data_cases_control_methy <- subset(data_cases_controls,select = selected_probes)
  data_cases_control_methy$CaseControl <- data_cases_controls$CaseControl

  nums <- unlist(lapply(data_cases_control_methy, is.numeric))
  data_cases_control_methy <- data_cases_control_methy[,nums]
  data_cases_control_methy <- data_cases_control_methy[ , colSums(is.na(data_cases_control_methy)) == 0]
  #write.csv(data_cases_control_methy,file = paste(data,".csv",sep = ""),row.names = F)
  data_control_methy <- data_cases_control_methy[data_cases_control_methy$CaseControl==0,]


  data_control_methy_trans <- as.data.frame(t(data_control_methy))
  data_control_methy_trans$CG <- rownames(data_control_methy_trans)

  data_control_methy_trans <- merge(data_control_methy_trans,probe_loc,by="CG")

  data_control_methy_trans$chr <- as.numeric(data_control_methy_trans$chr)
  data_control_methy_trans$pos <- as.numeric(data_control_methy_trans$pos)
  data_control_methy_trans <- data_control_methy_trans[with(data_control_methy_trans, order(chr, pos)),]

  probs_in_chr <- data_control_methy_trans[data_control_methy_trans$chr==chr,]
  probs_in_chr_pos <- subset(probs_in_chr,select = c(status,chr,CG,pos))

  positions <- data.frame(combinations(as.numeric(probs_in_chr_pos$pos),2,replace = F))
  positions$distance <- positions$X2 - positions$X1
  probe_combinations <- data.frame(combinations(probs_in_chr_pos$CG,2,replace = F))
  prob_position_combinations <- cbind(probe_combinations,positions)
  colnames(prob_position_combinations) <- c("CG1","CG2","pos1","pos2","dist")


  positions_all_permu <- data.frame(permutations(as.numeric(probs_in_chr_pos$pos),2,replace = T))
  positions_all_permu$distance <- abs(positions_all_permu$X2 - positions_all_permu$X1)
  distance_matrix <- matrix(positions_all_permu$distance,
                            nrow = dim(probs_in_chr_pos)[1],
                            ncol = dim(probs_in_chr_pos)[1])
  colnames(distance_matrix) <- as.character(probs_in_chr_pos$pos)
  row.names(distance_matrix) <- as.character(probs_in_chr_pos$pos)

  distance_matrix[distance_matrix < dist_cutoff] <- 0
  distance_matrix[distance_matrix >= dist_cutoff] <- 1


  probs_in_chr <- subset(probs_in_chr,select = -c(status,chr,pos))
  probs <- probs_in_chr$CG
  probs_in_chr$CG <- NULL
  probs_in_chr_df <- data.frame(t(probs_in_chr))
  colnames(probs_in_chr_df) <- probs
  probs_in_chr_df[] <- lapply(probs_in_chr_df, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  M_control<-cor(probs_in_chr_df,use="complete.obs")
  M_control_unfilterd <- M_control
  # mycor <- rcorr(as.matrix(probs_in_chr_df), type="pearson")
  # M_control_P <- mycor$P
  # M_control_R <- mycor$r
  # M_control_P[upper.tri(M_control_P,diag = T)] <-0
  # M_control_R[upper.tri(M_control_R,diag = T)] <-0
  #
  # hist(M_control_P[M_control_R>0.2])
  #
  # M_control_R <- M_control_R[upper.tri(M_control_R,diag = T)]
  # M_control_R <- M_control_R[M_control_R>0.2]

  # test <-data.frame(M_control_R[upper.tri(M_control_R,diag = F)],
  # M_control_P[upper.tri(M_control_P,diag = F)])



  control_correlation_df <- data.frame(row=rownames(M_control)[row(M_control)[upper.tri(M_control)]],
                                       col=colnames(M_control)[col(M_control)[upper.tri(M_control)]],
                                       corr_control=M_control[upper.tri(M_control)])

  #control_correlation_df <- control_correlation_df[control_correlation_df$corr_control > 0,]

  data_case_methy <- data_cases_control_methy[data_cases_control_methy$CaseControl==1,]

  data_case_methy_trans <- as.data.frame(t(data_case_methy))
  data_case_methy_trans$CG <- rownames(data_case_methy_trans)

  data_case_methy_trans <- merge(data_case_methy_trans,probe_loc,by="CG")

  data_case_methy_trans$chr <- as.numeric(data_case_methy_trans$chr)
  data_case_methy_trans$pos <- as.numeric(data_case_methy_trans$pos)
  data_case_methy_trans <- data_case_methy_trans[with(data_case_methy_trans, order(chr, pos)),]

  probs_in_chr <- data_case_methy_trans[data_case_methy_trans$chr==chr,]
  probs_in_chr <- subset(probs_in_chr,select = -c(status,chr,pos))
  probs <- probs_in_chr$CG
  probs_in_chr$CG <- NULL
  probs_in_chr_df <- data.frame(t(probs_in_chr))
  colnames(probs_in_chr_df) <- probs

  probs_in_chr_df[] <- lapply(probs_in_chr_df, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  M_case<-cor(probs_in_chr_df,use="complete.obs")
  M_case_unfilterd <- M_case
  case_correlation_df <- data.frame(row=rownames(M_case)[row(M_case)[upper.tri(M_case)]],
                                    col=colnames(M_case)[col(M_case)[upper.tri(M_case)]],
                                    corr_case=M_case[upper.tri(M_case)])

  merged_df <- merge(control_correlation_df, case_correlation_df, by.x=c("row", "col"), by.y=c("row", "col"))

  merged_df$diff <- merged_df$corr_control-merged_df$corr_case

  merged_df_pos <- merge(merged_df,prob_position_combinations,by.x=c("row", "col"), by.y=c("CG1", "CG2"))
  merged_df_pos$chr <- rep(chr,dim(merged_df_pos)[1])
  head(merged_df_pos)
  merged_df_pos$pos1 <- as.numeric(merged_df_pos$pos1)
  merged_df_pos$pos2 <- as.numeric(merged_df_pos$pos2)
  merged_df_pos <- merged_df_pos[
    with(merged_df_pos, order(pos1, pos2)),
    ]


  filtered_pairs <- !(r_low < abs(M_control) & abs(M_control) < r_high)
  M_case[filtered_pairs] <- 0
  M_control[filtered_pairs] <- 0
  M_control[upper.tri(M_control,diag = T)] <- 0
  M_case[lower.tri(M_case)] <-0
  merged <- M_control + M_case
  merged <- merged * distance_matrix

  return(list(merged,M_control_unfilterd,M_case_unfilterd,merged_df_pos))


}
