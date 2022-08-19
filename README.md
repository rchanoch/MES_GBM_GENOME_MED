# MES_GBM_GENOME_MEDICINE
"­­Elucidating the diversity of malignant mesenchymal states in glioblastoma by integrative analysis"

##### MES_sigs_makeup #### returns a dataframe with the percentage of genes from each core or function-specific MES signature that are included in each of the published MES signatures 
# MES_Sigs_List = Table S1 of manuscript
# MES_functional_list = Table S2 of manuscript

MES_sigs_makeup<- function(MES_Sigs_List,MES_functional_list)
{
intersect_programs<- lapply(MES_Sigs_List, function(x){lapply(MES_functional_list, function(y){length(intersect(y,x))/length(x) })  })

intersect_programs<- as.data.frame(do.call(cbind,intersect_programs))
intersect_programs<- apply(intersect_programs, 1, unlist)

sigs_makeup_df <- intersect_programs[rownames(intersect_programs)[order(intersect_programs[,"MES_CORE_GENES"],decreasing=T)],]

return(sigs_makeup_df)
}

##### get_score_dist #####  returns the distribution of scores + shuffled scores for a given signature
# tpm_matrix for each study downloaded from original publications and subset to include only malignant cells and 2000+ detected genes
# cell_type_df = metadata including cell name and patient of origin
# signature = name of signature from MES_functional_list in question


get_score_dist<- function(tpm_matrix, cell_type_df, signature=NULL, threshold=.99)
{
  log_expression_matrix <- log2(1+ tpm_matrix/10)
  
  # matrices were filtered to only include samples that had 50+ malignant cells 
  
  cells_per_tumor_in_study<- table(cell_type_df$Sample)
  cells_per_tumor_in_study<- cells_per_tumor_in_study[cells_per_tumor_in_study>50]
  log_expression_matrix<- log_expression_matrix[,cell_type_df$Cell[cell_type_df$Sample %in% names(cells_per_tumor_in_study)]]
  
  
  # matrices were shuffled to generate shuffled distribution 
  shuffled_log_tpm_matrix <- t(apply(log_expression_matrix, 1, function(x){sample(x,length(x))}))
  
  score_matrix<- scalop::sigScores(m = log_expression_matrix, sigs = MES_functional_list)
  score_matrix_shuffled<- scalop::sigScores(m = shuffled_log_tpm_matrix, sigs = MES_functional_list)
  

  #function get_threshold = calculates the threshold based on 99% percentile of shuffled dataset
  get_threshold<- function(score_matrix, cell_type_df,score_matrix_shuffled,signature=NULL, threshold=0.99)
  {
    
    cells_per_tumor_in_study<- table(cell_type_df$Sample)
    sample_number= 2* min(cells_per_tumor_in_study) ## number of cells to sample in patients with large number of cells
    avg_threshold<- c()
    
    for(tumor in names(cells_per_tumor_in_study)){
      cells_from_tumor<- cell_type_df[cell_type_df$Sample==tumor,]$Cell
      
      if(cells_per_tumor_in_study[tumor]> sample_number)
        cells_from_tumor<- sample(cells_from_tumor,size = sample_number, replace=T)
      
      shuffled_score_matrix_tumor<- score_matrix_shuffled[cells_from_tumor, ]
      tumor_threshold<- quantile(shuffled_score_matrix_tumor[,plot_factor], threshold)
      names(tumor_threshold)<- tumor
      avg_threshold<- c(avg_threshold, tumor_threshold)
    }
    return(mean(avg_threshold))
    
  }
  
  threshold<- get_threshold(score_matrix, cell_type_df,score_matrix_shuffled,signature, threshold)

  rownames(score_matrix)<- paste("normal",1:nrow(score_matrix))
  rownames(score_matrix_shuffled)<- paste("shuffled",1:nrow(score_matrix_shuffled))

  return(new_df)
}


#####  cor_expression #####  returns correlation of expression of functional MES programs (per tumor and average)

cor_expression<- function(score_matrix, cell_type_df)
{
  score_matrix[rownames(cell_type_df),"Sample"]<- cell_type_df$Sample
  score_matrix_split_by_sample<- split(score_matrix, score_matrix$Sample)
  
  pairwise_cor_calc <- function(DF) {
    n= which(colnames(DF)=="Sample")
    cor_calc<-cor(DF[,-n],method="spearman" )
    
    return(cor_calc)
  }
  
  cor_matrix_list<-lapply(score_matrix_split_by_sample,
                          function(x) pairwise_cor_calc(DF=x) )
  
  mean_cor_matrix<- Reduce(`+`, cor_matrix_list) / length(cor_matrix_list)
  
  return(list("cor_matrix_list"=cor_matrix_list, "mean_cor_matrix"= mean_cor_matrix))
  
}


##### TCGA_associations ##### returns the correlation of expression of MES CORE with T-cell/macrophage specific genes expression in TCGA bulk RNA seq
# TCGA_log_tpm = normalized log tpm expression matrix downloaded from TCGA RNA-seq dataset
# t_cell_specific_genes = genes that are at least 8-fold higher expressed in T-cells than in any other cell-type, as previously described in Hara et al
# macrophage_specific_genes = genes that are at least 8-fold higher expressed in T-cells than in any other cell-type, as previously described in Hara et al


TCGA_associations<- function(TCGA_log_tpm,t_cell_specific_genes,macrophage_specific_genes)
{
  TCGA_score_matrix = scalop::sigScores(m = TCGA_log_tpm, sigs = c(MES_functional_list,"T_Cell"=c("CD2", "CD3D", "CD3E", "CD3G"), "Macrophage"= c("CD14",   "AIF1" ,  "CD163" , "TYROBP", "CSF1R" )))
  corr_df<- data.frame(row.names = rownames(TCGA_log_tpm))
  corr_df[,"Correlation_T_Cell"]<- apply(TCGA_log_tpm,1,function(x){cor(x, TCGA_score_matrix$T_Cell)})
  corr_df[,"Correlation_Macrophage"]<- apply(TCGA_log_tpm,1,function(x){cor(x, TCGA_score_matrix$Macrophage)})
  corr_df[,"Correlation_MES"]<- apply(TCGA_log_tpm,1,function(x){cor(x, sigScores_tcga$MES_CORE_GENES)})
  
  return(list("T_cell_correlation"=cor.test(corr_df[t_cell_specific_genes,"Correlation_T_Cell"],corr_df[t_cell_specific_genes,"Correlation_MES"]), "Macrophage_correlation"=  cor.test(corr_df[t_cell_specific_genes,"Correlation_Macrophage"],corr_df[t_cell_specific_genes,"Correlation_MES"])))
  
}


##### TCGA_correlations ##### returns the correlation of each MES states with feature, and the partial correlation of each MES state with the MES CORE score for the feature in TCGA bulk RNA seq
# TCGA_log_tpm = normalized log tpm expression matrix downloaded from TCGA RNA-seq dataset
# feature = signature of T cell/macrophage state in question


TCGA_correlations<- function(TCGA_log_tpm, feature)
{
  TCGA_score_matrix = scalop::sigScores(TCGA_log_tpm, c(MES_functional_list,feature=feature))
  TCGA_score_matrix$MES_HYPOXIA= TCGA_score_matrix$MES_CORE_GENES- TCGA_score_matrix$HYPOXIA
  TCGA_score_matrix$MES_ASTROCYTE= TCGA_score_matrix$MES_CORE_GENES- TCGA_score_matrix$ASTROCYTE
  
  bar_corr<- data.frame(matrix(nrow=1,ncol=3)) ## bar_corr = pearson correlation table
  colnames(bar_corr)<- c("MES_CORE_GENES","MES_HYPOXIA","MES_ASTROCYTE")
  
  for(program in colnames(bar_corr))
  {
    x1 <- TCGA_score_matrix[,program]
    y1 <- TCGA_score_matrix[,feature]
    z1 <- replicate(1000, sample(1:length(x1),75, replace = T)) ### bootstrap
    bar_corr[1:1000,program]<-(unlist(apply(z1, 2, function(x){cor.test(x1[x], y1[x], method="pearson")$estimate})))
  }
  
  bar_prcor<- data.frame(row.names = 1:1000)  ## bar_prcor = partial correlation table
  bar_prcor[,"MES_CORE_GENES"]<- bar_corr[,"MES_CORE_GENES"]
  
  bar_prcor[1:1000,"MES_HYPOXIA_PCOR"]<-(unlist(apply(z1, 2, function(x){
    pcor(TCGA_score_matrix[x,c(feature,"HYPOXIA","MES_CORE_GENES")], method="pearson")$estimate["MES_CORE_GENES",feature]
  })))
  
  bar_prcor[1:1000,"MES_ASTROCYTE_PCOR"]<-(unlist(apply(z1, 2, function(x){
    pcor(TCGA_score_matrix[x,c(feature,"ASTROCYTE","MES_CORE_GENES")], method="pearson")$estimate["MES_CORE_GENES",feature]
  })))
  
  return(list("correlation_df"=bar_corr, "partial_correlation_df"=bar_prcor ))

}
