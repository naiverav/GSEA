library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpmisc)
library(gplots)



#Function to calculate correlation using t-test and rank order the list
cal_cor_tstat <- function(D, Sub_S, S_gene_set_list){
  ## Performing t.test for each gene and getting t-statistic
  cor_tstat <- apply(D, MARGIN = 1, function(x) {t.test(x[colnames(D)=='ALL'], x[colnames(D)=='AML'], alternative = 'two.sided')$statistic})
  L <- as.data.frame(cor_tstat)
  Sorted_L <- L[order(cor_tstat), , drop= FALSE]
  Sorted_L$rank <- rank(Sorted_L$cor_tstat, ties.method = 'average')
  
  return(Sorted_L)
  
}


#Function for calculating running sum statistic for each gene in expression set for each gene set
Cal_RunningSum <- function(S_gene_set_list,Sorted_L){
  
  L_cor_tstat_vector <- Sorted_L$cor_tstat
  L_gene_names_vector <- row.names(Sorted_L)

  
  for(i in 1:N_gs){
    run_sum <- 0
    run_sum_vector <- vector(mode = 'numeric')
    
    X_miss <- -(1/(N - Num_genes_in_Set[i]))
    for(j in 1:N){
      
      X_hit <- (abs(L_cor_tstat_vector[j])/Sum_of_cor_tstat_Set[i])
      
      
      ifelse((L_gene_names_vector[j] %in% (S_gene_set_list[[i]])), X<-X_hit, X<-X_miss )
      
      run_sum <- run_sum + X
      
      run_sum_vector[j] <- run_sum
    }
    Sorted_L[,i+2] <- run_sum_vector
  }
  colnames(Sorted_L)[3:ncol(Sorted_L)] <- rownames(GS_out)  
  
  return(Sorted_L)
}



#Function for Plotting graph fr running sum statistics

plot_lineGraph <- function(L_runningSum, gS_name){
  gs <- L_runningSum[ ,gS_name]
  max_dev <- GS_out[gS_name,]
  library(ggplot2)
  ggplot(L_runningSum, aes(x= rank, y = gs))+geom_line(size = 1, color='blue')+
    scale_x_continuous(limits = c(0, 2800), breaks = seq(0, 2800,200))+
    scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.2))+
    geom_hline(aes(yintercept = 0))+
    geom_point(aes(x=which.max(abs(gs)), y = max_dev), size = 3)+
    #geom_text(aes(x=which.max(abs(gs)), y = max_dev, label = "Maximum Deviation from Zero", hjust = 0.5, vjust = 1))+
    geom_step()+
    geom_segment(aes(x = which.max(abs(gs)), y = 0, xend = which.max(abs(gs)), yend = max_dev), color = 'red')+
    theme(axis.line.x = element_line(size = 2, colour = 'dark grey'))+
    xlab("Gene Ranking")+
    ylab("Running Sum Statistic")+
    ggtitle(paste0('Gene Set:\n', gS_name))
}


# Function to calculate ES = maximum deviation from zero, add it to GS_Out file
Cal_EnrichScore <- function(L_runningSum){
  ES_Score_List <- vector()
  for(i in 3:ncol(L_runningSum)){
    rs_min <- min(L_runningSum[,i])
    rs_max <- max(L_runningSum[,i])
    ES <- ifelse(abs(rs_min) > abs(rs_max), rs_min, rs_max)
    ES_Score_List <- append(ES_Score_List, ES)
  }
  return(ES_Score_List)
}




#Read in Lukemia.txt keeping 
D <- read.csv('leukemia.txt', header = TRUE, sep= '\t', row.names = 1, check.names = FALSE)
head(D, n=5)
#print(dim(D))

#N:number of Genes
N = nrow(D)
#k: number of samples
k=ncol(D)

#n_ALL: Number of ALL samples
n_ALL <- length(colnames(D)[(colnames(D)=="ALL")])
#n_AML: Number of AML samples
n_AML <- length(colnames(D)[(colnames(D)=="AML")])


#Reading pathways.txt
S <- read.table('pathways.txt', sep = '\t', fill = TRUE, quote = "", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
#head(S, n=5)
S <-  rename(S, gene_set_name = V1, Description = V2) 
names(S)[3:21] <- paste0("G", 1:19)


#Subset S to get gene sets with at least 15 member genes to focus on the robust signals
Sub_S <- subset(S,(rowSums(S != "")) >=17)
apply(Sub_S, 2, function(x) any(is.na(x)))

length(unique(Sub_S$gene_set_name)) #checking for uniuqe gene sets
Sub_S <- Sub_S[!duplicated(Sub_S$gene_set_name), ] # Removing duplicates


#Converting pathways dataframe S to list of lists
S_gene_set_list <- list()
for(row in 1:nrow(Sub_S)){
  gene_set_name <- Sub_S[row, 'gene_set_name']
  gene_set_vector <- Sub_S[row, 3:ncol(Sub_S)]
  clean_gene_vector <- gene_set_vector[(!is.na(gene_set_vector)) & (gene_set_vector != "")]
  gene_set_list <- list(clean_gene_vector)
  names(gene_set_list)  <- gene_set_name
  S_gene_set_list <- append(S_gene_set_list, gene_set_list)
}
N_gs <- length(S_gene_set_list)

#Calculated correlation to generate rank order gene list from original expressions data set. 
Sorted_L <- cal_cor_tstat(D)



# Create a new data frame with all gene_set_names and calculate sum of correlation statistic for each gene-set(for using in P-hits) 

GS_out <- data.frame(row.names = names(S_gene_set_list), stringsAsFactors = F)
Num_genes_in_Set <- vector(mode="numeric", length=0)
Sum_of_cor_tstat_Set <- vector(mode="numeric", length=0)

for(i in 1:length(S_gene_set_list)){
  
  Num_genes_in_Set[i] <- length(S_gene_set_list[[i]])
  Sum_of_cor_tstat_Set[i] <- sum(abs(Sorted_L[S_gene_set_list[[i]], "cor_tstat"]), na.rm = T)  
}

#Creating a dataframe with Gene set names as Rownames to keep all ES(perm)
ES_df_for_NES <- data.frame(row.names = names(S_gene_set_list)) 

n_perm <- 3 # Number of permutations
for(perm in 1:n_perm){              #First Enrichment Score calculated is bserved ES for the original expressions data set.Printing to output file GS_out 
  if(perm == 1){
    L_runningSum <- Cal_RunningSum(S_gene_set_list,Sorted_L) # Calling L_RunningSum for calculating runningSum for each permutation of expression set
    Observed_ES <- Cal_EnrichScore(L_runningSum)#Keeping a vector for Observed ES for later comparisons and calculations
    GS_out$Enrichment_Score <- signif(Observed_ES,4)

    #Plot graphs for Running Sum Statistic for randomly selected genes-sets from Original gene list
    pdf('GSEA_RunningSum.pdf')
    p1 <- plot_lineGraph(L_runningSum, 'NFKBIL2')
    p2 <- plot_lineGraph(L_runningSum, 'PGAM1')
    p3 <- plot_lineGraph(L_runningSum, 'CBF_LEUKEMIA_DOWNING_AML')
    p4 <- plot_lineGraph(L_runningSum, 'Electron_Transport_Chain')
    
    grid.arrange(p1,p2,p3,p4, ncol=2, top= "Plots for Running Enrichment Score for randomly selected gene-sets")
    dev.off()
    
  }
  else {            # Reassign the sample labels, and reorder the genes for each permutation and calculate ES
    
    names(D) <- names(D)[sample(ncol(D))]      #Randomly assign column names
    D <- D[sample(nrow(D)),]       #Reordering rows in D
    head(D, n=5)
    
    #Calling function cal_cor_tstat for calculating the correlation and get the ranked list
    
    Sorted_L <- cal_cor_tstat(D) #Rank order the gene list by correlation metric

    L_runningSum <- Cal_RunningSum(S_gene_set_list,Sorted_L) # Calling L_RunningSum for calculating runningSum for each permutation of expression set

    ES_df_for_NES[, paste0('ES_Permute_',perm, collapse = NULL)] <- signif(Cal_EnrichScore(L_runningSum), 4) #ES for each permutation
  }
}

### Calculating the Nominal P-Value: p = fraction of ES(S, perm???) values ??? ES(S)relative to null distribution
Nominal_P_vector <- vector(mode = 'numeric')
for(row in 1:length(Observed_ES)){
  
  ES_obs <- Observed_ES[row] 
  pos_ES_perm <- 0 #initializing vector to contain positive ES during permutations
  neg_ES_perm <- 0 #initializing vector to contain negative ES during permutations
  
  for(perm in 1:(n_perm-1)){
    ES_pi <- ES_df_for_NES[row,perm]
    ifelse((ES_pi >= 0), pos_ES_perm<- append(pos_ES_perm,ES_pi), neg_ES_perm <- append(neg_ES_perm, ES_pi))
  } 
  
  # Checking if Observed ES is positive or negative
  if(ES_obs >= 0){
    if(length(which(ES_df_for_NES[row,]>=ES_obs)) == 0){    #If there are no values greater than Observed ES, Nominal P =0
      Nominal_P_vector[row] <- 0
    } else {     #Nominal_P = sum of all ES(perm) which are greater than ES_obs/number of positive values
      Nominal_P_vector[row] <- sum(pos_ES_perm>=ES_obs)/length(pos_ES_perm)
    }
    ifelse(pos_ES_perm == 0, GS_out[row, 'NES'] <- NaN, GS_out[row, 'NES'] <-signif((ES_obs/mean(pos_ES_perm)),4))
  }
  else {
    if(length(which(ES_df_for_NES[row,]<ES_obs)) == 0){
      Nominal_P_vector[row] <- 0
    } else {
      Nominal_P_vector[row] <- sum(neg_ES_perm<ES_obs)/length(neg_ES_perm)
    }
    ifelse(neg_ES_perm == 0, GS_out[row, 'NES'] <- NaN, GS_out[row, 'NES'] <-signif((ES_obs/mean(neg_ES_perm)),4))
  }
}

GS_out$Nominal_P <- signif(Nominal_P_vector, digits = 4)   

### Sort the output file by NES 
GS_out <- GS_out[order(-GS_out$NES), ,] #arrange by descending order of NES
View(GS_out)


#Printing the top 20 to the output - Filename: 'GSEA_ES_TABLE.pdf'#
pdf('GSEA_ES_TABLE.pdf')
GS_top_20 <- head(GS_out, n=20) 
GS_top_20$Gene_Set_Names <- row.names(GS_top_20)
row.names(GS_top_20) <- NULL 
GS_top_20 <- GS_top_20[c(4,1,2,3)]
grid.table(GS_top_20)
dev.off()

View(GS_top_20)


