
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"                                                                # Setwd to location were output files have to be saved.

source(paste0(code_path,"R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.bioconductor.packages = c("GSVA","heatmap3", "gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/", download.method ,"/3.7.Proliferation_metagene_Enrichment/3.7.Proliferation_metagene_Enrichment_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("HML ICR clusters", "Bindea clusters")
Legend = c("ICR Low","ICR Med","ICR High", "Bindea Low", "Bindea High")
Legend_colors = c("blue","green","red", "pink", "purple")

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0(code_path, "Datalists/Proliferation.metagene.RData"))

# Create folders and log file
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folders to save Rdata.files
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/3.7.Proliferation_metagene_Enrichment"), showWarnings = FALSE)
cat("This is a log file for calculation of proliferation metagene scores using RNASeq data",   # Set-up logfile
    "_________________________________________________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Script Running Date :",
    capture.output(Sys.time()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),                                                          
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    "",
    "Calculating Deconvolution scores",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

start.time.process.all = Sys.time()
msg = paste0("Calculating ES scores for all Cancers", "\n")
cat(msg)

i =1
for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Calculating enrichment scores ",Cancer,"."))
  
  ## load RNASeq data
  if(Cancer == "LAML") 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  if(Cancer == "SKCM"){
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  } else{
    Cancer_path = paste0 ("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData")
    load(paste0(Cancer_path, "/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
  
  ## load cluster data
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  if(download.method == "TCGA_Assembler" | download.method == "Assembler_Panca_Normalized"){
    Expression.data = log(filtered.norm.RNAseqData +1, 2)
  }
  if(download.method == "Pancancer_matrix"){
    Expression.data = filtered.norm.RNAseqData
    load("./4_Analysis/Pancancer_matrix/Pancancer_matrix_QC/Genes_with_NA.Rdata")
    Expression.data = Expression.data[-which(rownames(Expression.data) %in% all_genes_that_have_NA_for_at_least_1_patient),]
  }
  available_genes = Proliferation_metagene_Miller[which(Proliferation_metagene_Miller %in% rownames(Expression.data))]
  unavailable_genes_RNAseq = Proliferation_metagene_Miller[-which(Proliferation_metagene_Miller %in% rownames(Expression.data))]
  
  
  ## Proliferation ssGSEA
  Proliferation.metagene.enrichment.score = gsva(Expression.data,list(Proliferation_metagene_Miller),method="ssgsea")
  Proliferation.metagene.enrichment.z.score = Proliferation.metagene.enrichment.score 
  for(j in 1: nrow(Proliferation.metagene.enrichment.z.score))  {
    Proliferation.metagene.enrichment.z.score[j,] = (Proliferation.metagene.enrichment.score[j,]-mean(Proliferation.metagene.enrichment.score[j,]))/sd(Proliferation.metagene.enrichment.score[j,]) # z-score the enrichment matrix
  }
  
  ## Bindea cluster (hierarchical)
  sHc = hclust(ddist <- dist(t(Proliferation.metagene.enrichment.score)), method = "ward.D2")
  
  plot(sHc,labels=FALSE)
  
  table_cluster_assignment$Proliferation_cluster = cutree(sHc,k = 2)[match(rownames(table_cluster_assignment),names(cutree(sHc,k = 2)))]
  table_cluster_assignment$Proliferation_score = Proliferation.metagene.enrichment.score[1,]
  
  proliferation_cluster_means = aggregate(Proliferation_score~Proliferation_cluster, data = table_cluster_assignment, FUN = mean)
  proliferation_cluster_means = proliferation_cluster_means[order(proliferation_cluster_means$Proliferation_score),]
  proliferation_cluster_means$Proliferation_cluster_name = c("Proliferation Low", "Proliferation High")
  table_cluster_assignment$Proliferation_cluster_name = proliferation_cluster_means$Proliferation_cluster_name[match(table_cluster_assignment$Proliferation_cluster,
                                                                                                                     proliferation_cluster_means$Proliferation_cluster)]
  
  table_cluster_assignment$Proliferation_cluster = NULL
  
  ## Save Scores
  dir.create(paste0("./4_Analysis/",download.method, "/", Cancer),showWarnings = FALSE)
  dir.create(paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment"))
  
  save(Proliferation.metagene.enrichment.score, 
       file = paste0("./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
                                                     "_Proliferation_metagene.Rdata"))
  save(optimal.calinsky, table_cluster_assignment, file = paste0("./4_Analysis/", download.method, "/", Cancer,"/Clustering/",
                                                                 Cancer, ".", download.method, ".EDASeq.ICR.reps5000/", Cancer, "_ICR_cluster_assignment_k2-6.Rdata"))

}

load("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata")
load(paste0("./4_Analysis/", download.method, "/ACC/Clustering/ACC.", download.method, ".EDASeq.ICR.reps5000/ACC_ICR_cluster_assignment_k2-6.Rdata"))
table_cluster_assignment$Cancer = "ACC"

ICR_cluster_assignment_allcancers = table_cluster_assignment

i=2
for (i in 2:N.sets){
  Cancer = CancerTYPES[i]
  if(Cancer == "LAML"){next}
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/", 
              Cancer, "_ICR_cluster_assignment_k2-6.Rdata"))
  table_cluster_assignment$Cancer = Cancer
  ICR_cluster_assignment_allcancers = rbind(ICR_cluster_assignment_allcancers, table_cluster_assignment)
}

save(ICR_cluster_assignment_allcancers, Cancer_color_table, file = paste0("./4_Analysis/", download.method,
                                                                           "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers.Rdata"))


#### Make pancancer file:

load(paste0("./4_Analysis/Assembler_Panca_Normalized/ACC/Signature_Enrichment/GSEA_ACC_Proliferation_metagene.Rdata"))
Proliferation.metagene.enrichment.score.all = Proliferation.metagene.enrichment.score

i=1
for (i in 2:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  cat (paste0 ("Calculating enrichment scores ",Cancer,"."))
  
  ## load RNASeq data
  if(Cancer == "LAML") 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer ,"_Proliferation_metagene.Rdata"))
  Proliferation.metagene.enrichment.score.all = cbind(Proliferation.metagene.enrichment.score.all, Proliferation.metagene.enrichment.score)
}

save(Proliferation.metagene.enrichment.score.all, file = "./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Signature_Enrichment/ssGSEA.Proliferation.metagene.Pancancer.Rdata")

### Append to oncogenic pathways

rm(list = ls())
load("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark.Pancancer.Rdata")
load("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Signature_Enrichment/ssGSEA.Proliferation.metagene.Pancancer.Rdata")

rownames(Proliferation.metagene.enrichment.score.all) = c("[LM] Proliferation")
Proliferation.metagene.enrichment.score.all = as.data.frame(Proliferation.metagene.enrichment.score.all, drop = FALSE)
Hallmark_enrichment_score_all = as.data.frame(Hallmark_enrichment_score_all)

Hallmark_enrichment_score_all = rbind(Hallmark_enrichment_score_all, Proliferation.metagene.enrichment.score.all[, order(colnames(Hallmark_enrichment_score_all))])
Hallmark_enrichment_score_all = as.matrix(Hallmark_enrichment_score_all)

save(Hallmark_enrichment_score_all ,file = "./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Signature_Enrichment/ssGSEA.Hallmark&Proliferation.Pancancer.Rdata")
