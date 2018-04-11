####################################################################
###
### This Script calculates 
### 
### Input data:
### ("./3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/")
### Output data are saved as Rdata file:
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("e1071", "doSNOW", "ipred", "xgboost", "caret")
install.packages(required.packages)
install.packages("ddalpha")
ipak(required.packages)

# Set parameters
download.method = "Assembler_Panca_Normalized"
pw_selection_version = "3.3"
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized"

# Load data
load(paste0("./4_Analysis/", 
            download.method,
            "/Pan_Cancer/Clustering/Hallmark_and_ICR_cluster_assignment_allcancers_ICRbenefit_cluster_v2.Rdata"))
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Select gene list
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/ICR_All_Pathway_High_Survival_analysis_High_vs_Low_Groups_Oncogenic_pathways.Rdata"))
All_survival_analysis_data_ONCO_High = All_survival_analysis_data_ONCO
rm(All_survival_analysis_data_ONCO)
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/ICR_All_Pathway_Low_Survival_analysis_High_vs_Low_Groups_Oncogenic_pathways.Rdata"))
All_survival_analysis_data_ONCO_Low = All_survival_analysis_data_ONCO
rm(All_survival_analysis_data_ONCO)

merged = merge(All_survival_analysis_data_ONCO_High, All_survival_analysis_data_ONCO_Low, by = "Oncogenic_Pathway",
               suffixes = c("_High", "_Low"))

merged$score_contribution = NA
merged$score_contribution[which(merged$HR_High >1 & merged$HR_Low <1 & merged$p_value_High < 0.05 & merged$p_value_Low < 0.05)] = "enabling"
merged$score_contribution[which(merged$HR_High <1 & merged$HR_Low >1 & merged$p_value_High < 0.05 & merged$p_value_Low < 0.05)] = "disabling"

enabling_pathways = gsub("_cluster_Pancancer", "", as.character(merged$Oncogenic_Pathway[which(merged$score_contribution == "enabling")]))
enabling_pathways[which(enabling_pathways == "Proliferation")] = c("[LM] Proliferation")
disabling_pathways = gsub("_cluster_Pancancer", "", as.character(merged$Oncogenic_Pathway[which(merged$score_contribution == "disabling")]))
pathways = c(enabling_pathways, disabling_pathways)

load(paste0(code_path,"/Datalists/Selected.pathways.", pw_selection_version,".RData"))
pathways = Selected.pathways[which(names(Selected.pathways) %in% pathways)]
genes = unique(unname(unlist(pathways)))

# Create subsetted expression matrix and log transform
load(paste0("./3_DataProcessing/", download.method, "/Pancancer/RNASeqData/Pancancer_gene_RNAseq_normalized_TissueType_Filtered.Rdata"))
filtered.norm.RNAseqData.Panca.subset = filtered.norm.RNAseqData.Panca[which(rownames(filtered.norm.RNAseqData.Panca) %in% genes),]
t_filtered.norm.RNAseqData.Panca.subset = t(filtered.norm.RNAseqData.Panca.subset)
RNAseqData.Panca.subset_log2 = log(t_filtered.norm.RNAseqData.Panca.subset +1, 2)
train = as.data.frame(RNAseqData.Panca.subset_log2)
train$Benefit_cluster = Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_cluster[match(rownames(train),
                                                                                                 rownames(Hallmark_and_ICR_cluster_assignment_allcancers))]

# Split data in training and test sets
set.seed(54321)
indexes = createDataPartition(train$Benefit_cluster,
                              times = 1,
                              p = 0.7,
                              list = FALSE)
pancancer.train = train[indexes,]
pancancer.test = train[-indexes,]

# Examine the proportions of the benefit clusters across the datasets
prop.table(table(pancancer.test$Benefit_cluster))
prop.table(table(pancancer.train$Benefit_cluster))

# Train model
## 10-fold cross validation (perform algorithm 10 times) repeated 3 times and use a
## grid search for optimal model hyperparameter values
train.control = trainControl(method = "repeatedcv",
                             number = 10,
                             repeats = 3,
                             search = "grid")

tune.grid <- expand.grid(eta = c(0.05, 0.075, 0.1),
                         nrounds = c(50, 75, 100),
                         max_depth = 6:8,
                         min_child_weight = c(2.0, 2.25, 2.5),
                         colsample_bytree = c(0.3, 0.4, 0.5),
                         gamma = 0,
                         subsample = 1)
View(tune.grid)

# Use the doSNOW package to enable caret to train in parallel.
# While there are many package options in this space, doSNOW
# has the advantage of working on both Windows and Mac OS X.
#
# Create a socket cluster using 10 processes. 
#
# NOTE - Tune this number based on the number of cores/threads 
# available on your machine!!!
#
cl <- makeCluster(10, type = "SOCK")

# Register cluster so that caret will know to train in parallel.
registerDoSNOW(cl)

caret.cv = train(Benefit_cluster ~ .,
                 data = pancancer.train,
                 method = "svmLinear",
                 tuneGrid = tune.grid,
                 trControl = train.control)

stopCluster(cl)

# Examine caret's processing results
caret.cv

# Make predictions on the test set using a xgboost model 
# trained on all rows of the training set using the 
# found optimal hyperparameter values.
preds = predict(caret.cv, pancancer.test)

# Use caret's confusionMatrix() function to estimate the 
# effectiveness of this model on unseen, new data.
confusionMatrix(preds, pancancer.test$Benefit_cluster)
