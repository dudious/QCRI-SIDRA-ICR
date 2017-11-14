
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                     # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                           # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R")) 

required.packages <- c("VennDiagram")
ipak(required.packages) 


load("/Users/jroelands/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/Selected Pathways/Selected.pathways.RData")

# Get the combinations of names of list elements
Pathway_combination = combn(names(Selected.pathways) , 2 , FUN = paste0 , collapse = "_vs_" , simplify = FALSE)

# Make the combinations of list elements
ll <- combn(Selected.pathways, 2 , simplify = FALSE )

# Intersect the list elements
Number_overlapping_genes = lapply( ll , function(x) length(intersect(x[[1]] , x[[2]])))
Overlapping_genes = lapply( ll , function(x) intersect(x[[1]] , x[[2]]))

# Display number of genes in each signature
Number_genes_in_signature1 = lapply(ll, function(x) length(x[[1]]))
Number_genes_in_signature2 = lapply(ll, function(x) length(x[[2]]))

df_pathway_analysis = data.frame(cbind(Pathway_combination, 
                                       Number_overlapping_genes,
                                       Overlapping_genes,
                                       Number_genes_in_signature1, 
                                       Number_genes_in_signature2))

df_pathway_analysis$Pathway_combination = as.character(df_pathway_analysis$Pathway_combination)
df_pathway_analysis$Number_genes_in_signature1 = as.numeric(df_pathway_analysis$Number_genes_in_signature1)
df_pathway_analysis$Number_genes_in_signature2 = as.numeric(df_pathway_analysis$Number_genes_in_signature2)
df_pathway_analysis$Number_overlapping_genes = as.numeric(df_pathway_analysis$Number_overlapping_genes)

df_pathway_analysis$Percentage_of_Sign1_also_in_Sign2 = round(( df_pathway_analysis$Number_overlapping_genes / 
                                                                  df_pathway_analysis$Number_genes_in_signature1) * 100, 2) 

df_pathway_analysis$Percentage_of_Sig2_also_in_Sign1 = round(( df_pathway_analysis$Number_overlapping_genes / 
                                                               df_pathway_analysis$Number_genes_in_signature2) * 100, 2)


df_pathway_analysis_subset = df_pathway_analysis[which(df_pathway_analysis$Percentage_of_Sign1_also_in_Sign2 > 30 | 
                                                         df_pathway_analysis$Percentage_of_Sig2_also_in_Sign1 > 30),]

percent_matrix = matrix(nrow = 87, ncol = 87, dimnames = list(col = names(Selected.pathways), row = names(Selected.pathways)))
mode(percent_matrix) = "numeric"


i = 1
for(i in 1:nrow(df_pathway_analysis)){
  x= strsplit(df_pathway_analysis$Pathway_combination[i], split = "_vs_")[[1]][1]
  y= strsplit(df_pathway_analysis$Pathway_combination[i], split = "_vs_")[[1]][2]
  percent_matrix[x,y] = df_pathway_analysis$Percentage_of_Sign1_also_in_Sign2[i]
  percent_matrix[y,x] = df_pathway_analysis$Percentage_of_Sig2_also_in_Sign1[i]
}


percent_matrix_filtered = percent_matrix
percent_matrix_filtered[percent_matrix_filtered < 30] = NA
percent_matrix_filtered = percent_matrix_filtered[-which(rowSums(is.na(percent_matrix_filtered)) == ncol(percent_matrix_filtered)),]
percent_matrix_filtered = percent_matrix_filtered[,-which(colSums(is.na(percent_matrix_filtered)) == nrow(percent_matrix_filtered))]

overlaps <- list(x = "a")
for (i in 1:nrow(percent_matrix_filtered)){
  combo = list(c(rownames(percent_matrix_filtered)[i],colnames(percent_matrix_filtered)[!is.na(percent_matrix_filtered[i,])]))
  overlaps = c(overlaps,combo)
}
overlaps = overlaps [-1]

for (i in 1:ncol(percent_matrix_filtered)){
  combo = list(c(colnames(percent_matrix_filtered)[i],rownames(percent_matrix_filtered)[!is.na(percent_matrix_filtered[,i])]))
  overlaps = c(overlaps,combo)
}

# delete complete overlapping groups
for (i in 1:length(overlaps)){
  pathways_to_match = overlaps[i][[1]]
  individual_cond = lapply(overlaps, function(x) pathways_to_match %in% x)
  cond = which(as.logical(lapply(individual_cond, function(x) all(x))))
  cond = cond[-which(cond == i)]
  if(length(cond)>0){overlaps = overlaps[-cond]}
}

# merge partially overlapping groups
overlaps.merged <- list(x = "a")
for (i in 1:length(overlaps)){
  pathways_to_match = overlaps[i][[1]]
  individual_cond = lapply(overlaps, function(x) pathways_to_match %in% x)
  cond = which(as.logical(lapply(individual_cond, function(x) sum(x))>1))
  cond = cond[-which(cond == i)]
  if(length(cond)>0){
    merged.pathway = list(unique(c(pathways_to_match,unlist(overlaps[cond]))))
    #overlap.index = c(overlap.index,cond)
  } else {merged.pathway = overlaps[i]}
  overlaps.merged = c(overlaps.merged,merged.pathway)
}
overlaps.merged = overlaps.merged [-1]

# delete complete overlapping groups
for (i in 1:length(overlaps.merged)){
  pathways_to_match = overlaps.merged[i][[1]]
  individual_cond = lapply(overlaps.merged, function(x) pathways_to_match %in% x)
  cond = which(as.logical(lapply(individual_cond, function(x) all(x))))
  cond = cond[-which(cond == i)]
  if(length(cond)>0){overlaps.merged = overlaps.merged[-cond]}
}

## Creating plot

percent_matrix_filtered_full = percent_matrix[rownames(percent_matrix_filtered),colnames(percent_matrix_filtered)]

percent_filtered_melt = melt(percent_matrix_filtered_full)
plot_matrix = ggplot(percent_filtered_melt, aes(x = col, y =  row)) +
  geom_tile(aes(fill = value)) +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 5))
dev.new()
print(plot_matrix)


## Creating Venn diagrams
dev.new()

dir.create("./5_Figures/Signature_Analysis", showWarnings = FALSE)
dir.create("./5_Figures/Signature_Analysis/VennDiagrams", showWarnings = FALSE)

overlaps.merged

  #number_venns = length(unlist(overlaps.merged[i]))
  
  #signature_list = list(x = "b")
  #for(j in 1:number_venns){
    #assign(paste0(signature[j]), unlist(overlaps.merged[j])[1])
    
    #signature_list = list(signature_list, signature[j])
  #}
  #signature_list = signature_list[-1]
  #unlist(overlaps.merged[i])[1]
  #unlist(overlaps.merged[i])[2]
  #signature_list = list(Selected.pathways$`[HM] HYPOXIA`, Selected.pathways$`[HM] GLYCOLYSIS`)
  
signature_list = list()
  venn.diagram(
        x = list(Selected.pathways$`[HM] INTERFERON ALPHA RESPONSE`,
               Selected.pathways$`[HM] INTERFERON GAMMA RESPONSE`, Selected.pathways$`[IPA] Interferon Signaling`),
    filename = paste0("./5_Figures/Signature_Analysis/VennDiagrams/HM_interferon_IPA_interferon.tiff"),
    col = "transparent",
    lty = "dotted",
    lwd = 4,
    fill = c("red", "blue", "green"),
    alpha = 0.60,
    units = "in",
    width = 15,
    height = 15,
    #label.col = c("darkred", "white", "darkblue"),
    cex = 5,
    category.names = c("[HM] INTERFERON ALPHA RESPONSE", "[HM] INTERFERON GAMMA RESPONSE",
                       "[IPA] Interferon Signaling"),
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkred", "darkblue", "darkgreen"),
    cat.cex = 2,
    cat.fontfamily = "serif",
    margin = 0.3
  )





