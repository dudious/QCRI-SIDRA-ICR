

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


df_pathway_analysis_subset = df_pathway_analysis[which(df_pathway_analysis$Percentage_of_Sign1_also_in_Sign2 > 20 | 
                                                         df_pathway_analysis$Percentage_of_Sig2_also_in_Sign1 > 20),]

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
percent_matrix_filtered[percent_matrix_filtered < 20] = NA
percent_matrix_filtered = percent_matrix_filtered[-which(rowSums(is.na(percent_matrix_filtered)) == ncol(percent_matrix_filtered)),]
percent_matrix_filtered = percent_matrix_filtered[,-which(colSums(is.na(percent_matrix_filtered)) == nrow(percent_matrix_filtered))]

percent_matrix_filtered_full = percent_matrix[rownames(percent_matrix_filtered),colnames(percent_matrix_filtered)]

percent_filtered_melt = melt(percent_matrix_filtered_full)
plot_matrix = ggplot(percent_filtered_melt, aes(x = col, y =  row)) +
  geom_tile(aes(fill = value)) +
  labs(x = "", y = "") + 
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 5))
dev.new()
print(plot_matrix)


