library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

#Load in DF - This should contain pvalues, degs, log2fcs
path <- "C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Proteomics\\Joyce_Proteomics\\CD4_M_vs_CD4_F.csv"
df <- read.csv(path)

#Extraction of only positive in DF (when comparing one group to another)
df<- subset(df, log2FoldChange > 0.01 & pvalue < 0.05)

#Extraction of only negative in DF (only one can be used per time obviously, # out one)
df<- subset(df, log2FoldChange < -0.01 & pvalue < 0.05)

#DataFrame is named 'df', use this to access files
gene_names_mouse <- df$genes
log2_fold_changes_mouse <- df$log2FoldChange
p_values_mouse <- df$pvalue

# Load gene annotations package for mice
gene_annot <- org.Mm.eg.db

# Perform pathway enrichment analysis
enrich_result_mouse <- enrichGO(gene = gene_names_mouse,
                                keyType = 'SYMBOL',
                                  OrgDb = org.Mm.eg.db, pvalueCutoff = 0.01)

# Visualize enriched pathways with a dot plot
dotplot(enrich_result_mouse, showCategory = 10) + ggtitle('CD4 Male vs CD4 Female')
