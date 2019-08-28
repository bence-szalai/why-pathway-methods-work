library(msigdbr)

m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
m_df= as.data.frame(m_df)[,c('gs_name','human_gene_symbol')]
write.csv(m_df,'../results/genesets/single//raw/KEGG.csv')

m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
m_df= as.data.frame(m_df)[,c('gs_name','human_gene_symbol')]
write.csv(m_df,'../results/genesets/single//raw/REACTOME.csv')

m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
m_df= as.data.frame(m_df)[,c('gs_name','human_gene_symbol')]
write.csv(m_df,'../results/genesets/single//raw/BIOCARTA.csv')

m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
m_df= as.data.frame(m_df)[,c('gs_name','human_gene_symbol')]
write.csv(m_df,'../results/genesets/single//raw/CGP.csv')