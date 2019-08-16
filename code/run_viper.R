library(viper)
ncores=11
set.seed(19890904)
regulon_type='msigdb'
regulon_name='KEGG'
load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'_set_gene.rdata'))
regulon=data
data=read.csv('../results/benchmark/gdsc/raw/gdsc_gex.csv',sep=',',header=T,row.names = 1,check.names=FALSE)
activities = viper(eset = data, regulon = regulon, nes = T, 
                   method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
write.csv()
