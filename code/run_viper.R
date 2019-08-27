library(viper)
ncores=6
set.seed(19890904)
#gdsc
gdsc=read.csv('../results/benchmark/gdsc/raw/gdsc_gex.csv',sep=',',header=T,row.names = 1,check.names=FALSE)
#msigdb
regulon_type='msigdb'
regulon_names=c('KEGG','BIOCARTA','REACTOME','CGP')

for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'_set_gene.rdata'))
  regulon=data
  
  activities = viper(eset = gdsc, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/gdsc/',regulon_name,'.csv'))
}
#dorothea
regulon_type='dorothea'
regulon_names=c('dorothea_A','dorothea_AB','dorothea_ABC','dorothea_ABCD','dorothea_ABCDE',
                'dorothea_B','dorothea_C','dorothea_D','dorothea_E','dorothea_BEST')
for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'_set_gene.rdata'))
  regulon=data
  
  activities = viper(eset = gdsc, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/gdsc/',regulon_name,'.csv'))
}

#progeny
progeny=read.csv('../results/benchmark/progeny/raw/progeny_data.csv',sep=',',header=T,row.names = 1,check.names=FALSE)
#msigdb
regulon_type='msigdb'
regulon_names=c('KEGG','BIOCARTA','REACTOME','CGP')

for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'_set_gene.rdata'))
  regulon=data
  
  activities = viper(eset = progeny, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/progeny/',regulon_name,'.csv'))
}
#dorothea
regulon_type='dorothea'
regulon_names=c('dorothea_A','dorothea_AB','dorothea_ABC','dorothea_ABCD','dorothea_ABCDE',
                'dorothea_B','dorothea_C','dorothea_D','dorothea_E','dorothea_BEST')
for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'_set_gene.rdata'))
  regulon=data
  
  activities = viper(eset = progeny, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/progeny/',regulon_name,'.csv'))
}

#overlaps
regulon_type='overlaps'
regulon_names=unlist(lapply(list.files('../results/genesets/overlaps/rdata/'),
                            function(x){substr(x,start = 1,stop = nchar(x)-6)}))
for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'.rdata'))
  regulon=data
  
  activities = viper(eset = progeny, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/progeny/',regulon_name,'.csv'))
}

#tcga
tcga=read.csv('../results/benchmark/tcga/raw/tcga_gex.csv',sep=',',header=T,row.names = 1,check.names=FALSE)
#msigdb
regulon_type='msigdb'
regulon_names=c('KEGG','BIOCARTA','REACTOME','CGP')

for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'_set_gene.rdata'))
  regulon=data
  
  activities = viper(eset = tcga, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/tcga/',regulon_name,'.csv'))
}
#dorothea
regulon_type='dorothea'
regulon_names=c('dorothea_A','dorothea_AB','dorothea_ABC','dorothea_ABCD','dorothea_ABCDE',
                'dorothea_B','dorothea_C','dorothea_D','dorothea_E','dorothea_BEST')
for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'_set_gene.rdata'))
  regulon=data
  
  activities = viper(eset = tcga, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/tcga/',regulon_name,'.csv'))
}

#overlaps
regulon_type='overlaps'
regulon_names=unlist(lapply(list.files('../results/genesets/overlaps/rdata/'),
                            function(x){substr(x,start = 1,stop = nchar(x)-6)}))
for (regulon_name in regulon_names){
  load(paste0('../results/genesets/',regulon_type,'/rdata/',regulon_name,'.rdata'))
  regulon=data
  
  activities = viper(eset = tcga, regulon = regulon, nes = T, 
                     method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
  write.csv(activities,paste0('../results/benchmark/scores/tcga/',regulon_name,'.csv'))
}
