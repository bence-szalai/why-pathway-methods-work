library(viper)
ncores=1
set.seed(19890904)

for (dname in c('tcga','progeny')){
  set.seed(19890904)
  gex=read.csv(paste0('../results/benchmark/datasets/',dname,'_data.csv'),sep=',',header=T,row.names = 1)
  genesets=list.files('../results/genesets/overlap/rdatas/')
  genesets=unlist(lapply(genesets, function(x){substr(x,1,nchar(x)-6)}))
  for (geneset in genesets){
    load(paste0('../results/genesets/overlap/rdatas/',geneset,'.rdata'))
    activities = viper(eset = gex, regulon = viper_geneset, nes = T, 
                       method = 'none' ,minsize = 4, eset.filter = F,cores = ncores)
    write.csv(activities,paste0('../results/benchmark/scores/',dname,'/',geneset,'.csv'))
  }
}
