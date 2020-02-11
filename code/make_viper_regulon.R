make_regulon=function(x){
  result=list()
  result$tfmode=rep(1,length(x))
  names(result$tfmode)=x
  result$likelihood=rep(1,length(x))
  result
}

dnames=c('single')
for (dname in dnames){
  fnames=c('BEST_dorothea_AB_filtered.csv','BEST_dorothea_CD_filtered.csv',
           'BIOCARTA_filtered.csv','CGP_filtered.csv','KEGG_filtered.csv',
           'REACTOME_filtered.csv')
  for (fname in fnames){
    data=read.csv(paste0('../results/genesets/',dname,'/csvs/',fname),sep = ',',header = T,row.names = 1)
    data=split(data$Gene,data$Set)
    viper_geneset=lapply(data,make_regulon)
    save(viper_geneset,file=paste0('../results/genesets/',dname,'/rdatas/',substr(fname,1,nchar(fname)-4),'.rdata'))
  }
}

