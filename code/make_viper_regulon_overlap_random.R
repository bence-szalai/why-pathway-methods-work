make_regulon=function(x){
  result=list()
  result$tfmode=rep(1,length(x))
  names(result$tfmode)=x
  result$likelihood=rep(1,length(x))
  result
}

dnames=c('overlap')
for (dname in dnames){
  fnames=list.files('../results/genesets/overlap/csvs/')
  fnames=fnames[grep('random',fnames)]
  for (fname in fnames){
    data=read.csv(paste0('../results/genesets/',dname,'/csvs/',fname),sep = ',',header = T,row.names = 1)
    data=split(data$Gene,data$Set)
    viper_geneset=lapply(data,make_regulon)
    save(viper_geneset,file=paste0('../results/genesets/',dname,'/rdatas/',substr(fname,1,nchar(fname)-4),'.rdata'))
  }
}