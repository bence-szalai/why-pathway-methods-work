make_regulon=function(x){
  temp=unique(x)
  fil=temp!='NONE'
  temp=temp[fil]
  result=list()
  result$tfmode=rep(1,length(temp))
  names(result$tfmode)=temp
  result$likelihood=rep(1,length(temp))
  result
}

fnames=c('A','B','C','D','E','AB','ABC','ABCD','ABCDE','BEST')

for (fname in fnames){
  data=read.csv(paste0('../results/genesets/dorothea/csvs/dorothea_',fname,'_set_gene.csv'),
                sep=',',check.names=FALSE,header = T,row.names = 1)
  data=as.matrix(data)
  data=apply(data,1,make_regulon)
  save(data,file=paste0('../results/genesets/dorothea/rdata/dorothea_',fname,'_set_gene.rdata'))
}

fnames=c('BIOCARTA','KEGG','CGP','REACTOME')

for (fname in fnames){
  data=read.csv(paste0('../results/genesets/msigdb/csvs/',fname,'_set_gene.csv'),
                sep=',',check.names=FALSE,header = T,row.names = 1)
  data=as.matrix(data)
  data=apply(data,1,make_regulon)
  save(data,file=paste0('../results/genesets/msigdb/rdata/',fname,'_set_gene.rdata'))
}

fnames=unlist(lapply(list.files('../results/genesets/overlaps/csvs/'),FUN=function(x) {substr(x,start = 1,stop = nchar(x)-4)}))

for (fname in fnames){
  data=read.csv(paste0('../results/genesets/overlaps/csvs/',fname,'.csv'),
                sep=',',check.names=FALSE,header = T,row.names = 1)
  data=as.matrix(data)
  data=apply(data,1,make_regulon)
  save(data,file=paste0('../results/genesets/overlaps/rdata/',fname,'.rdata'))
}
