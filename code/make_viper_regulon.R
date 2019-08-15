data=read.csv('../results/genesets/dorothea/csvs/dorothea_A_set_gene.csv',
              sep=',',check.names=FALSE,header = T,row.names = 1)
data=as.matrix(data)

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

data=apply(data,1,make_regulon)
save(data,file='../results/genesets/dorothea/rdata/dorothea_A_set_gene.rdata')
