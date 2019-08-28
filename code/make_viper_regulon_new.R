make_regulon=function(x){
  result=list()
  result$tfmode=rep(1,length(x))
  names(result$tfmode)=x
  result$likelihood=rep(1,length(x))
  result
}
