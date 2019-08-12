### saves wiper regulons as csv files

for (confindence in c('A','B','C','D','E','BEST')){
  load(paste0('../data/viper_rdata/',confindence,'_viperRegulon.rdata'))
  results=unlist(viper_regulon)
  write.csv(results,paste0('../results/genesets/dorothea/raw/',confindence,'_viperRegulon.csv'))
}