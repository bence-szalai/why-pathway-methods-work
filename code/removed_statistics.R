setwd("~/Documents/Projects/why-pathway-methods-work/code")
library(multcomp)

fnames = c('progeny_abs_remove', 'progeny_remove', 'gdsc_abs_remove', 'gdsc_remove')
stat_results = data.frame(row.names = c('anova_d', 'anova_r', 'anova_s',
                                        'pcd', 'pbio', 'pcgp', 'pkegg', 'preac',
                                        'cd', 'bio', 'cgp', 'kegg', 'reac'))
for (i in c(1,2,3,4)){
  stat_results[,i]=0
  
  data = read.csv(paste0('../results/benchmark/',fnames[i],'.csv'), sep=',', header=TRUE, row.names = 1)
  colnames(data) = c('database', 'score', 'size', 'random')
  data$database = as.factor(data$database)
  data$random = as.factor(data$random)
  model = lm('score ~ database + random + size', data=data)
  model = aov(model)
  anova_results = anova(model)
  stat_results[1:3,i] = anova_results$`Pr(>F)`[1:3]
  
  
  fil = data$random == 'Actual gene set'
  model = lm('score ~ database + size', data=data[fil,])
  model = aov(model)
  tukey_results = summary(glht(model, linfct = mcp(database = "Tukey")))
  stat_results[4:8, i] = tukey_results$test$pvalues[1:5]
  stat_results[9:13, i] = tukey_results$test$coefficients[1:5]
}
stat_results = t(stat_results)
rownames(stat_results) = fnames
stat_results = signif(stat_results, 3)
write.csv(stat_results, '../results/benchmark/remove_stats.csv')