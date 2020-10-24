data = read.csv('../results/benchmark/informative_similarity_stats.csv', sep=',',
                header=T)

data$r_part = signif(data$r_part, 3)
data$r_spearman = signif(data$r_spearman, 3)
data$p_part = signif(data$p_part, 3)
data$p_spearman = signif(data$p_spearman, 3)

write.csv(data, '../results/benchmark/informative_similarity_stats.csv')