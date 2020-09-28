setwd("~/Documents/Projects/why-pathway-methods-work/code")
library(multcomp)
data = read.csv('../results/benchmark/progeny_abs_stat.csv', sep=',', header=TRUE, row.names = 1)
data$database = as.factor(data$database)
data$random = as.factor(data$random)
model = lm('score ~ database + random + size', data=data)
model = aov(model)
anova(model)


fil = data$random == 'False'
model = lm('score ~ database', data=data[fil,])
model = aov(model)
summary(glht(model, linfct = mcp(database = "Tukey")))

model = lm('score ~ database + size', data=data[fil,])
model = aov(model)
summary(glht(model, linfct = mcp(database = "Tukey")))
