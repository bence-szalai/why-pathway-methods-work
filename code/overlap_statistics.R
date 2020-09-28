setwd("~/Documents/Projects/why-pathway-methods-work/code")
library(multcomp)
data = read.csv('../results/benchmark/progeny_abs_overlap.csv', sep=',', header=TRUE, row.names = 1)
data$Data = as.factor(data$Data)
data$Gene.set = as.factor(data$Gene.set)
model = lm('Delta ~ Data + Gene.set + Size', data=data)
model = aov(model)
anova(model)


fil = data$Gene.set == 'Actual gene set'
model = lm('Delta ~ Data', data=data[fil,])
model = aov(model)
summary(glht(model, linfct = mcp(Data = "Tukey")))

model = lm('Delta ~ Data + Size', data=data[fil,])
model = aov(model)
summary(glht(model, linfct = mcp(Data = "Tukey")))
