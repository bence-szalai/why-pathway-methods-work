setwd("~/Documents/Projects/why-pathway-methods-work/code")

data=read.csv('../results/benchmark/SFig4B.csv',sep=',',header = T)
data$Database=relevel(data$Database,'DoRothEA_AB')
model=lm('Score ~ Database + Random ',data=data)
summary(model)
summary(aov(model))


fil=data$Random=='Real'
model=lm(Score ~ Database,data=data[fil,])
summary(aov(model))
TukeyHSD(aov(model))
