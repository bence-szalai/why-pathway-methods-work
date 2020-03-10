setwd("~/Documents/Projects/why-pathway-methods-work/code")

data=read.csv('../results/statistics/similaritiy_informative/progeny_jaccard_abs_filtered.csv',
              sep=',',header=T,row.names = 1)
model=lm(r~to + dd_to,data=data)

summary(model)
summary(aov(model))

model=lm(r~to + of + dd_to + dd_from,data=data)
