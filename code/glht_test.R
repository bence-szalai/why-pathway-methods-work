data = data.frame(row.names = seq(1,120))
data$B = c(rep('A', 30), rep('B', 30), rep('C', 30), rep('D', 30))
data$C = rnorm(120)
data$A = (data$B=='A') * 2 + (data$B=='B') * 1 + data$C + rnorm(120)
data$B = as.factor(data$B)

model = lm('A ~ B + C', data=data)
summary(model)
anova(aov(model))
TukeyHSD(aov(model))

library(multcomp)
summary(glht(aov(model), linfct = mcp(B = "Tukey")))

model2 = lm('A ~ C', data=data)
data$R = residuals(model2)
model2 = lm('R ~ B', data=data)
TukeyHSD(aov(model2))
