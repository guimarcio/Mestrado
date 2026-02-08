# clean workspace
rm(list = ls())

# install required packages if needed
packages_needed <- c("ggplot2", "multicomp", "effects")
for (package_name in packages_needed) {      
  if (!(package_name %in% rownames(installed.packages()))){
    install.packages(package_name)
  }
}

# Leitura dos dados e manipulações

dados<-read.csv("todos_results_filt.csv")

# Filtrando o cojunto de dado e analisando

# Selecionar ambiente, ensaio e métrica

amb <- 'seminario'
met <- 'sdr'
ens <- 'P90'

dfilt<-dados[dados$ensaio==ens &
               dados$ambiente==amb,]

# selecionar apenas colunas relevantes
stay <- c("metodo", met)
dfilt<-dfilt[,names(dfilt) %in% stay]


dfilt[, 1]<-as.factor(dfilt[, 1])

summary(dfilt)

# Exploração

boxplot(sdr ~ metodo, data = dfilt, main = "Boxplot por Método")

# Gráfico qqplot
library(car)

# ANOVA
model <- aov(sdr ~ ., data=dfilt)
summary(model)
summary.lm(model)$r.squared

# Avaliar normalidade do resíduo do modelo

shapiro.test(model$residuals)

qp2<-qqPlot(model$residuals, 
            pch = 20, 
            las = 1,
            ylab = "Residuals")

# Avaliar homocedasticidade

fligner.test(sdr ~ metodo, data = dfilt)
plot(x = model$fitted.values, y = model$residuals)

# Comparações todos x todos

library(multcomp)
mc <- glht(model, linfct = mcp(metodo = "Tukey"))
summary(mc)
plot(mc, cex.axis = 0.5, cex = 1)

