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
dados


# tirar a primeira coluna e a saída indesejada 
drops <- c("X","sir","sar","sdr","pesq8")

dados<-dados[,!names(dados) %in% drops]

aggdados<-with(dados,
               aggregate(x   = pesq8,
                         by  = list(metodo,ensaio,ambiente),
                         FUN = mean))

names(aggdados) <- c("metodo", "ensaio", "ambiente", "pesq8")
for (i in 1:3){
  aggdados[, i] <- as.factor(aggdados[, i])
}
# mudando o primeiro metodo como o principal para a comparação por Dunnett

aggdados$metodo = factor(aggdados$metodo,
                         levels = c("IVAng", "FDICAnat", "fastFDICA", "FDICAaux","fastIVA","IVAaux"),ordered = TRUE)
# Análise Exploratória de Dados

summary(aggdados)


cores = colors = c("red3", "blue2", "orange2", 
                   "green2", "black", "cyan4",
                   "yellow3",'magenta3','cyan')
interaction.plot(dados$ensaio, dados$metodo, dados$pesq8, 
                 main = "Gráfico de Interação",
                 xlab = "Ensaio", ylab = "PESQ8",
                 trace.label = "Método:",
                 col=cores,
                 lty = 'solid',
                 lw=3)

interaction.plot(dados$ensaio, dados$ambiente, dados$pesq8, 
                 main = "Gráfico de Interação",
                 xlab = "Ensaio", ylab = "PESQ8",
                 trace.label = "Ambiente:",
                 col=cores,
                 lty = 'solid',
                 lw = 3)

interaction.plot(dados$metodo, dados$ambiente, dados$pesq8, 
                 main = "Gráfico de Interação",
                 xlab = "Método", ylab = "PESQ8",
                 trace.label = "Ambiente:",
                 col=cores,
                 lty = 'solid',
                 lw = 3)

interaction.plot(dados$metodo, dados$ensaio, dados$pesq8, 
                 main = "Gráfico de Interação",
                 xlab = "Método", ylab = "PESQ8",
                 trace.label = "Ensaio:",
                 col=cores,
                 lty = 'solid',
                 lw = 3)

# Análise do modelo completo

modelcompleto<-aov(pesq8 ~ .^3, data = aggdados)

summary(modelcompleto)

# Modelo de primeira ordem

model1<-aov(pesq8 ~ ., data = aggdados)

summary(model1)

summary.lm(model1)$r.squared

# Gráfico
library(car)

shapiro.test(model1$residuals)

qp1<-qqPlot(model1$residuals, 
            pch = 20, 
            las = 1,
            ylab = "Residuals")

# Análise do modelo alternativo 2 sem efeito combinado metodo:ambiente

model2<-aov(pesq8 ~ .^2, data = aggdados)
summary(model2)
summary.lm(model2)$r.squared


# Avaliar normalidade do resíduo do modelo 2

shapiro.test(model2$residuals)

qp2<-qqPlot(model2$residuals, 
            pch = 20, 
            las = 1,
            ylab = "Residuals")

# Homoscedasticidade fligner killeen
fligner.test(pesq8 ~ interaction(metodo,ensaio),
             data = aggdados)
plot(x = model2$fitted.values, y = model2$residuals)
par(mfrow = c(2, 2))
plot(model2, pch = 20, las = 1)


# CUIDADO! FATORES INTERAGINDO! AVALIAR 


# Análise multifatorial para o modelo 2

library(multcomp)

# Comparação para Algoritmo e o Gráfico de intervalos de confiança para modelo 1

algTM2 <- glht(model2, linfct = mcp(metodo = "Tukey"))
summary(algTM2)
plot(algTM2, cex.axis = 0.5, cex = 1)

# Gráfico dos efeitos para modelo 2

library(effects)
alg.effs2 <- allEffects(model2)
plot(alg.effs2, cex.lab=0.1,cex.axis=0.1, cex=1)
