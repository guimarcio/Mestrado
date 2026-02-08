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

amb <- c('estudio','seminario','cefala','2419')
met <- 'pesq8'
mtd <- 'IVAaux'
flag=0

for (x in 1:4) {
  dfilt_mpos<-dados[dados$metodo==mtd &
                 dados$ambiente==amb[x],]
  
  # selecionar apenas colunas relevantes
  stay <- c("ensaio", met)
  dfilt_mpos<-dfilt_mpos[,names(dfilt_mpos) %in% stay]
  dfilt_mdist<-dfilt_mpos

  # Agregar pela distância das fontes
  dfilt_mpos$ensaio[(dfilt_mpos$ensaio == 'P90') | (dfilt_mpos$ensaio == 'P57') | (dfilt_mpos$ensaio == 'P4')]<-paste0('P-',substr(amb[x],1,1))
  dfilt_mpos$ensaio[(dfilt_mpos$ensaio == 'M90') | (dfilt_mpos$ensaio == 'M57') | (dfilt_mpos$ensaio == 'M4')]<-paste0('M-',substr(amb[x],1,1))
  dfilt_mpos$ensaio[(dfilt_mpos$ensaio == 'L90') | (dfilt_mpos$ensaio == 'L57') | (dfilt_mpos$ensaio == 'L4')]<-paste0('L-',substr(amb[x],1,1))
  dfilt_mpos[, 1]<-as.factor(dfilt_mpos[, 1])
  
  
  # Agregar pela distância dos microfones
  dfilt_mdist$ensaio[(dfilt_mdist$ensaio == 'P90') | (dfilt_mdist$ensaio == 'M90') | (dfilt_mdist$ensaio == 'L90')]<-paste0('0-',substr(amb[x],1,1))
  dfilt_mdist$ensaio[(dfilt_mdist$ensaio == 'L4') | (dfilt_mdist$ensaio == 'M4') | (dfilt_mdist$ensaio == 'P4')]<-paste0('4-',substr(amb[x],1,1))
  dfilt_mdist$ensaio[(dfilt_mdist$ensaio == 'L57') | (dfilt_mdist$ensaio == 'M57') | (dfilt_mdist$ensaio == 'P57')]<-paste0('57-',substr(amb[x],1,1))
  dfilt_mdist[, 1]<-as.factor(dfilt_mdist[, 1])
  
  
  if (flag==0){
    dtot_mpos<-dfilt_mpos
    dtot_mdist<-dfilt_mdist
    flag=1
  }
  else{
    dtot_mpos<-rbind(dtot_mpos,dfilt_mpos)
    dtot_mdist<-rbind(dtot_mdist,dfilt_mdist)
  }
}

# Exploração

data<-dtot_mpos
data2<-dtot_mdist

 
mycolors<- ifelse(levels(data$ensaio)=="L-e"|levels(data$ensaio)=="M-e"|levels(data$ensaio)=="P-e", rgb(0.855,0.647,0.125,1) , 
           ifelse(levels(data$ensaio)=="L-s"|levels(data$ensaio)=="M-s"|levels(data$ensaio)=="P-s", rgb(0.18,0.545,0.341,1),
           ifelse(levels(data$ensaio)=="L-c"|levels(data$ensaio)=="M-c"|levels(data$ensaio)=="P-c", rgb(0.255,0.412,0.882,1) , 
           ifelse(levels(data$ensaio)=="L-2"|levels(data$ensaio)=="M-2"|levels(data$ensaio)=="P-2", rgb(0.698,0.133,0.133,1), 
                  "grey90"))))
mycolors2<- ifelse(levels(data2$ensaio)=="0-e"|levels(data2$ensaio)=="4-e"|levels(data2$ensaio)=="57-e", rgb(0.855,0.647,0.125,1) , 
           ifelse(levels(data2$ensaio)=="0-s"|levels(data2$ensaio)=="4-s"|levels(data2$ensaio)=="57-s", rgb(0.18,0.545,0.341,1),
           ifelse(levels(data2$ensaio)=="0-c"|levels(data2$ensaio)=="4-c"|levels(data2$ensaio)=="57-c", rgb(0.255,0.412,0.882,1) , 
           ifelse(levels(data2$ensaio)=="0-2"|levels(data2$ensaio)=="4-2"|levels(data2$ensaio)=="57-2", rgb(0.698,0.133,0.133,1), 
                  "grey90"))))
mysymbols <- c(16, 17, 15, 18)

par(mar=c(6.3, 5.1, 5.1, 9.1), xpd=TRUE)
boxplot(pesq8 ~ ensaio,col=mycolors, data = data, cex.lab=2,cex.axis=1.5,
        names= c("210","140","75","210","140","75","210","140","75","210","140","75"),
        xlab = 'Distância entre fonte e microfone (cm)',
        ylab = 'PESQ8 (ganho médio)')

legend("topright", legend = c("Estúdio","Seminário","Cefala","2419"), 
       col = c(rgb(0.855,0.647,0.125,1),rgb(0.18,0.545,0.341,1),rgb(0.255,0.412,0.882,1),rgb(0.698,0.133,0.133,1)), 
       bty = "n", pch=15 , pt.cex = 3, cex = 1.5, horiz = FALSE, inset=c(-0.15,0))

par(mar=c(6.3, 5.1, 5.1, 9.1), xpd=TRUE)
boxplot(pesq8 ~ ensaio,col=mycolors2, data = data2, cex.lab=2,cex.axis=1.5,
        names = c("0","4","57","0","4","57","0","4","57","0","4","57"),
        xlab = 'Distância entre microfones (cm)',
        ylab = 'PESQ8 (ganho médio)')
legend("topright", legend = c("Estúdio","Seminário","Cefala","2419"), 
       col = c(rgb(0.855,0.647,0.125,1),rgb(0.18,0.545,0.341,1),rgb(0.255,0.412,0.882,1),rgb(0.698,0.133,0.133,1)), 
       bty = "n", pch=15 , pt.cex = 3, cex = 1.5, horiz = FALSE, inset=c(-0.15,0))


# Gráfico qqplot
library(car)

# ANOVA
model <- aov(pesq8 ~ ., data=dtot_mdist)
summary(model)
summary.lm(model)$r.squared

# Avaliar normalidade do resíduo do modelo

shapiro.test(model$residuals)

qp2<-qqPlot(model$residuals, 
            pch = 20, 
            las = 1,
            ylab = "Residuals")

# Avaliar homocedasticidade

fligner.test(pesq8 ~ ensaio, data = dtot_mdist)
plot(x = model$fitted.values, y = model$residuals)

# Comparações todos x todos

library(multcomp)
mc <- glht(model, linfct = mcp(ensaio = "Tukey"))
summary(mc)
plot(mc, cex.axis = 0.5, cex = 1)

