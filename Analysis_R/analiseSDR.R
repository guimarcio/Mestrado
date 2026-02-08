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

install.packages("ggpattern")

# tirar a primeira coluna e a saída indesejada 
drops <- c("X","sir","sar","pesq8","pesq16")

dados<-dados[,!names(dados) %in% drops]

aggdados<-with(dados,
               aggregate(x   = sdr,
                         by  = list(metodo,ensaio,ambiente),
                         FUN = mean))

names(aggdados) <- c("metodo", "ensaio", "ambiente", "sdr")
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

line_ty <- c("solid", "dashed", "dotted", "dotdash", 
             "longdash", "twodash", "solid", "dashed", "dotted")

pch_ty <- c(1, 2, 3, 4, 8, 9, 15, 16, 17)

#start_add

mtd <- c("fastFDICA", "FDICAnat", "FDICAaux", "fastIVA", "IVAng", "IVAaux")
amb <- c("estudio","seminario","cefala","2419")
res <- matrix(0,ncol = 2, nrow = 24)
full <- 0

for (i in 1:6) {
  for (j in 1:4) {
    full <- full+1
    esp <- aggdados[aggdados$metodo == mtd[i] & aggdados$ambiente == amb[j], 'sdr']
    res[full,1] <- mean(esp)
    res[full,2] <- sd(esp)
  }
}
df <- data.frame(res)
colnames(df) <- c("Média","Desvio Padrão")
sd <- res[,2]

aggdados<-with(dados,
               aggregate(x   = sdr,
                         by  = list(metodo,ambiente),
                         FUN = mean))

colnames(aggdados) <- c("Métodos","Sala","SDR")

# library
library(ggplot2)
library(ggpattern)




# Grouped
p <- ggplot(aggdados, aes(fill=Sala, y=SDR, x=Métodos)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(
    aes(ymin = SDR - sd, ymax = SDR + sd),
    position = position_dodge(width = 0.9),
    width = 0.25,
  ) +
  
  geom_point(
    aes(shape = Sala),
    position = position_dodge(width = 0.9),
    size = 4,          # tamanho dos símbolos
    color = "black"    # imprime bem em preto e branco
  ) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  
  theme(
    axis.text.x = element_text(size = 16),      # tamanho da fonte dos labels do eixo X
    axis.text.y = element_text(size = 16),      # tamanho da fonte dos labels do eixo Y
    axis.title.x =element_blank(),              # tamanho do título do eixo X
    axis.title.y = element_text(size = 16),     # tamanho do título do eixo Y
    legend.text = element_text(size = 16),      # tamanho do texto da legenda
    legend.title = element_text(size = 16),     # tamanho do título da legenda
    axis.ticks = element_line(linewidth = 1),         # espessura dos ticks
    legend.position = "right",            # coloca a legenda no topo (fora do gráfico)
    legend.justification = "top",     # alinha a legenda à direita
    legend.box.just = "right",          # garante o alinhamento no box de legenda
    legend.box = "vertical"           # mantém os itens da legenda na horizontal
  )


p + ylab("SDR (dB)") + 
  scale_fill_manual(values = c(rgb(0.698,0.133,0.133,1),
                               rgb(0.255,0.412,0.882,1),
                               rgb(0.855,0.647,0.125,1),
                               rgb(0.18,0.545,0.341,1)),
                    labels = c("2419", "Cefala", "Estúdio", "Seminário"))

#end_add

interaction.plot(dados$ensaio, dados$metodo, dados$sdr, 
                 main = "Gráfico de Interação",
                 xlab = "Ensaio", ylab = "SDR (dB)",
                 trace.label = "Método:",
                 col=cores,
                 lty = 'solid',
                 lw=3)

interaction.plot(dados$ensaio, dados$ambiente, dados$sdr, 
                 main = "Gráfico de Interação",
                 xlab = "Ensaio", ylab = "SDR (dB)",
                 trace.label = "Ambiente:",
                 col=cores,
                 lty = 'solid',
                 lw = 3)

interaction.plot(dados$metodo, dados$ambiente, dados$sdr, 
                 main = "Gráfico de Interação",
                 xlab = "Método", ylab = "SDR (dB)",
                 trace.label = "Ambiente:",
                 col=cores,
                 lty = 'solid',
                 lw = 3)

par(mar = c(5, 4, 4, 10))   # aumenta a margem direita (último valor)

par(cex.axis = 1.2,   # números dos eixos
    cex.lab  = 1.2,   # labels xlab e ylab
    cex.main = 1.7,   # título
    cex      = 1.2)   # tamanho geral (linhas e símbolos)

interaction.plot(dados$metodo, dados$ensaio, dados$sdr, 
                 main = " ",
                 xlab = " ", ylab = "SDR (dB)",
                 trace.label = "Ensaio:",
                 col=cores,
                 type = 'b',
                 lty = line_ty,
                 pch = pch_ty,
                 legend=FALSE,
                 lw = 2)

par(xpd = FALSE)

grid(lwd = 1,col = 'gray50')

par(xpd = TRUE)  # permite legenda fora

legend("topright",
       inset = c(-0.12, 0),   # comece com -0.15 (fica visível!)
       legend = unique(dados$ensaio),
       col = cores,
       lty = line_ty,
       pch = pch_ty,
       bty = "n",
       y.intersp = 1.6,
       title = "Ensaio:")

# Análise do modelo completo

modelcompleto<-aov(sdr ~ .^3, data = aggdados)

summary(modelcompleto)

# Modelo de primeira ordem

model1<-aov(sdr ~ ., data = aggdados)

summary(model1)

summary.lm(model1)$r.squared

# Gráfico
library(car)

shapiro.test(model1$residuals)

qp1<-qqPlot(model1$residuals, 
            pch = 20, 
            las = 1,
            ylab = "Residuals")

# Análise do modelo alternativo 2 com interações de segunda ordem

model2<-aov(sdr ~ metodo, data = aggdados)
summary(model2)
summary.lm(model2)$r.squared


# Avaliar normalidade do resíduo do modelo 2

shapiro.test(model2$residuals)

qp2<-qqPlot(model2$residuals, 
            pch = 20, 
            las = 1,
            ylab = "Residuals")

# Homoscedasticidade fligner killeen
fligner.test(sdr ~ metodo,
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
