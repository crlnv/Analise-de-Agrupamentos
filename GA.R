library(readr)
dados = read.table("http://www.uio.no/studier/emner/matnat/math/STK4040/h09/computing/track-records-women.txt", header=T)
View(dados)
dados <- as.matrix(dados) 

library(ggplot2)
library(GGally)

ggpairs(dados)

#Detecção de variáveis atípicas: dados
# Calculando a distância de Mahalonobis 
p.cov <- var(scale(dados)) # standardize first
p.cov <- var(dados)
p.mean <- apply(dados,2,mean)
p.mah <- mahalanobis(dados, p.mean, p.cov)
sort(p.mah)
View(p.mah)


## Vemos o valor discrepante para a categoria CANNED_HA e EGGS do 
## do banco de dados promo
## Podemos remover usando:
remover <- c("COK", "KOR.N", "PNG", "SAM")
## Se tivermos mais variáveis para remover devemos separar cada um por 
## vírgula e colocar aspas em cada uma.
dados <-dados[!(row.names(dados) %in% remover), ]
## Como temos apenas uma variável, poderíamos entrar com ela diretamente
## no lugar de remover


library(mvShapiroTest)
dados_m <- as.matrix(dados)
mvShapiro.Test(dados_m)

#CONTINUAR TUDO DAQUI, RODAR TUDO DE NOVO

#======= LETRA A==================================
d <- dist(dados, method = "euclidean")
dist <- as.matrix(d)
View(dist)

tabela_d <- round(dist[1:54, 1:54], 1) #letra A
View(tabela_d)



#====================LETRA B====================================

#Método hierárquico de variância mínima de Ward ou distância média 
res.hcS <- hclust(d, method = "single")
res.hcS
res.hcC <- hclust(d, method = "complete")
res.hcC


# Calculando a matriz cofenética
## Compara as distâncias efetivamente observadas entre os objetos e
## as distâncias previstas a partir do processo de agrupamento.
res.cophS <- cophenetic(res.hcS)
res.cophC <- cophenetic(res.hcC)
# Correlação entre a distância cofenética e  a distância original
cor(d, res.cophS)
cor(d, res.cophC)

#Teste de Mantel
#H0:rho_cof = 0 H1: rho_cof > 0
#verificar
mantelTest(d,res.hcS)
mantelTest(d,res.hcC)

#Comparando os métodos
hc.S <- hclust(d, method = "single")
hc.S
hc.C <- hclust(d, method = "complete")
hc.C

# Correlação entre a distância cofenética e  a distância original
cor(d, cophenetic(hc.S))
cor(d, cophenetic(hc.C))

#Carregando o pacote factoextra
#install.packages("factoextra")
library("factoextra")
# Obtendo o dendograma
plot.new()
single <- fviz_dend(hc.S, 
          cex = 0.5, # tamanho do rótulo
          k_colors = c("#2E9FDF", "#c066c0", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # cores por grupo
          rect = TRUE, # Adicionar retângulo ao redor dos grupos
          rect_border = c("#2E9FDF", "#cc84cc", "#E7B800", "#FC4E07"),
          rect_fill = TRUE)

fviz_dend(hc.C, cex = 0.5)

#Pelo dendograma poderíamos ter de três a quatro grupos. 
#No estudo de Fader e Lodish (1990) utilizaram quatro agrupamentos. Sendo que os grupos
#foram classificados como:
# 1. Alta penetração com alta frequência de compras (ovos, leite e pão);
# 2. Alta penetração com baixa frequência de compras (condimentos e papelaria);
# 3. Baixa penetração com alta frequência de compras (ração, cigarros e produtos para bebês);
# 4. Baixa penetração com baixa frequência de compras;


# Pode-se utilizar alguns indicadores para auxiliar a escolha do número de agrupamentos. 
#Para calcular esses índices devemos instalar o pacote NbClust
install.packages("NbClust")
library("NbClust")


# Obs.: Podemos instalar mais de uma pacote por vez usando
pkgs <- c("factoextra", "NbClust")
install.packages(pkgs)

plot.new()
nbS <- NbClust(dados, distance = "euclidean", min.nc = 2,
               max.nc = 7, 
               method = "single" , index = "all")


nbC <- NbClust(dados, distance = "euclidean", min.nc = 2,
               max.nc = 7, 
               method = "complete" , index = "all")



#method = NULL deve ser substituído pelo algoritmo de agrupamento
#utilizado (“ward.D”, “ward.D2”, “single”, “complete”, “average”, “kmeans”, etc.)



nbS[["All.index"]]
nbS[["Best.partition"]]
nbS[["Best.nc"]]


#-----------------------------------------------------------------------
#5. Interpretação e validação dos agrupamentos.
#-----------------------------------------------------------------------
#Obtendo os agrupamentos
gS <- cutree(hc.S,k=3)
gC <- cutree(hc.C,k=3)

#Número de membros em cada agrupamento
table(gS)
table(gC)
#Obtendo o nome dos membros no agrupamento 1
rownames(dados)[gS == 1]
rownames(dados)[gS == 2]
rownames(dados)[gS == 3]
rownames(dados)[gS == 4]

rownames(dados)[gC == 1]
rownames(dados)[gC == 2]
rownames(dados)[gC == 3]
rownames(dados)[gC == 4]


#Podemos vizualidar o resultado do agrupamento no dendograma
dendo_complete <- fviz_dend(hc.C, k = 3, cex = 0.5,
                           k_colors = c("#2E9FDF", "#00BB0C", "#E7B800", "#FC4E07"),
                           color_labels_by_k = TRUE, rect = TRUE,
                           main = "Dendograma com Ligação Completa")
                           
          
plot.new()
                   
dendo_single <- fviz_dend(hc.S, k = 3, cex = 0.5,
                          k_colors = c("#2E9FDF", "#00BB0C", "#E7B800", "#FC4E07"),
                          color_labels_by_k = TRUE, rect = TRUE,
                          main = "Dendograma com Ligação Simples")

#==============LETRA C====================================
#Calculado a média em cada grupo
# Utilizando uma função para em g cada grupo g em i
clust.centroid = function(i, dat, gC) {
  ind = (gC == i)
  colMeans(dat[ind,])
}
#
sapply(unique(gC), clust.centroid, dados, gC)
mat<- t(as.matrix(sapply(unique(gC), clust.centroid, dados, gC)))
#Obtendo a média dos valores originais (não padronizados) por grupo
sapply(unique(gC), clust.centroid, dados, gC)
round(sapply(unique(gC), clust.centroid, dados, gC), 1)



#Análise se 3 é um número adequado
fviz_nbclust(dados, kmeans, method = "wss")+
  geom_vline(xintercept = 2, linetype = 3)

#Utilizando algoritmo de k-medias
#kmeans(x, centers, iter.max = 10, nstart = 1)
## Definindo uma semente. Isso permite que o resutlado seja 
## reproduzível, já que a semente interfere no resutlado final
set.seed(123)
km.res3 <- kmeans(dados, 3)
print(km.res)

?kmeans
# Utilizando os valores médios das variáveis em cada grupo
## para o método de k-medias
mat3<- t(as.matrix(sapply(unique(gC), clust.centroid, dados, gC)))

km.res3 <- kmeans(dados, mat3)

print(km.res3)

#-----------------------------------------------------------------
#Interpretação e validação dos agrupamentos
#----------------------------------------------------------------
#Acrescentando a coluna de clusters do k-media nos dados originais
dados.k <- cbind(dados, Grupos=km.res3$cluster)


#Calculando a média dos grupos para os dados originais
aggregate(dados, by=list(cluster=km.res3$cluster), mean)

aggregate(dados, by=list(cluster=km.res3$cluster), mean)
round(aggregate(dados, by=list(cluster=km.res3$cluster), mean),1)


fviz_cluster(km.res3, data = dados,
             palette = c("#2E9FDF", "#00BB0C", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal())

#=============================================================================
#Análise se 4 é um número adequado
fviz_nbclust(dados, kmeans, method = "wss")+
  geom_vline(xintercept = 4, linetype = 3)

#Utilizando algoritmo de k-medias
#kmeans(x, centers, iter.max = 10, nstart = 1)
## Definindo uma semente. Isso permite que o resutlado seja 
## reproduzível, já que a semente interfere no resutlado final
set.seed(123)
km.res4 <- kmeans(dados, 4)
?kmeans
print(km.res4)





#-----------------------------------------------------------------
#Interpretação e validação dos agrupamentos
#----------------------------------------------------------------
#Acrescentando a coluna de clusters do k-media nos dados originais
dados.k <- cbind(dados, Grupos=km.res4$cluster)


#Calculando a média dos grupos para os dados originais
aggregate(dados, by=list(cluster=km.res4$cluster), mean)

aggregate(dados, by=list(cluster=km.res4$cluster), mean)
round(aggregate(dados, by=list(cluster=km.res4$cluster), mean),1)

#
fviz_cluster(km.res4, data = dados,
             palette = c("#2E9FDF", "#00BB0C", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal())


#=============================================================================
#Análise se 2 é um número adequado
fviz_nbclust(dados, kmeans, method = "wss")+
  geom_vline(xintercept = 2, linetype = 3)

#Utilizando algoritmo de k-medias
#kmeans(x, centers, iter.max = 10, nstart = 1)
## Definindo uma semente. Isso permite que o resutlado seja 
## reproduzível, já que a semente interfere no resutlado final
set.seed(123)
km.res2 <- kmeans(dados, 2)
?kmeans
print(km.res2)

?kmeans
# Utilizando os valores médios das variáveis em cada grupo
## para o método de k-medias


#-----------------------------------------------------------------
#Interpretação e validação dos agrupamentos
#----------------------------------------------------------------
#Acrescentando a coluna de clusters do k-media nos dados originais
dados.k <- cbind(dados, Grupos=km.res2$cluster)


#Calculando a média dos grupos para os dados originais
aggregate(dados, by=list(cluster=km.res2$cluster), mean)

aggregate(dados, by=list(cluster=km.res2$cluster), mean)
round(aggregate(dados, by=list(cluster=km.res2$cluster), mean),1)


fviz_cluster(km.res2, data = dados,
             palette = c("#2E9FDF", "#00BB0C", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = TRUE, # Add segments from centroids to items
             repel = TRUE, # Avoid label overplotting (slow)
             ggtheme = theme_minimal())









