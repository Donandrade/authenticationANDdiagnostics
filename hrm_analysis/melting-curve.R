install.packages("tidyr")
if (!require("processx")) install.packages("processx")


library(ggplot2)
library(factoextra)
library(dendextend)
library(cluster)
library(ape)
library(circlize)
library(qpcR)
library(tidyverse)
library(plotly)

###################################################################################
###################################################################################
###################################################################################

data <- read.csv("dados_desconhecidos.txt", header = TRUE, sep = "\t", dec = ",")

data2 <- c()

for(i in 1:nrow(data)){
  if(data[i,1] > 73 & data[i,1] < 82){
    tmp <- c(data[i,])
    data2 <- as.data.frame(rbind(as.data.frame(data2), as.data.frame(tmp)))
  }
}

data <- data2

x <- meltcurve(data %>% select(1, 2), norm = TRUE)

c = x[[1]][1]

c = cbind(c, x[[1]][3])

d = x[[1]][1]
d = cbind(d, x[[1]][["Fluo"]])

for(i in 3:ncol(data)){
  x <- meltcurve(data %>% select(1, i), norm = TRUE, span.smooth = 0.05)
  c = cbind(c, x[[1]][3])
  d = cbind(d, x[[1]][["Fluo"]])
  
  print (i)
}

colnames(d) <- colnames(data)

normalized = d[,2:ncol(d)]


####################################################
####################################################
####################################################

#para calcular a dissimilaridade

matrix <- c()
tmpDiss <- c()
dissimilarity <- c()

for(a in 1:ncol(normalized)){
  for(j in 1:ncol(normalized)){
    tmp <- as.matrix((normalized[,a]-normalized[,j])^2)
    matrix <- as.matrix(cbind(matrix, tmp))
  }
  for(x in 1:ncol(matrix)){
    diss <- 1-(1.05^(sum(matrix[,x])*-0.02))
    diss <- diss*100
    tmpDiss <- as.matrix(cbind(tmpDiss, diss))
  }
  dissimilarity <- as.matrix(rbind(dissimilarity, tmpDiss))
  tmpDiss <- c()
  matrix <- c()
}


colnames(dissimilarity) <- colnames(normalized)
rownames(dissimilarity) <- colnames(normalized)


write.csv(dissimilarity, file = "dissimilaridade.txt", sep = "\t")
hc <- hclust(as.dist(dissimilarity), method = "ward.D2")

plot(hc)

hcd <- as.dendrogram(hc)



####################################################
####################################################
####################################################

#para calcular a Similaridade

matrix <- c()
tmpDiss <- c()
dissimilarity <- c()

for(a in 1:ncol(normalized)){
  for(j in 1:ncol(normalized)){
    tmp <- as.matrix((normalized[,a]-normalized[,j])^2)
    matrix <- as.matrix(cbind(matrix, tmp))
  }
  for(x in 1:ncol(matrix)){
    diss <- 1.05^(sum(matrix[,x])*-0.02)
    diss <- diss*100
    tmpDiss <- as.matrix(cbind(tmpDiss, diss))
  }
  dissimilarity <- as.matrix(rbind(dissimilarity, tmpDiss))
  tmpDiss <- c()
  matrix <- c()
}


colnames(dissimilarity) <- colnames(normalized)
rownames(dissimilarity) <- colnames(normalized)


write.csv(dissimilarity, file = "dissimilaridade.txt", sep = "\t")
hc <- hclust(as.dist(dissimilarity), method = "ward.D2")

plot(hc)

hcd <- as.dendrogram(hc)

####################################################
####################################################
####################################################


dend <- dissimilarity %>%
  as.dist() %>%
  hclust() %>%
  as.dendrogram()


#library(data.table)

#df <- as.data.frame(t(df))

dend <- dend %>%
  set("branches_lwd", c(10, 10, 10)) %>%
  set("branches_col", c(2)) %>%
  set("labels_colors") %>%
  set("labels_cex", c(5))

labels2 <- as.data.frame(labels(dend))

labels2
write.csv(labels2, file = "labels")
#setDT(df)

#setcolorder(df, as.character(labels2$`labels(hcd)`))

#df <- as.data.frame(t(df))

color_order <- read.csv("labels", header = TRUE, sep = "\t")

labels_colors(dend) <- color_order$V2

circlize_dendrogram(dend, bg.border = NA,
                    labels_track_height = 0.3,
                    dend_track_height = 0.5) 


###########################################################################
###########################################################################
###########################################################################
###########################################################################
data <- read.csv("data_media.txt", header = TRUE, sep = "\t", dec = ",")

data2 <- c()

for(i in 1:nrow(data)){
  if(data[i,1] > 70 & data[i,1] < 82){
    tmp <- c(data[i,])
    data2 <- as.data.frame(rbind(as.data.frame(data2), as.data.frame(tmp)))
  }
}

data <- data2

x <- meltcurve(data %>% select(1, 2), norm = TRUE)

c = x[[1]][1]

c = cbind(c, x[[1]][3])

d = x[[1]][1]
d = cbind(d, x[[1]][["Fluo"]])

for(i in 3:ncol(data)){
  x <- meltcurve(data %>% select(1, i), norm = TRUE, span.smooth = 0.05)
  c = cbind(c, x[[1]][3])
  d = cbind(d, x[[1]][["Fluo"]])
  
  print (i)
}

colnames(c) = colnames(data)

colnames(d) = colnames(data)

filter <- c()

for(i in 1:nrow(c)){
  if(c[i,1] > 70 & c[i,1] < 82){
    tmp <- c(c[i,])
    print(tmp)
    filter <- as.data.frame(rbind(as.data.frame(filter), as.data.frame(tmp)))
  }
}



rownames(filter) <- filter[,1]

colnames(filter) = colnames(data)

str(filter)

a <- plot_ly(x = filter$Temperature) %>%
  add_lines(y = filter$Robusta_1, color = I("#ff1a1a"), name ="Mp_1441_1") %>%
  add_lines(y = filter$Robusta_2, color = I("#ff1a1a"), name ="Mp_1441_2") %>%
  add_lines(y = filter$Robusta_3, color = I("#ff1a1a"), name ="Mp_1441_3") %>%
  add_lines(y = filter$Robusta_4, color = I("#ff1a1a"), name ="Mp_1441_4") %>%
  add_lines(y = filter$Robusta_5, color = I("#ff1a1a"), name ="Mp_4145_1") %>%
  add_lines(y = filter$Robusta_6, color = I("#ff1a1a"), name ="Mp_4145_2") %>%
  add_lines(y = filter$Robusta_7, color = I("#ff1a1a"), name ="Mp_4145_3") %>%
  add_lines(y = filter$Robusta_8, color = I("#ff1a1a"), name ="Mp_4145_4") %>%
  add_lines(y = filter$Robusta_9, color = I("#ff1a1a"), name ="Mp_60d_1") %>%
  add_lines(y = filter$Robusta_10, color = I("#ff1a1a"), name ="Mp_60d_2") %>%
  add_lines(y = filter$Robusta_11, color = I("#ff1a1a"), name ="Mp_F18_1") %>%
  add_lines(y = filter$Robusta_12, color = I("#ff1a1a"), name ="Mp_F18_2") %>%
  add_lines(y = filter$Robusta_13, color = I("#ff1a1a"), name ="Mp_F25_1") %>%
  add_lines(y = filter$Robusta_14, color = I("#ff1a1a"), name ="Mp_F25_2") %>%
  add_lines(y = filter$Robusta_15, color = I("#ff1a1a"), name ="Mp_F27_2") %>%
  add_lines(y = filter$Robusta_16, color = I("#ff1a1a"), name ="Mp_F28_1") %>%
  add_lines(y = filter$Robusta_17_1, color = I("#ff1a1a"), name ="Mp_F28_2") %>%
  add_lines(y = filter$Robusta_17_2, color = I("#ff1a1a"), name ="Mp_F28_2") %>%
  add_lines(y = filter$Conilon_1, color = I("#1a75ff"), name ="Mp_F30_1") %>%
  add_lines(y = filter$Conilon_2, color = I("#1a75ff"), name ="Mp_F30_2") %>%
  add_lines(y = filter$Conilon_3, color = I("#1a75ff"), name ="Mp_F31_1") %>%
  add_lines(y = filter$Conilon_4, color = I("#1a75ff"), name ="Mp_F31_2") %>%
  add_lines(y = filter$Conilon_5, color = I("#1a75ff"), name ="Mp_F36_1") %>%
  add_lines(y = filter$Conilon_6, color = I("#1a75ff"), name ="Mp_F36_2") %>%
  add_lines(y = filter$Conilon_5, color = I("#1a75ff"), name ="Mp_F40_1") %>%
  add_lines(y = filter$Conilon_6, color = I("#1a75ff"), name ="Mp_F40_2") %>%
  add_lines(y = filter$Conilon_9, color = I("#1a75ff"), name ="Mr_EqF1_1") %>%
  add_lines(y = filter$Conilon_10, color = I("#1a75ff"), name ="Mr_EqF1_2") %>%
  add_lines(y = filter$Conilon_11, color = I("#1a75ff"), name ="Mr_EqF2_1") %>%
  add_lines(y = filter$Conilon_12, color = I("#1a75ff"), name ="Mr_EqF2_2") %>%
  add_lines(y = filter$Conilon_13, color = I("#1a75ff"), name ="Mr_EqF3_1") %>%
  add_lines(y = filter$Conilon_14, color = I("#1a75ff"), name ="Mr_EqF3_2") %>%
  add_lines(y = filter$Conilon_15, color = I("#1a75ff"), name ="Mr_EqF5_1") %>%
  add_lines(y = filter$Conilon_16, color = I("#1a75ff"), name ="Mr_EqF5_2") %>%
  add_lines(y = filter$Conilon_17, color = I("#1a75ff"), name ="Mr_EqF6_1") %>%
  add_lines(y = filter$Hibrido_1, color = I("#ff66d9"), name ="Mr_EqF6_1") %>%
  add_lines(y = filter$Hibrido_2, color = I("#ff66d9"), name ="Mr_EqF6_1") %>%
  add_lines(y = filter$Hibrido_3, color = I("#ff66d9"), name ="Mr_EqF6_1") %>%
  add_lines(y = filter$Hibrido_6, color = I("#ff66d9"), name ="Mr_EqF6_1") %>%
  add_lines(y = filter$Hibrido, color = I("#ff66d9"), name ="Mr_EqF6_1") %>%
  add_lines(y = filter$Hibrido.1, color = I("#ff66d9"), name ="Mr_EqF6_1")
            

a

orca(a, "Fig_diffPlt_tmp.pdf")


