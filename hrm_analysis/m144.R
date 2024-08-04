install.packages("tidyr")
if (!require("processx")) install.packages("processx")

install.packages("MBmca")

library(ggplot2)
library(factoextra)
library(dendextend)
library(cluster)
library(ape)
library(circlize)
library(qpcR)
library(tidyverse)
library(plotly)
library(MBmca)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
<<<<<<< HEAD
data <- read.csv("../hrm_data/M144_Raw_Data.txt", header = TRUE, sep = "\t", dec = ".",
=======
data <- read.csv("../data/M144_Raw_Data.txt", header = TRUE, sep = "\t", dec = ".",
>>>>>>> e125b4f (merdeg)
                 check.names = FALSE)

data2 <- c()

for(i in 1:nrow(data)){
  if(data[i,1] > 77 & data[i,1] < 84){
    tmp <- c(data[i,])
    data2 <- as.data.frame(rbind(as.data.frame(data2), as.data.frame(tmp)))
  }
}

names(data2)[startsWith(names(data2), "X")] <- substring(names(data2)[startsWith(names(data2), "X")], 2)

data <- data2

x <- meltcurve(data %>% select(1, 2), norm = TRUE,  span.smooth = 0.05, 
               is.deriv = TRUE, peaklines = TRUE,
               calc.Area = TRUE)

c = x[[1]][1]

c = cbind(c, x[[1]][3])

d = x[[1]][1]
d = cbind(d, x[[1]][["Fluo"]])

for(i in 3:ncol(data)){
  x <- meltcurve(data %>% select(1, i), norm = TRUE,  span.smooth = 0.05, 
                 is.deriv = TRUE, peaklines = TRUE,
                 calc.Area = TRUE)
  c = cbind(c, x[[1]][3])
  d = cbind(d, x[[1]][["Fluo"]])
  
  print (i)
}

colnames(c) = colnames(data)

colnames(d) = colnames(data)

filter <- c()

for(i in 1:nrow(c)){
  if(c[i,1] > 77 & c[i,1] < 84){
    tmp <- c(c[i,])
    print(tmp)
    filter <- as.data.frame(rbind(as.data.frame(filter), as.data.frame(tmp)))
  }
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

hc <- hclust(as.dist(dissimilarity), method = "average")

pdf(file = "M144.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 6) # The height of the plot in inches

plot(hc, main="M144")



########################################
####define the cluster number
########################################

cluster <- rect.hclust(hc, k = 2, border = 2:5)

dev.off()

cluster_data <- matrix(nrow = 0, ncol = 2)  # Initialize an empty matrix

count =1

for (i in cluster) {
  for (x in 1:length(i)) {
    cluster_data <- rbind(cluster_data, c(count, names(i[x])))  # Add row to matrix
  }
  count <- count + 1
}

cluster_data

str(cluster)

write.table(cluster_data, file = "clusters_M144.txt",
            sep = "\t", quote = FALSE,
            col.names =  c("cluster", "indiv"))

#################################
#####color plant build option. Coloring based on clutser
#################################

library(RColorBrewer)

# Defina o número de cores para o número de clusters
num_clusters <- length(cluster)


#paletteM15 <- c("#ff4d4d", "#ffb3b3","#ff0000","#b3ccff", "#80aaff","#4d88ff","#9900FF","#ff944f","green")

#palette <- c( "#ffb3b3","#ff0000","#b3ccff", "#4d88ff")

palette <- c( "#ff0000", "#4d88ff")

##NAO mudar o ncol
sample_cluster_color <- matrix(nrow = 0, ncol = 4)

sample_cluster_color
count <- 1

for(i in cluster){
  hex_color <- palette[count]
  
  for(x in 1:length(cluster[[count]])){
    sample_cluster_color <- rbind(sample_cluster_color,c(count, names(cluster[[count]])[x], hex_color))
  }
  count <- count + 1  
}

sample_cluster_color
####################################################
####################################################
####################################################
colnames(data)
rownames(filter) <- filter[,1]
filter
colnames(filter) = colnames(data)

# Inicialize o gráfico

plot <- plot_ly(x = filter$Temp)

for(i in 1:nrow(sample_cluster_color)){
  plot <- plot %>% add_lines(y = filter[[sample_cluster_color[i, 2]]],
                             color = I(sample_cluster_color[i, 3]),
                             name = sample_cluster_color[i, 2],
                             line = list(width = 3))
}



plot

# Ajuste do layout para aumentar o tamanho das fontes
plot <- plot %>%
  layout(
    title = list(text = 'M144/k2', font = list(size = 20)),
    xaxis = list(title = 'Temperature (ºC)', titlefont = list(size = 20), tickfont = list(size = 20)),
    yaxis = list(title = 'Fluorescence (-dF/dT)', titlefont = list(size = 20), tickfont = list(size = 20)),
    legend = list(font = list(size = 14))
  )


plot

orca(plot, "df_M144.pdf", width=700, height=300)

