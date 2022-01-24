library(dplyr)

# Cargamos el archivo de datos
data <- read.csv("results/uniprot1.csv", header=T, sep=";")

# Filtramos las entradas y cambiamos el Entry name 
data2 <- dplyr::filter(data, grepl("HUMAN",Entry.Name))
data2$Entry.Name2 <- sapply(data2$Entry.Name, function(i) gsub("_HUMAN", "", i))

