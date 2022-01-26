library(igraph)
library(dplyr)
library(ggplot2)
library(zoo)
library(STRINGdb)

# Cargamos el archivo de datos
data <- read.csv("code/data/uniprot1.csv", header=T, sep=";")

# Filtramos las entradas y cambiamos el Entry name 
data2 <- dplyr::filter(data, grepl("HUMAN",Entry.Name))
data2$Entry.Name2 <- sapply(data2$Entry.Name, function(i) gsub("_HUMAN", "", i))

# Mapeamos con STRING
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )
data_mapped <- string_db$map(data2, "Entry", removeUnmappedRows = TRUE )
data_mapped_string_ids <- data_mapped$STRING_id

# Creamos un objeto igraph
string_network <- string_db$get_graph()

# Obtenemos la red para nuestros datos
covid_network <- string_db$get_subnetwork(data_mapped_string_ids)

# Guardamos los hits de string
hits <-data_mapped$STRING_id
png("results/01_string_hits.png")
string_db$plot_network(hits)
dev.off()

# Primera capa de la red
primer.vecino <- (neighbors(graph = string.network, v = V(hits.network)$name, mode = "all"))$name
hits.network <- string_db$get_subnetwork(unique(c(V(hits.network)$name, primer.vecino)))
cl <- components(hits.network)
borra.nodos <- names(cl$membership[cl$membership!=1]) 
hits.network <- delete_vertices(hits.network, borra.nodos)
names <- gsub("9606.ENSP00000", "", V(hits.network)$name)
V(hits.network)$name <- names

# Graficamos la red
png("results/02_hits.network.png")
plot(hits.network,
     vertex.color = "pink",
     vertex.size = degree(hits.network)/10,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.label.cex = 0.5,
     layout = layout.kamada.kawai
)
dev.off()


#### CLUSTERING

# Graficamos los clusters utilizando el paquete STRINGdb
clustersList <- string_db$get_clusters(data_mapped$STRING_id)
png("results/03_top4_clusters.png")
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}
dev.off()


#### ROBUSTEZ

robustness.random2 <- function(grafo, measure=degree){
  q = seq(from=0.01,to=1,by=0.01)
  g = grafo
  S = max(components(grafo)$csize)/vcount(grafo)
  contador = S
  removalset = NULL
  for(i in q){
    if(contador > 0.05 & length(removalset) < vcount(g)/2){
      removalset <- sample(x = V(g)$name, size = 10, replace = F)
      g <- delete.vertices(graph = g, v = removalset)
      S = c(S, max(components(g)$csize)/vcount(grafo))
      
      contador = max(components(g)$csize)/vcount(grafo)
    }
  }
  x <- as.numeric(q[1:length(S)])
  x <- x[!is.na(x)]
  y <- as.numeric(S)
  id <- order(x)
  return(sum(diff(x[id])*rollmean(y[id],2)))
} 

robustness.random2(covid_network)