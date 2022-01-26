library(igraph)
library(dplyr)
library(ggplot2)
library(zoo)
library(STRINGdb)

# Cargamos el archivo de datos
data <- read.csv("uniprot1.csv", header=T, sep=";")

# Filtramos las entradas y cambiamos el Entry name 
data2 <- dplyr::filter(data, grepl("HUMAN",Entry.Name))
data2$Entry.Name2 <- sapply(data2$Entry.Name, function(i) gsub("_HUMAN", "", i))

# Mapeamos con STRING
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )
data_mapped <- string_db$map(data2, "Entry", removeUnmappedRows = TRUE )
data_mapped_string_ids <- data_mapped$STRING_id

string_network <- string_db$get_graph()
#obtenemos la red para nuestros datps
covid_network <- string_db$get_subnetwork(data_mapped_string_ids )

# Guardamos la red de proteinas
hits <-data_mapped$STRING_id
png("01_string_hits.png")
string_db$plot_network(hits)
dev.off()



# ROBUSTEZ

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
