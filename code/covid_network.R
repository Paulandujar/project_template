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
png("results/string_hits.png")
string_db$plot_network(hits)
dev.off()

# Primera capa de la red
primer.vecino <- (neighbors(graph = string_network, v = V(covid_network)$name, mode = "all"))$name
hits.network <- string_db$get_subnetwork(unique(c(V(hits.network)$name, primer.vecino)))
cl <- components(hits.network)
borra.nodos <- names(cl$membership[cl$membership!=1]) 
hits.network <- delete_vertices(hits.network, borra.nodos)
names <- gsub("9606.ENSP00000", "", V(hits.network)$name)
V(hits.network)$name <- names

# Graficamos la red
png("results/hits.network.png")
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
png("results/top4_clusters.png")
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}
dev.off()


#### ROBUSTEZ
############### FUNCIONES ###############
robustness.random2 <- function(grafo, measure=degree){
  q <- seq(from=0.01,to=1,by=0.01)
  g <- grafo
  S <- max(components(grafo)$csize)/vcount(grafo)
  contador <- S
  removalset <- NULL
  for(i in q){
    if(contador > 0.05 & length(removalset) < vcount(g)/2){
      removalset <- sample(x = V(g)$name, size = 10, replace = F)
      g <- delete.vertices(graph = g, v = removalset)
      S <- c(S, max(components(g)$csize)/vcount(grafo))
      
      contador <- max(components(g)$csize)/vcount(grafo)
    }
  }
  x <- as.numeric(q[1:length(S)])
  x <- x[!is.na(x)]
  y <- as.numeric(S)
  id <- order(x)
  return(sum(diff(x[id])*rollmean(y[id],2)))
} 
heter <- function(network) {
  require(igraph)
  return(var(degree(network))/mean(degree(network))
  )
}

sequential.attacks.targeted <- function(grafo, measure=degree){
  q <- seq(from=0,to=1,by=0.01)
  g <- grafo
  S <- max(components(grafo)$csize)/vcount(grafo)
  contador <- S
  removalset <- NULL
  v <- vcount(grafo)
  s <- max(components(grafo)$csize)
  for(i in q){
    if(max(components(g)$csize)/vcount(grafo) >0.05){
      removalset <- names(sort(degree(g),decreasing = T)[1:(i*vcount(g))])
      g <- delete.vertices(graph = g, v = removalset)
      S <- c(S, max(components(g)$csize)/vcount(grafo))
      v <- c(v, vcount(g))
      s <- c(s, max(components(g)$csize))
      contador = max(components(g)$csize)/vcount(grafo)
    }
    
  }
  S.vs.q <- tbl_df(data.frame(cbind(q[1:length(S)],S,s,v)))
  names(S.vs.q) <- c("q", "S", "s", "v")
  return(S.vs.q)
}


sequential.attacks.random <- function(grafo, measure=degree){
  q <- seq(from=0,to=1,by=0.01)
  g <- grafo
  S <- max(components(grafo)$csize)/vcount(grafo)
  contador <- S
  removalset <- NULL
  v <- vcount(grafo)
  s <- max(components(grafo)$csize)
  for(i in q){
    if(contador > 0.05 & length(removalset) < vcount(g)/2){
      removalset <- sample(x = V(g)$name, size = 10, replace = F)
      g <- delete.vertices(graph = g, v = removalset)
      S <- c(S, max(components(g)$csize)/vcount(grafo))
      v <- c(v, vcount(g))
      s <- c(s, max(components(g)$csize))
      contador <- max(components(g)$csize)/vcount(grafo)
    }
    
  }
  S.vs.q <- tbl_df(data.frame(cbind(q[1:length(S)],S,s,v)))
  names(S.vs.q) <- c("q", "S", "s", "v")
  return(S.vs.q)
}



robustness.targeted2 <- function(grafo, measure=degree){
  q <- seq(from=0.01,to=1,by=0.01)
  g <- grafo
  S <- max(components(grafo)$csize)/vcount(grafo)
  contador <- S
  removalset <- NULL
  for(i in q){
    if(max(components(g)$csize)/vcount(grafo) >0.05){
      removalset <- names(sort(degree(g),decreasing = T)[1:(i*vcount(g))])
      g <- delete.vertices(graph = g, v = removalset)
      S <- c(S, max(components(g)$csize)/vcount(grafo))
      contador = max(components(g)$csize)/vcount(grafo)
    }
  }
  x <- as.numeric(q[1:length(S)])
  y <- as.numeric(S)
  id <- order(x)
  return(sum(diff(x[id])*rollmean(y[id],2)))
} 

# plot the graph
plot_graph <- function(graph){
  V(graph)$label <- NA
  V(graph)$name <- NA
  V(graph)$size <- degree(graph)/5
  E(graph)$edge.color <- "gray80"
  V(graph)$color <- "tomato"
  # graph_attr(graph, "layout") <- layout_with_lgl
  graph_attr(graph, "layout") <- layout.kamada.kawai
  plot(graph)
  
}

png(file="results/covid_network_graph.png")
plot_graph(covid_network)
dev.off()

# Calculo de la robustez de nuestra red
robusezRamdom <- robustness.random2(covid_network)
robustezDirigida <- robustness.targeted2(covid_network)
robustez_ambas <- rbind(robusezRamdom,robustezDirigida)
write.csv(robustez_ambas, file="results/RobustezRedCovid.csv")


# Ataques dirigidos
covid_AttackTargeted <- sequential.attacks.targeted(covid_network)
covid_AttackTargeted$attack <- rep("targeted")

# Ataques aleatorios

covid_AttackRandom <- sequential.attacks.random(covid_network)
covid_AttackRandom$attack <- rep("random")

attack <- rbind(covid_AttackTargeted,covid_AttackRandom)
write.csv(attack, file="results/AtaquesRedCovid.csv")
