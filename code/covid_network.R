library(igraph)
library(dplyr)
library(ggplot2)
library(zoo)
library(STRINGdb)
library(linkcomm)

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
hits.network <- string_db$get_subnetwork(unique(c(V(covid_network)$name, primer.vecino)))
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

###### GRADO DE DISTRIBUCION

# Obtengo el grado de la red
deg <- degree(hits.network, mode="all")

#Calculo el grado de distribucion y lo grafico
deg.dist <- degree_distribution(hits.network,cumulative = F, mode="all")
png("results/degree_distribution.png")
plot(deg.dist, pch=20, cex=1.2, col="orange", xlab="Degree", ylab="Frequency")
dev.off()


####coeficiente de agrupamiento
clusteringCoef<-transitivity(hits.network, type = "local", isolates = c("zero"))
cc_medio <- transitivity(hits.network, type = "average")
cat("El Coeficiente de Agrupaci�n medio es de:: ",cc_medio)
clusydegree <- data.frame(degree = deg, Clust_Coef = clusteringCoef, row.names = NULL)
png("results/coeficiente_agrupamiento.png")
plot(clusydegree,pch=20, cex=1.2, col="orange", xlab="Degree", ylab="Clustering Coeffient" )
dev.off()


###### distancia
distance <- distance_table(hits.network)
denominador <- 1:length(distance$res)
sumatorio <- sum(distance$res)
distance$res <- distance$res/s
cat("La distacnia media es de: ",mean_distance(hits.network))
png("results/distancia.png")
plot(distance$res,pch=20, cex=1.2, col="blue", xlab="Distancia", ylab="Pd")
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


###################################
# graficamos el grafo
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

png(file="results/hits.network_graph.png")
plot_graph(hits.network)
dev.off()

# Calculo de la robustez de nuestra red
robusezRamdom <- robustness.random2(hits.network)
robustezDirigida <- robustness.targeted2(hits.network)
robustez_ambas <- rbind(robusezRamdom,robustezDirigida)
write.csv(robustez_ambas, file="results/RobustezRedCovid.csv")


# Ataques dirigidos
covid_AttackTargeted <- sequential.attacks.targeted(hits.network)
covid_AttackTargeted$attack <- rep("targeted")

# Ataques aleatorios

covid_AttackRandom <- sequential.attacks.random(hits.network)
covid_AttackRandom$attack <- rep("random")

# Juntamos ambos ataques
attack <- rbind(covid_AttackTargeted,covid_AttackRandom)
write.csv(attack, file="results/AtaquesRedCovid.csv")


png(file = 'results/sequential_attacks.png')
ggplot(attack, aes(x=q, y=S, color=attack)) + geom_point(alpha=.4, size=2) +
  theme_bw() +
  theme(plot.background = element_blank(),  panel.grid.minor = element_blank(),plot.title=element_text(size=15)) +
  geom_line() +
  ggtitle("Random vs. Targeted sequential attacks") +
  ylab("S(q)") +
  xlab("q")
dev.off()

#### LINKED COMMUNITIES

# Se guarda la informaci?n en dataframe para usar linkcomm
covid_df = igraph::as_data_frame(hits.network, what="edges") 

# Obtenemos las comunidades vinculadas
covid_lc <- getLinkCommunities(covid_df, hcmethod = "single")

# Guardamos los resultamos
png("results/covid_lc_summary.png")
plot(covid_lc, type = "summary")
dev.off()
png("results/covid_lc_dend.png")
plot(covid_lc, type = "dend")
dev.off()

# Tama?o de los clusters
png("results/lc_larger_clusters.png")
par(mfrow = c(2,1))
pie(covid_lc$clustsizes[covid_lc$clustsizes > 8], radius = 1, main = "Tama?os de las comunidades m?s grandes")
barplot(covid_lc$clustsizes[covid_lc$clustsizes > 8], xlab = "Comunidades", ylab = "Tama?o (num genes)")
par(mfrow = c(1,1))
dev.off()

# Cluster por comunidad/modularidad
png("results/clusters_modularity.png")
plot(covid_lc, type = "commsumm", summmary = "mod")
dev.off()
cm <- getCommunityConnectedness(covid_lc, conn = "modularity")
print(head(sort(cm, decreasing = TRUE)))

# modularidad de las comunidades
community.connectedness <- getCommunityConnectedness(covid_lc, conn = "modularity") 
plot(covid_lc, type = "commsumm", summary = "modularity")

# Escogemos un cluster al azar (en este caso el 54) y lo graficamos
png(file="results/cluster54_graph.png")
plot(covid_lc, type = "graph", clusterids = 54, vlabel=FALSE)
dev.off()

png("results/6mejores_clusters_modularity.png")
plot(covid_lc, clusterids = c(78, 22, 54, 41, 21, 77), type = "commsumm", summary = "modularity")
dev.off()

# Cambio en el dise?o de los gr?ficos
png("results/lc_hits.network_layout_fruchterman.reingold.png")
plot(covid_lc, type = "graph", layout = layout.fruchterman.reingold, ewidth = 2, vlabel.cex = 0.5)
dev.off()
png("results/lc_hits.network_layout_spencer.circle.png")
plot(covid_lc, type = "graph", layout = "spencer.circle", ewidth = 2, vlabel.cex = 0.5)
dev.off()

# Mostramos solo los nodos que pertenecen a 10 o m?s comunidades
png("results/hits.network_layout_fruchterman.reingold_shownodesin_10.png")
plot(covid_lc, type = "graph", shownodesin = 10, node.pies = TRUE, ewidth = 2, vlabel.cex = 0.5)
dev.off()
png("results/lc_hits.network_layout_spencer.circle_shownodesin_10.png")
plot(covid_lc, type = "graph", layout = "spencer.circle", shownodesin = 10, ewidth = 2, vlabel.cex = 0.5)
dev.off()

# visualizamos la pertenencia a la comunidad de nodos para los nodos m?s conectados
png("results/members.png")
plot(covid_lc, type = "members")
dev.off() 


# obtener comunidades completamente anidadas dentro de una comunidad más grande de nodos
getAllNestedComm(covid_lc)

# comprobar comunidades anidadas
getAllNestedComm(covid_lc)[1]

#Elegimos las comunidades que son independientes de otras, en este caso 16 y 146

png("results/nested_comm.png")
plot(covid_lc, type = "graph", clusterids = c(16,146))
dev.off() 

#### ENRIQUECIMIENTO FUNCIONAL
enriquecimientoFuncional <- function(cluster) {
  # Observamos la representaci?n de los genes del cluster
  png(paste("results/cluster_", cluster, ".png", sep = ""))
  plot(covid_lc, type = "graph", clusterids = cluster)
  dev.off()
  
  #Extraemos los ids
  datos_ids <- paste0("9606.ENSP00000", linkcomm::getNodesIn(covid_lc, clusterids = cluster, type = "names"))
  
  # Enriquecimiento con GO
  enriGo <- string_db$get_enrichment(datos_ids, category = "Process")
  print(enriGo)
  write.csv(enriGo[, -c(1,7,8,9)], paste("results/funcionesbiologicas_GO_cluster", cluster, ".csv", sep = ""))
  
  # Enriquecimiento con KEGG
  enriKEGG <- string_db$get_enrichment(datos_ids, category = "KEGG")
  print(enriKEGG)
  write.csv(enriKEGG[, -c(1,7,8,9)], paste("results/funcionesbiologicas_KEGG_cluster", cluster, ".csv", sep = ""))
}

# Enriquecimiento funcional para cluster 112 por ser el de mayor tamaño
enriquecimientoFuncional(112)

# Los dos clusters con mayor modularidad son el 78 y el 22
enriquecimientoFuncional(78)
enriquecimientoFuncional(22)

