library(dplyr)

# Cargamos el archivo de datos
data <- read.csv("results/uniprot1.csv", header=T, sep=";")

# Filtramos las entradas y cambiamos el Entry name 
data2 <- dplyr::filter(data, grepl("HUMAN",Entry.Name))
data2$Entry.Name2 <- sapply(data2$Entry.Name, function(i) gsub("_HUMAN", "", i))

# Mapeamos con STRING
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )
data_mapped <- string_db$map(data2, "Entry", removeUnmappedRows = TRUE )

# Guardamos la red de proteinas
hits <-data_mapped$STRING_id
png("01_string_hits.png")
string_db$plot_network(hits)
dev.off()
