
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("igraph")
BiocManager::install("STRINGdb")
BiocManager::install("dplyr")
BiocManager::install("ggplot2")
BiocManager::install("zoo")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
