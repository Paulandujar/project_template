#Creacion del directorio para las librerias

chmod 777 setup.sh
mkdir software software/deps_r
chmod 755 software/deps_r


#Instalacion de paquetes necesarios

# brew install libffi
# sudo apt-get install curl libssl-dev libcurl4-openssl-dev libxml2-dev

#Instalacion de librerias R
#"'BiocManager''zoo' 'igraph' 'linkcomm' 'data.table' 'devtools' 'dplyr'"
PKGS_CRAN="'BiocManager' 'RCurl' 'zoo' 'igraph' 'linkcomm' 'data.table' 'devtools' 'dplyr'"
PKGS_BM="'STRINGdb'"

for PKG_CRAN in $PKGS_CRAN
    do
        Rscript -e 'install.packages('$PKG_CRAN', repos="https://cran.rstudio.com/", lib="software/deps_r", force = TRUE )'
        Rscript -e 'library('$PKG_CRAN', lib.loc="software/deps_r")'
    done

#Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.14",ask=FALSE)'

for PKG_BM in $PKGS_BM
    do
        Rscript -e 'library("BiocManager", lib.loc = "software/deps_r") ; install('$PKG_BM', lib="software/deps_r", force = TRUE )'
        Rscript -e 'library('$PKG_BM', lib.loc="software/deps_r")'
    done
    
#Añadir el path a la variable libPaths

Rscript -e '.libPaths(c(.libPaths(),paste(getwd(),"software/deps_r",sep="/"))))'
echo "*********************** INSTALLATION SUCCESSFULLY *****************************"
echo "Librerias cargadas en el directorio /software/deps_r"
