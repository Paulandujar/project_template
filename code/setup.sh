#Creacion del directorio para las librerias

mkdir software software/deps_r
chmod 755 software/deps_r


#Instalacion de paquetes necesarios

# brew install libffi openssl
# sudo apt-get install curl libssl-dev libcurl4-openssl-dev libxml2-dev

#Instalacion de librerias R

PKGS_CRAN="'BiocManager' 'zoo' 'igraph' 'linkcomm' 'data.table' 'devtools' 'dplyr' 'RCarl'"
PKGS_BM="'STRINGdb'"
R_INSTALL_DIR="'software/deps_r'"

for PKG_CRAN in $PKGS_CRAN
    do
        Rscript -e 'install.packages('$PKG_CRAN', repos="https://cran.rstudio.com/", lib='$R_INSTALL_DIR', force = TRUE )'
        Rscript -e 'library('$PKG_CRAN', lib.loc='$R_INSTALL_DIR')'
    done

Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.14",ask=FALSE)'

for PKG_BM in $PKGS_BM
    do
        Rscript -e 'library("BiocManager", lib.loc = '$R_INSTALL_DIR') ; install('$PKG_BM', lib='$R_INSTALL_DIR', force = TRUE )'
        Rscript -e 'library('$PKG_BM', lib.loc='$R_INSTALL_DIR')'
    done
    
#Añadir el path a la variable libPaths

Rscript -e '.libPaths(c(.libPaths(),paste(getwd(),'$R_INSTALL_DIR',sep="/")))'
echo "Librerias cargadas en el directorio $R_INSTALL_DIR"
