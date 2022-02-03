#Creacion del directorio para las librerias

mkdir software software/deps_r
chmod 755 software/deps_r

#Instalacion de paquetes necesarios

brew install curl libssl-dev libcurl4-openssl-dev libxml2-dev

#Instalacion de librerias R

PKGS_CRAN="'zoo' 'igraph' 'linkcomm' 'data.table' 'BiocManager' 'devtools' 'dplyr'"
PKGS_BM="'STRINGdb'"

for PKG_CRAN in $PKGS_CRAN
    do
        Rscript -e 'install.packages('$PKG_CRAN', lib="software/deps_r" )'
        Rscript -e 'library('$PKG_CRAN')'
    done

for PKG_BM in $PKGS_BM
    do
        Rscript -e 'BiocManager::install('$PKG_BM', lib="software/deps_r" )'
        Rscript -e 'library('$PKG_BM')'
    done
    
#Añadir el path a la variable libPaths

Rscript -e '.libPaths(c(.libPaths(),paste(getwd(),"software/deps_r",sep=""))))'
echo "*********************** INSTALLATION SUCCESSFULLY *****************************"
echo "Librerias cargadas en el directorio /software/deps_r"