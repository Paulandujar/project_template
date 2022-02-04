# Red de los targets de SARS-CoV2
## Participantes
  - Irene Romero Granados
  - Paula Andújar Zambrano
  - Rosario García Morales
  - Soledad del Castillo Carrera
  
## Como ejecutar el proyecto
### Instalar R
Para ejecutar el proyecto es necesario tener instalado R. En el caso de que ya lo tenga instalado, pase al siguiente paso. 
En el caso contrario, acceda al siguiente enlace https://cran.r-project.org/ y descargue la última versión disponible para su sistema operativo.
### Clonar el repositorio
Para poder ejecutar el proyecto, es necesario clonarlo en el equipo. Dírigase a la carpeta en la que lo quiera clonar y ejecute desde la terminal el siguiente comando:
``` 
git clone https://github.com/Paulandujar/project_template.git 
```
### Ejecución del proyecto
Una vez clonado el repositorio, desde la terminal es necesario ejecutar el archivo **launch.sh**, que a su vez ejecutará al **setup.sh** si es necesario.
> No debería de tener problemas de administrador a la hora de ejecutar dichos archivos. En el caso de que los tuviera, introduzca en la terminal la siguiente línea de comando: ``` chmod 755 setup.sh ```. Esto lo que hace es otorgarle permisos a la persona que lo va a ejecutar. En concreto, con chmod 755 se les está dando permiso de lectura y ejecución a todos los usuarios excepto al propietario que tiene además el permiso de edición.
> 
Para ello, introduzca el siguiente comando en la terminal: 
> ¡¡No olvide que debe ejecutarse desde el directorio en el que se encuentra el archivo!!
> 
``` 
./launch.sh
```
Desde el launch se instalarán las librerías necesarias para ejecutar el proyecto y, cuando estas estén instaladas, se ejecutará el archivo **covid_network.R**, en el cuál se encuentra el análisis realizado al conjunto de datos elegido.
