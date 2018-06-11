# skgen

Requerimientos:

*gsl:Libreria numerica para C. Se puede instalar desde repositorios de ubuntu.

*Libxc:Libreria de funcionales http://octopus-code.org/wiki/Libxc

*boost: 

*dftb+ y dptools. https://www.dftbplus.org/

Compilacion:

Ir a la carpeta scr/ y modificar el makefile, donde dice XCDIR escribir la direccion donde se encuentra la libreria LibXC.
gsl llama a BLAS, en mi configuracion uso openblas,hay que modificar eso en el makefile.

Correr:

Es necesario un archivo donde estan los parametros para las tablas sk y luego archivos individuales para cada atomo que aparece en las tablas.


Para correr :  skgen archivo


Ejemplos:
Crear archivo skf y comparar con bandas de DFT (QE).
Para graficar bandas es necesario dftb+ y dptools. Reemplazar en run.sh donde se encuentre el ejecutable de dftb+.

BN:

skgen test

Para graficar: bash run.sh

W:

skgen test

Para graficar: bash run.sh

