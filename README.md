# skgen

Requerimientos:

*gsl:Libreria numerica para C. 

*Libxc:Libreria de funcionales http://octopus-code.org/wiki/Libxc

*boost: (se incluye algunos headers en la carpeta include,los necesarios para compilar y correr)

*dftb+ y dptools. https://www.dftbplus.org/

Compilacion:

Ir a la carpeta scr/ y modificar el makefile, donde dice XCDIR escribir la direccion donde se encuentra la libreria LibXC.
gsl llama a BLAS, en mi configuracion uso openblas,hay que modificar eso en el makefile.

Correr:

Es necesario un archivo donde estan los parametros para las tablas sk (paramsk)y luego archivos individuales para cada atomo que aparece en las tablas.

paramsk:

Atoms

atomo1   basis={.. .. ..}   dens={.. .. ..}

atomo2   basis={.. .. ..}   dens={.. .. ..}
.
.
.
atomoN   basis={.. .. ..}   dens={.. .. ..}

/

Cada atomox es el nombre archivo (en general  el simbolo del elemento) ,el segmento " basis={W a r0} " son los parametros de confinamiento para la base (afecta S y H), "dens" para el confinamiento de la densidad (afecta solo H).Tanto basis como dens puede no estar presente, si dens no esta presente los parametros para la densidad son iguales a los de basis, si basis no esta presente el confinamiento es igual a {0 0 0}.
Importante no dejar lineas vacias entre atomoN y '/', ni entre atomos.

De forma optional se pueden especificar parametros para la construccion de las tablas, con el bloque

TableOption={

rmin=

rmax=

step=              

ngrid=             numero de puntos en la grilla (numero entero).

}

no es necesario que esten todas las opciones presentes.

     


atomox:


Relativistic=

GGA=

restart=
save= 

z=

sk_orb=                     Orbitales s p d en ese orden,con los que se calculara las tablas. En caso de que no se necesite algun momento angular reemplazar por 0.

config=

Son necesarias para poder correr z,sk_orb y config.


Para correr :  ./skgen paramsk


Ejemplos:
Para graficar bandas es necesario dftb+ y dptools.

BN:

../../scr/skgen test



W:

../../scr/skgen test

Para graficar: bash run.sh




  



