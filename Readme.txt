
%%%%%%%%%%%%%%%%%%%%%%%%%%
Name, version, description, and/or features of the program.
System requirements.
Install, uninstall, configuration, and operating instructions.
Files list.
Credit (authors), acknowledgments, contact information, and copyright.
Known bugs and a change log.
%%%%%%%%%%%%%%%%%%%%%%%%%%


WRail v0.1
-----

El programa WRail simula la propagacion de una onda ultrasonica a traves de un rail. Dado un rail especificado en un fichero externo con los parametros de fabricacion, se genera una malla (triangular) que discretiza la seccion transversal. A partir de esta malla se hallan las curvas de dispersion y escogiendo un punto sobre una de ellas se visualiza la propagacion de la onda correspondiente sobre el rail.


1. Cargar rail
--------------

1.1 En el menu "File" se carga un rail desde un fichero externo que contiene los parametros de fabricacion del rail (standard).

1.2 En el cuadro "Mesh Refinement" se especifica la densidad de la malla que discretiza la seccion transversal del rail con dos posibilidades: Medium o High. La malla se visualiza en la region inferior izquierda "Section".

1.3 En la region inferior derecha "Rail 3D" se visualiza el rail como la extension (como prisma) de una seccion transversal. En el panel adjunto "Modify Section" se seleccionan el numero de secciones (Sections Number) y el espaciamiento ("Spacing") entre ellas.

1.3.1 Con los botones "Zoom In" y "Zoom Out" se puede escalar el tamaño del rail, asi como con el scroll del raton. Con las flechas de movimiento del teclado se puede mover el rail a diferentes posiciones.


2. Hallar curvas de dispersion
------------------------------

2.1 En el panel superior derecho "Graph Properties" se seleccionan el numero maximo de ondas y el numero de curvas de dispersion. 

2.2 Con el boton "Compute" se muestran las curvas de dispersion en la region superior. 

2.3 Habilitando la opcion "Show Curves" en el menu "Display" (por defecto) se muestran las curvas de manera continua. Con dicha opcion deshabilitada se muestran los puntos sobre las curvas que corresponden a valores de frecuencia/velocidades de fase.


3. Animacion de la propagacion de una onda
------------------------------------------

3.1 Seleccionando un punto sobre una curva de dispersion en la region superior, se escoge una onda cuya propagacion se visualiza con el boton "Animate" en el cuadro "Rail 3D". 

3.2 Al animarse la propagacion en el rail, el boton cambia su nombre a "Stop Animation" y permite entonces detener la animacion.
