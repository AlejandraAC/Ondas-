# Ondas Acústicas-Electromagnéticas
Este código 2D simula la propagación de ondas viscoelásticas SH y de ondas electromagnéticas del modo TM. 
Para discretizar espacialmente las ecuaciones se utiliza el método de diferencias finitas y para la discretización en tiempo se programa el método de Runge Kutta, ambos de cuarto orden. 

El mismo código sirve para modelar ambos tipos de onda, únicamente se modifican los parámetros de la malla y tendido, de la fuente y del modelo geológico.
Este código fue basado y modificado del código realizado por J. Carcione (2015,p.548) en Fortran 77.
# Pre-requisitos
El código se realizó utilizando el software Matlab R2018b pero puede ser utilizado con cualquier otra versión de Matlab
Para obtener una prueba gratuita del software de 30 días ingrese a la URL "https://es.mathworks.com/campaigns/products/trials/targeted/dan.html"
# Ejecutando las pruebas
Una vez instalado Matlab, se revisan los 4 archivos. 

"SH.m" El main para el ejemplo viscoelástico

"TM.m" El main para el ejemplo electromagnético

"H.m" Subrutina para la ejecución del método de diferencias finitas

"wavelet.m" Subrutina para la ejecución del impulso de onda

Como ejemplo, se explica la estructura del ejemplo electromagnético "TM.m":

  * 1. Parámetros de la Malla y del Tendido
  
  El código se basa en datos reales de una adquisicón GPR (Georadar de Penetración Terrestre) utilizando el modo biestático con polarización H o modo TM.
  Para simular la adquisición se decide el número de dipolos receptores a utilizar y el largo del tendido que cubra el espaciamiento entre ellos. Después se decide el ancho de las fronteras absorbentes y se calcula la dimensión total de la malla numérica "nx(filas)", "nz(columnas)", "dx" y "dz" (intervalos entre nodos).
  
  * 2. Parámetros de la Fuente
  
  Se decide la frecuencia que se utilizará para simular la fuente impulso. En este ejemplo fue de 100 MHz.
  
  * 3. Modelo Heterogéneo
  
  El código funciona para un modelo geológico de una capa y un semi-espacio.
  Se decide cuál es el tipo de roca de cada capa y se ingresa el valor de las propiedades de conductividad[s/m], permeabilidad magnética[H/m], permitividad dieléctrica [F/m] y se obtiene el factor de calidad (Q). 
  
  * 4. Parámetros Absorbentes
  
  El código maneja fronteras absorbentes para evitar reflexiones de onda, para crearlas se utiliza un vector con factores de atenuación que disminuyen la amplitud de la onda conforme ésta hace contacto con los cuatro bordes de la malla.
  
  * 5. Subrutina Wavelet
  
  Se programa la fuente impulso, la cual es una función que genera una ondícula gaussiana (J. Carcione ,2006). Esta subrutina utiliza los datos de frecuencia y del muestreo en el tiempo "dt".
  
  * 6. Creación de los Snapshots
  
  Se programa un ciclo para generar cada uno de los snapshots. Dentro del mismo se generan las fronteras absorbentes y se programa el método de Runge Kutta de cuarto orden, el cual utiliza dentro de sus funciones de recurrencia, los operadores diferenciales del método de diferencias finitas de cuarto orden programados en la *subrutina H*.
  
  Con ello se resuelven los valores de campo magnético H_2[A/m], campo eléctrico E_3 y -E_1 [V/m]. 
  El código realiza 250 snapshots pero sólo desplega 25 de ellos para observar el comportamiento de la onda. 
 
 * 8. Gráfica de Snapshots
 
 Se generan las gráficas del modelo geológico simulado, se dibuja la línea del contacto entre capas, las fronteras absorbentes y la línea receptora. De este modo se visualiza la propagación de la onda y los valores resultantes en una escala de colores.
 
 * 9. Gráfica del Radiorgrama
 
  Para generar el radiograma resultante de la línea de receptores se genera una matriz que guarda los valores de los dipolos receptores en cada iteración. Esta matriz se grafica en forma de radiorgrama.
  
# Referencias
  Carcione, J. (2015). Wave Fields in Real Media Wave Propagation in Anisotropic,
Anelastic, Porous and Electromagnetic Media. Trieste, Italy:
Elsevier.

  Carcione, J. (2006). Geophysical Software and Algorithms. A spectral numerical
method for electromagnetic diffusion. Geophysics, 71, 11-19.
doi: 10.1190/1.2159050
  
  
