# Ondas Acústicas-Electromagnéticas
El código 2D simula la propagación de ondas viscoelásticas SH y de ondas electromagnéticas del modo TM. 
Para discretizar espacialmente las ecuaciones se utiliza el método de diferencias finitas y para la discretización en tiempo se programa el método de Runge Kutta, ambos de cuarto orden. 
# Pre-requisitos
El código se realizó utilizando el software Matlab R2018b pero puede ser utilizado con cualquier otra versión de Matlab
Para obtener una prueba gratuita del software de 30 días ingrese a la URL "https://es.mathworks.com/campaigns/products/trials/targeted/dan.html"
# Ejecutando las pruebas
Una vez instalado Matlab, en la carpeta Ondas Acústicas/Electromagnéticas vienen 4 archivos. 

"SH.m" El main para el ejemplo viscoelástico

"TM.m" El main para el ejemplo electromagnético

"H.m" Subrutina para la ejecución del método de diferencias finitas

"Wavelet.m" Subrutina para la ejecución del impulso de onda

Como ejemplo, se explican los datos que se ingresan para ejecutar el ejemplo electromagnético "SH.m" se siguen los siguientes pasos:
  *1. Parámetros de la malla y del tendido
  
  El código se basa en datos reales de una adquisicón GPR (Georadar de        Penetración Terrestre) utilizando el modo biestático con polarización H o modo TM.
  Para simular el modelo geológico, primero se decide el número de dipolos receptores a utilizar y el largo del tendido que cubra el espaciamiento entre ellos. Después se decide el ancho de las fronteras absorbentes y se calcula la dimensión total de la malla numérica "nx(filas)", "nz(columnas)", "dx" y "dz" (intervalos entre nodos).
  
  *2. Parámetros de la fuente
  
  Se decide la frecuencia que se utilizará para simular la fuente impulso, que en este ejemplo fue de 100 MHz.
  
  *3. Modelo Heterogéneo
  
  El código funciona para un modelo geológico de una capa y un semi-espacio.
  Se decide cuál es el tipo de roca de cada capa y se obtienen las propiedades de conductividad[s/m], permeabilidad magnética[H/m] y permitividad dieléctrica [F/m] y se obtiene su factor de calidad (Q) 
  
  *4.Parámetros absorbentes
  
  El código maneja fronteras absorbentes para evitar reflexiones de onda, para crearlas se utiliza un vector con factores de atenuación que disminuyen la amplitud de la onda conforme ésta hace contacto con los cuatro bordes de la malla.
  
  *5. Subrutina Wavelet
  
  Con esta subrutina se crea la fuente de impulso la cual es una función que genera una ondícula gaussiana, la cual utiliza los datos de frecuencia ingresados y los de muestreo en el tiempo "dt"
  
  *6. Creación de los snapshots
  
  El código realiza 250 snapshots pero sólo desplega 25 de ellos para observar el comportamiento de la onda. 
  Se crea un ciclo para generar cada uno de los snapshots, dentro del cual se crean las fronteras absorbentes, y se resuelven los valores de campo magnético H_2[A/m], campo eléctrico E_3 y -E_1  [V/m] a través del método de Runge Kutta de cuarto orden, el cual utiliza dentro de sus funciones de recurrencia, los operadores diferenciales del método de diferencias finitas de cuarto orden programado con la *subrutina H*.
  
 *7. Matriz de Radiograma
 
 Para generar el radiograma resultante de la línea de receptores se crea una matriz, la cual va guardando los valores de los dipolos receptores en cada iteración.
 
 *8. Gráfica de Snapshots
 
 Se generan las gráficas del modelo geológico simulado, junto con la línea del contacto entre capas, de las fronteras absorbentes, de la línea receptora y se visualiza a través de ellas la propagación de la onda y los valores resultantes en una escala de colores.
 
 *9. Gráfica del radiorgrama
 
 Por último se genera la gráfica del radiorgrama, la cual grafica las trazas de cada uno de los dipolos receptores.
  
  
  
