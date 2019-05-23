# Analogía de Ondas Acústicas-Electromagnéticas
Este código 2D simula la propagación de ondas viscoelásticas SH y de ondas electromagnéticas del modo TM, a partir de la correspondencia entre sus valores de campo y propiedades del medio. 

El programa fue elaborado en "Matlab R2018b" y fue basado y modificado del código realizado por J. Carcione (2015,p.548) en Fortran 77.

Alejandra Alvarado Contreras. 23 de Mayo de 2019. alejandra.alalcon@hotmail.com.

# Ecuaciones utilizadas
El prototipo computacional utiliza una malla escalonada y métodos numéricos para resolver las siguientes ecuaciones de onda:

 * Ecuaciones SH
 
 Para un medio viscoelástico e isótropo, las ecuaciones de onda SH con la formulación esfuerzo-velocidad son:
 
    δ_t v_y  =  1/ρ  (δ_x σ_xy + δ_z σ_yz + f_y), 
    δ_t σ_yz =  μ (δ_z v_y  - 1/η σ_yz), 
    δ_t σ_xy =  μ (δ_x v_y  - 1/η σ_xy). 
     
 donde "v_y" es la velocidad en la dirección de "y", "ρ" es la densidad, "σ_xy" y "σ_yz" son los esfuerzos en las direcciones correspondientes, "f_y" es la fuerza volumétrica en la dirección de "y", "μ" es la constante de rigidez y "η" es la viscosidad.
     
 * Ecuaciones TM
 
 Para un medio con pérdidas de energía e isótropo, las ecuaciones del modo TM son:
     
    δ_t H_y    =  1/Ϗ  ( δ_x E_z + δ_z (-E_x) - M_y ),
    δ_t (-E_x) =  1/ε ( δ_z H_y - α ̂(-E_x) ),	
    δ_t E_z    =  1/ε ̂(δ_x H_y - α ̂E_z ).
    
 donde "H_y" es el campo magnético en la dirección de "y", "Ϗ" es la permeabilidad magnética,  "E_z" y "E_x" es el campo eléctrico en la dirección corrrespondiente, "M_y" es la fuente magnética en la dirección de "y", "ε" es la permitividad relativa y "α" es la conductividad.
 
 Con la correspondencia
  
    v_y  ↔  H_y
    σ_yz ↔  (-E_x)
    σ_xy ↔  E_z
    η    ↔  1/α
    μ    ↔  1/ε
    ρ    ↔  Ϗ
    f_y  ↔  -M_y

Para discretizar espacialmente las ecuaciones se utiliza el método de diferencias finitas y para la discretización en tiempo se programa el método de Runge Kutta, ambos de cuarto orden. 
 
 # Pre-requisitos
El código se realizó utilizando el software Matlab R2018b. Para obtener una prueba gratuita del software de 30 días ingrese a la URL "https://es.mathworks.com/campaigns/products/trials/targeted/dan.html"
 
# Ejecutando las pruebas
Una vez instalado Matlab, se abre el archivo del caso deseado: 

  "SH.m" El main para el ejemplo viscoelástico

  "TM.m" El main para el ejemplo electromagnético

Una vez seleccionado el ejemplo se abren las dos subrutinas que utiliza el programa:

  "H.m" Subrutina para la ejecución del Método de Diferencias Finitas de cuarto orden

  "wavelet.m" Subrutina para la ejecución de la Fuente impulso.

Ambos casos se basan en datos geofísicos reales mencionados en la tesis "Alvarado, 2019". 
Para ejecutar el caso deseado únicamente se presiona "Run" desde el main en "Matlab" y se arrojarán los snapshots y gráficas resultantes  para observar el comportamiento de onda.

# Ejemplo caso electromagnético
Como ejemplo se explica la estructura del caso electromagnético, el cual se basa en datos reales de una adquisicón GPR (Georadar de Penetración Terrestre) utilizando el modo biestático con polarización H o modo TM. El código simula una adquisición similar al arreglo "common shot gathered" con un sólo dipolo transmisor y varios dipolos receptores.
 
 * Parámetros de la malla y del tendido
 
Para simular la adquisición se decide la cantidad de dipolos receptores a lo largo del tendido que cubra el espaciamiento entre ellos. Después se decide el ancho de las fronteras absorbentes y se calcula la dimensión total de la malla numérica "nx(filas)", "nz(columnas)"
 
 * Parámetros de la fuente
 
Se decide la frecuencia que se utilizará para simular la fuente impulso y la discretización en tiempo "dt". En este ejemplo la frecuencia fue de 100 (MHz) y el "dt" de 6 (ns).
  
 * Parámetros del modelo geológico
 
El código funciona para un modelo geológico de una capa y un semi-espacio. Se decide cuál es el tipo de roca de cada capa y se ingresa el valor de las propiedades de permeabilidad magnética[H/m], permitividad dieléctrica [F/m] y un valor de Factor de Calidad obtenido con la siguiente fórmula:
  
    Q = 2π f[Hz] ε[H/m] / α[S/m].
    
Con ello dichos parámetros se calcula también el intervalo de nodos de la malla "Δx = Δz":
    
    v_min = √(1/(Ϗε)),
    λ_min = ( v_min/f_max ),
    Δx = Δz = λ_min/n
 
donde "v_min" es la velocidad mínima de las dos capas, la "f_max" son 100 MHz y "n" es el número de nodos necesarios para discretizar la longitud de onda mínima, que por similitud al ejemplo de J. Carcione (2015) se toman 9 nodos. 
 
 * Parámetros absorbentes
 
El código maneja fronteras absorbentes para evitar reflexiones de onda, para crearlas se utiliza un vector con factores de atenuación que disminuyen la amplitud de la onda conforme ésta hace contacto con los cuatro bordes de la malla.
  
 * Creación de los snapshots 
 
Se genera un ciclo con 250 muestras. En cada iteración se crean las fronteras absorbentes y se resuelve el método de Runge Kutta, el cual utiliza dentro de sus funciones de recurrencia, los operadores diferenciales del método de diferencias finitas programados en la *subrutina H* y la implementación de la fuente en la velocidad utilizando la *subrutina wavelet* 
    
Con ello, se resuelven los valores de campo magnético H_2[A/m] y campo eléctrico E_3 y -E_1 [V/m] durante el tiempo de muestreo total, que en este caso es de 150 (ns). 
   
 # Resultados, Gráficas de Snapshots 
 
El código realiza 250 snapshots pero sólo desplega 25 de ellos para observar el comportamiento de la onda. 
Para obsverar las gráficas correspondientes a cada componentes de campo, se deben selecciónar las opciones correspondiente.
Por ejemplo para observar los resultados del campo magnético se seleccionan las siguientes opciones y las del campo eléctrico se dejan comentadas:
 
         %Se crea un vector con el tamaño de la matriz correspondiente
           [A,B]=size(transpose(v2)); %campo magnético v2 <-> H2
           ....
         %Se grafican los valores de la matriz 
           imagesc(x,y,transpose(v2)); %campo magnético v2 <-> H2
         ...
         %Escalas 
          caxis([-5e-12 5e-12]) %campo magnético
          ...
         %Titulos 
         % title(['Campo Magnético H_2 en t=' num2str(a) '[ns]'],'Fontsize',19,'FontName','Arial', 'FontWeight',     'bold','FontAngle','italic','HorizontalAlignment','center')
           
Se imprimen las gráficas del modelo geológico simulado, junto con la línea del contacto entre capas, las fronteras absorbentes y la línea receptora. De este modo se visualiza la propagación de la onda y los valores resultantes de campo debido a la fuente en una escala de colores.
 
Por último, se genera un radargrama donde se observan los valores de campo en los dipolos receptores.
  
# Referencias
  Alvarado, A. (2019). Analogía entre la propagación de ondas acústicas y electromagnéticas y su aplicación en un prototipo computacional 2D. Tesis de licenciatura, UNAM, Ciudad Universitaria, Cd. Mx.
  
  Carcione, J. (2015). Wave Fields in Real Media Wave Propagation in Anisotropic,
Anelastic, Porous and Electromagnetic Media. Trieste, Italy: Elsevier.

  Carcione, J. (2006). Geophysical Software and Algorithms. A spectral numerical
method for electromagnetic diffusion. Geophysics, 71, 11-19.
doi: 10.1190/1.2159050
  
  
