# Analogía de Ondas Viscoelásticas-Electromagnéticas
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
    
donde "v_y" es la velocidad en la dirección de "y", "ρ" es la densidad, "σ_xy" y "σ_yz" son los esfuerzos en las direcciones correspondientes, "f_y" es la fuerza volumétrica en la dirección de "y",  "μ" es la constante de rigidez y "η" es la viscosidad. 
     
 * Ecuaciones TM
 
 Para un medio con pérdidas de energía e isótropo, las ecuaciones del modo TM son:
     
    δ_t H_y    =  1/Ϗ  ( δ_x E_z + δ_z (-E_x) - M_y ),
    δ_t (-E_x) =  1/ε ( δ_z H_y - α ̂(-E_x) ),	
    δ_t E_z    =  1/ε ̂(δ_x H_y - α ̂E_z ),
    
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

Se debe tomar en cuenta que las siguientes subrutinas deben estar en la misma carpeta que el "main" para que se ejecute correctamente el programa:

  "H.m" Subrutina para la ejecución del Método de Diferencias Finitas de cuarto orden

  "wavelet.m" Subrutina para la ejecución de la Fuente impulso.

Ambos casos se basan en datos geofísicos reales mencionados en la tesis "Alvarado, 2019". 

# Ejemplo caso electromagnético
Como ejemplo se explica la estructura del caso electromagnético, el cual se basa en datos reales de una adquisicón GPR (Georadar de Penetración Terrestre) utilizando el modo biestático con polarización H o modo TM. El código simula una adquisición similar al arreglo "common shot gathered" con un sólo dipolo transmisor (fuente) y varios dipolos receptores.
 
 * Parámetros de la malla y del tendido
 
Para simular la adquisición se decide la cantidad de dipolos receptores y el largo del tendido que cubra el espaciamiento entre ellos. Después se decide el ancho de las fronteras absorbentes y se calcula la dimensión total de la malla rectangular numérica "nx(filas)", "nz(columnas)"

    LM = (LT + 2 * nab),
    
donde "LM" es la longitud de la malla, "LT" es la longitud del tendido y "nab" es el tamaño de las fronteras absorbentes.

Así mismo se debe posicionar el contacto entre capas y la línea de receptores donde se desee en el mallado.
 
 * Parámetros de la fuente
 
Se decide la frecuencia que se utilizará para simular la fuente impulso, la discretización en tiempo "dt" y las muestras en el tiempo. En este ejemplo la frecuencia fue de 100 (MHz), el "dt" de 6 (ns) y se realizan 250 muestras. 

También se posiciona la fuente en el mallado.
  
 * Parámetros del modelo geológico
 
El código funciona para un modelo geológico de una capa y un semi-espacio. Se decide cuál es el tipo de roca de cada capa y se ingresa el valor de las propiedades de permeabilidad magnética[H/m], permitividad dieléctrica [F/m] y un valor de Factor de Calidad obtenido con la siguiente fórmula:
  
    Q = 2π f[Hz] ε[F/m] / α[S/m].
    
A partir de estos parámetros, también se calcula el intervalo de nodos de la malla "Δx = Δz" de la siguiente manera:
    
    v_min = √(1/(Ϗε)),
    λ_min = v_min/f_max,
    Δx = Δz = λ_min/n,
    
donde "v_min" es la velocidad mínima de las dos capas, la "f_max" son 100 MHz y "n" es el número de nodos necesarios para discretizar la longitud de onda mínima, que por similitud al ejemplo de J. Carcione (2015) se toman 9 nodos. 
 
 * Parámetros absorbentes
 
El código maneja fronteras absorbentes para evitar reflexiones de onda, para crearlas se utiliza un vector con factores de atenuación que disminuye la amplitud de la onda conforme ésta hace contacto con los cuatro bordes de la malla.
  
    r=0.99;
    for i=1:nab
        ab(i)=r^i;
    end
Donde "r" es el factor de atenuación inicial y "ab" es el vector de atenuación.

 * Subrutina Wavelet 

Esta función crea una fuente impulso que genera una ondícula gaussiana establecida por J. Carcione, (2006):

    f(t) = exp⁡〖-(Δw^2 (t-t_0)^2))/4〗*cos⁡(ϖ(t-t_0)),
    donde
    t_0 = 6/5 F_s  es el tiempo de retraso con F_s como frecuencia máxima.
    n   = t_0/dt es el número de muestra del pulso,
    ϖ   = 2πF_s  es la frecuencia angular central,
    Δw  = 0.5ϖ es el ancho del pulso.

Para graficar la ondícula se debe seleccionar el caso en cuestión, en este ejemplo se activa el caso electromagnético y se deja comentado (%) el caso viscoelástico.

 * Subrutina H

Esta función resuelve las derivadas espaciales del sistema de ecuaciones TM a través del método de diferencias finitas de cuarto orden con la siguiente formulación:

    f'^(x_i) = 〖1/24 f(x_(i-3/2)) - 9/8 f(x_(i-1/2)) + 9/8 f(x_(i+1/2)) - 1/24 f(x_(i+3/2))〗/ h

donde "h" representa el intervalo "Δx o  Δz", f'^(x_i)" es la derivada espacial del campo magnético o eléctrico, y las constantes "∓1/24" y "∓9/8" representan los coeficientes de cuarto orden del método.

 * Creación de los snapshots 
 
Se genera un ciclo con 250 muestras y se decide mostrar sólo 25 de ellas para observar el comportamiento de la onda. 

En cada iteración se crean las fronteras absorbentes (2 horizontales y 2 verticales).

Así mismo se resuelve el método de Runge Kutta de cuarto orden para la discretización en el tiempo con la siguiente formulación:

    Vt^(n+1) = Vt^n + dt/6 (Δ_1 + 2Δ_2 + 2Δ_3 + Δ_4),
    donde
    Δ_1 = H (V^n) + f^n,
    Δ_2 = H (V^n + dt/2 Δ_1) + f^(n+1/2),
    Δ_3 = H (V^n + dt/2 Δ_2) + f^(n+1/2),
    Δ_4 = H (V^n + dt Δ_3) + f^(n+1) .
    
donde "Vt=(H_2t, E_3t, -E_1t)" es la dervidada temporal del vector "V=(H_2, E_3, -E_1), "Δ's" son las funciones de recurrencia, "H" representa el sistema de ecuaciones TM discretizadas espacialmente con el método de diferencias finitas en la *subrutina H* y "f" es la fuente implementada a través de la *subrutina wavelet*.
   
Con ello, se resuelven los valores de campo magnético H_2[A/m] y campo eléctrico E_3 y -E_1 [V/m] durante el tiempo de muestreo total, que en este caso es de 150 (ns).

 # Resultados, Gráficas de Snapshots y Radargrama

Para obtener los resultados del valor de campo correspondiente, se deben seleccionar distintas opciones al momento de graficar las snapshots.
Por ejemplo, para obtener los resultados del campo magnético se seleccionan las opciones de dicho campo y se dejan comentados (%) los de campo eléctrico:
 
         %Se crea un vector con el tamaño de la matriz correspondiente
           [A,B]=size(transpose(v2));  %campo magnético v2 <-> H2
            %  [A,B]=size(transpose(s12));  %s12 <-> E3
            %   [A,B]=size(transpose(s32));   %s32 <-> (-E2)
         %Se grafican los valores de la matriz 
           imagesc(x,y,transpose(v2)) %campo magnético v2 <-> H2
           ...
          %Escalas 
           caxis([-5e-12 5e-12]) %campo magnético v2 <-> H2
           ... 
         %Titulos 
         % title(['Campo Magnético H_2 en t=' num2str(a) '[ns]'],'Fontsize',19,'FontName','Arial', 'FontWeight','bold','FontAngle','italic','HorizontalAlignment','center')
           ...

 * Ejecución del programa 
 
Una vez ingresados los datos deseados, para ejecutar el programa se presiona "Run" desde el main en "Matlab" y se obtendrán los snapshots y gráficas resultantes para observar el comportamiento de onda.

Se imprimen las gráficas del modelo geológico simulado junto con la línea del contacto entre capas, las fronteras absorbentes y la línea receptora. De este modo se visualiza la propagación de la onda y los valores resultantes de campo debido a la fuente en una escala de colores.
 
Por último, se genera un radargrama donde se observan los valores de campo magnético en los dipolos receptores.
  
# Referencias
  Alvarado, A. (2019). Analogía entre la propagación de ondas viscoelásticas y electromagnéticas: Desarrollo de un prototipo computacional 2D. Tesis de licenciatura, UNAM, Ciudad Universitaria, Cd. Mx.
  
  Carcione, J. (2015). Wave Fields in Real Media Wave Propagation in Anisotropic,
Anelastic, Porous and Electromagnetic Media. Trieste, Italy: Elsevier.

  Carcione, J. (2006). Geophysical Software and Algorithms. A spectral numerical
method for electromagnetic diffusion. Geophysics, 71, 11-19.
doi: 10.1190/1.2159050
  
  
