                  %UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO. FACULTAD DE INGENIERÍA
                           %CÓDIGO BASADO Y MODIFICADO DE J.CARCIONE, 2015. 
               %ELABORADO POR ALVARADO CONTRERAS AEJANDRA PARA OBTENER EL TÍTULO DE
                                %INGENIERA GEOFÍSICA (Alvarado, 2019)
 
%Este código 2D simula la propagación de ondas viscoelásticas SH y de ondas electromagnéticas del modo TM, 
%a partir de la correspondencia entre sus valores de campo y propiedades del medio.

%-------------------------------------------------------------------------------------------------------
%CASO ELECTROMAGNÉTICO. MODELACIÓN DE LA PROPAGACIÓN DE ONDAS TM EN EL SUBSUELO DEL VALLE DE
%KONDRATOWA EN LAS MONTAÑAS TATRAS, POLONIA (Tomecka y col,(2017).
%-------------------------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------------------------
% Se limpia el workspace
clear all
% Se cierran todas las figuras
close all
%-------------------------------------------------------------------------------------------------------

%PARAMETROS DE LA MALLA Y DEL TENDIDO

rec = 10; %numero de dipolos receptores utilizados
tend = 308; %largo del tendido (nodos) %LT
srec = tend/(rec+1); %espaciamiento entre receptores (nodos)
nab = 50 ; %ancho de las fronteras absorbentes (nodos)

%Dimensiones de la malla rectangular
nx = tend + (2*nab); %numero de nodos en x [m] (filas) %LM
nz = tend + (2*nab); %numero de nodos en z [m] (columnas) %LM

%Ubicación de las lineas receptoras y del contacto en el modelo geológico
cont=180 ; %contacto entre capas (nodos)
lr=140;  %linea de receptores en el mallado (nodos)
rr=138; %ubicación de las letras de receptores (nodos)

%--------------------------------------------------------------------

%PARAMETROS DE LA FUENTE

freq=100000000; %Frecuencia central %100MHz
dt=0.6e-9; %muestreo en el tiempo(s),
nstep=250; %numero de muestras de tiempo %tiemp de muestreo 150 ns
nsp=nstep; % Cada nsp pasos se almacenara un snapshot

    %ubicacion de la fuente
     ix=204; %(nodos) R5 y R6
     %ix=78;  %R1 (nodos)
    % ix=330; %R10 %(nodos)
    iz=160; %(nodos)profundidad de la fuente

%--------------------------------------------------------------------

%PARÁMETROS DEL MODELO GEOLÓGICO

% Inicializacion de las propiedades
%PROPIEDADES VISCOELASTICAS y su correspondencia a PROPIEDADES ELECTROMAGNETICAS
mu = zeros(nx,nz); %mu(rigidez) <-------> 1/permitividad dielectrica 
rho = zeros(nx,nz); %rho(densidad) <--------> permeabilidad magnetica
eta = zeros(nx,nz); %eta(viscosidad) <---------> 1/conductividad

%Llenado de la malla con las propiedades correspondientes en las dos capas
for i=1:nx
    for j=1:nz

        %PRIMERA CAPA
         mu(i,j)=2.26E+10; %1/permitividad dieliectrica  (m/F)
         rho(i,j)=1.25664E-06; %permeabilidad magnetica mu (H/m)
         Q=50.55; %factor de calidad de la frecuencia central (Q = wpermitividad/conductividad)
         eta(i,j)=Q*mu(i,j)/(2*pi*freq); %1/conductividad
         v1=(mu(1,1)/rho(1,1))^(1/2);  
         
        %SEGUNDA CAPA
        if j>=cont
           mu(i,j)=1.25E+10; 
           rho(i,j)=1.25664E-06; 
           Q=25.02; 
           eta(i,j)=Q*mu(i,j)/(2*pi*freq);
           v3=(mu(cont,cont)/rho(cont,cont))^(1/2);
        end  
    end
end 

%Intervalo entre nodos en metros
vmin = min(v1,v3); %Se selecciona la velocidad mínima de las dos capas
lambdamin = vmin/freq; %Se calcula la longitud de onda mínima
dx = lambdamin/9; %intervalo entre nodos en x, (m)
dz = lambdamin/9; %intervalo entre nodos en z, (m)

%--------------------------------------------------------------------

%PARAMETROS ABSORBENTES

%Se crea el vector que tendra los factores absorbentes
r=0.99;
for i=1:nab
    ab(i)=r^i;
end 

%----------------------------------------------------------------------

%FUENTE

%Con la subrutina wavelet se crea la fuente de impulso
dt2=dt/2; %Se ingresa la mitad de dt para su posterior uso en el metodo de RG
[f,nw2]=wavelet(dt2,freq);
nw=nw2/2; %se divide a la mitad nw para su posterior uso en el metodo de RG

%-------------------------------------------------------------------------

%PARÁMETROS DE LOS MÉTODO NUMÉRICOS

%Pesos de cuarto orden referentes al metodo de diferencias finitas
x1 = 9/(8*dx);
x2 = -1/(24*dx);
z1 = 9/(8*dz);
z2 = -1/(24*dz);

% Se inicializan las variables de campo
%COMPONENTES DE CAMPO VISCOELASTICO <-----------> CAMPO ELECTROMAGNETICO
v2 = zeros(nx,nz); %velocidad <---------> Campo magnetico H2 
s12 = zeros(nx,nz); %Esfuerzo s12 <-------> Campo electrico -E1
s32 = zeros(nx,nz); %Esfuerzo s32 <------> Campo electrico E3
s = zeros(nstep,rec); %Matriz de radargramas 
 
%------------------------------------------------------------------------

%CICLO PARA CREAR LAS SNAPSHOTS

for n=1:250
    n1=(n*2)-1; % Numeros impares
    
    %Impresion del numero de snapshots realizados
    if (mod(n,10)==0) %Si el numero n es divisible entre 10, entonces imprime el numero n
        disp(n) %Se imprime en pantalla los numeros 10,20,30...250
    end
    
    %----------------------------------------------------------------------
    
    %CONDICIONES DE VARIABLE ABSORBENTES
    
    %Fronteras horizontales
    for j=1:nab
        j2=j+2;
        j3=nz-j-1;
        sab=ab(nab+1-j); %valor del coeficiente absorbente
        for i=3:nx-2
            v2(i,j2)=v2(i,j2)*sab; 
            v2(i,j3)=v2(i,j3)*sab;   
        end   
    end 

    %Fronteras verticales
    for i=1:nab
        i2=i+2;
        i3=nx-i-1;
        sab=ab(nab+1-i);
        for j=3:nz-2
            v2(i2,j)=v2(i2,j)*sab;
            v2(i3,j)=v2(i3,j)*sab;
        end 
    end
    
    %-------------------------------------------------------------------
    
    %METODO DE RUNGE-KUTTA
    
    %Este metodo resuelve las derivadas temporales de los valores de campo
    %magnético y eléctrico para ondas TM
    
    %vt^(n+1)= vt^n + dt/6 (D1 + 2D2 + 2D3 + D4)
    %D1= H(v^n)+ f^n
    %D2= H(v^n + dt/2 D1) + f^(n+1/2)
    %D3 = H(v^n + dt/2 D2) + f^(n+1/2)
    %D4 = H(v^n + dt D3) + f^(n+1)
    
    %Condiciones iniciales 
    for i=1:nx
        for j=1:nz
            %vt = (v2t,s12t,s32t) (Derivada del vector v)
            v2t(i,j)=v2(i,j);
            s12t(i,j)=s12(i,j);
            s32t(i,j)=s32(i,j);
            %v = (v2a, s12a, s32a) (Vector v) 
            v2a(i,j)=v2(i,j);
            s12a(i,j)=s12(i,j);
            s32a(i,j)=s32(i,j);
        end
    end
    
    %----------------------------------------------------------------------
    %D1 (Operador Delta1)
    %D1= H(v^n) 
 
    %Calculo de D1: H(v^n)= v =(v2a, s12a, s32a)
    [v2a,s12a,s32a] = H(v2a,s12a,s32a,mu,rho,eta,nx,nz,x1,x2,z1,z2);
    
    %Calculo parcial de vt: 
    %vt^(n+1)= vt^n + dt/6 (D1) 
    
    for i=1:nx
        for j=1:nz
            % Calculo parcial de vt 
            %vt^n = (v2t, s12t, s32t) (Vector vt inicial)
            %vt^(n+1)= vt^n + dt/6 D1 
            v2t(i,j)=v2t(i,j)+dt*v2a(i,j)/6; 
            s12t(i,j)=s12t(i,j)+dt*s12a(i,j)/6; 
            s32t(i,j)=s32t(i,j)+dt*s32a(i,j)/6; 
            
            %Se reescribe v para D2 
            %v^n= (v2,s12,s32) (Vector v de condiciones iniciales)
            %v = v^n + dt/2 D1 
            v2a(i,j)=v2(i,j)+0.5*dt*v2a(i,j);
            s12a(i,j)=s12(i,j)+0.5*dt*s12a(i,j);
            s32a(i,j)=s32(i,j)+0.5*dt*s32a(i,j);
        end
    end
    
    %D1= H(v^n) + f^n 
    if (n<=nw) 
           %Actualizacion de la velocidad mientras se inyecta la fuente 
           %Componentes impares de la fuente----> valores de dt 
         
           %Calculo parcial de vt
           %vt^n = vt^n + dt/6 [H(v^n)]
           %vt^(n+1)= vt^n + (dt/6) f^n
           v2t(ix,iz)=v2t(ix,iz)+dt*f(n1)/6; 
           
           %Se reescribe v para D2 (v^n + dt/2 D1)
           %v2a^n= v^n + dt/2 [H(v^n)]
           %v2a = v2a^n + (dt/2) f^n
           v2a(ix,iz)=v2a(ix,iz)+0.5*dt*f(n1);
    end
    
    %---------------------------------------------------------------------
    %D2 (Operador Delta2)
    %D2= H(v^n + dt/2 D1)  
   
    %Calculo de D2: H(v^n + dt/2*D1)= v =(v2a, s12a, s32a)
    [v2a,s12a,s32a] = H(v2a,s12a,s32a,mu,rho,eta,nx,nz,x1,x2,z1,z2);
   
    %Calculo parcial de vt: 
    %vt^(n+1)= vt^n + dt/6 (D1+ 2D2)
    
    for i=1:nx
        for j=1:nz     
            % Calculo parcial de vt 
            %vt^n= vt^n + (dt/6)D1 (Vector vt inicial)
            %vt^(n+1)= vt^n +(dt/3)D2 
            v2t(i,j)=v2t(i,j)+dt*v2a(i,j)/3;
            s12t(i,j)=s12t(i,j)+dt*s12a(i,j)/3;
            s32t(i,j)=s32t(i,j)+dt*s32a(i,j)/3;
            
            %Se reescribe v para D3 
            %v^n= (v2,s12,s32) (Vector v de condiciones iniciales)
            %v = v^n + dt/2 D2
            v2a(i,j)=v2(i,j)+0.5*dt*v2a(i,j);
            s12a(i,j)=s12(i,j)+0.5*dt*s12a(i,j);
            s32a(i,j)=s32(i,j)+0.5*dt*s32a(i,j);
        end
    end
    
    %D2= H(v^n + dt/2 D1) + f^(1/2) 
    if(n<=nw)
           %Actualizacion de la velocidad mientras se inyecta la fuente 
           %Componentes pares de la fuente----> valores de dt/2 
           
           %Calculo parcial de vt
           %vt^n= vt^n +(dt/3)H(v^n + dt/2*D1)
           %vt^(n+1)= vt^n +(dt/3)f^[n+(1/2)]
           v2t(ix,iz)=v2t(ix,iz)+dt*f(n1+1)/3;
           
           %Se reescribe v para D3 (v^n + dt/2 D2)
           %v2a^n= v^n + (dt/2) [H(v^n + dt/2*D1)]
           %v2a = v2a^n + (dt/2) f^[n+(1/2)]
           v2a(ix,iz)=v2a(ix,iz)+0.5*dt*f(n1+1);
    end
   
    %------------------------------------------------------------
    %D3 (Operador Delta3)
    %D3= H(v^n + dt/2 D2)  
   
    %Calculo de D3: H(v^n + dt/2*D2)= v =(v2a, s12a, s32a)
    [v2a,s12a,s32a] = H(v2a,s12a,s32a,mu,rho,eta,nx,nz,x1,x2,z1,z2);
   
    %Calculo parcial de vt: 
    %vt^(n+1)= vt^n + dt/6 (D1 + 2D2 + 2D3)
    
    for i=1:nx
        for j=1:nz
            %Calculo parcial de vt 
            %vt^n= vt^n + (dt/6) (D1+2D2) (Vector vt inicial)
            %vt^(n+1)= vt^n +(dt/3)D3 
            v2t(i,j)=v2t(i,j)+dt*v2a(i,j)/3;
            s12t(i,j)=s12t(i,j)+dt*s12a(i,j)/3;
            s32t(i,j)=s32t(i,j)+dt*s32a(i,j)/3;
            
            %Se reescribe v para D4 
            %v^n= (v2,s12,s32) (Vector v de condiciones iniciales)
            %v = v^n + dt D3
            v2a(i,j)=v2(i,j)+dt*v2a(i,j);
            s12a(i,j)=s12(i,j)+dt*s12a(i,j);
            s32a(i,j)=s32(i,j)+dt*s32a(i,j);
        end
    end
    
    %D3= H(v^n + dt/2 D2) + f^(1/2) 
    if(n<=nw)  
           %Actualizacion de la velocidad mientras se inyecta la fuente 
           %Componentes pares de la fuente----> valores de dt/2 
           
           %Calculo parcial de vt
           %vt^n= vt^n +(dt/3)H(v^n + dt/2*D2)
           %vt^(n+1)= vt^n +(dt/3) f^[n+(1/2)]
           v2t(ix,iz)=v2t(ix,iz)+dt*f(n1+1)/3;
           
           %Se reescribe v para D4 (v^n + dt D3)
           %v2a^n= v^n + dt [H(v^n + dt/2 *D2)]
           %v2a = v2a^n + dt f^[n+(1/2)]
           v2a(ix,iz)=v2a(ix,iz)+dt*f(n1+1);
    end

    %--------------------------------------------------------
    %D4 (Operador Delta4)
    %D4= H(v^n + dt*D3)  
   
    %Calculo de D4: H(v^n + dt*D3)= v =(v2a, s12a, s32a)
    [v2a,s12a,s32a] = H(v2a,s12a,s32a,mu,rho,eta,nx,nz,x1,x2,z1,z2);
   
    %Calculo de vt: 
    %vt^(n+1)= vt^n + dt/6 (D1 + 2D2 + 2D3 + D4)
   
    for i=1:nx
        for j=1:nz
            %Calculo de vt 
            %vt^n= vt^n + (dt/6) (D1+2D2+2D3) (Vector vt inicial)
            %vt^(n+1)= vt^n +(dt/6)D4
            v2t(i,j)=v2t(i,j)+dt*v2a(i,j)/6;
            s12t(i,j)=s12t(i,j)+dt*s12a(i,j)/6;
            s32t(i,j)=s32t(i,j)+dt*s32a(i,j)/6;
        end 
    end 
    
    %D4= H(v^n + dt*D3) + f^(1/2) 
    if(n<=nw)
           %Actualizacion de la velocidad mientras se inyecta la fuente 
           %Componentes pares de la fuente----> valores de dt/2 
           
           %Calculo de vt
           %vt^n= vt^n +(dt/6) H(v^n + dt*D3)
           %vt^(n+1)= vt^n +(dt/6) f^[n+(1/2)]
           v2t(ix,iz)=v2t(ix,iz)+dt*f(n1+1)/6;
    end
    
    %-----------------------------------------------------------------
    
    %Se renombran los vectores resultantes 
    for i=1:nx
        for j=1:nz
            v2(i,j)=v2t(i,j); %v2 <-> H2
            s12(i,j)=s12t(i,j); %s12 <-> E3
            s32(i,j)=s32t(i,j); %s32 <-> -E1
        end
    end
       
    %-----------------------------------------------------------------
    
    %MATRIZ DE RADARGRAMAS
               for j=1:rec
                s(n,j)=v2(nab+(j*srec),lr); %v2 <-> H2
               end 
                
    %-----------------------------------------------------------------
    
    %GRAFICAS DE SNAPSHOTS
    
    if (mod(n,10)==0)
        %valor de tiempo de muestreo
        a=(dt*n)*(10^9);
        %Escribir el numero de snapshot en curso
        disp('Snapshot'),n;
        
        %Se crea una ventana de figura
        n=figure; 
        %Se crea un vector con el tamaño de la matriz correspondiente
            [A,B]=size(transpose(v2));  %campo magnético v2 <-> H2
       %    [A,B]=size(transpose(s12));  %s12 <-> E3
       %    [A,B]=size(transpose(s32));   %s32 <-> (-E1)
         %Asignamos el rango de los ejes de la grafica en metros.
         x=(1:1:A)*dx;
         z=(1:1:B)*dz;
         %Se crea una cuadricula 2D basada en 'x' y 'z' 
         [X,Z]=meshgrid(x,z);
         
         %Se grafican los valores de la matriz 
          imagesc(x,z,transpose(v2)) %campo magnético v2 <-> H2
      %   imagesc(x,z,transpose(s12)) %s12 <-> E3
       %  imagesc(x,z,transpose(s32))  %s32 <-> (-E1)
         hold on;
         
         %Se grafica el contacto (linea horizontal)
         p0= plot([1*dx B*dx],[cont*dx cont*dx],'g','LineWidth',2);
         %Se dibujan las fronteras absorbentes
            %lineas horizontales
             p3= plot([1*dx B*dx], [nab*dx nab*dx], 'b','LineWidth',2);
             plot([1*dx B*dx], [(nx-nab)*dx (nx-nab)*dx], 'b','LineWidth',2);
            %lineas verticales
             plot([nab*dx nab*dx], [1*dx A*dx], 'b','LineWidth',2);
             plot([(nz-nab)*dx (nz-nab)*dx], [1*dx A*dx], 'b','LineWidth',2);
         %Se grafican y nombran los receptores
         for i=1:rec
            p4= plot([((nab+srec)*dx) (nab+(i*srec))*dx],[lr*dx lr*dx],'rv','MarkerFaceColor','r','MarkerSize',8);
            text((nab+(i*srec))*dx,rr*dx,['R' num2str(i)],'HorizontalAlignment','center', 'VerticalAlignment','bottom','fontsize',8,'FontWeight', 'bold')
         end
         
         %Formato de imagen
            axis('ij')  %Direccionn inversa. El eje 'z' es vertical y los valores aumentan de arriba a abajo.  
            %Barra de colores
            colorbar %mostrar barra de colores
           %Poner título a la barra de colores
            hcb=colorbar;
            title(hcb,'A / m','Fontsize',16,'FontName','Arial', 'FontWeight', 'bold'); %campo magnético
         %  title(hcb,'V / m','Fontsize',16,'FontName','Arial', 'FontWeight', 'bold'); %campo electrico
           
           %nombres y tamaños de los ejes
            xlabel('Distancia [m]','Fontsize',15,'FontWeight','bold' ) %nombre y tamaño de ejes
            ylabel('Profundidad [m]','Fontsize',15,'FontWeight', 'bold')
             %colocar legendas en la grafica
            legend([p0 p3 p4],'Contacto entre capas 19.9[m]','Fronteras Absorbentes','Receptores a 15.5[m]','Location','southeast','fontsize',10)
     
           %Escalas 
           caxis([-5e-12 5e-12]) %campo magnético v2 <-> H2
        %  caxis([-5e-10 5e-10]) %s12 <-> E3
         %  caxis([-5e-10 5e-10]) %s32 <-> (-E1)
          %Titulos 
          title(['Campo Magnético H_2 en t = ' num2str(a) ' [ns]' ],'Fontsize',17,'FontName','Arial', 'FontWeight', 'bold', 'HorizontalAlignment', 'Center')
        %  title(['Campo Eléctrico E_3 en t = ' num2str(a) ' [ns]'],'Fontsize',17,'FontName','Arial', 'FontWeight', 'bold','HorizontalAlignment','center')
         %  title(['Campo Eléctrico -E_1 en t=' num2str(a) '[ns]'],'Fontsize',17,'FontName','Arial', 'FontWeight', 'bold','HorizontalAlignment','center')
             
         %Imprimir y guardar imagenes
       print(n,'-dpng');
    end   
end

%Se grafica la matriz de radiorgramas
h=zeros(1,rec); 
hold off

for j=1:rec
%Posicionamos las figuras de cada traza 
    h(j)= subplot('position', [0.1, ((1/(rec+2)))*j, 0.7, 1/(rec+2)]); %[left bottom width high]  
    plot(h(j),s(:,j));
    %Titulo general de la figura 
    sgtitle('Radargrama','Fontsize',20,'FontName','Arial', 'FontWeight', 'bold','FontAngle','italic');
    %Nombramos a los ejes 'y' de cada figura
    ylabel(['R' num2str(j)],'Fontsize',20,'FontName','Arial','FontWeight', 'bold')
    %Escalar el eje y
    set(h(j),'YTick',-4e-12:4e-12:4e-12)
    %Quitamos los ticks y labels del eje 'x' de las figuras
    set(gca,'xlabel',[],'XTick',[],'Fontsize',10,'FontName','Arial')
    box off
    %Retomamos el tick y label de la primera figura
    set(h(1),'XTick',(25:25:250));
    xlabel(h(1),'Muestras','Fontsize',15,'FontName','Arial','FontWeight', 'bold')
end

%------------------------------------------------------------------------------------------------------------------- -         

%----------------------------------------------------------------------------------------
%REFERENCIAS

%Alvarado, A. (2019). Analogía entre la propagación de ondas viscoelásticas y electromagnéticas: 
%Desarrollo de un prototipo computacional 2D. Tesis de licenciatura, UNAM, Ciudad Universitaria, Cd. Mx.

%Carcione, J. (2015). Wave Fields in Real Media Wave Propagation in Anisotropic,
%"Anelastic, Porous and Electromagnetic Media. Trieste, Italy: Elsevier.

%Tomecka, S., Bogdan Zoga la, T. G., Grazyna Dzik, T. D. & Jochymczyk, K.
%(2017). "Application of electrical and electromagnetic methods to study
%sedimentary covers in high mountain areas." Acta Geophys, 65, 743-755.
%doi:10.1007/s11600-017-0068-z

%------------------------------------------------------------------------------------------------------------------
           
   




