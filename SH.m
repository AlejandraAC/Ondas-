%Este programa modela ondas SH
%-------------------------------------------------------------------

% Limpiar el workspace
clear all
% Cerrar todas las figuras
close all

%-------------------------------------------------------------------
%PARAMETROS DE LA MALLA
nx=121; %numero de nodos en x 
nz=121; %numero de nodos en z
dx = 1.38; %intervalo entre nodos en x, (m)
dz = 1.38; %intervalo entre nodos en z, (m)
dt=0.001; %muestreo en el tiempo(s), 1ms
nstep=250; %numero de muestras de tiempo (ms)
nsp=nstep; % Cada nsp pasos se almacenara un snapshot
cont=35; %contacto entre capas 
nab=15  ; %tamaño del vector absorbente
%------------------------------------------------------------
%Adquisición del tendido 
rec=12; %numero de receptores utilizados
tend=nx-(2*nab); %largo del tendido
srec=tend/(rec+1); %espaciamiento entre geófonos 
lr=25; %linea de receptores en el mallado 
%--------------------------------------------------------------------
%--------------------------------------------------------------------
%PARAMETROS DE LA FUENTE
ix=61; %(nodos)
iz=30; %(nodos)profundidad de la fuente
freq=25; %Frecuencia central %Hz
%--------------------------------------------------------------------
%--------------------------------------------------------------------
%MODELO HETEROGENEO
% InicializaciOn de las propiedades
%PROPIEDADES ACUSTICAS y su correspondencia a PROPIEDADES ELECTROMAGNETICAS
mu = zeros(nx,nz); %mu(rigidez) <-------> 1/permitividad dielectrica 
rho = zeros(nx,nz); %rho(densidad) <--------> permeabilidad magnetica
eta = zeros(nx,nz); %eta(viscosidad) <---------> 1/conductividad

% Se llena la malla de las propiedades correspondientes en las dos capas
for i=1:nx
    for j=1:nz

        %PRIMERA CAPA:
        mu(i,j)=273626056.8; %Pa
        rho(i,j)=2812; %kg/m3
        Q=93.582; %factor de calidad de la frecuencia central (Qs = 0.3*v_s = 0.3*311.94 = 93.582)
        eta(i,j)=Q*mu(i,j)/(2*pi*freq); 
        
        %SEGUNDA CAPA
         if j>=cont
           rho(i,j)=3017; %kg/m3
           mu(i,j)=1929046100; %Pa
            Q=159.924; % (Qs = 0.2*v_s = 0.2*799.62 = 159.924)
            eta(i,j)=Q*mu(i,j)/(2*pi*freq);
         end
    end
end 

%--------------------------------------------------------------------
%--------------------------------------------------------------------
%PARAMETROS ABSORBENTES
%Se crea el vector que tendra los factores absorbentes
%nab = 12; %tamaño del vector absorbente
r=0.99;
for i=1:nab
    ab(i)=r^i;
end 

%----------------------------------------------------------------------
%-----------------------------------------------------------------------
%FUENTE
%Con la subrutina wavelet se crea la fuente de impulso
dt2=dt/2; %Se ingresa la mitad de dt para su posterior uso en el metodo de RG
[f,nw2]=wavelet(dt2,freq);
nw=nw2/2; %se divide a la mitad nw para su posterior uso en el metodo de RG

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%Pesos de cuarto orden referentes al metodo de diferencias finitas
x1 = 9/(8*dx);
x2 = -1/(24*dx);
z1 = 9/(8*dz);
z2 = -1/(24*dz);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Se inicializan las variables 
%COMPONENTES DE CAMPO ACUSTICO <-----------> CAMPO ELECTROMAGNETICO
v2 = zeros(nx,nz); %velocidad <---------> Campo magnetico H2 
s12 = zeros(nx,nz); %Esfuerzo s12 <-------> Campo electrico -E1
s32 = zeros(nx,nz); %Esfuerzo s32 <------> Campo electrico E3
 s = zeros(nstep,rec); %Matriz de sismogramas 
 L = zeros(3,rec);
 

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%Creacion de los 250 snapshots

for n=1:nstep
    n1=(n*2)-1; % Numeros impares
    
    %Impresion del numero de snapshots realizados
    if (mod(n,10)==0) %Si el numero n es divisible entre 10, entonces imprime el numero n
        disp(n) %Se imprime en pantalla los numeros 10,20,30...250
    end
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %CONDICIONES DE VARIABLE ABSORBENTES:
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
    %--------------------------------------------------------------------
    %METODO DE RUNGE-KUTTA
    %Este metodo resuelve las derivadas temporales de la velocidad y los
    %esfuerzos para una onda SH
    %vt^(n+1)= vt^n + dt/6 (D1 + 2D2 + 2D3 + D4)
    %D1= H(v)+ f^n
    %D2= H(v^n + dt/2 D1) + f^(n+1/2)
    %D3 = H(v^n + dt/2 D2) + f^(n+1/2)
    %D4 = H(v^n + dt D3) + f^(n+1/2)
    
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
    %D1= H(v) 
 
    %Cï¿½lculo de D1: H(v)= v =(v2a, s12a, s32a)
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
    
    %D1= H(v) + f^n 
    if (n<=nw) 
           %Actualizacion de la velocidad mientras se inyecta la fuente 
           %Componentes impares de la fuente----> valores de dt 
         
           %Calculo parcial de vt
           %vt^n = vt^n + dt/6 [H(v)]
           %vt^(n+1)= vt^n + (dt/6) f^n
           v2t(ix,iz)=v2t(ix,iz)+dt*f(n1)/6; 
           
           %Se reescribe v para D2 (v^n + dt/2 D1)
           %v2a^n= v^n + dt/2 [H(v)]
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
            v2(i,j)=v2t(i,j); 
            s12(i,j)=s12t(i,j);
            s32(i,j)=s32t(i,j);
           % v2tr(i,j)=transpose(v2t(i,j)); 
        end
    end
    
    
    %-----------------------------------------------------------------
    %----------------------------------------------------------------- 
    %MATRIZ DE SISMOGRAMAS
                for j=1:rec
                s(n,j)=v2(nab+(j*srec),lr);
                end
    %-----------------------------------------------------------------
    %------------------------------------------------------------------
    %GRAFICAS DE SNAPSHOTS
    if (mod(n,10)==0)
        %valor de tiempo de muestreo
        a=dt*n;
        %Escribir el nï¿½mero de snapshot en curso
        disp('Snapshot'),n;
        %Se crea una ventana de figura
        n=figure; 
        %Se crea un vector con el tamaño de la matriz correspondiente
      %     [A,B]=size(transpose(v2)); 
           [A,B]=size(transpose(s12)); 
          % [A,B]=size(transpose(s32)); 
         %Asignamos el rango de los ejes de la grafica en metros.
         x=(1:1:A)*dx;
         y=(1:1:B)*dx;
         %Se crea una cuadricula 2D basada en 'x' y 'y' 
         [X,Y]=meshgrid(x,y);
         %Se grafican los valores de la matriz 
       %     imagesc(x,y,transpose(v2))
          imagesc(x,y,transpose(s12))
        %   imagesc(x,y,transpose(s32))
         hold on;
         %Se grafica el contacto (linea horizontal)
         p1= plot([1 B*dx],[cont*dx cont*dx],'g','LineWidth',2);
         %Se dibujan las fronteras absorbentes
            %lineas horizontales
             p2= plot([1 B*dx], [nab*dx nab*dx], 'b','LineWidth',2);
             plot([1 B*dx], [(nx-nab)*dx (nx-nab)*dx], 'b','LineWidth',2);
            %lineas verticales
             plot([nab*dx nab*dx], [1 A*dx], 'b','LineWidth',2);
             plot([(nx-nab)*dx (nx-nab)*dx], [1 A*dx], 'b','LineWidth',2);
         %Se grafican y nombran los geofonos
         for i=1:rec
            p3= plot([((nab+srec)*dx) (nab+(i*srec))*dx],[lr*dx lr*dx],'r^','MarkerFaceColor','r','MarkerSize',8);
            text((nab+(i*srec))*dx,31,['G' num2str(i)],'HorizontalAlignment','center', 'VerticalAlignment','bottom','fontsize',8,'FontWeight', 'bold')
         end
         %graficas de las lineas receptoras
       % for i=1:rec
         %  plot([((nab+(i*srec))*dx) ((nab+(i*srec)))*dx], [0 A*dx], 'r');
         %end  
         %Formato de imagen
            axis('ij')  %Direccionn inversa. El eje 'y' es vertical y los valores aumentan de arriba a abajo.  
            colorbar %mostrar barra de colores
            xlabel('Distancia [m]','Fontsize',15,'FontWeight','bold' ) %nombre y tamaño de ejes
            ylabel('Profundidad [m]','Fontsize',15,'FontWeight', 'bold')
           %colocar legendas en la grafica
            legend([p1 p2 p3],'Contacto entre capas a 48.3 [m]','Fronteras Absorbentes','Receptores a 34.5[m]','Location','southeast','fontsize',10)
           %Escalas 
          %  caxis([-20e-5 9e-5]) %velocidad
           caxis([-100 100]) %s12??
         %  caxis([-100 100]) %s32
          %Titulos 
         %   title(['Snapshot de Velocidad v_2 en t=' num2str(a) '[s]'],'Fontsize',19,'FontName','Arial', 'FontWeight', 'bold','FontAngle','italic','HorizontalAlignment','center')
           title(['Snapshot de Esfuerzos \sigma_{12} en t=' num2str(a) 's'],'Fontsize',19,'FontName','Arial', 'FontWeight', 'bold','FontAngle','italic','HorizontalAlignment','center')
         %  title(['Snapshot de Esfuerzos \sigma_{yz} en t=' num2str(a) 's'],'Fontsize',19,'FontName','Arial', 'FontWeight', 'bold','FontAngle','italic','HorizontalAlignment','center')
             
         %Imprimir y guardar imagenes
       %  print(n,'-dpng');
    end   
end

%Se grafica la matriz de sismogramas
h=zeros(1,rec); 
hold off
for j=1:rec
%Posicionamos las figuras de cada traza 
    h(j)= subplot('position', [0.1, ((1/(rec+2)))*j, 0.8, 1/(rec+2)]); %[left bottom width high]
    plot(h(j),s(:,j));
    %Titulo general de la figura 
    sgtitle('Sismograma','Fontsize',20,'FontName','Arial', 'FontWeight', 'bold','FontAngle','italic');
    %Nombramos a los ejes 'y' de cada figura
    ylabel(['G' num2str(j)],'Fontsize',20,'FontName','Arial','FontWeight', 'bold')
    %Escalar el eje y
    set(h(j),'YTick',-2e-4:1e-4:1e-4)
    %Quitamos los ticks y labels del eje 'x' de las figuras
    set(gca,'xlabel',[],'XTick',[],'Fontsize',10,'FontName','Arial')
    box off
    %Retomamos el tick y label de la primera figura
    set(h(1),'XTick',25:25:250)
    xlabel(h(1),'Tiempo [ms]','Fontsize',15,'FontName','Arial','FontWeight', 'bold')
end
            

           
   
