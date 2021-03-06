
function [f,nw2] = wavelet(dt2,freq)

%Esta función crea una fuente impulso que genera una ondícula gaussiana establecida por J. Carcione, (2006)

t0=(6/(5*freq)); %tiempo de retraso o fase
wb=2*pi*freq; %frecuencia angular central
nw2=(2*t0/dt2); %numero de muestras del pulso
Dw=0.5*wb; %ancho del pulso

% Se inicializa la fuente
f = zeros(1,nw2);
%Ciclo para obtener el historial de tiempo de la fuente
for n=1:nw2
    t=(n-1)*dt2;
    D=t-t0;
    f(n)=exp(-(Dw*Dw*D*D)/4)*cos(wb*D);
end

%CASO VISCOELÁSTICO
    %Se crea un vector con la escala en tiempo
 %   vect = [0:dt2:(nw2-1)*dt2];
    %Se grafica la fuente
  %   plot(vect,f)
   %  title(['Ondícula fuente con F = ' num2str(freq) ' [Hz] '],'Fontsize',19,'FontName','Arial', 'FontWeight', 'bold','HorizontalAlignment','center')
     %Nombre y tamaño de ejes
    % xlabel('Tiempo [s]','Fontsize',15,'FontWeight','bold' )
    % ylabel('Amplitud','Fontsize',15,'FontWeight', 'bold')


%CASO ELECTROMAGNÉTICO
    %Se reescribe la frecuencia en MHz 
    fh=freq*(10^-6);
    %Se crea un vector con la escala en tiempo
    vect = [0:dt2:(nw2-1)*dt2]*(10^9);
   % Se grafica la fuente
    plot(vect,f)
    title(['Ondícula fuente con F = ' num2str(fh) ' [MHz] '],'Fontsize',19,'FontName','Arial', 'FontWeight', 'bold','HorizontalAlignment','center')
    %Nombre y tamaño de ejes
    xlabel('Tiempo [ns]','Fontsize',15,'FontWeight','bold' ) 
    ylabel('Amplitud','Fontsize',15,'FontWeight', 'bold')
end


