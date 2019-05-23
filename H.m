
function [v2a,s12a,s32a] = H(v2a,s12a,s32a,mu,rho,eta,nx,nz,x1,x2,z1,z2)

%Ecuaciones SH/TM

%Esta función resuelve las derivadas espaciales del sistema de ecuaciones SH o TM 
%a través del método de diferencias finitas de cuarto orden 

%----------------------------------------------------------------------
%Se resuelve la velocidad de la ecuación de movimiento para una onda SH o
%el campo magnético de la ecuación de la Ley de Faraday para una onda TM

v2b = zeros(nx,nz);
for i=3:nx-2
    for j=3:nz-2
        v2b(i,j)=v2a(i,j);
        %cambios de malla
        % i-3/2---->i-1
        % i-1/2---->i
        % i+1/2---->i+1
        % i+3/2---->i+2
        %operadores diferenciales de esfuerzos
        ds4=z1*(s32a(i,j+1)-s32a(i,j))+z2*(s32a(i,j+2)-s32a(i,j-1));
        ds6=x1*(s12a(i+1,j)-s12a(i,j))+x2*(s12a(i+2,j)-s12a(i-1,j));
        %Resultado de velocidad o campo magnético. 
        v2a(i,j)=(ds4+ds6)/rho(i,j); %v2 <-> H2
    end 
end 
%-------------------------------------------------------------------------

%Se resuelven los esfuerzos de las ecuaciones correspondientes al modelo 
%viscoelastico de Maxwell para una onda SH o
%los valores de campo magnético de las ecuaciones correspondientes 
%a la Ley de Ampere para una onda TM

for i=3:nx-2
    for j=3:nz-2
    %cambios de malla
    % i-3/2---> i-2
    % i-1/2-----> i-1
    % i+1/2-----> i
    % i+3/2-----> i+1
        %operadores diferenciales de velocidad
        e4=z1*(v2b(i,j)-v2b(i,j-1))+z2*(v2b(i,j+1)-v2b(i,j-2));
        e6=x1*(v2b(i,j)-v2b(i-1,j))+x2*(v2b(i+1,j)-v2b(i-2,j));
        %Resultado de esfuerzos o valores de campo magnetico
        s32a(i,j)=mu(i,j)*(e4-s32a(i,j)/eta(i,j)); %s32 <-> (-E1) 
        s12a(i,j)=mu(i,j)*(e6-s12a(i,j)/eta(i,j)); %s12 <-> E3
    end 
end


end

