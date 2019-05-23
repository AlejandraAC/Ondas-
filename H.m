
function [v2a,s12a,s32a] = H(v2a,s12a,s32a,mu,rho,eta,nx,nz,x1,x2,z1,z2)

%Ecuaciones SH
%Esta función resuelve las derivadas espaciales de los esfuerzos y la
%velocidad utilizando el método de diferencias finitas de cuarto orden

%----------------------------------------------------------------------
%Se resuelve la velocidad de la ecuación de movimiento para una onda SH 
v2b = zeros(nx,nz);
for i=3:nx-2
    for j=3:nz-2
        v2b(i,j)=v2a(i,j);
        %momentum conservation 
        % i-3/2---->i-1
        % i-1/2---->i
        % i+1/2---->i+1
        % i+3/2---->i+2
        ds4=z1*(s32a(i,j+1)-s32a(i,j))+z2*(s32a(i,j+2)-s32a(i,j-1));
        ds6=x1*(s12a(i+1,j)-s12a(i,j))+x2*(s12a(i+2,j)-s12a(i-1,j));
        %velocidad 
        v2a(i,j)=(ds4+ds6)/rho(i,j);
    end 
end 
%-------------------------------------------------------------------------

%Se resuelven los esfuerzos de las ecuaciones correspondientes al modelo 
%viscoelastico de Maxwell para una onda SH

for i=3:nx-2
    for j=3:nz-2
    %esfuerzos y deformaciones
    % i-3/2---> i-2
    % i-1/2-----> i-1
    % i+1/2-----> i
    % i+3/2-----> i+1
        e4=z1*(v2b(i,j)-v2b(i,j-1))+z2*(v2b(i,j+1)-v2b(i,j-2));
        e6=x1*(v2b(i,j)-v2b(i-1,j))+x2*(v2b(i+1,j)-v2b(i-2,j));
        s32a(i,j)=mu(i,j)*(e4-s32a(i,j)/eta(i,j));
        s12a(i,j)=mu(i,j)*(e6-s12a(i,j)/eta(i,j));
    end 
end


end

