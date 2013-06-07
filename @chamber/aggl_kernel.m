% (c) Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0

%Agglomeraattien koagulaationopeus vapaamolekyylialueella
function [K] = aggl_kernel(Dp1,Dp2,dens,T,Df,r0)

if length(Dp1)>1,
    fprintf('koag_kernel error: The first argument MUST be scalar!! \n')
    return
end


k=1.3806488e-23;%Boltzmann
dens=1000*dens; %Hiukkasen tiheys

r1=Dp1./2;
r2=Dp2./2;


v1=4./3.*pi.*r1.^3;
v2=4./3.*pi.*r2.^3;


lambda=2/Df-1/2;

%Koagulaatiokerroin
K=sqrt(6.*k.*T./dens).*(3./(4.*pi)).^lambda.*r0^(2-6./Df).*sqrt(1./v1+1./v2).*(v1.^(1./Df)+v2.^(1./Df)).^2;
end
