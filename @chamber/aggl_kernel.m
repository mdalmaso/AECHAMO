function [K] = aggl_kernel(obj, Dp1, Dp2, dens, T, Df, r0)
% AGGL_KERNEL calculates the coagulation coefficients for agglomerates.
% [K] = chamber.aggl_kernel(Dp1, Dp2, dens, T, Df, r0)
% 
% Dp1 is a scalar volume-diameter of agglomerate
% Dp2 is a vector containing volume-diameters of all interacting
% agglomerates
% dens is the particulate matter density in g/cm3
% T is the temperature in Kelvins
% Df is the fractal dimension of agglomerates
% r0 is the primary particle diameter

% Based on Yu & Lin: "Taylor-expansion moment method for agglomerate 
% coagulation due to Brownian motion in the entire size regime." 2009

% (c) Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0 Free molecule regime
% 2013-08-05    0.2.0 Extended to work for all size regimes.

if length(Dp1)>1,
    error('The first argument MUST be scalar!!');
end

%%
% If both Dp1 and Dp2 are smaller than Dp0 (primary particle diameter), 
% the particles cannot be agglomerates, so use function koag_kernel
% instead.
if(Dp1 < 2*r0)
    if(isscalar(Dp2))
        if(Dp2 > 2*r0)
            aggl_begin = 1;
        else
            aggl_begin = 2;
        end
    else
        aggl_begin = find(Dp2 >= 2*r0);
    end
else
    aggl_begin = 1;
end

K = zeros(length(Dp2),1);

if(aggl_begin > 1)
    K(1:aggl_begin-1) = obj.koag_kernel(Dp1,Dp2(1:aggl_begin-1),dens,T);
end


%% Coagulation kernel for continuum regime
k=1.3806488e-23; %Boltzmann
% mu = 18.27e-6; % Gas viscosity
% lambda = 66e-9; % Mean free path
A = 1.591; % Constant
p = 1.0; %atm

dens = dens * 1000; % g/cm^3 to kg/m^3

lambda= (6.73e-8*T*(1+(110.4/T)))./(296*p*1.373); %[m]    
mu= (1.832e-5.*(T.^(1.5))*406.4)./(5093.*(T+110.4));

v1 = pi./6.*Dp1.^3;
v2 = pi./6.*Dp2.^3;


B2 = 2*k*T/(3*mu);
f=1/Df;
fii = lambda*A*(4*pi/3)^(1/3);
vp0 = 4/3*pi*r0^3;



K(aggl_begin:end) = B2.*(1./v1.^f + 1./v2(aggl_begin:end).^f).*(v1.^f+v2(aggl_begin:end).^f)+B2.*fii.*vp0.^(f-1/3).*(1./v1.^(2.*f)+1./v2(aggl_begin:end).^(2.*f)).*(v1.^f+v2(aggl_begin:end).^f);
 
%% Correction for free-molecule and transition regimes


% Dahneke's coefficients:
% B1=2/3;
% B2=4/3;
% B3=8/9;

% Calculate collision diameter using the power law:
kA= 1.0; % If the value of coefficient kA is changed, it should be taken
         % in account also when deriving the formula for coagulation kernel
         % in continuum range above.
Dp1 = 2.*r0.*(v1./(vp0.*kA)).^(1/Df);
Dp2 = 2.*r0.*(v2./(vp0.*kA)).^(1/Df);
Dp12 = Dp1 + Dp2;

D1 = diff_coeff(Dp1, T); % Diffusion coefficient. Diffusion is calculated
                         % assuming spherical particles. However, the
                         % collision diameter is used instead of volume
                         % diameter.
D2 = diff_coeff(Dp2, T);
D12 = D1 + D2;

m1 = v1.*dens; % mass
m2 = v2.*dens;

c1 = sqrt((8*k*T)/(m1*pi)); % mean thermal velocity
c2 = sqrt((8.*k.*T)./(m2.*pi));

c12 = sqrt(c1.^2 + c2.^2);

lambda_p = 3.*D12./c12;

Kn = (2.*lambda_p)./Dp12;

r1 = Dp1/2;
r2 = Dp2./2;
r12 = r1 + r2;

Kn_delta = 8.*D12./(pi.*c12.*r12);
delta12 = r12./Kn_delta.^2.*(1./5.*(1+Kn_delta).^5-1./3.*(1+Kn_delta).^3.*(1+Kn_delta.^2)+2/15.*(1+Kn_delta.^2).^(5./2))-r12;

B1 = delta12./lambda_p;
B2 = 4/3;
B3 = (4.*delta12)./(3.*lambda_p);

correction = (1+B1.*Kn)./(1+B2.*Kn+B3.*Kn.^2);

K(aggl_begin:end) = K(aggl_begin:end).*correction(aggl_begin:end);

end

%%
function [ coeff ] = diff_coeff(Dp,T)
%DIFF_COEFF Diffusion coefficient for spherical particles.
%   Calculation based on Friedlander 2000
p = 1.0; %atm

lambda= (6.73e-8*T*(1+(110.4/T)))./(296*p*1.373); %[m]    
mu= (1.832e-5.*(T.^(1.5))*406.4)./(5093.*(T+110.4));

r=Dp./2;
kn=lambda./r;
CC= 1. + (kn.*(1.142+(0.558.*exp((-.999)./kn))));


f=3.*pi.*mu.*Dp./CC; % friction coefficient

k=1.3806488e-23; %Boltzmann

coeff = k.*T./f;


end



%%
% function [K] = aggl_kernel(Dp1,Dp2,dens,T,Df,r0)
% % AGGL_KERNEL calculates the coagulation coefficients for agglomerates.
% % 
% % aggl_kernel(Dp1, Dp2, dens, T, Df, r0)
% 
% % (c) Pauli Simonen 2013
% %
% % Version history:
% % 2013-05-24    0.1.0
% 
% if length(Dp1)>1,
%     fprintf('koag_kernel error: The first argument MUST be scalar!! \n')
%     return
% end
% 
% 
% k=1.3806488e-23;%Boltzmann
% dens=1000*dens; %Hiukkasen tiheys
% 
% r1=Dp1./2;
% r2=Dp2./2;
% 
% 
% v1=4./3.*pi.*r1.^3;
% v2=4./3.*pi.*r2.^3;
% 
% 
% lambda=2/Df-1/2;
% 
% %Koagulaatiokerroin
% K=sqrt(6.*k.*T./dens).*(3./(4.*pi)).^lambda.*r0^(2-6./Df).*sqrt(1./v1+1./v2).*(v1.^(1./Df)+v2.^(1./Df)).^2;
% end
