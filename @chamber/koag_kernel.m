% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0

function[K] = koag_kernel(Dp1,Dp2,dens,T)

if length(Dp1)>1,
    fprintf('koag_kernel error: The first argument MUST be scalar!! \n')
    return
end


% setting the constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 	Mx=98.08;    %g/mol
% 	Mair=28.965; %g/mol
	p=1.0;   %atm
	kB=1.3807e-23; %J/K
    dens=dens*1e3;   % g/cm^3 -> kg/m^3
% 	dens=1000;
% 	Dair=19.7; %??
% 	Dx=51.96;   %??
% 	k = 8314.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = Dp1./2; 
r2 = Dp2./2; % big_R,

	lamda= (6.73e-8*T*(1+(110.4/T)))./(296*p*1.373); %[m]    
	myy= (1.832e-5.*(T.^(1.5))*406.4)./(5093.*(T+110.4));
	
    kn1=lamda./r1;               									%1*1
	kn2=lamda./r2;             									%28*1
	
    CC1= 1. + (kn1.*(1.142+(0.558.*exp((-.999)./kn1))));
    CC2= 1. + (kn2.*(1.142+(0.558.*exp((-.999)./kn2))));
	
    D1= (kB.*T.*CC1)./(6.*pi.*myy.*r1);                     		%1*1 koaguloituvan
	D2= (kB.*T.*CC2)./(6.*pi.*myy.*r2);                    %28*1 jakauman
    
    
    M2= 4./3.*pi.*(r2.^3).*dens;								%28*1
	M1= 4./3.*pi.*(r1.^3).*dens;											%1*1
    
   c2= sqrt((8.*kB.*T)./(pi.*M2));									%28*1
   c1= sqrt((8.*kB.*T)./(pi.*M1));								%1*1
   
   c12= sqrt((c2.^2)+(c1.^2));										%28*1
   
	r12= r2+r1;														%28*1
    
	D12= D2+D1;															%28*1
	
    CCONT= 4.*pi.*r12.*D12;											%28*1
    
	CFR= pi.*r12.*r12.*c12;											%28*1
    
	L2=(8.*D2)./(pi.*c2);												%28*1
	L1=(8.*D1)./(pi.*c1);												%1*1
    
	SIG2=(1./(3.*r12.*L2)).*((r12+L2).^3-(r12.*r12+L2.*L2).^1.5) - r12;    %28*1
	SIG1=(1./(3.*r12.*L1)).*((r12+L1).^3-(r12.*r12+L1.*L1).^1.5) - r12; 
    
	SIG12= sqrt((SIG2.^2)+(SIG1.^2));
    
   KO=CCONT./((r12./(r12+SIG12))+(CCONT./CFR));				%1*27  koagulaatiokerroin

   
   K = KO;