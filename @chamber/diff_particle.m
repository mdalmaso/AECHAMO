function [Diffcoeff] = diff_particle(Dp,T)
% function [CoagS,varargout] = CoagS_dR(Dp_small,Dp,N,T)
%
%
% Version 0.1, last update: 02.7.2001
%
% This function calculates the total coagulation sink
% value for a particle of size Dp_small
% Input:
% * Dp		: diameter vector. MUST INCREASE MONOTONICALLY!!
% * N 		: (absolute) concentration vector corrseponding to Dp
%             NOTE: NOT dN/dlogDp !!!!!!!!
% * T			: temperature in kelvins
% Output:
% * CoagsS		: the coagulation sink
%
% CoagsS calculated as per Pirjola code (hope it works :-)
%
% The sink is also calculated assuming that only particles with sizes
% > D_small contributing to the sink.
%

if T<100,
    error('Temperature should be in Kelvins (%iK  is VERY cold.)!',T)
    return
end

%% checking that the Dp vector increases monotonically
if any(diff(Dp)<0),
	fprintf('Caogulation Sink calculation error! Diameter vector not increasing monotonically!!');
	return;
end;

% setting the constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Mx=98.08;    %g/mol
	Mair=28.965; %g/mol
	p=1.0;   %atm
	kB=1.3807e-23; %J/K
	dens=1000;
	Dair=19.7; %??
	Dx=51.96;   %??
	k = 8314.7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N = zeros(size(Dp));
Dp_small = 1; 



big_R = Dp./2;   % diameter to radius
big_N = N;

%% find the nearest size bin for Dp_small
diff_vect = abs(Dp-Dp_small);
[m m_ix]  = min(diff_vect);               % m_ix is now the size for which the 													
														% coag. coeff. is multiplied by 0.5;


										
														
														
														
														
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    This part is from Liisa Pirjola. Beginning here...            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r0 = Dp_small./2;


	r1=r0;               %r0 on koaguloituvat, r on jakauman säteet
	lamda= (6.73e-8*T*(1+(110.4/T)))/(296*p*1.373);
	myy= (1.832e-5*(T^(1.5))*406.4)/(5093*(T+110.4));
	kn1=lamda/r1;               									%1*1
	kn=lamda./big_R;             									%28*1
	CC= 1. + (kn.*(1.142+(0.558.*exp((-.999)./kn))));     %28*1
	CC1= 1. + (kn1.*(1.142+(0.558.*exp((-.999)./kn1))));  %1*1
	D= (kB.*T.*CC)./(6.*pi.*myy.*big_R);                    %28*1 jakauman
	D1= (kB*T*CC1)/(6*pi*myy*r1);                     		%1*1 koaguloituvan
    Diffcoeff = D;
    
    return
    
    M= 4./3.*pi.*(big_R.^3).*dens;								%28*1
	M1= 4/3*pi*(r1^3)*dens;											%1*1
	c= sqrt((8.*kB.*T)./(pi.*M));									%28*1
   c1= sqrt((8.*kB.*T)/(pi.*M1));								%1*1
   c12= sqrt((c.^2)+(c1.^2));										%28*1
	r12= big_R+r1;														%28*1
	D12= D+D1;															%28*1
	CCONT= 4.*pi.*r12.*D12;											%28*1
	CFR= pi.*r12.*r12.*c12;											%28*1
	L= (8.*D)./(pi.*c);												%28*1
	L1=(8*D1)./(pi*c1);												%1*1
	SIG=(1./(3.*r12.*L)).*((r12+L).^3-(r12.*r12+L.*L).^1.5) - r12;    %28*1
	SIG1=(1./(3.*r12.*L1)).*((r12+L1).^3-(r12.*r12+L1.*L1).^1.5) - r12; %1*1
	SIG12= sqrt((SIG.^2)+(SIG1.^2));
   KO=CCONT./((r12./(r12+SIG12))+(CCONT./CFR));				%1*27  koagulaatiokerroin
   coags=0;
   for q=1:length(big_R);
      if ~isnan(big_N(q)),
         if q == m_ix,                  % if sizes equal, factor 0.5
            coags=coags + (0.5.*KO(q).*big_N(q).*1e6);	
         elseif q < m_ix,               % if smaller, no contribution
            coags = coags;
         else         	                % else full contribution
         	coags=coags + (KO(q).*big_N(q).*1e6);
         end;
      end;
   end;			

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    This part is from Liisa Pirjola. ...ending here               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 	
	
CoagS = coags;


	
