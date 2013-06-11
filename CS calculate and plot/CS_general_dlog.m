function [dCS, dCS_prime] = CS_general_dlog(Dp,dN,T,alpha)
%% function [dCS, dCS_prime] = CS_general(Dp,dN,T,alpha)
%%
%%
%% Version 0.5, last update: 04.02.2004
%%(c) Miikka Dal Maso
%%  calculates dCS/dlogDp if dN/dlogDp is known;
%%
%% Input:
%% * Dp		: diameter vector
%% * N 		: particle concentration vector as dN/dlogDp
%% * T		: temperature in KELVIN
%% * alpha	: sticking coefficient (usually using 1.0)
%%
%% Output:
%% * dCS		: the condensation sink [s^-1]
%% * sCS_prime: CS / (4*pi*Dif)	[cm^-2]
%%
%% CS calculated as per Pirjola: Effects of aerosol dynamics...
%% (1998), J. Aerosol. Sci
%%
%% notes: n_CS(log Dp) = 4.*pi.*D.*beta(Dp).*(Dp./2).*n_N(log Dp)
%% so in this formula is the RADIUS!! not diameter even if the 
%% size distribution is given as a function of diameter 

% setting the constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Mx=98.08;    %g/mol   % molar mass of H2SO4
	Mair=28.965; %g/mol   %  - " -        air
	Pr=1.0;   %atm        % pressure
	Dair=19.7; % ??       % diffusion volume of air
	Dx=51.96;  % ??       % diffusion volume of H2SO4
	k = 8314.7;           % actually R

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = Dp./2;   % diameter to radius


% The diffusion coefficient [cm^2/s]% from Reid et al.
Temp = T;	

Dif = (0.001 * (Temp.^1.75)*sqrt( (1./Mair)+(1./Mx))) / (Pr.*(Dair.^(1/3)+Dx.^(1/3)).^2);
	
% lambda [m]%
lam=3*(sqrt( (pi.*Mx)/(8.*k.*Temp) )) .* Dif .*1e-4;
	
% the knudsen number [dimensionless]%
knud=lam./R;
	
% beta, the correction coefficient [dimensionless]%
beta=(knud+1)./((0.377.*knud)+1+(4/(3.*alpha)).*(knud.^2)+(4/(3.*alpha)).*knud);

% the total CS calculation loop, cs in [1/s] %
% CS_prime is given in cm^-2;
% done this way to avoid total NaN:ing if only some values are NaN

for j=1:length(dN),
	if ~isnan(dN(j))
		dCS(j) = 4.*pi.*Dif.*dN(j).*beta(j).*R(j).*1e2; % for the size distribution
		dCS_prime = dN(j).*beta(j).*R(j).*1e2; % for the size distribution
	else
		dCS(j) = 0; % this is not very good; however, the error is not too bad
	end;
end;




	
