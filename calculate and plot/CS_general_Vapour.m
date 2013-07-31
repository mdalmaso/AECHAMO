function [CS, CS_prime,varargout] = CS_general_Vapour(Dp,N,T,alpha,yesmatrix)
% function [CS,CS_prime,varargout] = CS_general(Dp,N,T,alpha,yesmatrix)
%
%
% Version 0.1, last update: 02.7.2001
%
% This function calculates the total condensation sink
% value and also the CS size distribution if yesmatrix is
% nonzero.
% Input:
% * Dp		: diameter vector
% * N 		: (absolute) concentration vector corrseponding to Dp
%             NOTE: NOT dN/dlogDp !!!!!!!!
% * T			: temperature
% * alpha	: sticking coefficient
% * yesmatrix: [optional] output the cs size distribution
%
% Output:
% * CS		: the condensation sink [s^-1]
% * CS_prime: CS / (4*pi*Dif)	[cm^-2]
%
% CS calculated as per Pirjola: Effects of aerosol dynamics...
% (1998), J. Aerosol. Sci







% setting the constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Mx=300;    %g/mol   % molar mass of vapour
	Mair=28.965; %g/mol   %  - " -        air
	Pr=1.0;   %atm        % pressure
	Dair=19.7; % ??       % diffusion volume of air
	Dx=250;  % ??       % diffusion volume of vapour
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

CS=0;
CS_prime=0;
for j=1:length(N),
	if ~isnan(N(j))
		CS=CS + (4.*pi.*Dif).*N(j).*beta(j).*R(j).*1e2; % 1e2 comes for units
      CS_prime = CS_prime + N(j).*beta(j).*R(j).*1e2; % 1e2 comes for units
		
	   CS_m(j) = 4.*pi.*Dif.*N(j).*beta(j).*R(j).*1e2; % for the size distribution
      CS_m_prime = N(j).*beta(j).*R(j).*1e2; % for the size distribution
   end;
end;


% this outputs the size distribution if yesmatrix is given and nonzero

if nargin>4,
 	if yesmatrix~=0,
	 	varargout(1) = CS_m;
	   varargout(2) = CS_m_prime;
	end;
end;	 	
	






	