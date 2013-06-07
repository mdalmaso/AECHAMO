function[out] = distribution_info_Vtot(Dp,dN);
% DISTRIBUTION_INFO_VTOT calculates the total volume of distribution.
% 
% distribution_info_Vtot(Dp, dN)

% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0

<<<<<<< HEAD
=======
function[out] = distribution_info_Vtot(obj,Dp,dN);

>>>>>>> origin/poikkimaki
% surface and volume distributions
dS = dN.*(pi.*Dp.^2);
dV = dN.*(pi./6).*(Dp.^3);

% VOLUME
% concentrations
out = obj.integrate_distribution(Dp,dV,min(Dp),max(Dp));


