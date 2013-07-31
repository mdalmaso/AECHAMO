function[out] = M_dens(P,T)
% function[out] = M_dens(P,T)
% returns the amount of air molecules in a cubic centimeter of air 
% at given pressure and temperature
% P given in mbar (1013 for 1 atm)
% T given in Kelvin
% 
% Written by M. Dal Maso 2007

out = P./(T.*1.379e-19);