function [ Vtot, Ntot ] = Vtot_Ntot( Y, t, sections )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nSec = sections;
%initialize
Ntot = zeros(length(t),1);
Vtot = zeros(length(t),1);
Ni = zeros(1,nSec);
Dpi = zeros(1,nSec);

for i = 1:length(t),
    Ntot(i) = sum(Y(i,2:nSec+1)); % 1/cm3
    Ni = Y(i,2:nSec+1);
    Dpi = Y(i,nSec+2:(2*nSec+1));
    dV = (pi./6).*Dpi.^3.*Ni;
    Vtot(i) = sum(dV); 
end

end