function [ number ] = mass2number( Dp, mass, density )
% converts mass concentration to number consentration 
% Dp in m
% density in g/cm3
% mass in µg/m3
% number out #/cm3

rool = 1e3.*density; % now in kg/m3
M = 1e-9.*mass; % now in kg/m3

N = M./(rool.*pi./6.*Dp^3); % #/m3

number = 1e-6*N; % now in #/cm3

end

