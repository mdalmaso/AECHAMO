function [ mass ] = number2mass( Dp, number, density )
% converts mass concentration to number consentration 
% Dp in m
% density in g/cm3
% mass out µg/m3
% number in #/cm3

rool = 1e3.*density; % now in kg/m3
N = 1e6.*number; % now in #/m3

M = rool.*pi./6.*Dp^3*N; % kg/m3

mass = 1e9*M; % now in µg/m3

end

