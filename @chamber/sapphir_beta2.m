function[out] = sapphir_beta2(Dp,T)
% SAPPHIR_BETA2 calculates the wall losses for the sapphire chamber.
% 
% wall losses for the sapphire chamber: using the Anttila approach

% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0

out = 1.0.*sqrt(chamber.diff_particle(Dp,T)); % para from init
end