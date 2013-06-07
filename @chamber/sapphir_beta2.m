% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0


function[out] = sapphir_beta2(Dp,T)
% wall losses for the sapphire chamber: using the Anttila approach

out = 1.0.*sqrt(diff_particle(Dp,T)); % para from init
end