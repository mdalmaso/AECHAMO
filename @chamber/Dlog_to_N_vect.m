function [Dp, N]=Dlog_to_N_vect(Dp, dN)
% DLOG_TO_N_VECT calculates the concentration of particles in distribution.
% 
% [Dp, N] = Dlog_to_N_vect(Dp, dN)
% Gets the concentration of particles (N) in each section of
% distribution [Dp, dN].
% 
% N(i) is the concentration of particles in section i.

% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0



dN = dN(:)';




lkm=length(Dp);
vali=ones(lkm,1);

vpist=zeros(1,lkm); % Preallocate

for bi=2:lkm
   vpist(bi)=Dp(bi-1)+0.5*(Dp(bi)-Dp(bi-1));
end;

vpist(lkm+1)=Dp(lkm)+(Dp(lkm)-vpist(lkm));
vpist(1)=Dp(1)-(vpist(2)-Dp(1));
logpist=log10(vpist);
vali(:)=diff(logpist);           % dlogDp  jonkin verran approksimoiden

N = dN.*vali(:)';
end