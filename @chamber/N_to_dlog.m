function[dN] = N_to_dlog(obj,Dp,N);
% function[dN] = N_to_dlog(Dp,N);
% Dp = centers
dlogdp = diff(log10(Dp)); % spacing in log
lDp = log10(Dp); % centers in log
log_eDp = [(lDp(1)-dlogdp(1)./2) lDp(1:(end-1))+(dlogdp./2) (lDp(end)+dlogdp(end)./2)]; % edges in log
eDp = 10.^log_eDp; % edges

edlogDp = diff(log_eDp);

dN = N./edlogDp;

Ntot_sum = sum(N);
Ntot_int = obj.integrate_distribution(Dp,dN,min(Dp),max(Dp));

%Ntot_sum./Ntot_int;
