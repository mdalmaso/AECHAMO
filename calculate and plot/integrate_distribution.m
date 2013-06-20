function[out] = integrate_distribution(Dp,dN,dmin,dmax)

if dmin<min(Dp),
	dmin=min(Dp);
end;
if dmax>max(Dp),
	dmax=max(Dp);
end;

DpX = dmin:1e-10:dmax;
dNX = interp1(Dp,dN,DpX);          % data
lDpX = length(DpX);
out = trapz(log10(DpX),dNX);