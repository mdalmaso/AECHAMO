% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0

function[out] = log_normal(Dp_in,mu,sigma,N)
% function[out] = log_normal(Dp_in,mu,sigma,N)
% 
% defines a lognormal distribution over the diameter vector Dp_in

out = (N./(sqrt(2.*pi).*log10(sigma))).*exp(-((log10(Dp_in)-log10(mu)).^2)./(2.*(log10(sigma).^2)));
end