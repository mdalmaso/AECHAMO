function[out] = log_normal(Dp_in,mu,sigma,N)
% LOG_NORMAL creates a log-normal distribution.
% 
% function[out] = log_normal(Dp_in,mu,sigma,N)
% Creates a log-normal distribution based on mu, sigma and N. The
% distribution will have the same spacing as the diameter vector Dp_in.

% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0

out = (N./(sqrt(2.*pi).*log10(sigma))).*exp(-((log10(Dp_in)-log10(mu)).^2)./(2.*(log10(sigma).^2)));

end