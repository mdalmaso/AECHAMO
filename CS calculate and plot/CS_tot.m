function [ CS_tot, t ] = CS_tot( distr )
% Calculates CS vector from numberdistribution for each time-step 
% and plots CS(t)

alpha = 1.0;
dmin = 1e-9;
dmax = 1e-6;
T = 293;
Dp = distr(1, 3:end);
t = distr(2:end, 1); %time vector

sz = size(t); % size of time vector
CS_tot = zeros(sz(1),1);

for i = 1:sz(1)
    
    dN = distr(i+1, 3:end);
    
    dCSlogDp = CS_general_dlog(Dp,dN,T,alpha); % vector

    CS_tot(i) = integrate_distribution(Dp,dCSlogDp,dmin,dmax);

end

figure(1)
plot(t/(24*60*60),CS_tot,'b*-', 'MarkerSize', 5, 'LineWidth',1)
handle1 = xlabel('time (d)');
set(handle1,'Fontsize',9,'Fontname','Computermodern')
handle2 = ylabel('CS (s^{-1})','rotation',90);
set(handle2,'Fontsize',9,'Fontname','Computermodern')

end

