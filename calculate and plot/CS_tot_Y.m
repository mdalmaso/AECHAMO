function [ CS_tot, CS_prime ] = CS_tot_Y( Y,nSec, t )
% Calculates CS vector from Y matrix, sectionnumber and time vector for each time-step 
% and plots CS(t)

alpha = 1.0;
T = 293;
yesmatrix = 0;

CS_tot = zeros(length(t),1);
CS_prime = zeros(length(t),1);

for i = 1:length(t)
    
    Ni = Y(i,2:nSec+1)';
    Dpi = Y(i,nSec+2:(2*nSec+1))';
    
    [CS_tot(i), CS_prime(i)] = CS_general(Dpi,Ni,T,alpha,yesmatrix);    

end

% Get the number of open figures:
figs=findall(0,'type','figure');
num_figs = length(figs);

% Open a new figure:
figure(num_figs+1);

plot(t/(24*60*60),CS_tot,'b*-', 'MarkerSize', 5, 'LineWidth',1)
handle1 = xlabel('time (d)');
set(handle1,'Fontsize',9,'Fontname','Computermodern')
handle2 = ylabel('CS (s^{-1})','rotation',90);
set(handle2,'Fontsize',9,'Fontname','Computermodern')

end

