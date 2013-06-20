roo = 1.84; % Particle density (g/cm3)
M = 100; %g/mol
NA = 6.022e23; % 1/mol

% calculate deltaMoa
i = 1;
Vtot = chamb(i).output_data.Vtot; % m3
tim = chamb(i).output_data.tim;


deltaVtot = Vtot - Vtot(1);
%deltaMoa2(i2) = chamb(i).output_data.Mtot(i2) - chamb(i).output_data.Mtot(1); 
Moa = roo.*1e6*Vtot;
deltaMoa = roo*1e6*deltaVtot; % g/cm3 ilmaa

h1=figure(1);
% plot deltaMoa
plot(tim,deltaMoa,'b*')
%hold on;
%plot(tim,deltaMoa2,'m*')

%calculate deltaP
k = 9e-17*30e-9*2.6908e19;
P = 1.0*1e-9*2.6908e19;

deltaP = k.*P.*tim; % molkyyliä/(cm3)
deltaP_mass = deltaP.*M./NA; % g/(cm3 ilmaa)
% 
% hold on
% plot(tim,deltaP,'r*')
hold on;
plot(tim,deltaP_mass,'r')

Y = deltaMoa./deltaP_mass;

h2=figure(2);
hold on;
plot(tim, Y, 'c*')

h3=figure(3);
hold on;
plot(deltaMoa,Y,'m*')

saveas(h1,'deltaP_deltaMoa).fig')
saveas(h2,'Y(t).fig')
saveas(h2,'Y(deltaMoa).fig')