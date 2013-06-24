roo = 1.84; % Particle density (g/cm3)
M = 100; %g/mol
NA = 6.022e23; % 1/mol
alfa = 0.3;

%% calculate deltaMoa
i = 1;
Vtot = chamb_temp(i).output_data.Vtot; % m3
tim = chamb_temp(i).output_data.tim;


deltaVtot = Vtot - Vtot(1);
%deltaMoa2(i2) = chamb(i).output_data.Mtot(i2) - chamb(i).output_data.Mtot(1); 
Moa = roo.*1e6*Vtot;
deltaMoa = roo*1e6*deltaVtot + chamb_temp(i).output_data.Mdilu; % syntynyt aerosoli g/cm3 ilmaa 

h1=figure(1);
hold on;
% plot deltaMoa
plot(tim/(24*3600),deltaMoa,'b*')
%hold on;
%plot(tim,deltaMoa2,'m*')

%% calculate deltaP
k = 2e-16*60e-9*2.6908e19;
P = chamb_temp(i).initials.gas_source(1:end,2)/(k*alfa); %vector
kP = k.*P;

deltaP(1:length(tim),1) = 0;
deltaP_mass(1:length(tim),1) = 0;
for i2 = 2:length(tim)
deltaP(i2) = trapz(tim(1:i2),kP(1:i2)); % molkyyliä/(cm3)
deltaP_mass(i2) = deltaP(i2).*M./NA; % g/(cm3 ilmaa)
end

% hold on
% plot(tim,deltaP,'r*')
hold on;
plot(tim/(24*3600),deltaP_mass,'r')

%%
Y = deltaMoa./deltaP_mass;

h2=figure(2);
hold on;
plot(tim/(24*3600), Y, 'c*')

h3=figure(3);
hold on;
plot(deltaMoa/(24*3600),Y,'m*')

saveas(h1,'deltaP_deltaMoa).fig')
saveas(h2,'Y(t).fig')
saveas(h3,'Y(deltaMoa).fig')

plot(tim/(24*3600),chamb_temp(i).output_data.Mdilu,'c*')