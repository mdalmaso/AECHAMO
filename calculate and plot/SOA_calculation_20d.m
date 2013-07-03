roo = 1.84; % Particle density (g/cm3)
M = 100; %g/mol
NA = 6.022e23; % 1/mol
alfa = 0.3;

%% calculate and plot deltaMoa
i = 18;
Vtot = chamb(i).output_data.Vtot; % m3
tim = chamb(i).output_data.tim;


deltaVtot = Vtot - Vtot(1);
%deltaMoa2(i2) = chamb(i).output_data.Mtot(i2) - chamb(i).output_data.Mtot(1); 
Moa = roo.*1e6*Vtot;
deltaMoa = roo*1e6*deltaVtot + chamb(i).output_data.Mdilu; % syntynyt aerosoli g/cm3 ilmaa 

h1=figure(1);
hold on;
% plot deltaMoa
plot(tim/(24*3600),deltaMoa,'b*')
%hold on;
%plot(tim,deltaMoa2,'m*')

%% calculate and plot deltaP
if isscalar(chamb(i).initials.gas_source) == 0 
    
    kP = chamb(i).initials.gas_source(1:end,2)/alfa; %vector   

    deltaP(1:length(tim),1) = 0;
    deltaP_mass(1:length(tim),1) = 0;
    for i2 = 2:length(tim)
        deltaP(i2) = trapz(tim(1:i2),kP(1:i2)); % molekyyliä/(cm3)
        deltaP_mass(i2) = deltaP(i2).*M./NA; % g/(cm3 ilmaa)
    end

else
    kP = chamb(i).initials.gas_source/alfa;
    deltaP = kP.*tim; % molkyyliä/(cm3)
    deltaP_mass = deltaP.*M./NA; % g/(cm3 ilmaa)    
end

% hold on
% plot(tim,deltaP,'r*')
hold on;
plot(tim/(24*3600),deltaP_mass,'r')

%% calculate Y

for i3 = 1:length(deltaP_mass)    
    % if deltaP_mass = 0, Y = inf or NaN
    if deltaP_mass(i3) ~= 0
        Y(i3) = deltaMoa(i3)/deltaP_mass(i3);
    else
        Y(i3) = 0;
    end    
end

%% plot Y(t) and Y(Moa)
h2=figure(2);
hold on;
plot(tim/(24*3600), Y, 'c*')

h3=figure(3);
hold on;
plot(deltaMoa/(24*3600),Y,'m*')

% h4=figure(4);
% hold on;
% plot(tim/(24*3600),chamb(i).output_data.Mdilu,'c*')

%% calculate CS, Yend and plot Yend

CS = CS_tot_Y( chamb(i).output_data.Y, chamb(i).initials.sections, tim );

if chamb(i).initials.vap_wallsink_on ~= 0
    Yend = alfa./(1+(chamb(i).initials.vap_wallsink)./CS);
else
    Yend = alfa;
end

figure(2);
hold on;
plot(tim/(24*3600), Yend, 'r.')

%% calculate (Yend-Y)/Y and plot it

Yend = Yend';
error_of_Yend = (Yend-Y)./Y;
h7 = figure(7);
plot(tim/(24*3600),error_of_Yend,'*')

%% calculate and plot fraction of ELVOC forming aerosol = Y and fraction of ELVOC lost to wall Mvwall/deltaP

for i4 = 1:length(deltaP_mass)    
    % if deltaP_mass = 0, Y = inf or NaN
    if deltaP_mass(i4) ~= 0
        Wall_loss(i4) = chamb(i).output_data.Mvwall(i4)/deltaP_mass(i4);
    else
        Wall_loss(i4) = 0;
    end    
end

h5=figure(5);
plot(tim/(24*3600), Y, 'g*')
hold on;
plot(tim/(24*3600), Wall_loss, 'r*')

h6=figure(6);
plot(tim/(24*3600), Y/alfa, 'g*')
hold on;
plot(tim/(24*3600), Wall_loss/alfa, 'r*')

%% create folder and save pictures into it

name = 'Run0';
str = num2str(i);
new_name = strrep(name, '0', str);
mkdir(new_name)

saveas(h1,[new_name '/deltaP_deltaMoa.fig'])
saveas(h2,[new_name '/Y(t).fig'])
saveas(h3,[new_name '/Y(deltaMoa).fig'])
saveas(h5,[new_name '/YandMvwall.fig'])
saveas(h6,[new_name '/YandMvwall_alfa.fig'])
saveas(h7,[new_name '/error_of_Yend.fig'])

% save CS
h4 = figure(4);
saveas(h4,[new_name '/CS(t).fig'])