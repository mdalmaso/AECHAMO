clear all;
close all;

% constants
roo = 1.84; % Particle density (g/cm3)
M = 100; %g/mol
NA = 6.022e23; % 1/mol
alfa = 0.3;

% initiate variables 
count1 = 1;
count2 = 1;

Yend_real_save_90s = zeros(1,12);
CSsave_90s = zeros(1,12);
Vtot_end_save_90s = zeros(1,12);

Yend_real_save_900s = zeros(1,10);
CSsave_900s = zeros(1,10);
Vtot_end_save_900s = zeros(1,10);

h_CS = zeros(1,13);
h2_CS =zeros(1,11);
h_Vtot = zeros(1,13);
h2_Vtot = zeros(1,11);

for file = 1:3
% data files and runs for calc
if file == 1
    runs = [5:1:8, 17,18,19,20];
    load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\run_20130628T122218\run_20130628T122218.mat')
elseif file == 2
    runs = [5:1:8, 17,18,19,20];
    load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\run_20130627T171225\run_20130627T171225.mat')
elseif file == 3
    runs = [5:1:8, 17,18];
    load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\run_20130629T085050\run_20130629T085050.mat')
end
%% i is index of run
for i = runs    
%% calculate deltaMoa
Vtot = chamb(i).output_data.Vtot; % m3
tim = chamb(i).output_data.tim;
deltaVtot = Vtot - Vtot(1);
%deltaMoa2(i2) = chamb(i).output_data.Mtot(i2) - chamb(i).output_data.Mtot(1); 
Moa = roo.*1e6*Vtot;
deltaMoa = roo*1e6*deltaVtot + chamb(i).output_data.Mdilu; % syntynyt aerosoli g/cm3 ilmaa 

%% calculate deltaP
if isscalar(chamb(i).initials.gas_source) == 0 
    
    kP = chamb(i).initials.gas_source(1:end,2)/alfa; %vector   
   
    deltaP = zeros(length(tim),1);
    deltaP_mass = zeros(length(tim),1);
    for i2 = 2:length(tim)
        deltaP(i2) = trapz(tim(1:i2),kP(1:i2)); % molekyyliä/(cm3)
        deltaP_mass(i2) = deltaP(i2).*M./NA; % g/(cm3 ilmaa)
    end

else
    kP = chamb(i).initials.gas_source/alfa;
    deltaP = kP.*tim; % molkyyliä/(cm3)
    deltaP_mass = deltaP.*M./NA; % g/(cm3 ilmaa)    
end

%% calculate Y
Y = zeros(1,length(deltaP_mass));
for i3 = 1:length(deltaP_mass)    
    % if deltaP_mass = 0, Y = inf or NaN
    if deltaP_mass(i3) ~= 0
        Y(i3) = deltaMoa(i3)/deltaP_mass(i3);
    else
        Y(i3) = 0;
    end    
end

%% calculate CS and Yend
CS = CS_tot_Y( chamb(i).output_data.Y, chamb(i).initials.sections, tim, 0);

if chamb(i).initials.vap_wallsink_on ~= 0
    Yend = alfa./(1+(chamb(i).initials.vap_wallsink)./CS);
else
    Yend = alfa;
end

%% calculate (Yend-Y)/Y
Yend = Yend';
error_of_Yend = (Yend-Y)./Y;

%% calculate and loglog fraction of ELVOC forming aerosol = Y and fraction of ELVOC lost to wall Mvwall/deltaP
Wall_loss = zeros(1,length(deltaP_mass));
for i4 = 1:length(deltaP_mass)    
    % if deltaP_mass = 0, Y = inf or NaN
    if deltaP_mass(i4) ~= 0
        Wall_loss(i4) = chamb(i).output_data.Mvwall(i4)/deltaP_mass(i4);
    else
        Wall_loss(i4) = 0;
    end    
end

%% calculate Yend_real(CSend) and Yend_real(Vtot) loglog
% file is current file, i is current run and tau is time place of Ymax
% correct times tau for different runs and files
if file == 1
    if i < 10
        tau = 81; % 0.5556days
    elseif i > 10
        tau = 23; % 0.1528days
    end
elseif file == 2
    tau = 81; % 0.5556days 
elseif file == 3
    tau = 88; % 0.6042days
end

% what is loglog symbol (mark) its different for runs and files and same
% for diff gamma but same other values
if file == 1
    if (i == 5) || (i == 7)
        mark = '.';        
    elseif (i == 6) || (i == 8)
        mark = 'o';
    elseif (i == 17) || (i == 19)
        mark = 'x';
    elseif (i == 18) || (i == 20)  
        mark = '+';
    end
elseif file == 2
    if (i == 5) || (i == 7)
        mark = '*';
    elseif (i == 6) || (i == 8)
        mark = 's';
    elseif (i == 17) || (i == 19)
        mark = 'd';
    elseif (i == 18) || (i == 20)  
        mark = 'v';
    end 
elseif file == 3
    if (i == 5) || (i == 7)
        mark = 'p';
    elseif (i == 6) || (i == 8)
        mark = 'h';
    elseif i == 17
        mark = '^';
    elseif i == 18  
        mark = '<';
    end
end % mark

% calc Yend and Vtot if gamma = 1/90s
if (5<=i && i<=6) || (17<=i && i<=18)     
    % Yend at right time and save every value to vector
    Yend_real_90s = Y(tau);
    Yend_real_save_90s(count1) = Yend_real_90s;
    % same for CSend
    CSend_90s = CS(tau);  
    CSsave_90s(count1) = CSend_90s;     
    % same for Vtot_end
    Vtot_end_90s = Vtot(tau);
    Vtot_end_save_90s(count1) = Vtot_end_90s;
    
    % loglog Yend(CSend)(one point each time)
    h8 = figure(8);
    h_CS(count1) = loglog(CSend_90s, Yend_real_90s, mark);      
    hold on;
    % loglog Yend(Vtot_end)(one point each time)
    h10 = figure(10);
    h_Vtot(count1) = loglog(Vtot_end_90s, Yend_real_90s, mark);      
    hold on;
    
    count1 = count1 + 1;     
    
% calc Yend and Vtot if gamma = 1/900s
elseif (7<=i && i<=8) || (19<=i && i<=20)
    % Yend at right time and save every value to vector
    Yend_real_900s = Y(tau);
    Yend_real_save_900s(count2) = Yend_real_900s;
    % same for CS
    CSend_900s = CS(tau);    
    CSsave_900s(count2) = CSend_900s;  
    % same for Vtot_end
    Vtot_end_900s = Vtot(tau);
    Vtot_end_save_900s(count2) = Vtot_end_900s;
    
    % loglog Yend(CSend) (one point each time)
    h9 = figure(9);     
    h2_CS(count2) = loglog(CSend_900s, Yend_real_900s, mark);       
    hold on;
    % loglog Yend(Vtot_end)(one point each time)
    h11 = figure(11);
    h2_Vtot(count2) = loglog(Vtot_end_900s, Yend_real_900s, mark);      
    hold on;
    
    count2 = count2 + 1;
end

%% loglog
% % loglog deltaMoa and deltaP
% h1=figure(2);
% hold on;
% loglog(tim/(24*3600),deltaMoa,'b*')
% %hold on;
% %loglog(tim,deltaMoa2,'m*')
% % hold on
% % loglog(tim,deltaP,'r*')
% hold on;
% loglog(tim/(24*3600),deltaP_mass,'r')
% 
% % loglog Y(t) and Yend
% h2=figure(3);
% hold on;
% loglog(tim/(24*3600), Y, 'c*')
% hold on;
% loglog(tim/(24*3600), Yend, 'r.')
% 
% % loglog Y(Moa)
% h3=figure(4);
% hold on;
% loglog(deltaMoa/(24*3600),Y,'m*')
% 
% % h4=figure(4);
% % hold on;
% % loglog(tim/(24*3600),chamb(i).output_data.Mdilu,'c*')
% 
% h5=figure(5);
% loglog(tim/(24*3600), Y, 'g*')
% hold on;
% loglog(tim/(24*3600), Wall_loss, 'r*')
% 
% h6=figure(6);
% loglog(tim/(24*3600), Y/alfa, 'g*')
% hold on;
% loglog(tim/(24*3600), Wall_loss/alfa, 'r*')
% 
% % loglog (Yend-Y)/Y
% h7 = figure(7);
% loglog(tim/(24*3600),error_of_Yend,'*')

%% create folder and save pictures into it
% % make directory
% name = 'Run0';
% str = num2str(i);
% new_name = strrep(name, '0', str);
% mkdir(new_name)
% 
% saveas(h1,[new_name '/deltaP_deltaMoa.fig'])
% saveas(h2,[new_name '/Y(t).fig'])
% saveas(h3,[new_name '/Y(deltaMoa).fig'])
% saveas(h5,[new_name '/YandMvwall.fig'])
% saveas(h6,[new_name '/YandMvwall_alfa.fig'])
% saveas(h7,[new_name '/error_of_Yend.fig'])
% 
% % save CS
% h4 = figure(1);
% saveas(h4,[new_name '/CS(t).fig'])
end % i
end % file

%% edit figures Yend_real(CSend)

% add theoretical loglog to fig8
figure(8);
hold on;
CSarea_90s = 1e-5:max(CSsave_90s)/100:max(CSsave_90s)*1.1;
Yend_kaava_90s = alfa./(1+(1/90)./CSarea_90s);
h_CS(count1) = loglog(CSarea_90s, Yend_kaava_90s, 'r');
% fit data
fitted_90s = fit_formula_mp(CSsave_90s',Yend_real_save_90s',max(CSsave_90s),0);
% edit legend 
hleg1 = legend([h_CS,(fitted_90s.pict_fit)'],'1','2','3','4','5','6','7','8','9','10','11','12','\alpha /[1+( \gamma /CSend)]',fitted_90s.leg_name1,fitted_90s.leg_name2);
set(hleg1,'Location','EastOutside')
% add labels
xhandle90 = xlabel('CS_{end} (s^{-1})');
yhandle90 = ylabel('Y_{end}','rotation',90); 
% add title
title('\alpha = 0.3 and \gamma = 1/90s');

% add theoretical loglog to fig9
figure(9);
hold on;
CSarea_900s = 1e-5:max(CSsave_900s)/100:max(CSsave_900s)*1.1;
Yend_kaava_900s = alfa./(1+(1/900)./CSarea_900s);
h2_CS(count2) = loglog(CSarea_900s, Yend_kaava_900s, 'r');
% fit data
fitted_900s = fit_formula_mp(CSsave_900s',Yend_real_save_900s',max(CSsave_900s),0);
% edit legend 
hleg2 = legend([h2_CS,(fitted_90s.pict_fit)'],'1','2','3','4','5','6','7','8','9','10','\alpha /[1+( \gamma /CSend)]',fitted_900s.leg_name1,fitted_900s.leg_name2);
set(hleg2,'Location','EastOutside')
% add labels
xhandle900 = xlabel('CS_{end} (s^{-1})');
yhandle900 = ylabel('Y_{end}','rotation',90); 
% add title
title('\alpha = 0.3 and \gamma = 1/900s');

% add theoretical loglog to fig10
figure(10);
hold on;
Vtot_area_90s = 1e-22:max(Vtot_end_save_90s)/100:max(Vtot_end_save_90s)*1.1;
Yend_kaava_Vtot_90s = alfa./(1+(1/90)./((roo.*1e6.*Vtot_area_90s).^0.63));
h_Vtot(count1) = loglog(Vtot_area_90s, Yend_kaava_Vtot_90s, 'r');
% fit data
fitted_Vtot_90s = fit_formula_mp(Vtot_end_save_90s',Yend_real_save_90s',max(Vtot_end_save_90s),1);
% edit legend 
hleg3 = legend([h_Vtot,(fitted_Vtot_90s.pict_fit)'],'1','2','3','4','5','6','7','8','9','10','11','12','\alpha /[1+( \gamma / M^{0.63})]',fitted_Vtot_90s.leg_name1,fitted_Vtot_90s.leg_name2);
set(hleg3,'Location','EastOutside')
% add labels
xhandle90_V = xlabel('Vtot_{end} (m^{3})');
yhandle90_V = ylabel('Y_{end}','rotation',90); 
% add title
title('\alpha = 0.3 and \gamma = 1/90s');

% add theoretical loglog to fig11
figure(11);
hold on;
Vtot_area_900s = 1e-22:max(Vtot_end_save_900s)/100:max(Vtot_end_save_900s)*1.1;
Yend_kaava_Vtot_900s = alfa./(1+(1/900)./((roo.*1e6.*Vtot_area_900s).^0.63));
h2_Vtot(count2) = loglog(Vtot_area_900s, Yend_kaava_Vtot_900s, 'r');
% fit data
fitted_Vtot_900s = fit_formula_mp(Vtot_end_save_900s',Yend_real_save_900s',max(Vtot_end_save_900s),1);
% edit legend 
hleg4 = legend([h2_Vtot,(fitted_Vtot_900s.pict_fit)'],'1','2','3','4','5','6','7','8','9','10','\alpha /[1+( \gamma / M^{0.63})]',fitted_Vtot_900s.leg_name1,fitted_Vtot_900s.leg_name2);
set(hleg4,'Location','EastOutside')
% add labels
xhandle900_V = xlabel('Vtot_{end} (m^{3})');
yhandle900_V = ylabel('Y_{end}','rotation',90); 
% add title
title('\alpha = 0.3 and \gamma = 1/900s');

%% save pictures
% saveas(h8,'Yend(CSend)_90s.fig')
% saveas(h9,'Yend(CSend)_900s.fig')
% saveas(h10,'Yend(Vtotend)_90s.fig')
% saveas(h11,'Yend(Vtotend)_900s.fig')