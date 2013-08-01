clear all;
close all;

%% constants
roo = 1.4; % Particle density (g/cm3)
M = 300; %g/mol precursor molar mass
NA = 6.022e23; % 1/mol
alfa = 0.3;
ms = 8;

%% initiate variables 
count1 = 1;
count2 = 1;
c1 = 1;
c2 = 1;
c3 = 1;
c4 = 1;
c5 = 1;
c6 = 1;
c7 = 1;
c8 = 1;
c10 = 1;
c20 = 1;
c30 = 1;
c40 = 1;
c50 = 1;
c60 = 1;
c70 = 1;
c80 = 1;

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

Vtot_3ng_50s = zeros(1,4);
Yend_3ng_50s = zeros(1,4);
Vtot_15ng_50s = zeros(1,4);
Yend_15ng_50s = zeros(1,4);
Vtot_30ng_50s = zeros(1,4);
Yend_30ng_50s = zeros(1,4);
Vtot_100ng_50s = zeros(1,4);
Yend_100ng_50s = zeros(1,4);
Vtot_1000ng_50s = zeros(1,4);
Yend_1000ng_50s = zeros(1,4);
Vtot_5000ng_50s = zeros(1,4);
Yend_5000ng_50s = zeros(1,4);
Vtot_10000ng_50s = zeros(1,4);
Yend_10000ng_50s = zeros(1,4);
Vtot_30000ng_50s = zeros(1,4);
Yend_30000ng_50s = zeros(1,4);
Vtot_3ng_500s = zeros(1,4);
Yend_3ng_500s = zeros(1,4);
Vtot_15ng_500s = zeros(1,4);
Yend_15ng_500s = zeros(1,4);
Vtot_30ng_500s = zeros(1,4);
Yend_30ng_500s = zeros(1,4);
Vtot_100ng_500s = zeros(1,4);
Yend_100ng_500s = zeros(1,4);
Vtot_1000ng_500s = zeros(1,4);
Yend_1000ng_500s = zeros(1,4);
Vtot_5000ng_500s = zeros(1,4);
Yend_5000ng_500s = zeros(1,4);
Vtot_10000ng_500s = zeros(1,4);
Yend_10000ng_500s = zeros(1,4);
Vtot_30000ng_500s = zeros(1,4);
Yend_30000ng_500s = zeros(1,4);

Ntot_3ng_50s = zeros(1,4);
Ntot_15ng_50s = zeros(1,4);
Ntot_30ng_50s = zeros(1,4);
Ntot_100ng_50s = zeros(1,4);
Ntot_1000ng_50s = zeros(1,4);
Ntot_5000ng_50s = zeros(1,4);
Ntot_10000ng_50s = zeros(1,4);
Ntot_30000ng_50s = zeros(1,4);
Ntot_3ng_500s = zeros(1,4);
Ntot_15ng_500s = zeros(1,4);
Ntot_30ng_500s = zeros(1,4);
Ntot_100ng_500s = zeros(1,4);
Ntot_1000ng_500s = zeros(1,4);
Ntot_5000ng_500s = zeros(1,4);
Ntot_10000ng_500s = zeros(1,4);
Ntot_30000ng_500s = zeros(1,4);

%% calc and plot starts
for file = 1:2
% data files and runs for calc
if file == 1
    runs = 1:32;
    load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\Batch\Fixed_diff_molemass\10nm\fixed correct lambda\fixed correct density\run_20130726T203251.mat')
elseif file == 2
    runs = 1:32;
    load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\Batch\Fixed_diff_molemass\80nm\fixed correct lambda\fixed correct density\run_20130726T200039.mat')
end
%% i is index of run
for i = runs    
%% calculate deltaMoa
Ntot = chamb(i).output_data.Ntot; % 1/cm3
Vtot = chamb(i).output_data.Vtot; % m3/cm3
tim = chamb(i).output_data.tim;
deltaVtot = Vtot - Vtot(1);
%deltaMoa2(i2) = chamb(i).output_data.Mtot(i2) - chamb(i).output_data.Mtot(1); 
Moa = roo.*1e6*Vtot; % g/cm3
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
CS = CS_tot_Vapour( chamb(i).output_data.Y, chamb(i).initials.sections, tim, 0);

if chamb(i).initials.vap_wallsink_on ~= 0
    Yend = alfa./(1+(chamb(i).initials.vap_wallsink)./CS);
else
    Yend = alfa;
end

%% calculate (Yend-Y)/Y
Yend = Yend';
error_of_Yend = (Yend-Y)./Y;

%% calculate fraction of ELVOC forming aerosol = Y and fraction of ELVOC lost to wall Mvwall/deltaP
Wall_loss = zeros(1,length(deltaP_mass));
for i4 = 1:length(deltaP_mass)    
    % if deltaP_mass = 0, Y = inf or NaN -> now Y = 0
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
    tau = 189; % 0.1306days    
elseif file == 2
    tau = 117; % 0.08056days 
end

% what is loglog symbol (mark) its different for runs and files and same
% for diff gamma but same other values

if file == 1
    if (i == 1) || (i == 2)
        mark = '.';        
    elseif (i == 3) || (i == 4)
        mark = 'o';
    elseif (i == 5) || (i == 6)
        mark = 'x';
    elseif (i == 7) || (i == 8)  
        mark = '+';
    elseif (i == 9) || (i == 10)
        mark = '*';
    elseif (i == 11) || (i == 12)
        mark = 's';
    elseif (i == 13) || (i == 14)
        mark = 'd';
    elseif (i == 15) || (i == 16)  
        mark = 'v';
    elseif (i == 17) || (i == 18)
        mark = 'k.';
    elseif (i == 19) || (i == 20)
        mark = 'ko';
    elseif (i == 21) || (i == 22)
        mark = 'kx';
    elseif (i == 23) || (i == 24)  
        mark = 'k+';
    elseif (i == 25) || (i == 26)  
        mark = 'k*';
    elseif (i == 27) || (i == 28)
        mark = 'ks';        
    elseif (i == 29) || (i == 30)
        mark = 'kd';
    elseif (i == 31) || (i == 32)
        mark = 'kv';
    end
elseif file == 2
    if (i == 1) || (i == 2)
        mark = 'r.';        
    elseif (i == 3) || (i == 4)
        mark = 'ro';
    elseif (i == 5) || (i == 6)
        mark = 'rx';
    elseif (i == 7) || (i == 8)  
        mark = 'r+';
    elseif (i == 9) || (i == 10)
        mark = 'r*';
    elseif (i == 11) || (i == 12)
        mark = 'rs';
    elseif (i == 13) || (i == 14)
        mark = 'rd';
    elseif (i == 15) || (i == 16)  
        mark = 'rv';
    elseif (i == 17) || (i == 18)
        mark = 'm.';
    elseif (i == 19) || (i == 20)
        mark = 'mo';
    elseif (i == 21) || (i == 22)
        mark = 'mx';
    elseif (i == 23) || (i == 24)  
        mark = 'm+';
    elseif (i == 25) || (i == 26)  
        mark = 'm*';
    elseif (i == 27) || (i == 28)
        mark = 'ms';        
    elseif (i == 29) || (i == 30)
        mark = 'md';
    elseif (i == 31) || (i == 32)
        mark = 'mv';
    end
end % mark

% calc Yend and Vtot if gamma = 1/500s
if mod(i,2) == 0
    mark2 = 'k+';
%i is even     
    % Yend at right time and save every value to vector
    Yend_real_900s = Y(tau);
    Yend_real_save_900s(count2) = Yend_real_900s;
    % same for CS
    CSend_900s = CS(tau);    
    CSsave_900s(count2) = CSend_900s;  
    % same for Vtot_end
    Vtot_end_900s = Vtot(tau);
    Vtot_end_save_900s(count2) = Vtot_end_900s;
    
    % different mass lines
    if file == 1
        if  (i == 2) || (i == 4) || (i == 6) || (i == 8)  
            mark_3ng_500s = 'k+';
            Vtot_3ng_500s(c1) = Vtot(tau);
            Yend_3ng_500s(c1) = Y(tau);
            Ntot_3ng_500s(c1) = Ntot(tau);
            c1 = c1 + 1;
        elseif (i == 10) || (i == 12) || (i == 14) || (i == 16)  
            mark_15ng_500s = 'kx';
            Vtot_15ng_500s(c2) = Vtot(tau);
            Yend_15ng_500s(c2) = Y(tau);
            Ntot_15ng_500s(c2) = Ntot(tau);
            c2 = c2 + 1;
        elseif (i == 18) || (i == 20) || (i == 22) || (i == 24)  
            mark_30ng_500s = 'kd';
            Vtot_30ng_500s(c3) = Vtot(tau);
            Yend_30ng_500s(c3) = Y(tau);
            Ntot_30ng_500s(c3) = Ntot(tau);
            c3 = c3 + 1;
        elseif (i == 26) || (i == 28) || (i == 30) || (i == 32)
            mark_100ng_500s = 'kv';
            Vtot_100ng_500s(c4) = Vtot(tau);
            Yend_100ng_500s(c4) = Y(tau);
            Ntot_100ng_500s(c4) = Ntot(tau);
            c4 = c4 + 1;
        end
    elseif file == 2
        if  (i == 2) || (i == 4) || (i == 6) || (i == 8)  
            mark_1000ng_500s = 'k+';
            Vtot_1000ng_500s(c5) = Vtot(tau);
            Yend_1000ng_500s(c5) = Y(tau);
            Ntot_1000ng_500s(c5) = Ntot(tau);
            c5 = c5 + 1;
        elseif (i == 10) || (i == 12) || (i == 14) || (i == 16)  
            mark_5000ng_500s = 'kx';
            Vtot_5000ng_500s(c6) = Vtot(tau);
            Yend_5000ng_500s(c6) = Y(tau);
            Ntot_5000ng_500s(c6) = Ntot(tau);
            c6 = c6 + 1;
        elseif (i == 18) || (i == 20) || (i == 22) || (i == 24)  
            mark_10000ng_500s = 'kd';
            Vtot_10000ng_500s(c7) = Vtot(tau);
            Yend_10000ng_500s(c7) = Y(tau);
            Ntot_10000ng_500s(c7) = Ntot(tau);
            c7 = c7 + 1;
        elseif (i == 26) || (i == 28) || (i == 30) || (i == 32)
            mark_30000ng_500s = 'kv';
            Vtot_30000ng_500s(c8) = Vtot(tau);
            Yend_30000ng_500s(c8) = Y(tau);
            Ntot_30000ng_500s(c8) = Ntot(tau);
            c8 = c8 + 1;
        end
    end
    
    
    % loglog Yend(CSend) (one point each time)
    h9 = figure(9);     
    h2_CS(count2) = semilogx(CSend_900s, Yend_real_900s, mark2,'MarkerSize',ms);       
    hold on;
%     % loglog Yend(Vtot_end)(one point each time)
%     h11 = figure(11);
%     h2_Vtot(count2) = semilogx(Vtot_end_900s, Yend_real_900s, mark);      
%     hold on;
%     
%     % plot CS
%     h13=figure(13);
%     loglog(tim,CS,mark);
%     handle1 = xlabel('time (s)');
%     set(handle1,'Fontsize',9,'Fontname','Computermodern')
%     handle2 = ylabel('CS (s^{-1})','rotation',90);
%     set(handle2,'Fontsize',9,'Fontname','Computermodern')
%     hold on;
%     % add title
%     %title('\alpha = 0.3 and \gamma = 1/500s');
%     
%     % plot Moa
%     h15=figure(15);
%     loglog(tim,Moa.*1e12,mark);
%     handle1 = xlabel('time (s)');
%     set(handle1,'Fontsize',9,'Fontname','Computermodern')
%     handle2 = ylabel('M_{OA} (\mugm^{-3})','rotation',90);
%     set(handle2,'Fontsize',9,'Fontname','Computermodern')
%     hold on;
%     % add title
%     %title('\alpha = 0.3 and \gamma = 1/500s');

    count2 = count2 + 1;
    
% calc Yend and Vtot if gamma = 1/50s
else 
    mark1 = 'kx';
  %i is odd  
    
    % Yend at right time and save every value to vector
    Yend_real_90s = Y(tau);
    Yend_real_save_90s(count1) = Yend_real_90s;
    % same for CSend
    CSend_90s = CS(tau);  
    CSsave_90s(count1) = CSend_90s;     
    % same for Vtot_end
    Vtot_end_90s = Vtot(tau);
    Vtot_end_save_90s(count1) = Vtot_end_90s;
    
    % different mass lines
    if file == 1
        if  (i == 1) || (i == 3) || (i == 5) || (i == 7)  
            mark_3ng_50s = 'k.';
            Vtot_3ng_50s(c10) = Vtot(tau);
            Yend_3ng_50s(c10) = Y(tau);
            Ntot_3ng_50s(c10) = Ntot(tau);
            c10 = c10 + 1;
        elseif (i == 9) || (i == 11) || (i == 13) || (i == 15)  
            mark_15ng_50s = 'k*';
            Vtot_15ng_50s(c20) = Vtot(tau);
            Yend_15ng_50s(c20) = Y(tau);
            Ntot_15ng_50s(c20) = Ntot(tau);
            c20 = c20 + 1;
        elseif (i == 17) || (i == 19) || (i == 21) || (i == 23)  
            mark_30ng_50s = 'ko';
            Vtot_30ng_50s(c30) = Vtot(tau);
            Yend_30ng_50s(c30) = Y(tau);
            Ntot_30ng_50s(c30) = Ntot(tau);
            c30 = c30 + 1;
        elseif (i == 25) || (i == 27) || (i == 29) || (i == 31)
            mark_100ng_50s = 'ks';
            Vtot_100ng_50s(c40) = Vtot(tau);
            Yend_100ng_50s(c40) = Y(tau);
            Ntot_100ng_50s(c40) = Ntot(tau);
            c40 = c40 + 1;
        end
    elseif file == 2
        if  (i == 1) || (i == 3) || (i == 5) || (i == 7)   
            mark_1000ng_50s = 'k.';
            Vtot_1000ng_50s(c50) = Vtot(tau);
            Yend_1000ng_50s(c50) = Y(tau);
            Ntot_1000ng_50s(c50) = Ntot(tau);
            c50 = c50 + 1;
        elseif (i == 9) || (i == 11) || (i == 13) || (i == 15)   
            mark_5000ng_50s = 'k*';
            Vtot_5000ng_50s(c60) = Vtot(tau);
            Yend_5000ng_50s(c60) = Y(tau);
            Ntot_5000ng_50s(c60) = Ntot(tau);
            c60 = c60 + 1;
        elseif (i == 17) || (i == 19) || (i == 21) || (i == 23)   
            mark_10000ng_50s = 'ko';
            Vtot_10000ng_50s(c70) = Vtot(tau);
            Yend_10000ng_50s(c70) = Y(tau);
            Ntot_10000ng_50s(c70) = Ntot(tau);
            c70 = c70 + 1;
        elseif (i == 25) || (i == 27) || (i == 29) || (i == 31)
            mark_30000ng_50s = 'ks';
            Vtot_30000ng_50s(c80) = Vtot(tau);
            Yend_30000ng_50s(c80) = Y(tau);
            Ntot_30000ng_50s(c80) = Ntot(tau);
            c80 = c80 + 1;
        end
    end
    
    % loglog Yend(CSend)(one point each time)
    h8 = figure(9);
    h_CS(count1) = semilogx(CSend_90s, Yend_real_90s, mark1,'MarkerSize',ms);      
    hold on;
%     % loglog Yend(Vtot_end)(one point each time)
%     h10 = figure(11);
%     h_Vtot(count1) = semilogx(Vtot_end_90s, Yend_real_90s, mark);      
%     hold on;
    
%     % plot CS
%     h12=figure(12);
%     loglog(tim,CS,mark);
%     handle1 = xlabel('time (s)');
%     set(handle1,'Fontsize',9,'Fontname','Computermodern')
%     handle2 = ylabel('CS (s^{-1})','rotation',90);
%     set(handle2,'Fontsize',9,'Fontname','Computermodern')
%     hold on;
%     % add title
%     %title('\alpha = 0.3 and \gamma = 1/50s');
%     
%     % plot Moa
%     h14=figure(14);
%     loglog(tim,Moa.*1e12,mark);
%     handle1 = xlabel('time (s)');
%     set(handle1,'Fontsize',9,'Fontname','Computermodern')
%     handle2 = ylabel('M_{OA} (\mugm^{-3})','rotation',90);
%     set(handle2,'Fontsize',9,'Fontname','Computermodern')
%     hold on;
%     % add title
%     %title('\alpha = 0.3 and \gamma = 1/50s');
    
    count1 = count1 + 1;     
end

%% plot final CS and Moa examples
lw = 1.5;
if file == 1
    if (i == 3)        
        % plot CS
        h30=figure(30);
        CS3 = loglog(tim,CS,'k-','LineWidth',lw);
        hold on;  
        % plot Moa
        h31=figure(31);
        M3 = loglog(tim,Moa.*1e12,'k-','LineWidth',lw);
        hold on;    
    elseif (i == 4)
        % plot CS
        h30=figure(30);
        CS4 = loglog(tim,CS,'k--','LineWidth',lw);        
        hold on;  
        % plot Moa
        h31=figure(31);
        M4 = loglog(tim,Moa.*1e12,'k--','LineWidth',lw);        
        hold on;
    end
elseif file == 2
    if (i == 23) 
        % plot CS
        h30=figure(30);
        CS31 = loglog(tim,CS,'k-.','LineWidth',lw);
        hold on;  
        % plot Moa
        h31=figure(31);
        M31 = loglog(tim,Moa.*1e12,'k-.','LineWidth',lw);
        hold on;    
    elseif (i == 24)
        % plot CS
        h30=figure(30);
        CS32 = loglog(tim,CS,'k:','LineWidth',lw);
        handl1 = xlabel('time (s)');
        handl2 = ylabel('CS (s^{-1})','rotation',90);
        hold on;  
        axis([5e1 4e4 8e-5 1.5e-1])
        matlab2tikz('CS_time.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );
        % plot Moa
        h31=figure(31);
        M32 = loglog(tim,Moa.*1e12,'k:','LineWidth',lw);
        handle1 = xlabel('time (s)');
        handle2 = ylabel('M_{OA} (\mugm^{-3})','rotation',90);    
        hold on;
        axis([5e1 4e4 4e-3 1e3]) 
        matlab2tikz('Moa_time.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );
    end
end 
% leg_name_CS = [CS3 CS4 CS31 CS32];
% leg_name_M = [M3 M4 M31 M32];
% figure(30);
% legend(leg_name_CS,'
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
Loc = 'SouthEast';
Orient = 'vertical';
% add theoretical loglog to fig8
figure(9);
hold on;
CSarea_90s = min(CSsave_90s)*0.9:max(CSsave_90s)/100:max(CSsave_90s)*1.1;
Yend_kaava_90s = alfa./(1+(1/50)./CSarea_90s);
h_CS(count1) = loglog(CSarea_90s, Yend_kaava_90s, 'r','LineWidth',lw);
% fit data
fitted_90s = fit_formula_mp(CSsave_90s',Yend_real_save_90s',0,0);
% edit legend 
%hleg1 = legend([h_CS(end),(fitted_90s.pict_fit)'],'\alpha /[1+(\gamma /CS_{end})]',fitted_90s.leg_name1,fitted_90s.leg_name2);
%         '10nm,3ng/m^3,10ppb',      '                     50ppb','                     100ppb','                     200ppb',...
%         '       15ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '       30ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '      100ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '80nm,1\mug/m^3,10ppb',    '                     50ppb','                     100ppb','                     200ppb',... 
%         '         5\mug/m^3,10ppb','                     50ppb','                     100ppb','                     200ppb',...
%         '       10\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
%         '       30\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
        
%set(hleg1,'Location',Loc,'Orientation',Orient)
% add labels
xhandle90 = xlabel('CS_{end} (s^{-1})');
yhandle90 = ylabel('Y_{end}','rotation',90); 
% add title
%title('\alpha = 0.3 and \gamma = 1/50s');

% add theoretical loglog to fig9
figure(9);
hold on;
CSarea_900s = min(CSsave_900s)*0.9:max(CSsave_900s)/100:max(CSsave_900s)*1.1;
Yend_kaava_900s = alfa./(1+(1/500)./CSarea_900s);
h2_CS(count2) = loglog(CSarea_900s, Yend_kaava_900s, 'r--','LineWidth',lw);
% fit data
fitted_900s = fit_formula_mp(CSsave_900s',Yend_real_save_900s',0,1);
% edit legend 
hleg2 = legend([h2_CS(1),h2_CS(end),(fitted_900s.pict_fit)',h_CS(1),h_CS(end),(fitted_90s.pict_fit)'],'data \gamma = 0.002','fit \alpha = 0.3     \gamma = 0.002',fitted_900s.leg_name1,fitted_900s.leg_name2,'data \gamma = 0.02','fit \alpha = 0.3     \gamma = 0.02',fitted_90s.leg_name1,fitted_90s.leg_name2);
%         '10nm,3ng/m^3,10ppb',      '                     50ppb','                     100ppb','                     200ppb',...
%         '       15ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '       30ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '      100ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '80nm,1\mug/m^3,10ppb',    '                     50ppb','                     100ppb','                     200ppb',... 
%         '         5\mug/m^3,10ppb','                     50ppb','                     100ppb','                     200ppb',...
%         '       10\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
%         '       30\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
        
%hleg2 = legend([fitted_90s.pict_fit(1),fitted_900s.pict_fit(1)],'\gamma = 1/50s','\gamma = 1/500s');
set(hleg2,'Location',Loc,'Fontsize',8)
%legend(hleg2,'hide')
legend(hleg2,'boxoff')
% add labels
xhandle900 = xlabel('CS_{end} (s^{-1})');
yhandle900 = ylabel('Y_{end}','rotation',90); 
%set(yhandle900,'ylim',[1 4])
% add title
%title('\alpha = 0.3 and \gamma = 1/500s');
%fix axis
axis([1e-3 5e-1 0 0.30000001]);
set(gca,'YTick',[0,0.1,0.2,0.3])

%TextBox('\gamma = 1/50s',[20 20 5 0],figure(9))
matlab2tikz('Yend_CSend.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

% % add theoretical loglog to fig10
% figure(11);
% hold on;
% Vtot_area_90s = min(Vtot_end_save_90s)*0.9:max(Vtot_end_save_90s)/100:max(Vtot_end_save_90s)*1.1;
% Yend_kaava_Vtot_90s = alfa./(1+(1/50)./(2e-4.*(1.0e4.^0.37).*(1e6.*1e6.*roo.*1e6.*Vtot_area_90s).^0.63)); % mass in kg
% h_Vtot(count1) = loglog(Vtot_area_90s, Yend_kaava_Vtot_90s, 'r');
% % fit data
% fitted_Vtot_90s = fit_formula_mp(Vtot_end_save_90s',Yend_real_save_90s',1,1);
% % edit legend 
% hleg3 = legend([h_Vtot,(fitted_Vtot_90s.pict_fit)'],...
%         '10nm,3ng/m^3,10ppb',      '                     50ppb','                     100ppb','                     200ppb',...
%         '       15ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '       30ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '      100ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '80nm,1\mug/m^3,10ppb',    '                     50ppb','                     100ppb','                     200ppb',... 
%         '         5\mug/m^3,10ppb','                     50ppb','                     100ppb','                     200ppb',...
%         '       10\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
%         '       30\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
%         '\alpha /[1+( \gamma / M^{0.63})]',fitted_Vtot_90s.leg_name1,fitted_Vtot_90s.leg_name2);
% set(hleg3,'Location',Loc)
% % add labels
% xhandle90_V = xlabel('Vtot_{end} (m^{3})');
% yhandle90_V = ylabel('Y_{end}','rotation',90); 
% % add title
% %title('\alpha = 0.3 and \gamma = 1/50s');
% 
% % add theoretical loglog to fig11
% figure(11);
% hold on;
% Vtot_area_900s = min(Vtot_end_save_900s)*0.9:max(Vtot_end_save_900s)/100:max(Vtot_end_save_900s)*1.1;
% Yend_kaava_Vtot_900s = alfa./(1+(1/500)./(2e-4.*(1.0e4.^0.37).*(1e6.*1e6.*roo.*1e6.*Vtot_area_900s).^0.63)); % mass in µg/m3
% h2_Vtot(count2) = loglog(Vtot_area_900s, Yend_kaava_Vtot_900s, 'r');
% % fit data
% fitted_Vtot_900s = fit_formula_mp(Vtot_end_save_900s',Yend_real_save_900s',1,0);
% % edit legend 
% hleg4 = legend([h2_Vtot,(fitted_Vtot_900s.pict_fit)'],...
%         '10nm,3ng/m^3,10ppb',      '                     50ppb','                     100ppb','                     200ppb',...
%         '       15ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '       30ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '      100ng/m^3,10ppb',   '                     50ppb','                     100ppb','                     200ppb',...
%         '80nm,1\mug/m^3,10ppb',    '                     50ppb','                     100ppb','                     200ppb',... 
%         '         5\mug/m^3,10ppb','                     50ppb','                     100ppb','                     200ppb',...
%         '       10\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
%         '       30\mug/m^3,10ppb', '                     50ppb','                     100ppb','                     200ppb',...
%         '\alpha /[1+( \gamma / M^{0.63})]',fitted_Vtot_900s.leg_name1,fitted_Vtot_900s.leg_name2);
% set(hleg4,'Location',Loc)
% % add labels
% xhandle900_V = xlabel('Vtot_{end} (m^{3})');
% yhandle900_V = ylabel('Y_{end}','rotation',90); 
% % add title
% %title('\alpha = 0.3 and \gamma = 1/500s');
% %fix axis
% %axis([1e-18 1e-15 1e-1 1])

%% mean values for Nend

mean1 = mean([Ntot_3ng_50s Ntot_1000ng_50s]);
mean2 = mean([Ntot_15ng_50s Ntot_5000ng_50s]);
mean3 = mean([Ntot_30ng_50s Ntot_10000ng_50s]);
mean4 = mean([Ntot_100ng_50s Ntot_30000ng_50s]);

mean5 = mean([Ntot_3ng_500s Ntot_1000ng_500s]);
mean6 = mean([Ntot_15ng_500s Ntot_5000ng_500s]);
mean7 = mean([Ntot_30ng_500s Ntot_10000ng_500s]);
mean8 = mean([Ntot_100ng_500s Ntot_30000ng_500s]);

double_mean1 = mean([mean1 mean5]);
double_mean2 = mean([mean2 mean6]);
double_mean3 = mean([mean3 mean7]);
double_mean4 = mean([mean4 mean8]);

%% semilogx different N, Vtot

figure(32);

fig32 = semilogx(Vtot_1000ng_50s, Yend_1000ng_50s, mark_1000ng_50s, ...
                Vtot_5000ng_50s, Yend_5000ng_50s, mark_5000ng_50s, ...
                Vtot_10000ng_50s, Yend_10000ng_50s, mark_10000ng_50s, ...
                Vtot_30000ng_50s, Yend_30000ng_50s, mark_30000ng_50s, ...
                Vtot_1000ng_500s, Yend_1000ng_500s, mark_1000ng_500s, ...
                Vtot_5000ng_500s, Yend_5000ng_500s, mark_5000ng_500s, ...
                Vtot_10000ng_500s, Yend_10000ng_500s, mark_10000ng_500s, ...
                Vtot_30000ng_500s, Yend_30000ng_500s, mark_30000ng_500s,'MarkerSize',ms);
            
% fig32 = semilogx(Vtot_3ng_50s, Yend_3ng_50s, mark_3ng_50s, ... 
%                 Vtot_15ng_50s, Yend_15ng_50s, mark_15ng_50s, ...
%                 Vtot_30ng_50s, Yend_30ng_50s, mark_30ng_50s, ...
%                 Vtot_100ng_50s, Yend_100ng_50s, mark_100ng_50s, ...
%                 Vtot_1000ng_50s, Yend_1000ng_50s, mark_1000ng_50s, ...
%                 Vtot_5000ng_50s, Yend_5000ng_50s, mark_5000ng_50s, ...
%                 Vtot_10000ng_50s, Yend_10000ng_50s, mark_10000ng_50s, ...
%                 Vtot_30000ng_50s, Yend_30000ng_50s, mark_30000ng_50s, ...
%                 Vtot_3ng_500s, Yend_3ng_500s, mark_3ng_500s, ...
%                 Vtot_15ng_500s, Yend_15ng_500s, mark_15ng_500s, ...
%                 Vtot_30ng_500s, Yend_30ng_500s, mark_30ng_500s, ...
%                 Vtot_100ng_500s, Yend_100ng_500s, mark_100ng_500s, ...
%                 Vtot_1000ng_500s, Yend_1000ng_500s, mark_1000ng_500s, ...
%                 Vtot_5000ng_500s, Yend_5000ng_500s, mark_5000ng_500s, ...
%                 Vtot_10000ng_500s, Yend_10000ng_500s, mark_10000ng_500s, ...
%                 Vtot_30000ng_500s, Yend_30000ng_500s, mark_30000ng_500s);
            
% fits and legends

% fitted to be vectors
Vtot1_50s = Vtot_1000ng_50s;
Vtot2_50s = Vtot_5000ng_50s;
Vtot3_50s = Vtot_10000ng_50s;
Vtot4_50s = Vtot_30000ng_50s;

Yend1_50s = Yend_1000ng_50s;
Yend2_50s = Yend_5000ng_50s;
Yend3_50s = Yend_10000ng_50s;
Yend4_50s = Yend_30000ng_50s;

Vtot1_500s = Vtot_1000ng_500s;
Vtot2_500s = Vtot_5000ng_500s;
Vtot3_500s = Vtot_10000ng_500s;
Vtot4_500s = Vtot_30000ng_500s;

Yend1_500s = Yend_1000ng_500s;
Yend2_500s = Yend_5000ng_500s;
Yend3_500s = Yend_10000ng_500s;
Yend4_500s = Yend_30000ng_500s;

% % fitted to be vectors
% Vtot1_50s = [Vtot_3ng_50s Vtot_1000ng_50s];
% Vtot2_50s = [Vtot_15ng_50s Vtot_5000ng_50s];
% Vtot3_50s = [Vtot_30ng_50s Vtot_10000ng_50s];
% Vtot4_50s = [Vtot_100ng_50s Vtot_30000ng_50s];
% 
% Yend1_50s = [Yend_3ng_50s Yend_1000ng_50s];
% Yend2_50s = [Yend_15ng_50s Yend_5000ng_50s];
% Yend3_50s = [Yend_30ng_50s Yend_10000ng_50s];
% Yend4_50s = [Yend_100ng_50s Yend_30000ng_50s];
% 
% Vtot1_500s = [Vtot_3ng_500s Vtot_1000ng_500s];
% Vtot2_500s = [Vtot_15ng_500s Vtot_5000ng_500s];
% Vtot3_500s = [Vtot_30ng_500s Vtot_10000ng_500s];
% Vtot4_500s = [Vtot_100ng_500s Vtot_30000ng_500s];
% 
% Yend1_500s = [Yend_3ng_500s Yend_1000ng_500s];
% Yend2_500s = [Yend_15ng_500s Yend_5000ng_500s];
% Yend3_500s = [Yend_30ng_500s Yend_10000ng_500s];
% Yend4_500s = [Yend_100ng_500s Yend_30000ng_500s];

% fit
hold on;
fitted_Vtot1_50s = fit_formula_mp(Vtot1_50s',Yend1_50s',1,0,mean(Ntot_1000ng_50s));
fitted_Vtot2_50s = fit_formula_mp(Vtot2_50s',Yend2_50s',1,1,mean(Ntot_5000ng_50s));
fitted_Vtot3_50s = fit_formula_mp(Vtot3_50s',Yend3_50s',1,2,mean(Ntot_10000ng_50s));
fitted_Vtot4_50s = fit_formula_mp(Vtot4_50s',Yend4_50s',1,3,mean(Ntot_30000ng_50s));

fitted_Vtot1_500s = fit_formula_mp(Vtot1_500s',Yend1_500s',1,4,mean(Ntot_1000ng_500s));
fitted_Vtot2_500s = fit_formula_mp(Vtot2_500s',Yend2_500s',1,5,mean(Ntot_5000ng_500s));
fitted_Vtot3_500s = fit_formula_mp(Vtot3_500s',Yend3_500s',1,6,mean(Ntot_10000ng_500s));
fitted_Vtot4_500s = fit_formula_mp(Vtot4_500s',Yend4_500s',1,7,mean(Ntot_30000ng_500s));

% % % theoretical fit
% % Vtot1_area_50s = min(Vtot1_50s)*0.9:max(Vtot1_50s)/100:max(Vtot1_50s)*1.1;
% % Yend_kaava_Vtot1_50s = alfa./(1+(1/500)./(2e-4.*(double_mean1.^0.37).*(1e6.*1e6.*roo.*1e6.*Vtot1_area_50s).^0.63)); % mass in µg
% % h_Vtot1 = loglog(Vtot1_area_50s, Yend_kaava_Vtot1_50s, 'r');

% legends and axes

hleg2 = legend([(fig32)',fitted_Vtot1_50s.pict_fit,fitted_Vtot2_50s.pict_fit,fitted_Vtot3_50s.pict_fit,fitted_Vtot4_50s.pict_fit,fitted_Vtot1_500s.pict_fit,fitted_Vtot2_500s.pict_fit,fitted_Vtot3_500s.pict_fit,fitted_Vtot4_500s.pict_fit],'data \gamma = 0.002, M_0 = 1','data \gamma = 0.002, M_0 = 5','data \gamma = 0.002, M_0 = 10','data \gamma = 0.002, M_0 = 30','data \gamma = 0.02, M_0 = 1','data \gamma = 0.02, M_0 = 5','data \gamma = 0.02, M_0 = 10','data \gamma = 0.02, M_0 = 30',fitted_Vtot1_50s.leg_name1,fitted_Vtot2_50s.leg_name1,fitted_Vtot3_50s.leg_name1,fitted_Vtot4_50s.leg_name1,fitted_Vtot1_500s.leg_name1,fitted_Vtot2_500s.leg_name1,fitted_Vtot3_500s.leg_name1,fitted_Vtot4_500s.leg_name1);
Loc = 'NorthOutside';
set(hleg2,'Location',Loc,'Fontsize',8)
%legend(hleg2,'hide')
legend(hleg2,'boxoff')
% add labels
xhandle9 = xlabel('V_{end} (m^3/cm^3)');
yhandle9 = ylabel('Y_{end}','rotation',90); 
%set(yhandle900,'ylim',[1 4])
% add title
%title('\alpha = 0.3 and \gamma = 1/500s');
%fix axis
axis([2e-18 5e-16 0 0.30000001])
set(gca,'YTick',[0,0.1,0.2,0.3])

%TextBox('\gamma = 1/50s',[20 20 5 0],figure(9))
matlab2tikz('Yend_Vend.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

%% save pictures
% saveas(h8,'Yend(CSend)_50s.fig')
% saveas(h9,'Yend(CSend)_500s.fig')
% saveas(h10,'Yend(Vtotend)_50s.fig')
% saveas(h11,'Yend(Vtotend)_500s.fig')
% saveas(h12,'CS(t)_50s.fig')
% saveas(h13,'CS(t)_500s.fig')
% saveas(h14,'Moa(t)_50s.fig')
% saveas(h15,'Moa(t)_500s.fig')