ajoja = 24;

for i = [1:8]
    chamb(i).plot('original')
     chamb(i).mass_conserv_check
end


% jakaumaa eri kohdista
for ajo = 1:32
%chamb(ajo).plot('original','dist')
sect = chamb(ajo).initials.sections;
semilogx(chamb(ajo).output_data.distr_original(10,sect+3:2*sect+2),chamb(1).output_data.distr_original(10,3:sect+2),'r.-')
hold on
semilogx(chamb(ajo).output_data.distr_original(2,sect+3:2*sect+2),chamb(1).output_data.distr_original(2,3:sect+2),'.-')
hold on
semilogx(chamb(ajo).output_data.distr_original(100,sect+3:2*sect+2),chamb(1).output_data.distr_original(100,3:sect+2),'g.-')
hold on
semilogx(chamb(ajo).output_data.distr_original(600,sect+3:2*sect+2),chamb(1).output_data.distr_original(600,3:sect+2),'m.-')
hold on
end

CS = zeros(2881,ajoja);
for i = 1:ajoja
    CS(1:end,i) = CS_tot(chamb(i).output_data.distr);
end

%Ntot end
sect = [10, 20, 30, 40, 60, 90];
for i = 1:ajoja
   Ntot = chamb(i).output_data.Ntot(end);
   figure(1);
   plot(i ,Ntot, 'm*')   
   hold on;
   axis([0,9, 0, 6000])
end

%Vtot end
sect = [10, 20, 30, 40, 60, 90];
for i = 1:ajoja   
   Vtot =chamb(i).output_data.Vtot(end);
   figure(2);
   plot(i, Vtot, 'm*')
   hold on;
   axis([0,9, 0, 8e-18])
end

% particle source
h20 = figure(20);
plot(tvect/3600,part_source(:,2),'k-','LineWidth',1.5,'MarkerSize',8)
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
set(gca,'YTick',[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], 'YTickLabel',{'0' ' ' ' ' ' ' ' ' '0.5' ' ' ' ' ' ' ' ' '1'})
ylabel('J (cm^{-3}s^{-1})')
xlabel('time')
axis([0 96 -0.1 1.1])
saveas(h20,'nucleation.fig') 
matlab2tikz('nucleation.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

