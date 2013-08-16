function[out] = MT_kinetics_night_and_day(plotting)

% run variable values
o3_OH_MT_T_Variable_values


ppb = M_dens(1013,290).*1e-9;


odepara.tim = time.*3600;
odepara.O3 = O3.*ppb;
odepara.MT = E;
odepara.ppb = ppb;
odepara.OH = OH;
odepara.kO3 = kO3;
odepara.kOH = kOH;

% % Used by Dal Maso
% odepara.kO3 = 15e-17;
% odepara.kOH = 6e-11;

tspan = 0:60:time_max*3600;
y0 = 0.1975.*ppb;

[tim, Y] = ode45(@(t,y) kinefun(t,y,odepara),tspan,y0);

out.meas.time = time;
out.meas.E_MT = E;
out.meas.O3 = O3;
out.meas.OH = OH./ppb*1e3; 
out.meas.T = T;
out.ppb = ppb;
out.model.MT = Y(:)./ppb;
out.model.time= tim;

% calculate source of Cvap

iO3 = interp1(odepara.tim,odepara.O3,tim);
iOH = interp1(odepara.tim,odepara.OH,tim);
out.model.Q_Cvap = (kOH.*iOH + kO3.*iO3).*Y(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 1.5;
ms = 8;
if plotting == 0
    h1=figure(1);
    plot(out.model.time/3600, out.model.MT, 'k-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2,0.205,0.21], 'YTickLabel',{' ' '0.17' ' ' '0.18' ' ' '0.19' ' ' '0.2' ' ' '0.21'})
    ylabel('MT (ppb)')
    xlabel('time')
    axis([0 96 0.1649999 0.21])
    saveas(h1,'concentration_MT.jpg') 
    saveas(h1,'concentration_MT.fig') 
    matlab2tikz('concentration_MT.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

    h2=figure(2);
    plot(time, out.meas.T, 'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[9,10,11,12,13,14,15,16,17,18], 'YTickLabel',{'9' ' ' ' ' '12' ' ' ' ' '15' ' ' ' ' '18'})
    ylabel('T (^{\circ}C)')
    xlabel('time')
    axis([0 96 9 18])
    saveas(h2,'temperature.jpg') 
    saveas(h2,'temperature.fig') 
    matlab2tikz('temperature.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

    h3=figure(3);
    plot(time, out.meas.O3, 'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40], 'YTickLabel',{'25' ' ' ' ' ' ' ' ' '30' ' ' ' ' ' ' ' ' '35' ' ' ' ' ' ' ' ' '40'})
    ylabel('O_3 (ppb)')
    xlabel('time')
    axis([0 96 25 40])
    saveas(h3,'concentration_O3.jpg') 
    saveas(h3,'concentration_O3.fig') 
    matlab2tikz('concentration_O3.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

    h4=figure(4);
    semilogy(time, out.meas.OH, 'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    ylabel('OH (ppt)')
    xlabel('time')
    axis([0 96 7e-4 3.1e-2])
    saveas(h4,'concentration_OH.jpg') 
    saveas(h4,'concentration_OH.fig') 
    matlab2tikz('concentration_OH.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

    h5=figure(5);
    plot(time,out.meas.E_MT,'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[2.5e5,3e5,3.5e5,4e5,4.5e5,5e5,5.5e5], 'YTickLabel',{'2.5' ' ' '3.5' ' ' '4.5' ' ' '5.5'})
    ylabel('E (cm^{-3}s^{-1})')
    xlabel('time')
    axis([0 96 2.5e5 5.5e5])
    saveas(h5,'emission_MT.jpg') 
    saveas(h5,'emission_MT.fig') 
    matlab2tikz('emission_MT.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

    h6=figure(6);
    plot(out.model.time/3600, out.model.Q_Cvap, 'k-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[2.5e5,3e5,3.5e5,4e5,4.5e5,5e5,5.5e5], 'YTickLabel',{'2.5' ' ' '3.5' ' ' '4.5' ' ' '5.5'})
    ylabel('Q_{vap} (cm^{-3}s^{-1})')
    xlabel('time')
    axis([0 96 2.5e5 5.5e5])
    saveas(h6,'source_Cvap.jpg') 
    saveas(h6,'source_Cvap.fig') 
    matlab2tikz('source_Cvap.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );

elseif plotting == 1

    figure(1);
    subplot(5,3,5)
    plot(out.model.time/3600, out.model.MT, 'k-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2,0.205,0.21], 'YTickLabel',{' ' '0.17' ' ' '0.18' ' ' '0.19' ' ' '0.2' ' ' '0.21'})
    ylabel('MT (ppb)')
    xlabel('time')
    axis([0 96 0.1649999 0.21])

    subplot(5,3,1)
    plot(time, out.meas.T, 'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[9,10,11,12,13,14,15,16,17,18], 'YTickLabel',{'9' ' ' ' ' '12' ' ' ' ' '15' ' ' ' ' '18'})
    ylabel('T (^{\circ}C)')
    xlabel('time')
    axis([0 96 9 18])

    subplot(5,3,3)
    plot(time, out.meas.O3, 'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40], 'YTickLabel',{'25' ' ' ' ' ' ' ' ' '30' ' ' ' ' ' ' ' ' '35' ' ' ' ' ' ' ' ' '40'})
    ylabel('O_3 (ppb)')
    xlabel('time')
    axis([0 96 25 40])

    subplot(5,3,4)
    semilogy(time, out.meas.OH, 'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    ylabel('OH (ppt)')
    xlabel('time')
    axis([0 96 7e-4 3.1e-2])

    subplot(5,3,2)
    plot(time,out.meas.E_MT,'kx-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[2.5e5,3e5,3.5e5,4e5,4.5e5,5e5,5.5e5], 'YTickLabel',{'2.5' ' ' '3.5' ' ' '4.5' ' ' '5.5'})
    ylabel('E (cm^{-3}s^{-1})')
    xlabel('time')
    axis([0 96 2.5e5 5.5e5])

    subplot(5,3,6)
    plot(out.model.time/3600, out.model.Q_Cvap, 'k-','LineWidth',lw,'MarkerSize',ms)
    set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
    set(gca,'YTick',[2.5e5,3e5,3.5e5,4e5,4.5e5,5e5,5.5e5], 'YTickLabel',{'2.5' ' ' '3.5' ' ' '4.5' ' ' '5.5'})
    ylabel('Q_{vap} (cm^{-3}s^{-1})')
    xlabel('time')
    axis([0 96 2.5e5 5.5e5])

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = kinefun(t,y,odepara)

t
iO3 = interp1(odepara.tim,odepara.O3,t);
iMT = interp1(odepara.tim,odepara.MT,t);
iOH = interp1(odepara.tim,odepara.OH,t);


% inflow of MT
E = iMT; % molecules/second*cm3;

% MT + O3 reaction rate
rO3 = iO3.*odepara.kO3;

% MT + OH reaction rate
rOH = iOH.*odepara.kOH;


dydt = E - rOH.*y - rO3.*y;

%dcdt = rOH.*y + rO3.*y;

end














