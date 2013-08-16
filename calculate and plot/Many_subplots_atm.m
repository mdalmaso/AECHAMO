% Atm modelling, subplot of many parametres variation suring 4 days
close all
clearvars -except chamb

figure(1)
% Plot the distribution:
subplot(5,3,9)
chamb.plot('dist')

MT_kinetics_night_and_day(1)

load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\ATM modelling\v4\run_20130814T114856.mat')

lw = 1.5;
ms = 8;

figure(1)
% Plot the total particle volume in aerosol:
subplot(5,3,7)
plot(chamb.output_data.tim/(60*60),chamb.output_data.Vtot,'k-','LineWidth',lw,'MarkerSize',ms)
xhandle = xlabel('time');
yhandle = ylabel('V_{tot}(m^3/cm^{3})','rotation',90);
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
%set(gca,'YTick',[0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2,0.205,0.21], 'YTickLabel',{' ' '0.17' ' ' '0.18' ' ' '0.19' ' ' '0.2' ' ' '0.21'})
axis([0 96 0 2.5e-18])

% Plot the particle concentration:
subplot(5,3,8)
plot(chamb.output_data.tim/(60*60),chamb.output_data.Ntot,'k-','LineWidth',lw,'MarkerSize',ms)
xhandle = xlabel('time');
yhandle = ylabel('N_{tot}(cm^{-3})','rotation',90);
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
axis([0 96 0 10000])

% Plot the vapor concentration:
subplot(5,3,10)
plot(chamb.output_data.tim/(60*60),chamb.output_data.vap,'k-','LineWidth',lw,'MarkerSize',ms)
xhandle = xlabel('time');
yhandle = ylabel('C_{vap}(cm^{-3})','rotation',90);
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
axis([0 96 0 6e7])

figure(1)
% particle source
subplot(5,3,11)
tvect = 0:60:4*24*3600-1*3600;
part_source = create_part_source(tvect,1.0, 1*3*3600, 39600, 3e-9);
plot(tvect/3600,part_source(:,2),'k-','LineWidth',1.5,'MarkerSize',8)
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
set(gca,'YTick',[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], 'YTickLabel',{'0' ' ' ' ' ' ' ' ' '0.5' ' ' ' ' ' ' ' ' '1'})
ylabel('J (cm^{-3}s^{-1})')
xlabel('time')
axis([0 96 -0.1 1.1])

% copy CS, yield(M) and Y(t)
h1 = openfig('CS(t).fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
h2 = openfig('Y(Moa).fig','reuse');
ax2 = gca;
h3 = openfig('YandMvwall.fig','reuse');
ax3 = gca;
h4 = openfig('Moa.fig','reuse');
ax4 = gca;
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots

figure(1)
s1 = subplot(5,3,12); %create and get handle to the subplot axes
s2 = subplot(5,3,13);
s3 = subplot(5,3,14);
s4 = subplot(5,3,15);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
fig3 = get(ax3,'children');
fig4 = get(ax4,'children');

copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
copyobj(fig3,s3);
copyobj(fig4,s4);

subplot(5,3,15)
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
%set(gca,'YTick',[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], 'YTickLabel',{'0' ' ' ' ' ' ' ' ' '0.5' ' ' ' ' ' ' ' ' '1'})
ylabel('M_{OA} (\mug/m^3)')
xlabel('time')
axis([0 96 0 3])

subplot(5,3,12)
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
%set(gca,'YTick',[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], 'YTickLabel',{'0' ' ' ' ' ' ' ' ' '0.5' ' ' ' ' ' ' ' ' '1'})
ylabel('CS (s^{-1})')
xlabel('time')
axis([0 96 0 4e-3])

subplot(5,3,13)
%set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
set(gca,'YTick',[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1], 'YTickLabel',{'0' ' ' ' ' ' ' ' ' '0.05' ' ' ' ' ' ' ' ' '0.1'})
ylabel('Yield')
xlabel('M_{OA} (\mug/m^3)')
axis([0.4 3 0 0.1])

subplot(5,3,14)
set(gca,'XTick',[12,24,36,48,60,72,84], 'XTickLabel',{'12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00' ' ' '12:00'})
set(gca,'YTick',[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1], 'YTickLabel',{'0' ' ' ' ' ' ' ' ' '0.05' ' ' ' ' ' ' ' ' '0.1'})
ylabel('Yield')
xlabel('time')
axis([0 96 0 0.1])


%matlab2tikz('conclusion.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );



