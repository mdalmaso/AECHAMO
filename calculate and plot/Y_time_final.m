
figure(30)
figfig30= semilogx(tim,Y,'k-'); % run 11 80nm
hold on;
figure(30)
figfig302 =semilogx(tim,Y,'k--'); % run 12 80nm
hold on;

% leg=legend([figfig30 figfig302],'\tau = 1/\gamma = 50 s','\tau = 1/\gamma = 500 s');
% set(leg,'Location','SouthEast')
% legend(leg,'boxoff')
axis([50 40000 0.04 0.30001])
set(gca,'YTick',[0,0.1,0.2,0.3])
xlabel('time (s)');
ylabel('Yield','rotation',90);

matlab2tikz('Y_time2.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth');