
figure(30)
figfig30= semilogx(tim,Y,'k-');
hold on;
figure(30)
figfig302 =semilogx(tim,Y,'k--');
hold on;

leg=legend([figfig30 figfig302],'\tau = 1/\gamma = 50 s','\tau = 1/\gamma = 500 s');
set(leg,'Location','SouthEast')
legend(leg,'boxoff')
axis([50 40000 0.05 0.3])
xlabel('time (s)');
ylabel('Yield','rotation',90);

matlab2tikz('Y_time.tikz','checkForUpdates',false,'showInfo', false);