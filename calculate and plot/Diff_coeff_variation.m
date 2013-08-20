clear all;
close all;

C = [10,10,10,10,20];
H = [15,14,16,16,32];
O = [11,9,10,9,12];
N = [1,0,0,0,0];


M = C*12.01 + H*1.008 + O*16.00 + N*14.01;
V = C*15.9 + H*2.31 + O*6.11 + N*4.54 - 18.3;

[Dif, lambda, lambda2] = Diff_coeff(M,V);

p = polyfit(M,Dif,1);
x = min(M)*0.985:max(M)/100:max(M)*1.015;
f = polyval(p,x);

figure(1)
plot(M,Dif,'k.')
hold on;
f1 = plot(x,f,'k','LineWidth',1.5);
xhandle = xlabel('M_{vap} (g/mol)');
yhandle = ylabel('D_{AB} (cm^2/s)','rotation',90); 
legend(f1,'linear fit -0.0001 x + 0.0736')
legend('boxoff')
axis([250 500 0.035 0.052])
set(gca,'XTick',[250,300,350,400,450,500])
set(gca,'YTick',[0.036, 0.04, 0.044, 0.048, 0.052])

matlab2tikz('DAB.tikz','checkForUpdates',false,'showInfo', false, 'height', '\fheight', 'width', '\fwidth' );
% 
% figure(2)
% plot(M,lambda,'r.',M,lambda2,'k.')
% xhandle = xlabel('M_{vap} (g/mol)');
% yhandle = ylabel('\lambda (m)','rotation',90); 
