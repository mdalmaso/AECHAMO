C = [10,10,10,10,20];
H = [15,14,16,16,32];
O = [11,9,10,9,12];
N = [1,0,0,0,0];


M = C*12.01 + H*1.008 + O*16.00 + N*14.01;
V = C*15.9 + H*2.31 + O*6.11 + N*4.54 - 18.3;

[Dif, lambda, lambda2] = Diff_coeff(M,V);

p = polyfit(M,Dif,1);
f = polyval(p,M);

figure(1)
plot(M,Dif,'.')
hold on;
plot(M,f,'r')
xhandle = xlabel('M (g/mol)');
yhandle = ylabel('D (cm^2/s)','rotation',90); 

figure(2)
plot(M,lambda,'r.',M,lambda2,'k.')
xhandle = xlabel('M (g/mol)');
yhandle = ylabel('\lambda (m)','rotation',90); 
