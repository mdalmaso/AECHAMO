function fit
% fit ap data using 1 product model

ap = load('ap_data.txt');


y = ap(:,2)./ap(:,1);




x = ap(:,2);


s=fitoptions('Method','NonLinearLeastSquares','Startpoint',[0.5 1000.^0.37 1e-4],'upper',[0.9 10000.^0.37 1/20],'lower',[0.1 100.^ 0.37 1./(5.*3600)])

f = fittype('a*(1-1/(1+6e-5*b*x^0.63/c))','options',s)
x0 = 1:100;
r = fit(x,y,f)


plot(x,y,'ks')
hold on
plot(x0,feval(r,x0),'r-')


end
