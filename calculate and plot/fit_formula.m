function fit_formula(x,y)
% fit ap data using 1 product model


% a vapaa
s=fitoptions('Method','NonLinearLeastSquares','Startpoint',[0.1  100 ],'upper',[0.40 1e3],'lower',[0.03 10])
f = fittype('a*(1/(1+b./x^0.63))','options',s)

% a kiinnitetty
a = 0.36;
s2=fitoptions('Method','NonLinearLeastSquares','Startpoint',[ 100 ],'upper',[1e3],'lower',[10])
f2 = fittype('0.36*(1/(1+b./x^0.63))','options',s)


x0 = 1:100;



r = fit(x,y,f)
r2 = fit(x,y,f2)

plot(x,y,'ks')
hold on
plot(x0,feval(r,x0),'r-')
plot(x0,feval(r2,x0),'b-')


end
