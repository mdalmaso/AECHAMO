function [out] = fit_formula_mp(x,y,xmax)
% fit ap data using 1 product model


% a vapaa
s=fitoptions('Method','NonLinearLeastSquares','Startpoint',[0.1  0.1],'upper',[0.40 1],'lower',[0.03 0.001])
custom = 'a.*(1/(1+b./x))';
f = fittype(custom,'options',s)

% a kiinnitetty
a = 0.3;
s2=fitoptions('Method','NonLinearLeastSquares','Startpoint',[ 0.1 ],'upper',[1],'lower',[0.001])
custom2 = '0.3.*(1/(1+b./x))';
f2 = fittype(custom2 ,'options',s2)


x0 = 0:0.001:(xmax*1.1);



r = fit(x,y,f)
r2 = fit(x,y,f2)
out.coeff_r = coeffvalues(r);
out.coeff_r2 = coeffvalues(r2);

%plot(x,y,'ks')
hold on
out.pict_fit = plot(x0,feval(r,x0),'k-',x0,feval(r2,x0),'b-');
a_string = num2str(out.coeff_r(1));
b_string = num2str(out.coeff_r(2));
b2_string = num2str(out.coeff_r2(1));
out.leg_name1 = ['fit ' custom ' with a = ' a_string ' ja b = ' b_string];
out.leg_name2 = ['fit ' custom2 ' with b = ' b2_string];

end
