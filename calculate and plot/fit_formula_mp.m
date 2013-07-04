function [out] = fit_formula_mp(x,y,xmax,mass_on)
% fit ap data

% if x = CS
if mass_on == 0
    % a vapaa
    s=fitoptions('Method','NonLinearLeastSquares','Startpoint',[0.1  0.1],'upper',[0.40 1],'lower',[0.03 0.001]);
    custom = 'a.*(1/(1+b./x))';
    f = fittype(custom,'options',s);

    % a kiinnitetty a = 0.3;
    s2=fitoptions('Method','NonLinearLeastSquares','Startpoint', 0.1 ,'upper',1,'lower',0.001);
    custom2 = '0.3.*(1/(1+b./x))';
    f2 = fittype(custom2 ,'options',s2);

    log_xmax = log10(xmax*1.1);
    x0 = logspace(-5,log_xmax);

    r = fit(x,y,f);
    r2 = fit(x,y,f2);
    out.coeff_r = coeffvalues(r);
    out.coeff_r2 = coeffvalues(r2);

    %plot(x,y,'ks')
    hold on;
    out.pict_fit = plot(x0,feval(r,x0),'k-',x0,feval(r2,x0),'b-');
    a_string = num2str(out.coeff_r(1));
    b_string = num2str(out.coeff_r(2));
    b2_string = num2str(out.coeff_r2(1));
    out.leg_name1 = ['fit ' custom ' with a = ' a_string ' ja b = ' b_string];
    out.leg_name2 = ['fit ' custom2 ' with b = ' b2_string];

% if x = Vtot
elseif mass_on == 1
    % a vapaa
    s=fitoptions('Method','NonLinearLeastSquares','Startpoint',[0.1  100 ],'upper',[0.40 1e3],'lower',[0.03 10]);
    custom = 'a*(1/(1+b./(1e6.*1.84.*1e6.*x)^0.63))';
    f = fittype(custom,'options',s);

    % a kiinnitetty a = 0.3;
    s2=fitoptions('Method','NonLinearLeastSquares','Startpoint', 100 ,'upper',1e3,'lower',10);
    custom2 = '0.3.*(1/(1+b./(1e6.*1.84.*1e6.*x)^0.63))';
    f2 = fittype(custom2 ,'options',s2);

    log_xmax = log10(xmax*1.1);
    x0 = logspace(-22,log_xmax);

    r = fit(x,y,f);
    r2 = fit(x,y,f2);
    out.coeff_r = coeffvalues(r);
    out.coeff_r2 = coeffvalues(r2);

    %plot(x,y,'ks')
    hold on;
    out.pict_fit = plot(x0,feval(r,x0),'k-',x0,feval(r2,x0),'b-');
    a_string = num2str(out.coeff_r(1));
    b_string = num2str(out.coeff_r(2));
    b2_string = num2str(out.coeff_r2(1));
    out.leg_name1 = ['fit ' custom ' with a = ' a_string ' ja b = ' b_string];
    out.leg_name2 = ['fit ' custom2 ' with b = ' b2_string];

end

end
