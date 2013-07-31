function [out] = fit_formula_mp(x,y,mass_on,dashed)
%linewidth
lw=1.3;
% fit ap data
xmin = min(x);
xmax = max(x);
% if x = CS
if mass_on == 0
    % a vapaa
    s=fitoptions('Method','NonLinearLeastSquares','Startpoint',[0.1  0.1],'upper',[1 10],'lower',[0.003 0.0000001]);
    custom = 'a.*(1/(1+b./x))';
    f = fittype(custom,'options',s);

    % a kiinnitetty 
    a = 0.3;
    s2=fitoptions('Method','NonLinearLeastSquares','Startpoint', 0.1 ,'upper',10,'lower',0.0000001);
    custom2 = '0.3.*(1/(1+b./x))';
    f2 = fittype(custom2 ,'options',s2);
    
    log_xmin = log10(xmin*0.9);
    log_xmax = log10(xmax*1.1);
    x0 = logspace(log_xmin,log_xmax);

    r = fit(x,y,f);
    r2 = fit(x,y,f2);
    out.coeff_r = coeffvalues(r);
    out.coeff_r2 = coeffvalues(r2);

    %plot(x,y,'ks')
    hold on;
    if dashed ~= 0
        out.pict_fit = plot(x0,feval(r,x0),'k--',x0,feval(r2,x0),'b--','LineWidth',lw);
    else
        out.pict_fit = plot(x0,feval(r,x0),'k-',x0,feval(r2,x0),'b-','LineWidth',lw);
    end
    a_string = num2str(out.coeff_r(1));
    b_string = num2str(out.coeff_r(2));
    b2_string = num2str(out.coeff_r2(1));
    a2_string = num2str(a);
    out.leg_name1 =  ['fit \alpha = ' a_string  ' \gamma = ' b_string];
    out.leg_name2 =  ['fit \alpha =      ' a2_string ' \gamma = ' b2_string];
%     out.leg_name1 = ['fit ' custom ' with a = ' a_string ' ja b = ' b_string];
%     out.leg_name2 = ['fit ' custom2 ' with b = ' b2_string];

% if x = Vtot
elseif mass_on == 1
    % a vapaa
    s=fitoptions('Method','NonLinearLeastSquares','Startpoint',[0.1  100 ],'upper',[1 1e4],'lower',[0.0003 0.00000000010]);
    custom = 'a*(1/(1+b./(1e6.*1e6.*1.4.*1e6.*x)^0.63))'; % mass in µg/m3
    f = fittype(custom,'options',s);

    % a kiinnitetty a = 0.3;
    s2=fitoptions('Method','NonLinearLeastSquares','Startpoint', 100 ,'upper',1e4,'lower',0.00000000010);
    custom2 = '0.3.*(1/(1+b./(1e6.*1e6.*1.4.*1e6.*x)^0.63))'; % mass in µg/m3
    f2 = fittype(custom2 ,'options',s2);

    log_xmin = log10(xmin*0.9);
    log_xmax = log10(xmax*1.1);
    x0 = logspace(log_xmin,log_xmax);

    r = fit(x,y,f);
    r2 = fit(x,y,f2);
    out.coeff_r = coeffvalues(r);
    out.coeff_r2 = coeffvalues(r2);

    %plot(x,y,'ks')
    hold on;
    if dashed == 0
        out.pict_fit = plot(x0,feval(r,x0),'k-','LineWidth',lw);
    elseif dashed == 1
        out.pict_fit = plot(x0,feval(r,x0),'k--','LineWidth',lw);
    elseif dashed == 2
        out.pict_fit = plot(x0,feval(r,x0),'k:','LineWidth',lw);
    elseif dashed == 3
        out.pict_fit = plot(x0,feval(r,x0),'k-.','LineWidth',lw);        
    end   
    a_string = num2str(out.coeff_r(1));
    b_string = num2str(out.coeff_r(2));
    b2_string = num2str(out.coeff_r2(1));
    out.leg_name1 =  sprintf(['fit a = ' a_string ' and \n b = ' b_string]);
    out.leg_name2 =  sprintf(['fit b = ' b2_string]);
%     out.leg_name1 =  sprintf(['fit ' custom ' with \n a = ' a_string ' ja b = ' b_string]);
%     out.leg_name2 =  sprintf(['fit ' custom2 ' with \n b = ' b2_string]);

end

end
