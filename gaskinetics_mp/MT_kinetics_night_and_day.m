function[out] = MT_kinetics_night_and_day()
% test new theory
close all
clear all

% run variable values
o3_OH_MT_T_Variable_values


ppb = M_dens(1013,290).*1e-9;


odepara.tim = time.*3600;
odepara.O3 = O3.*ppb;
odepara.MT = E;
odepara.ppb = ppb;
odepara.OH = OH;
odepara.kO3 = kO3;
odepara.kOH = kOH;

% % Used by Dal Maso
% odepara.kO3 = 15e-17;
% odepara.kOH = 6e-11;

tspan = 0:60:time_max*3600;
y0 = 0.1975.*ppb;

[tim, Y] = ode45(@(t,y) kinefun(t,y,odepara),tspan,y0);

out.meas.time = time;
out.meas.E_MT = E;
out.meas.O3 = O3;
out.meas.OH = OH./ppb*1e3; 
out.meas.T = T;
out.ppb = ppb;
out.model.MT = Y(:)./ppb;
out.model.time= tim;

% calculate source of Cvap

iO3 = interp1(odepara.tim,odepara.O3,tim);
iOH = interp1(odepara.tim,odepara.OH,tim);
out.model.Q_Cvap = (kOH.*iOH + kO3.*iO3).*Y(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotting = 1;

if plotting ~= 0
h1=figure(1);
plot(out.model.time/3600, out.model.MT, 'k.-')
set(gca,'XTick',[12,36,60,84], 'XTickLabel',{'12' '12' '12' '12' '12'})
ylabel('MT (ppb)')
saveas(h1,'concentration_MT.jpg') 
saveas(h1,'concentration_MT.fig') 

h2=figure(2);
plot(time, out.meas.T, 'b.-')
set(gca,'XTick',[12,36,60,84], 'XTickLabel',{'12' '12' '12' '12' '12'})
ylabel('T (\circC)')
saveas(h2,'temperature.jpg') 
saveas(h2,'temperature.fig') 

h3=figure(3);
plot(time, out.meas.O3, 'k.-')
set(gca,'XTick',[12,36,60,84], 'XTickLabel',{'12' '12' '12' '12' '12'})
ylabel('O_3 (ppb)')
saveas(h3,'concentration_O3.jpg') 
saveas(h3,'concentration_O3.fig') 

h4=figure(4);
semilogy(time, out.meas.OH, 'r.-')
set(gca,'XTick',[12,36,60,84], 'XTickLabel',{'12' '12' '12' '12' '12'})
ylabel('OH (ppt)')
saveas(h4,'concentration_OH.jpg') 
saveas(h4,'concentration_OH.fig') 

h5=figure(5);
plot(time,out.meas.E_MT,'b.-')
set(gca,'XTick',[12,36,60,84], 'XTickLabel',{'12' '12' '12' '12' '12'})
ylabel('E_{MT} (cm^{-3}s^{-1})')
saveas(h5,'emission_MT.jpg') 
saveas(h5,'emission_MT.fig') 

h6=figure(6);
plot(out.model.time/3600, out.model.Q_Cvap, 'b.-')
set(gca,'XTick',[12,36,60,84], 'XTickLabel',{'12' '12' '12' '12' '12'})
ylabel('Q_{C_{vap}} (cm^{-3}s^{-1})')
saveas(h6,'concentration_Cvap.jpg') 
saveas(h6,'concentration_Cvap.fig') 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = kinefun(t,y,odepara)

t
iO3 = interp1(odepara.tim,odepara.O3,t);
iMT = interp1(odepara.tim,odepara.MT,t);
iOH = interp1(odepara.tim,odepara.OH,t);


% inflow of MT
E = iMT; % molecules/second*cm3;

% MT + O3 reaction rate
rO3 = iO3.*odepara.kO3;

% MT + OH reaction rate
rOH = iOH.*odepara.kOH;


dydt = E - rOH.*y - rO3.*y;

%dcdt = rOH.*y + rO3.*y;

end














