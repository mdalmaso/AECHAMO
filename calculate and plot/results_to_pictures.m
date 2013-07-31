ajoja = 24;

for i = [1:32]
    chamb(i).plot('original')
     chamb(i).mass_conserv_check
end


% jakaumaa eri kohdista
for ajo = 1:32
%chamb(ajo).plot('original','dist')
sect = chamb(ajo).initials.sections;
semilogx(chamb(ajo).output_data.distr_original(10,sect+3:2*sect+2),chamb(1).output_data.distr_original(10,3:sect+2),'r.-')
hold on
semilogx(chamb(ajo).output_data.distr_original(2,sect+3:2*sect+2),chamb(1).output_data.distr_original(2,3:sect+2),'.-')
hold on
semilogx(chamb(ajo).output_data.distr_original(100,sect+3:2*sect+2),chamb(1).output_data.distr_original(100,3:sect+2),'g.-')
hold on
semilogx(chamb(ajo).output_data.distr_original(600,sect+3:2*sect+2),chamb(1).output_data.distr_original(600,3:sect+2),'m.-')
hold on
end

CS = zeros(2881,ajoja);
for i = 1:ajoja
    CS(1:end,i) = CS_tot(chamb(i).output_data.distr);
end

%Ntot end
sect = [10, 20, 30, 40, 60, 90];
for i = 1:ajoja
   Ntot = chamb(i).output_data.Ntot(end);
   figure(1);
   plot(i ,Ntot, 'm*')   
   hold on;
   axis([0,9, 0, 6000])
end

%Vtot end
sect = [10, 20, 30, 40, 60, 90];
for i = 1:ajoja   
   Vtot =chamb(i).output_data.Vtot(end);
   figure(2);
   plot(i, Vtot, 'm*')
   hold on;
   axis([0,9, 0, 8e-18])
end

