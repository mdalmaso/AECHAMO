ajoja = 24;

for i = [1:14,15, 16:24]
    chamb(i).plot('original')
     chamb(i).mass_conserv_check
end


chamb(1).plot('dist')

CS = zeros(2881,ajoja);
for i = 1:ajoja
    CS(1:end,i) = CS_tot(chamb(i).output_data.distr);
end

%Vtot end
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

