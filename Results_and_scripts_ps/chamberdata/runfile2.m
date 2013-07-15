kam=chamber;

clear tvect;
clear part_source;
clear gas_source;
clear dilu_coeff;

tvect = 0:60:24*3600;

kam.initialize('fixed_sections', 1,...
               'sedi_on', 1,...
               'coag_on', 1,...
               'dilu_on',1,...
               'vap_wallsink_on', 1,...
               'vap_wallsink', 1/300,...
               'Dp_min', -9,...
               'Dp_max', -6,...
               'sections', 25,...
               'output_sections', 250,...
               'tvect', tvect,...
               'Cvap0', 0,...
               'N', 0,...
               'mu', 37e-9,...
               'sigma', 1.6);

part_source(:,:,1) = [tvect' [zeros(866,1);5.5.*ones(28,1);zeros(length(tvect)-866-28,1)], 10e-9.*ones(length(tvect),1)];
part_source(:,:,2) = [tvect' [zeros(869,1);zeros(26,1);1.7.*ones(108,1);2.6.*ones(170,1);zeros(length(tvect)-869-26-108-170,1)], 18e-9.*ones(length(tvect),1)];



gas_source = [tvect' [interp1(time_a.*3600,Q,tvect(1:1241)),zeros(1,length(tvect)-1241)]'];

dilu_coeff = [tvect' [interp1(time_a.*3600,dilu,tvect)]'];

gas_source(find(isnan(gas_source)))=nanmean(Q);
dilu_coeff(find(isnan(dilu_coeff)))=nanmean(dilu);


kam.initialize('part_source',part_source,...
               'gas_source', gas_source,...
               'dilu_coeff', dilu_coeff);
kam.run;

kam.plot('original','dist');
ylim([Dps(1) Dps(end)]);
caxis([1 4]);
gcc=findall(gcf,'tag','Colorbar');
delete(gcc);
% 
% figure;
% subplot_dmps(dist_fixed,111);
figure('Color',[1 1 1]);
p1=plot(time,Ntot,'r--');
set(p1,'LineWidth',2);
hold on;
p1=plot(kam.output_data.tim./3600,kam.output_data.Ntot);
set(p1,'LineWidth',2);
xlim([14 24]);
xlabel('time (h)');
ylabel('Ntot (1/cm^3)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);

figure('Color',[1 1 1]);
p2=plot(time,Vtot,'r--');
set(p2,'LineWidth',2);
hold on;
p2=plot(kam.output_data.tim./3600,kam.output_data.Vtot);
set(p2,'LineWidth',2);
xlim([14 24]);
xlabel('time (h)');
ylabel('Vtot (m^3/cm^3)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);

figure('Color',[1 1 1]);
p3=plot(time,Vtot./Ntot','r--');
set(p3,'LineWidth',2);
hold on;
p3=plot(kam.output_data.tim./3600,kam.output_data.Vtot./kam.output_data.Ntot);
set(p3,'LineWidth',2);
xlim([14 24]);
ylim([0 18e-23]);
xlabel('time (h)');
ylabel('Vtot / Ntot (m^3)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);
hold off;
