kam=chamber;

clear tvect;
clear part_source;
clear gas_source;
clear dilu_coeff;

tvect = 51480:60:24*3600;

kam.initialize('fixed_sections', 1,...
               'sedi_on', 1,...
               'coag_on', 1,...
               'dilu_on',1,...
               'vap_wallsink_on', 1,...
               'vap_wallsink', 1/120,...
               'Dp_min', -9,...
               'Dp_max', -6,...
               'sections', 25,...
               'output_sections', 250,...
               'tvect', tvect,...
               'Cvap0', 0,...
               'N', 0,...
               'mu', 37e-9,...
               'sigma', 1.6,...
               'vap_molmass', 300,...
               'diff_coeff', 0.0489,...
               'particle_dens', 1.4,...
               'lambda', 1.0265e-7,...               
               'T', 273.15+16);

part_source(:,:,1) = [tvect' [24.0.*ones(21,1);zeros(length(tvect)-21,1)], 10e-9.*ones(length(tvect),1)];
part_source(:,:,2) = [tvect' [zeros(21,1);10.*ones(11,1);1.15.*ones(98,1);3.97.*ones(64,1);3.4.*ones(133,1);zeros(length(tvect)-21-11-98-64-133,1)], 20e-9.*ones(length(tvect),1)];

timevector = datevec(time_a);

begin = 1;
last = 1679;

timehours = timevector(begin:last,4)+timevector(begin:last,5)./60+timevector(begin:last,6)./3600;
clear timevector;


% gas_source = [tvect' [interp1(timehours.*3600,Q(begin:last),tvect(1:327),'linear',0),zeros(1,length(tvect)-327)]'];
gas_source = [tvect' [interp1(source_time,0.38.*Q_MTOH+0.09.*Q_MTO3,tvect(1:end))]'];

dilu_coeff = [tvect' [interp1(timehours.*3600,dilu(begin:last),tvect)]'];
clear begin;


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

timevector=datevec(time);
begin = 1;
last = 121;
timevector = timevector(begin:last,:);

timeax = timevector(:,4)+timevector(:,5)./60+timevector(:,6)./3600;


figure('Color',[1 1 1]);
p1=plot(timeax,Ntot(begin:last),'r--');
set(p1,'LineWidth',2);
hold on;
p1=plot(kam.output_data.tim./3600,kam.output_data.Ntot);
set(p1,'LineWidth',2);
xlim([14 24]);
% ylim([0 5000]);
xlabel('time (h)');
ylabel('Ntot (1/cm^3)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);

figure('Color',[1 1 1]);
p2=plot(timeax,Vtot(begin:last),'r--');
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
p3=plot(timeax,Vtot(begin:last)./Ntot(begin:last)','r--');
set(p3,'LineWidth',2);
hold on;
p3=plot(kam.output_data.tim./3600,kam.output_data.Vtot./kam.output_data.Ntot);
set(p3,'LineWidth',2);
xlim([14 24]);
ylim([0 3.5e-22]);
xlabel('time (h)');
ylabel('Vtot / Ntot (m^3)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);
hold off;

figure('Color',[1 1 1]);
p4=plot(kam.initials.tvect./3600,kam.initials.part_source(:,2,1)+kam.initials.part_source(:,2,2));
set(p4,'LineWidth',2);
xlim([14 24]);
ylim([0 4.5]);
xlabel('time (h)');
ylabel('Nucleation rate (1/cm^3s)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);
hold off;

% clear begin;