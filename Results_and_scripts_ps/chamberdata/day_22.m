kam=chamber;

clear tvect;
clear part_source;
clear gas_source;
clear dilu_coeff;

tvect = 53100:60:24*3600;

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
               'vap_molmass', 300,...
               'diff_coeff', 0.0489,...
               'particle_dens', 1.4,...
               'lambda', 1.0265e-7,...
               'N', 0,...
               'mu', 37e-9,...
               'sigma', 1.6);

part_source(:,:,1) = [tvect' [0.3.*ones(30,1);0.7.*ones(35,1);zeros(length(tvect)-30-35,1)], 10e-9.*ones(length(tvect),1)];
part_source(:,:,2) = [tvect' [zeros(30,1);zeros(35,1);1.92.*ones(161,1);zeros(length(tvect)-30-35-161,1)], 18e-9.*ones(length(tvect),1)];

timevector = datevec(time_a);

A=find(timevector(:,4) == 0);
begin = A(1);
clear A;

timehours = timevector(begin:end,4)+timevector(begin:end,5)./60+timevector(begin:end,6)./3600;
clear timevector;


gas_source = [tvect' [interp1(source_time,0.38.*Q_MTOH+0.09.*Q_MTO3,tvect(1:end))]'];

dilu_coeff = [tvect' [interp1(timehours.*3600,dilu(begin:end),tvect)]'];
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


timevector=datevec(time);
A=find(timevector(:,4) == 0);
begin = A(1);
clear A;
timevector = timevector(begin:end,:);

timeax = timevector(:,4)+timevector(:,5)./60+timevector(:,6)./3600;


figure('Color',[1 1 1]);
p1=plot(timeax,Ntot(begin:end),'r--');
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
ylim([0 5000]);

figure('Color',[1 1 1]);
p2=plot(timeax,Vtot(begin:end),'r--');
set(p2,'LineWidth',2);
hold on;
p2=plot(kam.output_data.tim./3600,kam.output_data.Vtot);
set(p2,'LineWidth',2);
xlim([14 24]);
xlabel('time (h)');
ylabel('Vtot (m^3/cm^3)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);
ylim([0 6e-19]);

figure('Color',[1 1 1]);
p3=plot(timeax,Vtot(begin:end)./Ntot(begin:end)','r--');
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

figure('Color',[1 1 1]);
% p4=plot(kam.initials.tvect./3600,kam.initials.part_source(:,2,1)+kam.initials.part_source(:,2,2));
p4=plot(kam.initials.tvect./3600,kam.initials.part_source(:,2,1));
set(p4,'LineWidth',2);
xlim([14 24]);
ylim([0 4.5]);
xlabel('time (h)');
ylabel('Nucleation rate (1/cm^3s)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);
hold off;