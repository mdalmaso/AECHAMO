kam=chamber;

clear tvect;
clear part_source;
clear gas_source;
clear dilu_coeff;


tvect = 50520:60:24*3600;

kam.initialize('fixed_sections', 1,...
               'sedi_on', 0,...
               'coag_on', 1,...
               'dilu_on',1,...
               'vap_wallsink_on', 1,...
               'vap_wallsink', 1/250,...
               'Dp_min', -9,...
               'Dp_max', -6,...
               'sections', 25,...
               'output_sections', 250,...
               'tvect', tvect,...
               'Cvap0', 0,...
               'N', 0,...
               'mu', 37e-9,...
               'sigma', 1.6,...
               'vap_molmass',250,...
               'T', 273.15+16);

part_source(:,:,1) = [tvect' [36.3.*ones(22,1);zeros(length(tvect)-22,1)], 3e-9.*ones(length(tvect),1)];
part_source(:,:,2) = [tvect' [zeros(22,1);6.5.*ones(77,1);11.5.*ones(230,1);zeros(length(tvect)-22-77-230,1)], 3e-9.*ones(length(tvect),1)];


timevector = datevec(time_a);

A=find(timevector(:,4) == 0);
begin1 = A(1);
clear A;

timehours = timevector(begin1:end,4)+timevector(begin1:end,5)./60+timevector(begin1:end,6)./3600;
clear timevector;


gas_source = [tvect' [interp1(timehours.*3600,Q(begin1:end),tvect(1:328)),zeros(1,length(tvect)-328)]'];

dilu_coeff = [tvect' [interp1(timehours.*3600,dilu(begin1:end),tvect)]'];
% clear begin;


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
A=find(timevector(:,4) == 0);
begin = A(1);
clear A;
timevector = timevector(begin:end,:);

timeax = timevector(:,4)+timevector(:,5)./60+timevector(:,6)./3600;

Y = kam.output_data.Y;
sections = kam.initials.sections;
% Y_smps = [Y(:,indices),Y(:,sections+indices:end)];
Ntot_smps = 0;
for i = 1:length(Y(:,1))
    indices = find(Y(i,sections+2:2*sections+1) > 15e-9);
    Ntot_smps(i) = sum(Y(i,indices+1));
    Vtot_smps(i) = sum(Y(i,indices+1).*(pi/6.*Y(i,sections+indices+1).^3));
end

figure('Color',[1 1 1]);
% p1=plot(timeax,Ntot(begin:end),'r--');
p1=plot(timehours,data.mcpc(begin1:end),'r--');
set(p1,'LineWidth',2);
hold on;
p1=plot(kam.output_data.tim./3600,kam.output_data.Ntot);
set(p1,'LineWidth',2);
p1=plot(kam.output_data.tim./3600,Ntot_smps,'--');
set(p1,'LineWidth',2);
p1=plot(timeax,Ntot(begin:end),'r');
set(p1,'LineWidth',2);
xlim([14 24]);
% ylim([0 5000]);
xlabel('time (h)');
ylabel('Ntot (1/cm^3)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);

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
p4=plot(kam.initials.tvect./3600,kam.initials.part_source(:,2,1)+kam.initials.part_source(:,2,2));
set(p4,'LineWidth',2);
xlim([14 24]);
ylim([0 4.5]);
xlabel('time (h)');
ylabel('Nucleation rate (1/cm^3s)');
set(gca,'FontSize',18);
set(findall(gcf,'type','text'),'FontSize',18);
hold off;