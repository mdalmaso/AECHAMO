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

part_source(:,:,1) = [tvect' [zeros(885,1);0.3.*ones(30,1);0.7.*ones(35,1);zeros(length(tvect)-885-30-35,1)], 10e-9.*ones(length(tvect),1)];
part_source(:,:,2) = [tvect' [zeros(885,1);zeros(30,1);zeros(35,1);1.92.*ones(161,1);zeros(length(tvect)-885-30-35-161,1)], 18e-9.*ones(length(tvect),1)];

timevector = datevec(time_a);

A=find(timevector(:,4) == 0);
begin = A(1);
clear A;

timehours = timevector(begin:end,4)+timevector(begin:end,5)./60+timevector(begin:end,6)./3600;
clear timevector;


gas_source = [tvect' [interp1(timehours.*3600,Q(begin:end),tvect(1:1111)),zeros(1,length(tvect)-1111)]'];

dilu_coeff = [tvect' [interp1(timehours.*3600,dilu(begin:end),tvect)]'];
clear begin;


gas_source(find(isnan(gas_source)))=nanmean(Q);
dilu_coeff(find(isnan(dilu_coeff)))=nanmean(dilu);


kam.initialize('part_source',part_source,...
               'gas_source', gas_source,...
               'dilu_coeff', dilu_coeff);
kam.run;