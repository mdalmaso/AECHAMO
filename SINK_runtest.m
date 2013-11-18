
% run the sink tests

% method = 'moving_center';
% sedi_on = 0;
% coag_on = 1;
% dilu_on = 1;
% dilu_coeff = 1/(3.*24.*3600);
% 
% mu=100e-9;
% sigma = 1.4;
% N = 1000;
% 
% Dp_min = -9;
% Dp_max = -6;
% sections = 30;
% 
% tvect = 0:600:15*24*3600;
% Cvap0 = 0;
% 
% gas_source = create_gas_source(tvect,6*3600,5e4.*(6*3600),10*3600);
% 
% part_source(:,:,1) = create_part_source(tvect,1000/(3.*24.*3600),1000./(3.*3600),0,3e-9);
% part_source(:,:,2) = create_part_source(tvect,1,2*3.6e3,10*3600,3e-9);
% run a load of 5e5 runs with sinks



close all
clear classes

NN = 300;
DpSINK = 100e-9;
filen = 'Q00_1e5_300_2.mat';

a = chamber;
a.initialize('method','moving_center','sedi_on',0,'coag_on',1,'dilu_on',1,'dilu_coeff',1/(2.*24.*3600));
a.initialize('mu',100e-9,'N',NN,'sigma',1.4);
a.initialize('Dp_min',-9,'Dp_max',-6,'sections',30);
timev = 0:600:15*24*3600;
a.initialize('tvect',timev);
a.initialize('Cvap0',0);

gas_q = create_gas_source(timev,6*3600,1e5.*(6*3600),10*3600); 
big_rate = NN./(3.*24.*3600);
source_length = 24.*3600; % = Ntot/rate_max
Ntot_big = big_rate.*source_length;
a.initialize('gas_source',gas_q)

part_q= create_part_source(timev,1,2*3.6e3,10*3600,3e-9);

a.initialize('part_source',part_q)

a.run

save(filen,'a')
%*******


close all
clear classes

NN = 300;
DpSINK = 100e-9;
filen = 'Q00_1e5_300_1.mat';

a = chamber;
a.initialize('method','moving_center','sedi_on',0,'coag_on',1,'dilu_on',1,'dilu_coeff',1/(1.*24.*3600));
a.initialize('mu',100e-9,'N',NN,'sigma',1.4);
a.initialize('Dp_min',-9,'Dp_max',-6,'sections',30);
timev = 0:600:15*24*3600;
a.initialize('tvect',timev);
a.initialize('Cvap0',0);

gas_q = create_gas_source(timev,6*3600,1e5.*(6*3600),10*3600); 
big_rate = NN./(3.*24.*3600);
source_length = 24.*3600; % = Ntot/rate_max
Ntot_big = big_rate.*source_length;
a.initialize('gas_source',gas_q)

part_q= create_part_source(timev,1,2*3.6e3,10*3600,3e-9);

a.initialize('part_source',part_q)

a.run

save(filen,'a')
%*******
