% Example of using functions create_gas_source and create_part_source.
% The simulation consists of two days, so that both days vapor source
% begins at 9:00 and particle source at 10:00. Initially there are 1000
% particles/cc.

#
method = 'moving_center';
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(3.*24.*3600);

mu=100e-9;
sigma = 1.4;
N = 1000;

Dp_min = -9;
Dp_max = -6;
sections = 30;

tvect = 0:600:15*24*3600;
Cvap0 = 0;

gas_source = create_gas_source(tvect,6*3600,5e4.*(6*3600),10*3600);

part_source(:,:,1) = create_part_source(tvect,1000/(3.*24.*3600),1000./(3.*3600),0,3e-9);
part_source(:,:,2) = create_part_source(tvect,1,2*3.6e3,10*3600,3e-9);

# 