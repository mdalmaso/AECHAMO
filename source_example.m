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

Dp_min = -9;
Dp_max = -6;
sections = 30;

tvect = 0:600:2*24*3600;
Cvap0 = 0;

mu=[180e-9];
sigma = 1.4;
N = [1000];

gas_source = create_gas_source(tvect,2*3600,7.2e8,9*3600);
part_source = create_part_source(tvect,5,3.6e4,10*3600,3e-9);


# 