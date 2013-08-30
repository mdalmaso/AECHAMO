% Example script for chamber_runfile. 3 different runs.

#
method = 'moving_center';
sedi_on = 0;
coag_on = 1;
dilu_on = 1;

Dp_min = -9;
Dp_max = -6;
sections = 30;

tvect = 0:600:2*24*3600;
gas_source = 1e5;
Cvap0 = 0;

mu=[10e-9];
sigma = 1.3;
N = [4000];

# 

method = 'moving_center';
sedi_on = 0;
coag_on = 1;
dilu_on = 1;

Dp_min = -9;
Dp_max = -6;
sections = 30;

tvect = 0:600:2*24*3600;
gas_source = 1e6;
Cvap0 = 0;

mu=[10e-9];
sigma = 1.3;
N = [4000];

#
method = 'moving_sectional';
sedi_on = 0;
coag_on = 1;
dilu_on = 1;

Dp_min = -9;
Dp_max = -6;
sections = 30;

tvect = 0:600:2*24*3600;
gas_source = 1e6;
Cvap0 = 0;

mu=[10e-9];
sigma = 1.3;
N = [4000];
#
