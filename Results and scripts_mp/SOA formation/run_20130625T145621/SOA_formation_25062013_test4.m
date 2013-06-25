% 20d test, SOA formation

#
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;
dilu_coeff = 1/(24*3600);

vap_wallsink_on = 1;
vap_wallsink  = 1/9000;

Dp_min = -9;
Dp_max = -5;
sections = 20;

tvect = 0:600:1728000; % time starts and stops at 00:00

% part_source(:,:,1) = [tvect' (1.0.*ones(1,length(tvect)))' (3e-9.*ones(1,length(tvect)))'];

gas_source(:,:,1) = 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19;

mu=100e-9;
sigma = 1.3;
N = 1000;

# 