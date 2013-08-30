% Nucleation test

# 1, 30 sections, nucleation at 0-30 min
method = 'moving_center';
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 30;

tvect = 0:60:7200;

part_source(:,:,1) = [tvect' [1.0.*ones(1,31), zeros(1,length(tvect)-31)]' (3e-9.*ones(1,length(tvect)))'];

Cvap_const = 1;
Cvap0 = 5e7;

mu=[80e-9];
sigma = [1.4];
N = [3e3];
#