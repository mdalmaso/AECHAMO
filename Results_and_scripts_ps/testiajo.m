% Nucleation test

# 1, 30 sections, 0-30 min
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;

sections = 90;


output_sections = 90;

tvect = 0:60:18000;

% part_source(:,:,1) = [tvect' [1.0.*ones(1,31), zeros(1,length(tvect)-31)]' (3e-9.*ones(1,length(tvect)))'];

Cvap_const = 1;
Cvap0 = 5e7;

mu=[10e-9];
sigma = [1.4];
N = [1e5];
#