#
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:3600;
sigma = [1.6];
N = [1e6];
mu=[50e-9];
sections = 25;
output_sections = 10*sections;
Cvap_const = 0;
Cvap0 = 1e7;
dilu_on = 1;
fixed_sections = 1;

part_source(:,:,1) = [tvect' [20.*ones(50,1);zeros(length(tvect)-50,1)], 3e-9.*ones(length(tvect),1)];
gas_source(:,:,1) = 1e9;
#