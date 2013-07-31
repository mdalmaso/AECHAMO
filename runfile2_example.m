#
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.6];
N = [1e4];
mu=[50e-9];
sections = 25;
output_sections = 10*sections;
Cvap_const = 0;
Cvap0 = 1e7;
dilu_on = 1;


dilu_coeff = [tvect' 10.*tvect'];
#