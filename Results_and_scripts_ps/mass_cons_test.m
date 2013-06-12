% Lets particles grow over the grid limits

#
coag_on = 0;
dilu_on = 0;
fixed_sections = 1;
Cvap_const = 0;

Cvap0 = 6e7;
gas_source = 1e5;

tvect = 0:60:10800;

sections = 30;
output_sections = 90;

Dp_min = -9;
Dp_max = -8;

N = 1e5;
mu = 7e-9;
sigma = 1.05;

#