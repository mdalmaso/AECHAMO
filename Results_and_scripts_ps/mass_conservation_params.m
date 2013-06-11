% Mass conservation test with fixed sections.
% Cvap0 = 1e5
% gas_source is a vector
% sedimentation is on
% dilution is a vector
% part_source is a vector

#
dilu_on = 1;
coag_on = 1;
sedi_on = 1;
fixed_sections = 1;
Cvap_const = 0;

Cvap0 = 1e5;

tvect = 0:60:7200;

part_source = [tvect', [1.0.*ones(1,31), zeros(1,length(tvect)-31)]', (3e-9.*ones(1,length(tvect)))'];

gas_source = [tvect', (1e2.*tvect)'];

dilu_coeff = [tvect', (1.00.*370/1460000.*tvect./6800)'];

N = 1e5;
mu = 15e-9;
sigma = 1.4;

sections = 40;
output_sections = 400;

#