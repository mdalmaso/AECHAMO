# % 1, 30 sections, Nucleation event everyday between 11:00-13:00
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 30;

tvect = 0:600:3600*24; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 0;
Cvap0 = 0;
gas_source = 5e4;

mu=100e-9;
sigma = 1.3;
N = 0;

#