% 20d test, timestep 10min

# % 1, dilu_on, 1day
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(24*3600);

Dp_min = -9;
Dp_max = -6;
sections = 40;

tvect = 0:600:1728000; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66),zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

mu=100e-9;
sigma = 1.3;
N = 10;

# % 2, Cvap
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 40;

tvect = 0:600:1728000; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66),zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e6;

mu=100e-9;
sigma = 1.3;
N = 10;

# % 3, Cvap
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 40;

tvect = 0:600:1728000; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66),zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 6e7;

mu=100e-9;
sigma = 1.3;
N = 10;

# % 4, Cvap, Dil_on, 1/2day
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(12*3600);

Dp_min = -9;
Dp_max = -6;
sections = 40;

tvect = 0:600:1728000; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66),zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 6e7;

mu=100e-9;
sigma = 1.3;
N = 10;

# % 5, Dil_on, 1/2day
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(12*3600);

Dp_min = -9;
Dp_max = -6;
sections = 40;

tvect = 0:600:1728000; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66),zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

mu=100e-9;
sigma = 1.3;
N = 10;

#