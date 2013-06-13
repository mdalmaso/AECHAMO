% Validation 1day, timestep 10min

# % 1, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [5000, 2500];

# % 2, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [1000, 500];

# % 3, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 6e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [5000, 2500];

# % 4, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 1e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [5000, 2500];

# % 5, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 2.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [5000, 2500];

# % 6, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 3.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [5000, 2500];

# % 7, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dill_coeff = 1/(24*60*60);

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [5000, 2500];

# % 8, 10sect
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dill_coeff = 1/(48*60*60);

Dp_min = -9;
Dp_max = -6;
sections = 10;

tvect = 0:600:86400; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

Cvap_const = 1;
Cvap0 = 3e7;

% aithken and accumulation
mu=[30e-9, 110e-9];
sigma = [1.3, 1.3];
N = [5000, 2500];

#