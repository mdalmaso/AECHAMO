% SOA formation 1day, timestep 10min, N 1000 1/cm3

# % 1, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(24*60*60);

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.3*9e-17*30e-9*2.6908e19*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

# % 2, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.3*9e-17*30e-9*2.6908e19*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];
#

% 3, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(24*60*60);

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.6*9e-17*30e-9*2.6908e19*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

# % 4, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.6*9e-17*30e-9*2.6908e19*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

#
% 5, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(24*60*60);

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.9*9e-17*30e-9*2.6908e19*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

# % 6, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.9*9e-17*30e-9*2.6908e19*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

# % 7, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(24*60*60);

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

# % 8, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

# % 9, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 1/(24*60*60);

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.3*9e-17*30e-9*2.6908e19*1.5e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];

# % 10, 
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -9;
Dp_max = -6;
sections = 60;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.3*9e-17*30e-9*2.6908e19*1.5e-9*2.6908e19;

mu=[10e-9];
sigma = [1.3];
N = [1000];
#