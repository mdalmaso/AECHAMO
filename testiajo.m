# % 3,
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 0;

Dp_min = -8;
Dp_max = -6;
sections = 30;

tvect = 0:600:86400; % time starts and stops at 00:00

% Q = alfa*k*P, k = k_O_3*[O_3], P
gas_source = 0.6*2e-16*60e-9*2.6908e19*100e-9*2.6908e19;

mu=[100e-9];
sigma = [1.3];
N = [30000]; 

#