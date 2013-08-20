% 4day test3, SOA formation

#
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 5e-6;

vap_wallsink_on = 1;
vap_wallsink  = 5e-5;

Dp_min = -9;
Dp_max = -6;
sections = 60;

out = load('MT_night_and_day.mat');
tvect = 0:60:4*24*3600-1*3600;
gas_source(:,:,1) = [tvect' 0.08*out.out.model.Q_Cvap];

part_source = create_part_source(tvect,1.0, 1*3*3600, 39600, 3e-9);

Cvap0 = 0;

diff_coeff = 0.0489; % The condensing vapor diffusion coefficient (cm^2/s)
vap_molmass = 300; % Molecular mass of condensing vapor (g/mol)
particle_dens = 1.4; % Density of particle matter. (g/cm^3)
lambda =  1.0265e-07; % The condensing vapor mean free path (m).

mu= 80e-9;
sigma = 1.6;
N = 500; % 'aithken' from dalmaso2008

# 
