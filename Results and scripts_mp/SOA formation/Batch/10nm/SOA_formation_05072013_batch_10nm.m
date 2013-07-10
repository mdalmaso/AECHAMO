% 10h test, SOA formation

#
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 5e-5;

vap_wallsink_on = 1;
vap_wallsink  = [1/50; 1/500];

Dp_min = -9;
Dp_max = -5;
sections = 40;

tvect = 0:60:10*3600;

gas_source(:,:,1) = [tvect' [0.3*2e-16*60e-9*2.6908e19*10e-9*2.6908e19*exp(-2e-16.*60e-9.*2.6908e19.*tvect)]'];
gas_source(:,:,2) = [tvect' [0.3*2e-16*60e-9*2.6908e19*50e-9*2.6908e19*exp(-2e-16.*60e-9.*2.6908e19.*tvect)]'];
gas_source(:,:,3) = [tvect' [0.3*2e-16*60e-9*2.6908e19*100e-9*2.6908e19*exp(-2e-16.*60e-9.*2.6908e19.*tvect)]'];
gas_source(:,:,4) = [tvect' [0.3*2e-16*60e-9*2.6908e19*200e-9*2.6908e19*exp(-2e-16.*60e-9.*2.6908e19.*tvect)]'];

Cvap0 = 0;

mu=[10e-9];
sigma = 1.3;
N = mass2number(mu,[0.003; 0.015; 0.03; 0.1],1.84); % [...] mass in µg/m3

% diff_coeff = 1.093392158075694e-1; % The condensing vapor diffusion coefficient (cm^2/s)
% vap_molmass = 300; % Molecular mass of condensing vapor (g/mol)
% particle_dens = 1.4; % Density of particle matter. (g/cm^3)
% lambda = 129.29e-9; % The condensing vapor mean free path (m).

# 
