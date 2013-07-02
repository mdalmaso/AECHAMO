#
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = [0];

vap_wallsink_on = 1;
vap_wallsink  = [0];

Dp_min = -9;
Dp_max = -6;
sections = 40;

tvect = 0:600:3*24*3600; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event
% part_source(:,:,1) = create_part_source(tvect, 1.0, 2*1*3600, 39600, 3e-9);

% gas_source(:,:,1) = create_source(tvect, 1.4662e05, 1.0557e9,39600);
% gas_source(:,:,1) = [tvect' [zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,67)]'];
gas_source(:,:,1) = [tvect' [zeros(1,66), 0.3*2e-16*60e-9*2.6908e19*100e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*2e-16*60e-9*2.6908e19*100e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*2e-16*60e-9*2.6908e19*100e-9*2.6908e19.*ones(1,12), zeros(1,67)]'];
Cvap0 = 0;

mu=5e-9;
sigma = 1.3;
N = 1e-5;

# 