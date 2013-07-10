#
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = 1;
dilu_coeff = 3.5e-4;

vap_wallsink_on = 0;
vap_wallsink  = [0];

Dp_min = -9;
Dp_max = -6;
sections = 30;
output_sections = 300;

tvect = 0:420:10*3600; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [4.7.*ones(9,1);zeros(7,1);3.0.*ones(15,1);zeros(length(tvect)-31,1)], 10e-9.*ones(length(tvect),1)];

% gas_source(:,:,1) = [tvect' [zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,67)]']; % 11h no event, 2h event, 11h no event

gas_source = 8.4e5;
Cvap0 = 8.1e6;

mu=34e-9;
sigma = 1.65;
N = 0;

# 