% 20d test, SOA formation

#
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
dilu_on = [1;0];
dilu_coeff = [1/(24*3600); 1/(5*24*3600)];

vap_wallsink_on = [1;0];
vap_wallsink  = [1/9; 1/90; 1/900; 1/9000; 1/90000];

Dp_min = -9;
Dp_max = -5;
sections = 40;

tvect = 0:600:1728000; % time starts and stops at 00:00

part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66), zeros(1,66), 1.0.*ones(1,12), zeros(1,66),zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))']; % 11h no event, 2h event, 11h no event

gas_source(:,:,1) = [tvect' [zeros(1,66), 0.3*2e-16*60e-9*2.6908e19*100e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66), zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,66),zeros(1,66), 2.606545670400000e+08.*ones(1,12), zeros(1,67)]']; % 11h no event, 2h event, 11h no event
gas_source(:,:,2) = [tvect' [zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66), zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,66),zeros(1,66), 0.3*9e-17*30e-9*2.6908e19*0.25*1e-9*2.6908e19.*ones(1,12), zeros(1,67)]']; % 11h no event, 2h event, 11h no event


mu=100e-9;
sigma = 1.3;
N = 0;

# 