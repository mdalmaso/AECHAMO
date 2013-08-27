% Example script for chamber_runfile2. 
% Four different runs: dilution coefficient and initial N is varied.

#
tvect = 0:600:24*3600;
gas_source = 1e6;
mu= [5e-9, 180e-9]
sigma = [1.2, 1.7];

dilu_coeff(:,:,1) = 2e-4;
dilu_coeff(:,:,2) = 0;
N = [4e3 0; 4e3 1000];
#
