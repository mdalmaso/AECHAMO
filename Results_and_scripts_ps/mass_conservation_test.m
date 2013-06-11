% Mass conservation test.

kam=chamber_runfile('mass_conservation_params.m');

initials = kam.initials;
output = kam.output_data;

NA = 6.022e23;

% Mass at the beginning = mass of vapor + mass of particles:
M0 = initials.Cvap0*initials.vap_molmass/NA + output.Mtot(1);

% Addition of mass by gas_source:
dM_gas = trapz(initials.tvect, initials.gas_source(:,2).*initials.vap_molmass./NA);

% Addition of mass by part_source:
Dp_nucl_particle = initials.part_source(1,3);
M_nucl_particle = initials.particle_dens*Dp_nucl_particle^3*pi/6;

dM_particles = trapz(initials.tvect,initials.part_source(:,2).*M_nucl_particle);

% Total mass at the end = initial mass + dM_gas + dM_particles:
M_final_1 = M0 + dM_gas + dM_particles;


% On the other hand, the final mass can be calculated:
% M_final = M_aerosol + M_diluted,
% where M_aerosol is the mass of vapor and particles at the end that are
% not diluted, and M_diluted is the mass that has diluted away from aerosol
% during the simulation.

M_aerosol = output.Mtot(end);

M_diluted = output.Mwall(end) + output.Mdilu(end) + output.Mvdilu(end);

M_final_2 = M_aerosol + M_diluted;

diff_M = abs(M_final_2 - M_final_1);

display('M_final_1: %d', M_final_1);
display('M_final_2: %d', M_final_2);
display('Difference: %d', diff_M);
