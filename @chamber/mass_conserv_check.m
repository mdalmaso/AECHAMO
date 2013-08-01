function [] = mass_conserv_check( chamber )
%MASS_CONSERV_CHECK Finds out if the mass is conserved during simulation.


initials = chamber.initials;
output = chamber.output_data;

NA = 6.022e23;

% Mass at the beginning = mass of vapor + mass of particles:
M0 = initials.Cvap0*initials.vap_molmass/NA + output.Mtot(1);

% Addition of mass by gas_source:
if(initials.Cvap_const ~= 1)
    if(~isscalar(initials.gas_source)) 
        dM_gas = trapz(initials.tvect, initials.gas_source(:,2).*initials.vap_molmass./NA);
    else
        dM_gas = initials.gas_source*initials.vap_molmass/NA*(initials.tvect(end)-initials.tvect(1));
    end
else
    error('The mass conservation check is not possible if vapor concentration is constant.');
end


% Addition of mass by part_source:
if(~isscalar(initials.part_source))
    Dp_nucl_particle = initials.part_source(1,3);
    

    M_nucl_particle = initials.particle_dens*Dp_nucl_particle^3*pi/6*1e6;

    dM_particles = trapz(initials.tvect,initials.part_source(:,2).*M_nucl_particle);
else
    dM_particles = 0;
end

% Total mass at the end = initial mass + dM_gas + dM_particles:
M_final_1 = M0 + dM_gas + dM_particles;


% On the other hand, the final mass can be calculated:
% M_final = M_aerosol + M_diluted,
% where M_aerosol is the mass of vapor and particles at the end that are
% not diluted, and M_diluted is the mass that has diluted away from aerosol
% during the simulation.

% First, the mass of aerosol at the end is the mass of particles + mass of
% vapor:
M_aerosol = output.Mtot(end) + output.vap(end)*initials.vap_molmass/NA;

% Diluted mass is:
%   mass deposited on walls (particles) 
% + mass deposited on walls (vapor)
% + mass diluted as particles
% + mass diluted as vapor:
M_diluted = output.Mwall(end) + output.Mvwall(end) + output.Mdilu(end) + output.Mvdilu(end);

M_final_2 = M_aerosol + M_diluted;

diff_M = abs(M_final_2 - M_final_1);

fprintf('M_final_1: %d\n', M_final_1);
fprintf('M_final_2: %d\n', M_final_2);
fprintf('Difference: %d\n', diff_M);
fprintf('Error (difference/M_final_1): %d %%\n',diff_M/M_final_1*100);

end

