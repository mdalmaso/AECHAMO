function[out_struct] = model_convert(obj,t,Y)
% MODEL_CONVERT creates a structure of results calculated by simulation.
% 
% [out_struct] = model_convert_moving_center(obj, t, Y)
% t is the time vector and Y is the matrix calculated by simulation.
% The data from Y is exported to out_struct, so that out_struct contains
% following fields:
% 
% CMD - Count Median Diameter
% Ntot - Total particle concentration
% Vtot - Total particle volume
% Dpmean - Mean diameter
% Mtot - Total mass of particles
% Mwall - Mass lost to wall as particles
% Mdilu - Mass diluted as aerosols
% Mvdilu - Mass diluted as vapor
% Mvwall - Mass condensed to wall as vapor
% distr - Distribution as a function of time
% distr_smoothed - Smoothed distribution as a function of time. This is
%                  created only if sections are fixed. Otherwise there is
%                  no need for smoothed distribution.

% (c) Miikka Dal Maso & Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-06-06    0.1.1 Dp0 will be defined differently if fixed sections is
%                     on.
% 2013-06-10    0.1.2 Added distribution smoothing if fixed sections are
%                     used.
% 2013-06-20    0.1.3 Added condensation sink to chamber walls. The mass of
%                     vapor condensed to walls is saved into Mvwall.

initials = obj.initials;

% converts schlosser model results and plots them
% translate back 
mv     = initials.vap_molmass;
rool   = initials.particle_dens;

% If initials.sections is scalar, it tells the number of sections.
% Otherwise it is an array of all sections, and the number of secs is the
% length of this array.
if(isscalar(initials.sections))
    nSec   = initials.sections;
else
    nSec = length(unique(initials.sections));
end


% Create a new lognormal vector containing the diameters. The number of
% sections in this vector should be much more than the number of sections
% in original distribution. This is because interpolation may cause a big
% error if this grid is too sparse and the distribution is simultaniously
% very narrow.
Dp0 = logspace(initials.Dp_min, initials.Dp_max, initials.output_sections);

% Preallocate the variables for the loop.
Ntot=zeros(1,length(t));
dN=zeros(length(t),length(Dp0));

dist_original=zeros(length(t)+1,2*length(Y(1,nSec+2:(2*nSec+1)))+2);

dist_original(2:end,1) = t;

Vtot=zeros(1,length(t));
out_struct.CMD = zeros(1,length(t));

for i = 1:length(t),
    Ntot(i) = sum(Y(i,2:nSec+1));
    Ni = Y(i,2:nSec+1);
    Dpi = Y(i,nSec+2:(2*nSec+1));
    
%     if(initials.fixed_sections ~= 0)
%         for j=1:floor(length(Dpi)/2)
%             Ni_smoothed(j)=Ni(2*j) + Ni(2*j-1);
%             if(Ni_smoothed(j) > 0)
%                 Dpi_smoothed(j) = ((Ni(2*j)*Dpi(2*j)^3 + Ni(2*j-1)*Dpi(2*j-1)^3)/Ni_smoothed(j))^(1/3);
%             else
%                 Dpi_smoothed(j) = (Dpi(2*j)^3 + Dpi(2*j-1)^3)^(1/3);
%             end
%         end
%         dNi_smoothed = obj.N_to_dlog(Dpi_smoothed,Ni_smoothed);
%         dN_smoothed(i,:) = interp1(log10(Dpi_smoothed),dNi_smoothed,log10(Dp0),'linear',0);
%     end

    dNi = obj.N_to_dlog(Dpi,Ni);
    
    dist_original(i+1, 3:end)=[dNi, Dpi];
    dist_original(i+1, 2) = Ntot(i);
%     dN(i,:) = interp1(Dpi,dNi,Dp0,'linear',0);
    
    dN(i,:) = interp1(log10(Dpi),dNi,log10(Dp0),'linear',0);
    
%     Vtot(i) = obj.distribution_info_Vtot(Dp0,dN(i,:));
    
%     Vtot(i) = obj.distribution_info_Vtot(Dpi,dNi);
%     dV = (pi./6).*Dpi.^3.*dNi;
%     Vtot(i) = trapz(log10(Dpi),dV);
    % Vtot = N(Dp)*V(Dp)
    dV = (pi./6).*Dpi.^3;
    Vtot(i) = sum(Ni.*dV);
    
    % Calculate CMD. Interpolate Ni to Dp0 and find such a point in Dp0
    % that the amount of particles that have smaller diameter equal the
    % amount of particles that have bigger diameter.
    Ni2=interp1(Dpi,Ni,Dp0,'linear',0); %First interpolate Ni to denser grid.
    B=0;  % Preallocate
    j=1;  % Preallocate
    while(B < Ntot(i)/2 && j < length(Dp0)) % Run the loop as long as the sum of particles is smaller than the total number of particles.
        B=sum(Ni2(1:j));  % Sum the number of particles from the beginning of the distribution to j:s index of distribution.
        j = j+1;
    end
    % Now we know that the half of total number of particles is in sections
    % Dp0(1:j). In other words, CMD = Dp0(j).
    out_struct.CMD(i)=Dp0(j);
    
end

out_struct.distr_original = dist_original;

% Get the distribution and save it in format:
% 
% (time) (Ntot)     Dp0(1)         Dp0(2)        ...  Dp0(n)
%  t1    Ntot(t1)   N(Dp0(1), t1)  N(Dp0(2), t1)      N(Dp0(n), t1)
%  .      .         .              .                  .
%  .      .         .              .                  .
%  .      .         .              .                  .

[ro, co] = size(dN);
dist = zeros(ro+1,co+2); % Preallocate the distribution array.
dist(2:end,1) = t(:);    % Insert the time vector to the first column of dist.
dist(2:end,2) = Ntot(:); % Insert Ntot to the second column of dist.
dist(1,3:end) = Dp0;     % The first row of dist tells the Dp0s of distribution.
dist(2:end,3:end) = dN;  % Then each column tells the particle concentration
                         % of corresponding section (Dp0) for all time
                         % points.
                         
out_struct.distr = dist; % Save the distribution to output.

% if(initials.fixed_sections ~= 0)
%     % Make the smoothed distribution:
%     [ro, co] = size(dN_smoothed);
%     dist_smoothed = zeros(ro+1,co+2); % Preallocate the distribution array.
%     dist_smoothed(2:end,1) = t(:);    % Insert the time vector to the first column of dist.
%     dist_smoothed(2:end,2) = Ntot(:); % Insert Ntot to the second column of dist.
%     dist_smoothed(1,3:end) = Dp0;     % The first row of dist tells the Dp0s of distribution.
%     dist_smoothed(2:end,3:end) = dN_smoothed;  % Then each column tells the particle concentration
%                              % of corresponding section (Dp0) for all time
%                              % points.
% 
%     out_struct.distr_smoothed = dist_smoothed; % Save the distribution to output.
% else
%     out_struct.distr_smoothed = NaN;
% end


out_struct.Y = Y; % Save the original Y:

out_struct.vap = Y(:,1); % Save the vapor concentration.
out_struct.tim = t(:);   % Save the time vector.


% get a mean volume per particle:
vpart = Vtot./Ntot;
% get a mean diameter of the particles:
Dpmean = (6.*vpart./pi).^(1./3);

out_struct.Ntot = Ntot(:); % Save Ntot values.
out_struct.Vtot = Vtot(:); % Save Vtot values
out_struct.Dpmean= Dpmean(:); % Save Dpmean values

% mass balance stuff
NA=6.022e23;
% all [M] = g/cm^3 !!!
out_struct.Mtot = Vtot.*rool*1e6; % Total mass of particles ([Vtot]=m^3/cm^3, [rool]=g/cm^3)
out_struct.Mwall= Y(:,2*nSec+2).*mv./NA; % Mass lost to wall
out_struct.Mdilu= Y(:,2*nSec+3).*mv./NA; % Mass diluted as aerosols
out_struct.Mvdilu=Y(:,2*nSec+4).*mv./NA; % Mass diluted in gas phase
if(length(Y(1,:)) == 2*nSec+5)
    out_struct.Mvwall = Y(:,2*nSec+5).*mv./NA; % Mass lost to wall as vapor
else out_struct.Mvwall = 0;
end

end
