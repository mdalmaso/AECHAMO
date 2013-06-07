function plot(obj,varargin)
% PLOT plots data calculated by simulation.
% 
% Examples:
% 
% chamber.plot;  % Plots the distribution, total volume of particles,
%                % number of particles and vapor concentration as a
%                % function of time.
% 
% chamber.plot('dist'); % Plots only the distribution.

% (c) Miikka Dal Maso & Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-06-03    0.1.1 Finds now the number of open figures and makes new
%                     figs so that the old ones are not overwritten.

% Get the number of open figures:
figs=findall(0,'type','figure');
num_figs = length(figs);

% Open a new figure:
figure(num_figs+1);

% If user has typed obj.plot('dist'), plot only the distribution.
if(nargin > 1)
    if(strcmp(varargin,'dist'))
        obj.subplot_dmps(111)
        return;
    elseif(nargin > 2) 
        error('Plot: Too many arguments.');
    else
        error('Plot: Invalid argument ''%s''.',varargin{1});
    end
end

% Plot the distribution:
obj.subplot_dmps(311)

% Plot the total particle volume in aerosol:
subplot(312)
plot(obj.output_data.tim,obj.output_data.Vtot,'b*-')

% Plot the particle concentration:
subplot(313)
plot(obj.output_data.tim,obj.output_data.Ntot,'b*-')

% Plot the vapor concentration to a new figure:
figure(num_figs + 2)
plot(obj.output_data.tim,obj.output_data.vap,'b*-')

end