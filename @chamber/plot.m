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
% 
% chamber.plot('original'); % Plots the figures using the original
% distribution, not the interpolated one.
% 
% chamber.plot('dist', 'original') % Plots only the original distribution.

% (c) Miikka Dal Maso & Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-06-03    0.1.1 Finds now the number of open figures and makes new
%                     figs so that the old ones are not overwritten.
% 2013-06-10    0.1.2 Takes now also argument 'smoothed' to plot the
%                     smoothed distribution instead of the original.

% If user has typed obj.plot('dist'), plot only the distribution.
% In case of obj.plot('smoothed'), plot the smoothed distribution.

plot_distr = 0;
% plot_smoothed = 0;
plot_original = 0;

% Check if there is additional arguments:
if(nargin > 1)
    for i=2:nargin
        switch(varargin{i-1})
            case('dist')
                plot_distr = 1;
%             case('smoothed')
%                 plot_smoothed = 1;
            case('original')
                plot_original = 1;
            otherwise
                error('Invalid argument: ''%s''.', varargin{i-1});
        end
    end
end

% Get the number of open figures:
figs=findall(0,'type','figure');
num_figs = length(figs);

% Open a new figure:
figure(num_figs+1);

% If one of the arguments is 'dist', plot only the distribution.
if(plot_distr == 1)
%     if(plot_smoothed == 1)
%         % If there is argument 'smoothed' in addition, plot only the
%         % smoothed distribution:
%         obj.subplot_dmps(111,'smoothed');
%         
%         % And terminate the function:
%         return;
    if(plot_original == 1)
        obj.subplot_dmps(111,'original');
        return;
    else
        % Otherwise plot the interpolated distribution:
        obj.subplot_dmps(111);
        
        % And terminate the function:
        return;
    end
end

% Plot the distribution:
% if(plot_smoothed == 1)
%     % If there was argument 'smoothed', plot the smoothed distribution:
%     obj.subplot_dmps(311,'smoothed');
if(plot_original == 1)
    obj.subplot_dmps(311,'original');
else
    % Otherwise plot the original distribution:
    obj.subplot_dmps(311);
end

% Plot the total particle volume in aerosol:
subplot(312)
plot(obj.output_data.tim/(24*60*60),obj.output_data.Vtot,'b*-')
xhandle = xlabel('time (d)');
yhandle = ylabel('V_{tot}(cm^{3})','rotation',90);

% Plot the particle concentration:
subplot(313)
plot(obj.output_data.tim/(24*60*60),obj.output_data.Ntot,'b*-')
xhandle = xlabel('time (d)');
yhandle = ylabel('N_{tot}(cm^{-3})','rotation',90);

% Plot the vapor concentration to a new figure:
figure(num_figs + 2)
plot(obj.output_data.tim/(24*60*60),obj.output_data.vap,'b*-')
xhandle = xlabel('time (d)');
yhandle = ylabel('C_{vap}(cm^{-3})','rotation',90);

% Plot the distribution:
figure(num_figs + 3)
% if(plot_smoothed == 1)
%     % If there was argument 'smoothed', plot the smoothed distribution:
%     obj.subplot_dmps(111,'smoothed');
if(plot_original == 1)
    obj.subplot_dmps(111,'original');
else
    % Otherwise plot the interpolated distribution:
    obj.subplot_dmps(111);
end

end