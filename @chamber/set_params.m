function set_params(obj, varargin)
% SET_PARAMS Sets the parameters for chamber object.
%
% Sets the initial parameters of the chamber object (obj) but does not
% check the correctness of them. In addition, this function does not apply
% the parameters Dp, Dplims, nubmer_distr or center_diameters to the
% chamber object. These are applied by calling function form_distribution
% after setting the parameters.
% 
% For initializing the values, use method chamber.initialize instead!
%
% obj.set_params('field_name_1', value_1, ...
%               'field_name_2', value_2, ...
%               'field_name_n', value_n)
% sets the obj's initial parameters. Fields with names field_name_1, 
% field_name_2 and field_name_n are replaced by values value_1, value_2 
% and value_n, respectively. Other fields will not be changed.
% 
% This function is used by chamber.initialize.
% 
% THE PARAMETERS OF THE CHAMBER MODEL
% 
% For more information, see the documentation.
% 
% SWITCHES:
% dilu_on           Defines whether the dilution is on or not. If
%                   dilu_on = 0, there will not happen dilution in the
%                   chamber model. If dilu_on = 1, the aerosol will dilute
%                   with rate defined by parameter dilu_coeff.
% 
% coag_on           Defines whether the coagulation is on or not. If
%                   coag_on = 0, the particles won't coagulate. If
%                   coag_on = 1, the coagulation is set on.
% 
% sedi_on           Defines whether the sedimentation is on or not. If
%                   sedi_on = 0, sedimentation is turned off. If
%                   sedi_on = 1, sedimentation will occur. Only usable for
%                   sedimentation in SAPPHIR chamber!
% 
% coag_mode         Defines whether the particles coagulate 'normally' or
%                   agglomerate. Value should be either 'coag' for normal
%                   coagulation or 'aggl' for agglomeration. Agglomeration
%                   works only for particles in the free-molecule region.
%                   Condensation and deposition might not work correctly
%                   for agglomerates.
% 
% method            Defines which of the sectional methods will be used.
%                   The alternatives are 'moving_sectional',
%                   'moving_center' and 'moving_center_beta'.
% 
% max_timestep      Defines ode's MaxStep value. If max_timestep == 0,
%                   ode's MaxStep option will not be defined, so the step
%                   size will not be restricted.
% 
% Cvap_const        Defines whether the vapor concentration is constant or
%                   not. If Cvap_const == 1, the vapor concetration stays
%                   at value Cvap0 during the whole simulation time, so
%                   that the value gas_source has no effect on vapor
%                   concentration. Otherwise the vapor concentration is not
%                   kept constant, but its value depends on Cvap0 and
%                   gas_source.
% 
% vap_wallsink_on   Defines whether the vapor deposits on walls or not.
% 
% BASIC VALUES:
% gas_source        The condensing vapor source rate (1/cm^3/s). 
%                   Can be defined as a scalar or two column array. When 
%                   defined as scalar, the source rate will be constant 
%                   during the simulation.
%                   If gas_source is an array, it must have two columns;
%                   the first column is a time vector and the second one
%                   tells the gas_source value at respective time.
%                   The time vector doesn't need to have same length as
%                   parameter tvect. However, the first and last element
%                   of gas_source's first column must have the same values
%                   as respective elements of tvect. If the length of the
%                   array is different than tvect's length, gas_source will 
%                   be interpolated to same length.
% 
% part_source       The particle source rate (1/cm^3/s). Must be a
%                   three-colum array. The first column is the time vector
%                   and the second column the particle source rates at
%                   corresponding time points in a similar way as in
%                   gas_source. The first cell of third column defines the
%                   size of particles (in meters). The rest of the cells in
%                   the third column are not used.
% 
% dilu_coeff        Dilution coefficient (1/s). Dilution affects particle
%                   concentration in following way:
%                   dN/dt = -dilu_coeff*N   (N = particle concentration)
%                   dilu_coeff can be either a scalar or an array. When
%                   defined as a scalar, dilution coefficient will be 
%                   constant during the simulation.
%                   If dilu_coeff is an array, it must have two columns;
%                   the first column is a time vector and the second one
%                   tells the dilu_coeff value at respective time.
%                   The time vector doesn't need to have same length as
%                   parameter tvect. However, the first and last element
%                   of dilu_coeff's first column must have the same values
%                   as respective elements of tvect. If the length of the
%                   array is different than tvect's length, dilu_coeff will 
%                   be interpolated to same length.
% 
% vap_wallsink      The flux of vapor molecules that condense on walls
%                   (1/s).
% 
% satu_conc         The condensing vapor saturation concentration (1/cm^3).
% 
% lambda            The condensing vapor mean free path (m).
% 
% diff_coeff        The condensing vapor diffusion coefficient (cm^2/s).
% 
% vap_molmass       Molecular mass of condensing vapor (g/mol).
% 
% particle_dens     Density of particle matter. (g/cm^3).
% 
% stick_coeff       Sticking coefficient. The probability that vapor
%                   molecules will stick to aerosol particles.
% 
% Cvap0             Initial condensing vapor concentration (1/cm^3)
% 
% T                 The temperature (K)
% 
% AGGLOMERATE-RELATED PARAMETERS:
% Df                The fractal dimension of agglomerates. 1.0 < Df < 3.0.
% 
% r0                The radius of agglomerate primary particles.
% 
% TIME VECTOR:
% tvect             The time vector (seconds). Define as row vector.
% 
% DISTRIBUTION PARAMETERS:
% N                 Initial total particle concentration (1/cm^3). This can
%                   be a vector, too. If N is a vector, the distribution 
%                   will be a superposition of several distributions. For
%                   example if N is [1e5 1e6], the final distribution will
%                   consist of two distributions: one contains 1e5
%                   particles and the other 1e6 particles. Each
%                   distribution has respective values of mu and sigma, so
%                   these values have to be vectors of same length as N.
%
% Dp_min            The exponent of minimum diameter of size distribution.
%                   The program will form a lognormal distribution from
%                   10^(Dp_min) to 10^(Dp_max).
% 
% Dp_max            The exponent of maximum diameter of size distribution.
%                   The program will form a lognormal distribution from
%                   10^(Dp_min) to 10^(Dp_max).
%                   
% mu                The mean of the lognormal size distribution. If this is
%                   defined as a vector, the distribution will be a
%                   superposition of several distributions. For example if
%                   mu is [10e-9 50e-9], the final distribution will
%                   consist of two distributions: one has mean diameter of
%                   10 nm and the other 50 nm. If mu is a vector, sigma and 
%                   N has to be vectors of same length.
% 
% sigma             Sigma (standard deviation) of lognormal size 
%                   distribution. Can be defined also as vector, see N and
%                   mu.
% 
% sections          The number of sections in the size distribution.
% 
% output_sections   Defines the number of sections in output size grid. The
%                   dN/dlogDp distribution is interpolated to a denser grid
%                   after the simulation is finished. The number of grid
%                   points is defined by this parameter. The more output
%                   sections, the smoother the plot of size distribution
%                   will be (except when using fixed sectional model).
% 
% Dp                Defines the diameters of sections. The number of
%                   sections will be the same as Dp's length. If this is
%                   used, the parameters Dp_min, Dp_max and sections are
%                   not in use.
% 
% Dplims            Defines the limits between sections. Used only for
%                   fixed sectional model. The length must be one less than
%                   the number of sections because the last section does
%                   not have upper limit and the first section does not
%                   have lower limit. Usually this vector is calculated by
%                   the model, but can be also defined by user.
% 
% number_distr      A vector that defines the particle concentration in
%                   each section. The length of number_distr must equal the
%                   number of sections. If this parameter is used, the
%                   parameters N, mu and sigma are not in use.
% 
% center_diameters  Defines the initial center diameters of sections. Used
%                   only for fixed sectional model. The length must equal
%                   the number of sections. Usually this vector is
%                   calculated by the model, but can be also defined by
%                   user.
% 
% TOLERANCE PARAMETERS:
% Define the tolerance settings for ode45
% 
% Cvap_tol          Vapor concentration tolerance.
% 
% N_tol             Particle concentartion tolerance.
% 
% Dp_tol            Particle diameter tolerance.
%
% *****
% Example:
% 
% obj=chamber;
% obj.initialize('sigma',1.6); % Change the parameter obj.sigma from
%                              % default 1.33 to 1.6.

% (c) Pauli Simonen 2013
%
% Version history:
% 2013-06-04    0.1.0 Separated from function initialize.
% 2013-06-20    0.1.1 Added vap_wallsink_on and vap_wallsink for vapor
%                     condensation onto chamber walls.


% Initialize with default values only if the object has not been
% initialized befor OR if user has not input any parameters.
if(isempty(obj.initials) || nargin == 1)

    obj.initials.dilu_on = 1;             % 1 if dilution is switched on
    obj.initials.coag_on = 1;             % 1 if coagulation is switched on
    obj.initials.sedi_on = 0;             % 1 if sedimentation is on
    obj.initials.coag_mode = 'coag';      % 'coag' if particles coagulate and 'aggl' if they agglomerate
    obj.initials.Cvap_const = 0;          % 1 if vapor concentration is kept constant.
                             % If Cvap_const == 1, the program will set
                             % gas_source to 0.
    obj.initials.vap_wallsink_on = 0;     % 1 if vapor deposits on walls
    obj.initials.method = 'moving_sectional';
    obj.initials.retracking = 0;

    % Basic values
    obj.initials.part_source   = 0;
    obj.initials.gas_source    = 1e5;                  % The condensing vapor source rate. Can be defined as scalar or two-column array.
    obj.initials.dilu_coeff    = 1.00.*370/1460000;    % Dilution coefficient (1/(cm^3s)). Can be defined as scalar or two-column array.
    obj.initials.vap_wallsink  = 0;                    % Vapor deposition velocity to walls (1/(cm^3s)).
    obj.initials.satu_conc     = 0;                    % The condensing vapor saturation concentration (1/cm^3).
    obj.initials.lambda        = 129.29e-9;            % The condensing vapor mean free path (m)
    obj.initials.diff_coeff    = 1.093392158075694e-1; % The condensing vapor diffusion coefficient (cm2/s)
    obj.initials.vap_molmass   = 100;                  % Molecular mass of condensing vapor (g/mol)
    obj.initials.particle_dens = 1.84;                  % Particle density (g/cm3)
    obj.initials.stick_coeff   = 1.0;                  % Sticking coefficient
    obj.initials.Cvap0         = 1e7;                  % Initial condensing vapor concentration
    obj.initials.T             = 290;                  % temperature (K)

    % Time vector
    obj.initials.tvect = 0:60:10800;  % The time vector (seconds)

    % Distribution parameters
    obj.initials.number_distr  = 0;
    obj.initials.Dp            = 0;
    obj.initials.Dplims        = 0;
    obj.initials.center_diameters = 0;
    obj.initials.N             = 1e5;                  % Initial particle concentration
    obj.initials.Dp_min        = -9;     % Exponent of the minimum diameter of size distribution. Diameter will be 10^(Dp_min).
    obj.initials.Dp_max        = -6;     % Exponent of the maximum diameter of size distribution. Diameter will be 10^(Dp_max).
    obj.initials.mu            = 15e-9;  % The mean of lognormal size distribution.
    obj.initials.sigma         = 1.33;   % Sigma of lognormal size distribution.
    obj.initials.sections      = 30;     % The number of size sections in the size distribution.
    obj.initials.output_sections = 300;  % The number of size sections in the final
                            % interpolated distribution. This should be
                            % large enough to reduce interpolating
                            % error. Increase of this number doesn't
                            % slow down the program in the same way as
                            % increase of chamber.sections.

    % Agglomeration-related parameters. Only for free-molecule range.
    obj.initials.Df = 1.8;   % Fractal dimension of agglomerates.
    obj.initials.r0 = 3e-9;  % Radius of agglomerate primary particles.

    % Ode45 tolerance parameters
    obj.initials.Cvap_tol = 1e2; % Vapor concentration tolerance.
    obj.initials.N_tol = 0.01;   % Particle concentration tolerance.
    obj.initials.Dp_tol = 1e-11;% Particle diameter tolerance.
    obj.initials.max_timestep = 0;
end 
% End of default values.


param_names=fieldnames(obj.initials);      % Get the names of parameter fields
total_params=length(param_names);   % Get the number of parameter fields


% If user has set parameters in the input, replace default values with
% user-defined values.


% Check that the input is even (e.g. initialize('gas_source',1e5)):
if(rem((nargin-1),2) ~= 0)
    error('set_initials: The number of arguments must be even.');
end

set_parameters=(nargin-1)/2;    % Get the number of arguments set in the input

% Check for duplicate definitions:
for i=1:2:set_parameters*2
    summary = 0;
    for j=1:2:set_parameters*2
        summary=summary+strcmp(varargin{i},varargin{j});
        if(summary > 1)
            error('set_initials: Parameter %s is defined twice.', varargin{i});
        end
    end
end



for i=1:2:set_parameters*2  % Go through all parameter names in the input
      parameter_name=varargin{i};
      for j=1:total_params % Go through all param fields
          % Compare the input name with the names of param fields:
          if(strcmp(parameter_name,param_names{j}))
              % If the names match, replace the default param value with 
              % input value:
              obj.initials.(parameter_name)=varargin{i+1};
              break;
          elseif(j==total_params)
              % If we are in the end of loop and no matching name is found,
              % such a field doesn't exist.
              warning('set_initials: Invalid argument: %s', parameter_name);              
          end
      end
end


% Set the coag_num. If coag_mode = 'coag' => coag_num = 1.
%                   If coag_mode = 'aggl' => coag_num = 0.
% It is later faster to compare this numerical value than string coag_mode.
if(strcmp(obj.initials.coag_mode,'coag'))
    obj.initials.coag_num = 1;
elseif(strcmp(obj.initials.coag_mode,'aggl'))
    obj.initials.coag_num = 0;
else
    error('chamber.initialize: coag_mode must be either ''coag'' or ''aggl''.');
end
% 
% % Make a logarithmically spaced diameter vector between 10^(Dp_min) and 
% % 10^(Dpmax). Number of cells is params.sections.
% % If params.Dp_vector is not a scalar, the user has already defined the
% % Dp vector, so Dp = params.Dp_vector.
% if(isscalar(obj.initials.Dp))    
%     % Make logarithmically spaced diameter vector:
% %     obj.initials.Dp = logspace(obj.initials.Dp_min, obj.initials.Dp_max, obj.initials.sections);
%     obj.Dps = logspace(obj.initials.Dp_min, obj.initials.Dp_max, obj.initials.sections);
%     obj.sections = obj.initials.sections;
% else
%     % Make sure that all diameter values are unique:
% %     obj.initials.Dp = unique(obj.initials.Dp);
% %     obj.initials.sections = length(obj.initials.Dp);
%     obj.Dps = unique(obj.initials.Dp);
%     obj.sections = length(obj.Dps);
% end
% 
% % Make the particle number distribution:
% if(isscalar(obj.initials.number_distr))
%    % If initials.number_distr is scalar, user has not defined it. In that
%    % case, form a lognormal distribution dN/dlogDp where the minimum 
%    % diameter is vector Dp's first value and maximum is Dp's last value.
%    % Number of cells is the same as Dp's length (i.e. params.sections).
%    dNdlogDp = zeros(size(obj.Dps));
%    [rows, cols] = 
%     for i=1:length(obj.initials.mu)
%         dNdlogDp=dNdlogDp + obj.log_normal(obj.Dps,obj.initials.mu(i),obj.initials.sigma(i),obj.initials.N(i));
%     end
%     
%     % Get the concentration of particles (N) in each cell (=section) of the
%     % distribution (Dp,dNdlogDp). N(i) is the concentration of particles in
%     % section i.
%     [~, obj.number_distribution]  = obj.Dlog_to_N_vect(obj.Dps,dNdlogDp);
%     clear dNdlogDp;
% 
%     obj.number_distribution = abs(obj.number_distribution); % Make sure that concentrations are positive 
%                 % (Dlog_to_N_vect may make near-zero concentrations negative).
%                 
% else
%     obj.number_distribution = obj.initials.number_distr;
% end
% 
% 
% % If fixed sectional model is used, find the limits for each section and
% % save them to initials. If initials.Dplims is a vector, don't set the
% % limits because they are set by user.
% if(isscalar(obj.initials.Dplims))
%     Dp = obj.Dps;
%     
%     % The length of Dplims is (length(Dp)-1) because the first section
%     % does not have lower limit and the last section does not have upper
%     % limit.
%     obj.Dplims=zeros(1,length(Dp)-1);
%     for i = 1:length(obj.Dplims)
%         exp1 = log10(Dp(i))/log10(10);  % Find the exponents of current and
%         exp2 = log10(Dp(i+1))/log10(10);% next Dp.
%         % Then create a logarithmically spaced 3-element vector between these Dp:s.
%         temp = logspace(exp1,exp2,3);
% 
%         % Now the second value of the vector is the logarithmic center between
%         % Dp:s and will be the upper limit of section i.
%         obj.Dplims(i)=temp(2);
%     end
%     clear temp exp1 exp2 Dp;
% else
%     obj.Dplims = obj.initials.Dplims;
% end
% 
% if(isscalar(obj.initials.center_diameters))
%     obj.center_diameters = obj.Dps;
% else
%     obj.center_diameters = obj.initials.center_diameters;
% end


% Set the vector switches on if necessary
if(isscalar(obj.initials.gas_source)==0)
    obj.initials.gas_source_is_vect = 1;
else obj.initials.gas_source_is_vect = 0;
end

if(isscalar(obj.initials.dilu_coeff)==0)
    obj.initials.dilu_vect_on  = 1;
else obj.initials.dilu_vect_on  = 0;
end

if(isscalar(obj.initials.part_source)==0)
    obj.initials.part_source_is_vect = 1;
else obj.initials.part_source_is_vect = 0;
end

% Set gas_source to 0 if Cvap_const == 1:
if(obj.initials.Cvap_const == 1)
    obj.initials.gas_source = 0;
end


end

