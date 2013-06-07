% (c) Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-06-04    0.1.1 Separated function initialize to functions set_params
%                     and check_initials.


function initialize(obj,varargin)
% Sets the initial parameters of the chamber object (obj).
%
% obj.set_params('field_name_1', value_1, ...
%               'field_name_2', value_2, ...
%               'field_name_n', value_n)
% sets the obj's initial parameters. Fields with names field_name_1, 
% field_name_2 and field_name_n are replaced by values value_1, value_2 
% and value_n, respectively. Other fields will not be changed.
%
%
% THE PARAMETERS OF THE CHAMBER MODEL
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
% gas_source_vect_on  Defines whether the parameter gas_source is an array
%                     or scalar. This cannot be set by user. Instead, the
%                     program checks if gas_source is an array or not and
%                     sets the value to 1 or 0 respectively.
% 
% dilu_vect_on      Defines whether the parameter dilu_coeff is an array or
%                   scalar. This cannot be set by user. Instead, the
%                   program checks if dilu_coeff is an array or scalar and
%                   sets the value to 1 or 0 respectively.
% 
% coag_mode         Defines whether the particles coagulate 'normally' or
%                   agglomerate. Value should be either 'coag' for normal
%                   coagulation or 'aggl' for agglomeration. Agglomeration
%                   works only for particles in the free-molecule region.
%                   Condensation and deposition might not work correctly
%                   for agglomerates.
%
% coag_num          Numerical representative of coag_mode. 
%                   If coag_mode == 'coag' => coag_num = 1.
%                   If coag_mode == 'aggl' => coag_num = 0.
%                   This cannot be set directly by user, but the program
%                   sets it based on the value of coag_mode.
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
% *****
% Example:
% 
% obj=chamber;
% obj.initialize('sigma',1.6); % Change the parameter obj.sigma from
%                              % default 1.33 to 1.6.


obj.set_params(varargin{:});
obj.check_initials;



end