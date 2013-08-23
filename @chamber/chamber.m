classdef chamber < handle

% Class chamber 
% 
% A class for running the chamber model. The class contains the needed
% initial parameters, such as chamber.Ntot0 (initial particle
% concentration) and chamber.tvect (the time vector). The calculated data
% will be in chamber.output_data.
% 
% The class has also functions related to chamber model.
% 
% Function chamber.initialize(varargin) defines the default initial values.
% User can change default values through input. For example
% chamber.initialize('Ntot0',1e6) sets the value chamber.Ntot0 = 1e6. Other
% initial values will be default. For more details, type 'help
% chamber.initialize'.
% 
% Function chamber.run() runs the chamber model. Before using this
% function, the initial values must be set by chamber.initialize().
% Function chamber.run calculates the data and saves it to
% chamber.output_data. 
% 
% Function chamber.plot(varargin) plots some calculated data. Note that
% you must execute function chamber.run first to get the data.
% chamber.plot() plots several graphs whereas chamber.plot('dist') plots
% only the particle distribution as a function of time.
% 
% Function values = chamber.export_initials outputs a structure that 
% contains all the initial values used in calculations.
% 
% *****
% 
% Example:
% 
% chamb=chamber %Create a chamber object named chamb
% 
% % Initialize chamb with default values, except for chamb.Ntot and chamb.T
% % that are defined as chamb.Ntot = 1e6 and chamb.T = 300.
% chamb.initialize('N',1e6,'T',300);
% 
% chamb.run; %Run the simulation.
% 
% chamb.plot; %Plot the results.
% 
% % Get the calculated data and save it to dat.
% dat=chamb.output_data;
% 
% % Plot the total particle concentration as a function of time:
% plot(dat.tim, dat.Ntot);
% 

% (c) Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-05-31    0.1.1 Added method run_moving_center.
% 2013-06-04    0.1.2 Splitted method initialize to methods set_params and
%                     check_initials.

    properties (Access = public)
        error_messages  % Contains warnings and errors during simulation.
        %CALCULATED DATA:
        output_data  % Calculated data in a handy structure.
    end
    % Properties that user can see, but not modify:
    properties (GetAccess = public, SetAccess = private)
        initials; % Contains all initial parameters. Defined in initialize.m.
    end
    
    properties (Access = private)
        Dps;
        Dplims;
        center_diameters;
        number_distribution;
    end
    % Public methods:
    methods
        % Creates the object.
        function obj = chamber(varargin)
            obj.initialize;
        end
        
        % Runs the simulation with fixed sections, but the diameter inside
        % these sections moves.
        [out_t, out_Y] = run_moving_center(obj)       
        
        [out_t, out_Y] = run_retracking(obj)   
        
        % Initializes object with user input values:
        initialize(obj,varargin)
        
        % Sets the initial parameters. Function initialize uses this
        % function.
        set_params(obj,varargin)
        
        % Checks that the initial values are correct. Function initialize
        % uses this function.
        check_initials(obj)
        
        % Runs the chamber simulation with initialized values:
        run(obj)
        
        % Converts the calculated data to nice form
        out_struct = model_convert(obj, t, Y)
        
        % Plots some data:
        plot(obj,varargin)
        
        % Makes a mass conservation check:
        mass_conserv_check(chamber);
        
        % Makes a copy of a chamber object:
        function [new_object] = copy(object)
            %Copies all the data of a chamber object and makes a new 
            %chamber object using the copied values.
            
            % Make a new object of the same class.
            new_object = feval(class(object));
 
            % Copy all properties.
            p = fieldnames(object);
            for i = 1:length(p)
                new_object.(p{i}) = object.(p{i});
            end
        end

    end
    
    % Private methods:
    methods (Access = private)       
        % Plots the distribution, used in public method chamber.plot.
        subplot_dmps(obj,sub,varargin);
        
        % Runs the simulation with moving sections.
        [t, Y] = run_movsec(obj)
       
        % Converts a N-Dp distribution to dN/dlogDp distribution.
        [dN] = N_to_dlog(obj, Dp,N);
        
        
%         [out] = distribution_info_Vtot(obj,Dp,dN);
        
        % Adds nucleated particles to dy:
        [dy] = add_nucleation(obj,dy,y,t,part_source);
        
        % Makes the agglomeration kernel
        [K] = aggl_kernel(obj, Dp1,Dp2,dens,T,Df,r0)
        
    end
    
    % Private static methods:
    methods (Static, Access = private)
        % Creates a lognormal distribution over the diameter vector Dp_in
        [out] = log_normal(Dp_in,mu,sigma,N)
        
        % Gets the concentration of particles in each section of
        % dN/dlogDp distribution
        [Dp, N] = Dlog_to_N_vect(Dp, dN)
        
        % Makes the coagulation kernel
        [K] = koag_kernel(Dp1,Dp2,dens,T)
        
        % Makes the coagulation matrix
        [out] = coagulationMatrix(Dp,ind);
        
        % Calculates particle deposition into chamber walls:
        [out] = sapphir_beta2(Dp,T);
        
        % Calculates the coagulation sink for particles
        [Diffcoeff] = diff_particle(Dp,T);
        
        
        [out] = integrate_distribution(Dp,dN,dmin,dmax);
        
        % Adds the effect of condensation to dy:
        [dy] = add_condensation(dy, y, initials, index);
    end
    
    % Static methods:
    methods (Static)
        function hymyile()
            fprintf(':-)\n');
        end
    end
end