% (c) Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-05-31    0.1.1 Added method run_moving_center.
% 2013-06-04    0.1.2 Splitted method initialize to methods set_params and
%                     check_initials.

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
% chamb.initialize('Ntot0',1e6,'T',300);
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

    % Properties that user can see, but not modify:
    properties (GetAccess = public, SetAccess = private)
        initials; % Contains all initial parameters. Defined in initialize.m.
                
        %CALCULATED DATA:
        output_data  % Calculated data in a handy structure.
    end
    
    % Public methods:
    methods
        % Creates the object.
        function obj = chamber(varargin)
            obj.initialize;
        end
        
        % Initializes object with user input values:
        initialize(obj,varargin) % Defined in initialize.m
        
        set_params(obj,varargin)
        
        check_initials(obj)
        
        % Runs the chamber simulation with initialized values:
        run(obj) % Defined in run.m
                
        % Plots some data:
        plot(obj,varargin) % Defined in plot.m       

    end
    
    % Private methods:
    methods (Access = private)       
        % Converts the calculated data to nice form
        out_struct = model_convert(obj, t, Y) % Defined in model_convert.m
               
        % Plots the distribution, used in public method chamber.plot.
        subplot_dmps(obj,sub);
        
        % Runs the simulation with moving sections.
        run_movsec(obj)
       
        % Runs the simulation with fixed sections, but the diameter inside
        % these sections moves.
        run_moving_center(obj)
        

    end
    
    % Private static methods:
    methods (Static, Access = private)
        % Creates a lognormal distribution over the diameter vector Dp_in
        [out] = log_normal(Dp_in,mu,sigma,N)
        
        % Gets the concentration of particles in each section of
        % distribution [Dp, dN]
        [Dp, N] = Dlog_to_N_vect(Dp, dN)
        
        % Makes the coagulation kernel
        [K] = koag_kernel(Dp1,Dp2,dens,T)
        
        % Makes the agglomeration kernel (Free-molecule range)
        [K] = aggl_kernel(Dp1,Dp2,dens,T,Df,r0)
        
        % Makes the coagulation matrix
        [out] = coagulationMatrix(Dp,ind);
        
        [out] = sapphir_beta2(Dp,T)
        
        [out] = distribution_info_Vtot(Dp,dN);
        
        [Diffcoeff] = diff_particle(Dp,T)
    end
    
    % Static methods:
    methods (Static)
        function hymyile()
            fprintf(':-)\n');
        end
    end
end