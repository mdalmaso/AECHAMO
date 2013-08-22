function run(obj)
% RUN Runs the simulation based on initial parameters and saves the results
% to output_data.
% 
% Example:
% chamber.run;  % Runs the simulation. Results are saved in
%               % chamber.output_data

% (c) Pauli Simonen 2013
% Version history:
% 2013-05-31    0.1.0
tic
if(obj.initials.fixed_sections == 0)
    [t, Y] = obj.run_movsec;
else
    [t, Y] = obj.run_moving_center;
end

display('Ode45 finished, processing data...');

% This makes a handy structure of the results:
obj.output_data = obj.model_convert(t,Y);
toc
end