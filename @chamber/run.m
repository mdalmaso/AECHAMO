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

switch obj.initials.method
    case 'moving_sectional'
        [t,Y] = obj.run_movsec;
        display('Ode45 finished, processing data...');
        obj.output_data = obj.model_convert(t,Y);
    case 'moving_center'
        [t,Y] = obj.run_moving_center;
        display('Ode45 finished, processing data...');
        obj.output_data = obj.model_convert(t,Y);
    case 'moving_center_beta'
        [t,Y] = obj.run_movcent2;
        display('Ode45 finished, processing data...');
        obj.output_data = obj.model_convert2(t,Y);
end

toc

end