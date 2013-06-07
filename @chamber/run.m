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

if(obj.initials.fixed_sections == 0)
    obj.run_movsec;
else
    obj.run_moving_center;
end

end