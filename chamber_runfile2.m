function [chamb,elapsed] = chamber_runfile2(filename)
% CHAMBER_RUNFILE Runs the simulations initialized in file.
% 
% [chamb, elapsed] = chamber_runfile(filename).
% 
% Reads file 'filename' and calculates the number of runs (n) defined in
% the file. Creates n chamber objects as arrays chamb(1), chamb(2) ... 
% and runs the simulations.
% 
% Returns the chamber objects and the time needed to run the simulation.
% In addition, the results are saved to file 'run_timestamp.mat'.
% If the function is interrupted before the final file is saved, the ready
% simulations are saved in files 'temp_timestamp.mat'.
% 
% The settings file must be in following form:
% The definitions of parameters begin with character # and end with the 
% same character #. The code between these characters is run and
% the variable names in this code are saved as well as their values. These
% parameters will be forwarded to chamber.initialize to apply the values to
% chamber objects. 

% A parameter can have multiple definitions, and these
% definitions are usually separated by defining the parameter as a row
% vector. Example:
% coag_on = [1; 0];
% 
% If there are several multi-definitions, all their combinations are
% calculated.
% 
% Multiple definitions of gas_source and part_source are not separated by
% adding rows. Instead, these parameters can be defined as 3d-matrices.
% Example: two definitions for gas_source
% gas_source(:,:,1) = [tvect', 1.4*tvect'];
% gas_source(:,:,2) = [tvect', 0.7*tvect'];
% 
% Example of settings file:
% % Defines some constants and 2 definitions for coag_on, mu and
% % gas_source. Total number of chamber objects will be 2*2*2=8.
% 
% #
% dilu_on = 0;
% sections = 30;
% output_sections = 10*sections;
% Cvap0 = 2e7;
% tvect = 0:60:32400;
% N = [1e3 1e2];
% 
% coag_on = [1; 0];
% mu=[50e-9, 140e-9; 100e-9, 200e-9];
% gas_source(:,:,1) = [tvect', 1.337.*tvect'];
% gas_source(:,:,2) = [tvect', 3.14159.*tvect'];
% #

% (c) Pauli Simonen 2013
%
% Version history:
% 2013-06-19    0.1.0


% Read the settings from a file:
settings = read_file(filename); % Function defined in the end of this file.

% Get the amount of chamber objects defined in the file:
[rows, cols] = size(settings);

if(rows > 1)
    error('Multiple definitions in %s. Define multiple values for parameters as a vector or use function chamber_runfile.',filename);
end

%% Remove duplicate definitions, save the last definition.
% Example 1: 
% gas_source = 1e5;
% gas_source = 2e5;
% % The final value of gas_source will be 2e5. 

% Example 2: 
% 
% gas_source(:,:,1) = 1e5;
% gas_source(:,:,2) = 2e5;
% % The final value of gas_source will be (:,:,1)=1e5 and
% % (:,:,2)=2e5, because the last definition does not replace the first
% % definition.
index = 1;
for i=cols:-1:1
    sum = 0;
    for j=cols:-1:1
        sum=sum+strcmp(settings(i).param_name,settings(j).param_name);
        if(sum == 1 && i == j)
            settings_temp(index).param_name=settings(j).param_name;
            settings_temp(index).param_value=settings(j).param_value;
            index = index + 1;
            break;
        elseif(sum == 1 && i < j)
            warning('Duplicate definition of parameter ''%s''.',settings(i).param_name);
            break;
        end
    end
end
settings=settings_temp;
[rows, cols] = size(settings);
clear index;
clear settings_temp;
%%


j=1;
k=1;
total_num = 1;
variables = [];
constants = [];
ranges = [];
for i=1:cols
    [ro co dep] = size(settings(i).param_value);
    if((strcmp(settings(i).param_name,'gas_source') ~= 1) && (strcmp(settings(i).param_name,'part_source') ~= 1) && (strcmp(settings(i).param_name,'dilu_coeff') ~= 1))
        if(ro > 1)
            variables(j).param_name=settings(i).param_name;
            % Param_value is a vector of multiple rows:
            variables(j).param_value=settings(i).param_value;
            total_num = total_num * ro;
            ranges(j) = ro;
            j=j+1;
        else
            constants(k).param_name=settings(i).param_name;
            constants(k).param_value = settings(i).param_value;
            k = k+1;
        end
    elseif(dep > 1)
            variables(j).param_name=settings(i).param_name;
            % Param value is a 3d-matrix:
            variables(j).param_value=settings(i).param_value;
            total_num = total_num * dep;
            ranges(j) = dep;
            j=j+1;
    else
            constants(k).param_name=settings(i).param_name;
            constants(k).param_value = settings(i).param_value;
            k=k+1;
    end
end
clear j;
clear k;

chamb_temp(total_num) = chamber;


[chamb_temp, ~] = recfor([],ranges,1,1,constants, variables,chamb_temp);



% Run the simulations now as we are sure that all the initializations went
% OK.
s = warning('error', 'MATLAB:ode45:IntegrationTolNotMet');
for i=1:total_num
    tic
    try 
        chamb_temp(i).run;
    catch
        err_information=lasterror;
        chamb_temp(i).error_messages = err_information.message;
        clear err_information;
    end
    elapsed_temp(i) = toc;
    out_filename_temp{i} = strcat('temp_', datestr(now,30),'.mat');
    save(out_filename_temp{i},'chamb_temp','elapsed_temp');
    delete(chamb_temp(i));
end



chamb(total_num) = chamber;

out_filename_final = strcat('run_', datestr(now,30));

for i=1:total_num
    load(out_filename_temp{i});
    chamb(i) = chamb_temp(i).copy;
    delete(chamb_temp(i));
    elapsed(i) = elapsed_temp(i);
    save(out_filename_final,'chamb','elapsed');
    delete(out_filename_temp{i});
end
delete(chamb_temp);

end

% Recursive for-loop
function [chamb, obj_index] = recfor(indices,ranges,n,obj_index,constants, variables,chamb)
if(n <= length(ranges))
    for i=1:ranges(n)
        indices(end+1)=i;
        [chamb, obj_index]=recfor(indices,ranges,n+1,obj_index,constants,variables,chamb);
        indices(end) = [];
    end
else
    for i=1:length(constants)
        chamb(obj_index).set_params(constants(i).param_name,constants(i).param_value);
    end
    for i=1:length(indices)
        if((strcmp(variables(i).param_name,'gas_source') || strcmp(variables(i).param_name,'part_source') || (strcmp(variables(i).param_name,'dilu_coeff'))))
            chamb(obj_index).set_params(variables(i).param_name,variables(i).param_value(:,:,indices(i)));
        else
            chamb(obj_index).set_params(variables(i).param_name,variables(i).param_value(indices(i),:));
        end
    end
    chamb(obj_index).check_initials;
    obj_index = obj_index + 1;
end
end

function [settings] = read_file(filename)
% function [settings] = read_file(filename)

% Reads all variables between characters # and # and stores their names and
% values to a matrix. The next values between # and # will be added to the
% matrix as a next row. Function is created to read the settings for
% chamber objects.
%
% The settings file should be in following form:
% The definition of first chamber object begins with character # and ends
% with the same character #. The code between these characters is run and
% the variable names in this code are saved as well as their values. These
% parameters will be forwarded to chamber.initialize to apply the values to
% chamber objects. When the program finds the next character #, it knows
% that the definition of first chamber ends and the next begins.


% Open the file:
file = fopen(filename);

% Read the first line:
line = fgetl(file);

matches = [];

% Search for character # to know where the definitions begin.
while (ischar(line) && isempty(matches) == 1)
    % Check if there is a # character
    matches = strfind(line,'#');
    
    % And read the next line
    line = fgetl(file);
end


i=1;
j=1;

% Loop until the end of the file is reached:
while ischar(line)
    % Look for character #:
    matches = strfind(line,'#');
    
    % Remove the leading and trailing whitespace from line:
    line = strtrim(line);
    
    % If # was not found, read the parameter name from that line:
    if(isempty(matches) == 1)
        % If the line is not empty or comment line,
        % read the parameter name:
        if(~isempty(line)  && line(1) ~= '%')            
            % Parameter name is on the left side of character '='.
            param_names{j} = strtok(line,'=');
            % If there is for example a vector definition, the parameter
            % name is before parenthesis:
            param_names{j} = strtok(param_names{j},'(');
            
            % Remove all whitespaces from param_name:
            param_names{j}(param_names{j} == ' ') = '';
            
            % Run the line with eval, so the variable will be saved in the
            % workspace.
            eval(line);
            
            % Begin the next column:
            j = j+1;
        end
    % If # was found and j > 0, store all parameter names 
    % and their values to settings array. If j = 0, there are no
    % parameters.
    elseif(j > 0)
        for k=1:j-1
            % Store first the param name
            settings(i,k).param_name=param_names{k};
            
            % And then the value. Because the rows were run with eval, all
            % variables are in workspace. So writing a parameter name
            % Matlab returns its value.
            settings(i,k).param_value=eval(param_names{k});
        end
        % Begin the next row in settings array:
        i = i+1;
        % And reset the column index to one:
        j = 1;
        
        clear param_names;
    end
    
    % Read the next line:
    line = fgetl(file);
end

% TODO: 
%       -Better error messaging.


% Close the file:
fclose(file);

end