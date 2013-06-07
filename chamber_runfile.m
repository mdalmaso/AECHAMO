% (c) Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-06-04    0.1.1 Now all parameters are first initialized and then
%                     checked all at once.

% TODO:
% -Muistin säästäminen: aja yksi ajo kerrallaan, tallenna tulokset
%  tiedostoon, puhdista workspace ja aja seuraava ajo.
% 


function [chamb,elapsed] = chamber_runfile(filename)
% function chamber_runfile(filename)
% 
% Reads file 'filename' and calculates the number of runs (n) defined in
% the file. Returns n chamber objects as arrays chamb(1), chamb(2) ... 
% and runs the simulations.
%
% The settings file must be in following form:
% The definition of first chamber object begins with character # and ends
% with the same character #. The code between these characters is run and
% the variable names in this code are saved as well as their values. These
% parameters will be forwarded to chamber.initialize to apply the values to
% chamber objects. When the program finds the next character #, it knows
% that the definition of first chamber ends and the next begins.
%

% Read the settings from a file:
settings = read_file(filename); % Function defined in the end of this file.

% Get the amount of chamber objects defined in the file:
[rows, cols] = size(settings);

% Create the needed amount of chamber objects:
chamb(rows) = chamber;

% Initialize first all chamber objects to see if there are any problems
% with initializing before running the simulations
vars = '';
for i=1:rows
    for j=1:cols
        if(~isempty(settings(i,j).param_name))
            chamb(i).set_params(settings(i,j).param_name,settings(i,j).param_value);
        end
    end
    chamb(i).check_initials;
end


% Run the simulations now as we are sure that all the initializations went
% OK.
for i=1:rows
    tic
    chamb(i).run;
    elapsed = toc;
    save(datestr(now,30),'chamb','elapsed');
    delete(chamb(i));
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