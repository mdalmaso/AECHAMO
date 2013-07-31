function [] = check_initials(obj)
%CHECK_INITIALS Checks the initial parameters and reports errors.

% (c) Pauli Simonen 2013
%
% Version history:
% 2013-06-04    0.1.0 Separated from function initialize.


initials = obj.initials;

% Check that coag_mode is either 'coag' or 'aggl'.
if(strcmp(initials.coag_mode,'coag') == 0)
    if(strcmp(initials.coag_mode,'aggl') == 0)
        error('set_initials: Argument coag_mode must be either ''coag'' or ''aggl''.');
    end
end

% Check that Dp_min and Dp_max are exponents:
if(abs(initials.Dp_min) < 1) 
    error('set_initials: Argument ''Dp_min'' must be an exponent. The minimum diameter will be 10^(Dp_min).');
end

if(abs(initials.Dp_max) < 1)
    error('set_initials: Argument ''Dp_max'' must be an exponent. The maximum diameter will be 10^(Dp_max).');
end

% Then check that Dp_max > Dp_min:
if(initials.Dp_min >= initials.Dp_max)
    error('set_initials: ''Dp_min'' must be smaller than ''Dp_max''.');
end

% Make sure that sigma(s) of distribution is > 1.0:
for i=1:length(initials.sigma)
    if(initials.sigma(i) <= 1.0)
        error('set_initials: ''sigma'' must be bigger than 1.0.');
    end

% Check that the time vector is a vector.
if(length(initials.tvect) < 2)
    error('set_initials: Argument ''tvect'' must be a vector.');
end

% If 'dilu_coeff' has been defined as an array, check that it has two
% columns and the same last and first point as obj.tvect.
if(isscalar(initials.dilu_coeff) == 0)  % Is it an array?
    % Check that there are two columns in dilu_coeff:
    [rows cols]=size(initials.dilu_coeff);
    if(cols ~=2)
        error('set_initials: The argument ''dilu_coeff'' must consist of two columns.');
    end
    % Check that the first element of dilu_coeff(:,1) is the same as the
    % first element of tvect:
    if(initials.dilu_coeff(1,1) ~= initials.tvect(1))
        error('set_initials: The first column of argument ''dilu_coeff'' must have the same first value as ''tvect''.');
    end
    % Check that the last element of dilu_coeff(:,1) is the same as the last
    % element of tvect:
    if(initials.dilu_coeff(rows,1) ~= initials.tvect(length(initials.tvect)))
        error('set_initials: The first column of argument ''dilu_coeff'' must have the same last value as ''tvect''.');
    end
end

% The same check for gas_source if it is a vector.
if(isscalar(initials.gas_source) == 0)  % Is it an array?
    % Check that there are two columns in gas_source:
    [rows cols]=size(initials.gas_source);
    if(cols ~=2)
        error('set_initials: The argument ''gas_source'' must consist of two columns.');
    end
    % Check that the first element of gas_source(:,1) is the same as the
    % first element of tvect:
    if(initials.gas_source(1,1) ~= initials.tvect(1))
        error('set_initials: The first column of argument ''gas_source'' must have the same first value as ''tvect''.');
    end
    % Check that the last element of gas_source(:,1) is the same as the last
    % element of tvect:
    if(initials.gas_source(rows,1) ~= initials.tvect(length(initials.tvect)))
        error('set_initials: The first column of argument ''gas_source'' must have the same last value as ''tvect''.');
    end
end

% Check the particle source as well.
if(isscalar(initials.part_source) == 0)  % Is it an array?
    % Check that there are three columns in part_source:
    [rows, cols, depth]=size(initials.part_source);
    if(cols ~=3)
        error('set_initials: The argument ''part_source'' must consist of three columns.');
    end
    % Check that the first element of part_source(:,1) is the same as the
    % first element of tvect:
    if(initials.part_source(1,1) ~= initials.tvect(1))
        error('set_initials: The first column of argument ''part_source'' must have the same first value as ''tvect''.');
    end
    % Check that the last element of part_source(:,1) is the same as the last
    % element of tvect:
    if(initials.part_source(rows,1) ~= initials.tvect(length(initials.tvect)))
        error('set_initials: The first column of argument ''part_source'' must have the same last value as ''tvect''.');
    end
end


% Check the value of Df. Should be 1.0 < Df < 3.0.
if(initials.Df < 1.0)
    error('chamber.initialize: Parameter Df must be larger than 1.0.');
elseif(initials.Df > 3.0)
    error('chamber.initialize: Parameter Df must be smaller than 3.0.');
end

% If mu, N and sigma are vectors, make sure that all these have same
% length, because the distribution is then sum 
% D(mu(1),sigma(1),N(1)) + ... + D(mu(n),sigma(n),N(n)), where D denotes
% lognormal distribution.
num_of_distr = length(initials.mu);
if(length(initials.N) ~= num_of_distr || length(initials.sigma) ~= num_of_distr)
    error('chamber.initialize: Parameters ''mu'', ''N'' and ''sigma'' must have equal lengths.');
end

end

