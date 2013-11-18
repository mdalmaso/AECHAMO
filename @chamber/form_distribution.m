function form_distribution(obj)
%FORM_DISTRIBUTION Defines vectors Dp, Dplims, center_diameters and N
%   Use the user set values for diameter vector, limits, center diameters
%   and number distribution, or calculate them if they are not defined.


% Make a logarithmically spaced diameter vector between 10^(Dp_min) and 
% 10^(Dpmax). Number of cells is params.sections.
% If params.Dp_vector is not a scalar, the user has already defined the
% Dp vector, so Dp = params.Dp_vector.
if(isscalar(obj.initials.Dp))    
    % Make logarithmically spaced diameter vector:
    obj.Dps = logspace(obj.initials.Dp_min, obj.initials.Dp_max, obj.initials.sections);
    obj.sections = obj.initials.sections;
else
    % User has defined the diameter vector, so use it, but make sure that
    % all diameter values are unique:
    obj.Dps = unique(obj.initials.Dp);
    obj.sections = length(obj.Dps);
end

% Make the particle number distribution:
if(isscalar(obj.initials.number_distr))
   % If initials.number_distr is scalar, user has not defined it. In that
   % case, form a lognormal distribution dN/dlogDp where the minimum 
   % diameter is vector Dp's first value and maximum is Dp's last value.
   % Number of cells is the same as Dp's length (i.e. params.sections).
   dNdlogDp = zeros(size(obj.Dps));
    for i=1:length(obj.initials.mu)
        dNdlogDp=dNdlogDp + obj.log_normal(obj.Dps,obj.initials.mu(i),obj.initials.sigma(i),obj.initials.N(i));
    end
    
    % Get the concentration of particles (N) in each cell (=section) of the
    % distribution (Dp,dNdlogDp). N(i) is the concentration of particles in
    % section i.
    [~, obj.number_distribution]  = obj.Dlog_to_N_vect(obj.Dps,dNdlogDp);
    clear dNdlogDp;
    % Make sure that concentrations are positive
    % (Dlog_to_N_vect may make near-zero concentrations negative).
    obj.number_distribution = abs(obj.number_distribution);               
else
    % User has defined the distribution, so use it.
    obj.number_distribution = obj.initials.number_distr;
end


% Find the limits for each section for fixed sectional model and
% save them. If initials.Dplims is a vector, don't set the
% limits because they are set by user.
if(isscalar(obj.initials.Dplims))
    Dp = obj.Dps;
    % The length of Dplims is (length(Dp)-1) because the first section
    % does not have lower limit and the last section does not have upper
    % limit.
    obj.Dplims=zeros(1,length(Dp)-1);
    for i = 1:length(obj.Dplims)
        % Now the second value of the vector is the geometric mean, that
        % is, logarithmic center between Dp:s and will be the upper limit
        % of section i and lower limit of section (i+1).
        obj.Dplims(i)=geomean([Dp(i),Dp(i+1)]);
    end
    clear temp exp1 exp2 Dp;
else
    % User has set the limits, so use them:
    obj.Dplims = obj.initials.Dplims;
end

if(isscalar(obj.initials.center_diameters))
    % If user has not set the center diameters, use the Dp values as
    % initial center diameters:
    obj.center_diameters = obj.Dps;
else
    % Otherwise, use the user set center diameters:
    obj.center_diameters = obj.initials.center_diameters;
end


end

