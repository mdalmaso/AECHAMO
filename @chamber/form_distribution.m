function form_distribution(obj)
%FORM_DISTRIBUTION Defines vectors Dp, Dplims, center_diameters and N
%   Detailed explanation goes here


% Make a logarithmically spaced diameter vector between 10^(Dp_min) and 
% 10^(Dpmax). Number of cells is params.sections.
% If params.Dp_vector is not a scalar, the user has already defined the
% Dp vector, so Dp = params.Dp_vector.
if(isscalar(obj.initials.Dp))    
    % Make logarithmically spaced diameter vector:
%     obj.initials.Dp = logspace(obj.initials.Dp_min, obj.initials.Dp_max, obj.initials.sections);
    obj.Dps = logspace(obj.initials.Dp_min, obj.initials.Dp_max, obj.initials.sections);
    obj.sections = obj.initials.sections;
else
    % Make sure that all diameter values are unique:
%     obj.initials.Dp = unique(obj.initials.Dp);
%     obj.initials.sections = length(obj.initials.Dp);
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

    obj.number_distribution = abs(obj.number_distribution); % Make sure that concentrations are positive 
                % (Dlog_to_N_vect may make near-zero concentrations negative).
                
else
    obj.number_distribution = obj.initials.number_distr;
end


% If fixed sectional model is used, find the limits for each section and
% save them to initials. If initials.Dplims is a vector, don't set the
% limits because they are set by user.
if(isscalar(obj.initials.Dplims))
    Dp = obj.Dps;
    
    % The length of Dplims is (length(Dp)-1) because the first section
    % does not have lower limit and the last section does not have upper
    % limit.
    obj.Dplims=zeros(1,length(Dp)-1);
    for i = 1:length(obj.Dplims)
        exp1 = log10(Dp(i))/log10(10);  % Find the exponents of current and
        exp2 = log10(Dp(i+1))/log10(10);% next Dp.
        % Then create a logarithmically spaced 3-element vector between these Dp:s.
        temp = logspace(exp1,exp2,3);

        % Now the second value of the vector is the logarithmic center between
        % Dp:s and will be the upper limit of section i.
        obj.Dplims(i)=temp(2);
    end
    clear temp exp1 exp2 Dp;
else
    obj.Dplims = obj.initials.Dplims;
end

if(isscalar(obj.initials.center_diameters))
    obj.center_diameters = obj.Dps;
else
    obj.center_diameters = obj.initials.center_diameters;
end

end

