function [ out ] = create_gas_source(tvect, source_length, vap_tot, t_start)
%GAS_SOURCE Creates a vapor source vector
%   Creates a sigmoid-form vapor source. The total amount of vapor will be
%   vap_tot and the duration of source will be source_length. 
%   This source is created for every 24 hours, so the total amount will be
%   number_of_days * vap_tot. t_start defines the beginning time of source.
%   If the vapor source begins for example at 12:00, t_start should be
%   12*3600 = 43200. tvect should equal the time vector used in the
%   simulation. The smoothness of sigmoid function can be adjusted by
%   changing the variable 'smoothness'.
%
%   With a small modification, this function can be used also so that the
%   maximum vapor source rate will be defined instead of source duration.
%   See function create_part_source to see how.

smoothness = 0.003;
num_of_days = ceil(tvect(end)/(3600*24));

vap_max = vap_tot/source_length;

out = [tvect', zeros(length(tvect),1)];

for i=1:num_of_days
    out(:,2) = out(:,2) + sigmoid(tvect, t_start, t_start + source_length, vap_max, smoothness)';
    t_start = 24*3600 + t_start;
end

end


function [out] = sigmoid(x, x1, x2, max_val, coeff)

out = max_val./(1+exp(-coeff*(x-x1)))-max_val./(1+exp(-coeff*(x-x2)));

tolerance = 1;
out(out<tolerance) = 0;
end