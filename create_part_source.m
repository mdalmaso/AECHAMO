function [ out ] = create_part_source(tvect, rate_max, N_tot, t_start, particle_size)
%CREATE_PART_SOURCE Creates particle source vector
%   Creates a sigmoid-form particle source. The total amount of particles
%   will be N_tot and the maximum particle source rate will be rate_max.
%   This source is created for every 24 hours, so the total amount will be
%   number_of_days * N_tot. t_start defines the beginning time of source.
%   If the nucleation happens for example at 12:00, t_start should be
%   12*3600 = 43200. tvect should equal the time vector used in the
%   simulation. The smoothness of sigmoid function can be adjusted by
%   changing the variable 'smoothness'.
%
%   With a small modification, this function can be used also so that the
%   duration of particle source will be defined instead of maximum rate.
%   See function create_gas_source to see how.

smoothness = 0.004;
num_of_days = ceil(tvect(end)/(3600*24));
source_length = N_tot/rate_max;

out = [tvect', zeros(length(tvect),1), particle_size.*ones(length(tvect),1)];

for i=1:num_of_days
    out(:,2) = out(:,2) + sigmoid(tvect, t_start, t_start + source_length, rate_max, smoothness)';
    t_start = 24*3600 + t_start;
end

end


function [out] = sigmoid(x, x1, x2, max_val, coeff)

out = max_val./(1+exp(-coeff*(x-x1)))-max_val./(1+exp(-coeff*(x-x2)));

tolerance = 0.005;
out(out<tolerance) = 0;
end



