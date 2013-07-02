function [ out ] = create_source(tvect, source_length, vap_tot, t_start)
%GAS_SOURCE Summary of this function goes here
%   Detailed explanation goes here

num_of_days = ceil(tvect(end)/(3600*24));

vap_max = vap_tot/source_length;

out = [tvect', zeros(length(tvect),1)];

for i=1:num_of_days
    out(:,2) = out(:,2) + sigmoid(tvect, t_start, t_start + source_length, vap_max, 0.003)';
    t_start = 24*3600 + t_start;
end

end


function [out] = sigmoid(x, x1, x2, max_val, coeff)

out = max_val./(1+exp(-coeff*(x-x1)))-max_val./(1+exp(-coeff*(x-x2)));

tolerance = 1;
out(out<tolerance) = 0;
end