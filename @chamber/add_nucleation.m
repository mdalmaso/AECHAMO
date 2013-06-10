function [ dy ] = add_nucleation(obj,dy,t, part_source)
%ADD_NUCLEATION Calculates the nucleation and adds it to dy.
% 
% [dy] = add_nucleation(obj,dy,t,part_source)
% Calculates the dN/dt for each section and adds it to dy.

% (c) Pauli Simonen 2013
% Version history:
% 2013-06-10    0.1.0

initials = obj.initials;

if(initials.part_source_is_vect)
    for i=1:length(part_source(1,1,:))
        part_source_temp = interp1(part_source(:,1,i), part_source(:,2:3,i), t,'linear',0);
        index = part_source_temp(2);
        dy(1+index) = dy(1+index) + part_source_temp(1);
    end
%     else % Korjaa: index yms, part_source on aina vektori tai nolla
%         dy(1+index) = dy(1+index)+part_source(1);
end

end

