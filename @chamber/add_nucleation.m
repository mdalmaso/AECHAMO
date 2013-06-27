function [ dy ] = add_nucleation(obj,dy,y,t, part_source)
%ADD_NUCLEATION Calculates the nucleation and adds it to dy.
% 
% [dy] = add_nucleation(obj,dy,t,part_source)
% Calculates the dN/dt for each section and adds it to dy.

% (c) Pauli Simonen 2013
% Version history:
% 2013-06-10    0.1.0
% 
% TODO:
% -Get the part_source from obj

initials = obj.initials;
nSec = initials.sections;

if(initials.part_source_is_vect)
    for i=1:length(part_source(1,1,:))
        part_source_temp = interp1(part_source(:,1,i), part_source(:,2:4,i), t,'linear',0);
        index = part_source_temp(2);
        diam = part_source_temp(3);
        dy(1+index) = dy(1+index) + part_source_temp(1);
        
%         if(abs(y(1+index)) < eps)
%             y(2*nSec+5+index) = diam;
%             
% %             Ntot = eps;
%         else
%             Ntot = y(1+index);
% %         end
%         % Move the center diameter of corresponding section:
%             dy(2*nSec+5+index) = part_source_temp(1)/(3*Ntot*y(2*nSec+5+index)^2)*(diam^3-y(2*nSec+5+index)^3);
%         end
%     end
%     else % Korjaa: index yms, part_source on aina vektori tai nolla
%         dy(1+index) = dy(1+index)+part_source(1);
end

end
