function [ dy ] = add_nucleation(obj,dy,y,t, part_source)
%ADD_NUCLEATION Calculates the nucleation and adds it to dy.
% 
% [dy] = add_nucleation(obj, dy, y, t, part_source)
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
        
        if(part_source_temp(1) > 0 && y(2*nSec+5+index) ~= diam)
            % If there are particles that will nucleate and the diameter of
            % the corresponding section is not the same as the diameter of
            % the nucleating particles, the diameter of the section must be
            % changed.
            if(abs(y(1+index)) < eps || y(3*nSec+6) < 0)
                % If the section to which particles will nucleate is empty,
                % try to set y(3*nSec+6) negative so ode will stop and the
                % diameter of section 'index' will be changed to the
                % diameter of nucleating particles.
                dy(3*nSec+6) = -10*eps;
            elseif(y(3*nSec+6) > 0)
                Ntot = y(1+index);
                % If there are already particles in the section, just move 
                % the center diameter of corresponding section:
                dy(2*nSec+5+index) = part_source_temp(1)/(3*Ntot*y(2*nSec+5+index)^2)*(diam^3-y(2*nSec+5+index)^3);
            end
        end
    end
%     else % Korjaa: index yms, part_source on aina vektori tai nolla
%         dy(1+index) = dy(1+index)+part_source(1);
end

end
