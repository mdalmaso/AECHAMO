function [ dy ] = add_nucleation(obj,dy,y,t, part_source)
%ADD_NUCLEATION Calculates the nucleation and adds it to dy.
% 
% [dy] = add_nucleation(obj, dy, y, t, part_source)
% Calculates the dN/dt for each section and adds it to dy.

% (c) Pauli Simonen 2013
% Version history:
% 2013-06-10    0.1.0
% 2013-07-03    0.1.1 Added cell to y and dy to indicate if particles
%                     nucleate to a section that has no particles and which
%                     diameter is different to the diameter of nucleating
%                     particles. If this new cell goes below zero, ode is
%                     stopped and the diameter of the section is changed.
% 
% TODO:
% -Get the part_source from obj

initials = obj.initials;
nSec = initials.sections;





if(initials.part_source_is_vect)
    for i=1:length(part_source(1,1,:))
        part_source_temp = interp1(part_source(:,1,i), part_source(:,2:4,i), t,'linear',0);
        
        
        %% Miikka: fixing the particle source so that it depends on the CS
        if (part_source_temp(1)>0 && part_source_temp(3)<5e-9),
        load('J_ave_CS.mat');
        CSx = [1e-4 5e-4 1e-3 5e-3 1e-2];
       
        CS = CS_general(y((nSec+2):(2*nSec+1)),y(2:nSec+1),288,1.0)
        if CS>1e-2,
            CS = 1e-2;
        elseif CS<1e-4
            CS = 1e-4;
        end
          
        Jx = interp1(CSx,mat(:,3),CS)

        part_source_temp(1) = part_source_temp(1).*Jx;
        if part_source_temp(1)<1e-3,
            part_source_temp(1) = 0;
        end
        end
        
        index = part_source_temp(2);
        diam = part_source_temp(3);
        dy(1+index) = dy(1+index) + part_source_temp(1);
     
        if part_source_temp(1)>0 & part_source_temp(3)>5e-9,
            disp('adding big')
        end
        
        if(part_source_temp(1) > 0 && y(nSec+1+index) ~= diam)
            % If there are particles that will nucleate and the diameter of
            % the corresponding section is not the same as the diameter of
            % the nucleating particles, the diameter of the section must be
            % changed.
            if(abs(y(1+index)) < eps || y(2*nSec+6) < 0)
                % If the section to which particles will nucleate is empty
                % (or nearly empty), try to set y(3*nSec+6) negative so ode
                % will stop and the diameter of section 'index' will be
                % changed to the diameter of nucleating particles.
                dy(2*nSec+6) = -10*eps;
            else
                Ntot = y(1+index);
                % If there are already particles in the section, just move 
                % the center diameter of corresponding section:
                dy(nSec+1+index) = part_source_temp(1)/(3*Ntot*y(nSec+1+index)^2)*(diam^3-y(nSec+1+index)^3);
            end
        end
    end
%     else % Korjaa: index yms, part_source on aina vektori tai nolla
%         dy(1+index) = dy(1+index)+part_source(1);
end

end
