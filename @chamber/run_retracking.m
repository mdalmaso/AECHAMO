function [out_t, out_Y] = run_retracking(obj)
%RUN_RETRACKING Summary of this function goes here
%   Detailed explanation goes here
% RUN_MOVING_CENTER Runs the simulation with fixed sections using ode45.
% Saves the results to chamber.output_data.
%  
% A private function of class chamber, used by public method chamber.run.

% Moving center solver
% (c) Pauli Simonen 2013
%
% Version history:
% 2013-05-31    0.1.0
% 2013-06-03    0.1.1 Commented and tidied.
% 2013-06-04    0.1.2 Added the option to create particle source(s). In
%                     addition, user can now define sections by hand.
% 2013-06-06    0.1.3 Modified the main while loop so that only the points
%                     of ode's output t that equal user defined time vector
%                     are saved to the output data. So there is no need to
%                     interpolate the output data to match user defined
%                     time vector anymore.
% 2013-06-12    0.1.4 Now the condensation is calculated only for sections
%                     that contain particles. Huge improvement of
%                     efficiency and makes also the distribution look more
%                     natural.
% 2013-06-20    0.1.5 Added condensation sink to chamber walls. 
% 2013-07-03    0.1.6 Added a new cell to y to indicate if particles
%                     nucleate to a section that has no particles and which
%                     diameter is different to the diameter of nucleating
%                     particles. If this new cell goes below zero, ode is
%                     stopped and the diameter of the section is changed.

% TODO:
% -Constant particle sources.

% Get the initial parameters:
params=obj.initials;

% And diameter vector:
Dp = obj.Dps;

% And the particle number distribution. N(i) represents the particle
% concentration in section i, and the initial diameter of section i is
% Dp(i).
N = obj.number_distribution;

nSec = params.sections;

% Dp_variable keeps track of the diameter inside a section. In the
% beginning it is the same as Dp unless user has set the center diameters.
Dp_variable = obj.center_diameters;

% If the particle source is defined, search the index that corresponds the
% diameter defined in part_source. Replace the diameter in part_source with
% the index.
if(params.part_source_is_vect)
    % If part_source is a 3d-matrix, length(part_source(1,1,:)) > 1, and
    % the loop will go through all particle sources.
    for i=1:length(params.part_source(1,1,:))
        ind = 1;
        while(obj.Dplims(ind) < params.part_source(1,3,i))
            ind = ind+1;
        end
        params.part_source(:,4,i) = params.part_source(1,3,i);
        params.part_source(:,3,i) = ind;
    end
    clear ind;
end

AE_Wall  = 0; % aerosol lost to wall (molecules)
AE_dilu  = 0; % aerosol lost to dilution (molecules)
Vap_dilu = 0; % vapor lost to dilution (molecules)
Vap_Wall = 0; % vapor lost to walls (molecules)

% Initial conditions:
y0 = [params.Cvap0, N, Dp_variable, AE_Wall, AE_dilu, Vap_dilu, Vap_Wall, eps];

% Error tolerance options for ode45:
absTol = ones(size(y0)).*1e-6; % Preallocate the tolerance vector.
absTol(1) = params.Cvap_tol; % Vapor concentration tolerance
absTol(2:params.sections+1) = params.N_tol; % Particle concentration tolerance
absTol(params.sections+2:(2*params.sections+1)) = params.Dp_tol; % Particle diameter tolerance
% absTol(3*params.sections+6) = 20*eps;

% Apply the tolerance settings to ode45 options:
options = odeset('absTol',absTol, 'NonNegative', 1:params.sections+1); 

% Set the function 'events' as Event-function:
options = odeset(options,'Events',@events);

% Summary:
% concentrations: Y(2:nSec+1)
% diameters:      Y((nSec+2):(2*nSec+1))
% lost:           Y(2*nSec+2 2*nSec+3 2*nSec+4 2*nSec+5)
% variable diams: Y((2*nSec+6 : 3*nSec+5))
% yout = [];
% tout = [];

ye = 0;

tvect = params.tvect;

% Preallocate yout and tout:
yout = zeros(length(tvect),length(y0));
tout = zeros(length(tvect),1);

tout(1)=tvect(1);
yout(1,:)=y0;

t_start=tvect(1);
delta_t = tvect(2)-tvect(1);
t_end = tvect(end);

% Set the MaxStep to delta_t so ode will not miss for example short
% nucleation events. This slows down the simulation remarkably and is
% unnecessary in most cases. Consider commenting this out or handling the
% problem in some other way.
if(params.max_timestep)
    options = odeset(options, 'MaxStep', params.max_timestep);
end

% Define the time span for ode so that it equals the user defined tvect.
t_span = t_start:delta_t:t_end;

run_ind = 2;

% Create a visual waitbar (slows the program down a little, a few seconds)
h = waitbar(0,'0 %','Name','Running simulation...');


retrack_interval = 60;
retrack_moment = retrack_interval;

% Run the simulation as long as the beginning of the ode's time span vector
% is smaller than the end of the user defined tvect.
while(t_span(1) < tvect(end))
    % Run ode45 until one or more Dp:s grow over limit (i.e. the event
    % occurs) OR until the end of t_span is reached.

    [t,y,te,ye,ie]=ode45(@chamberODE,t_span,y0,options);

    % Now ode45 has stopped, so some Dp has grown over limit OR ode45 has 
    % reached the end of t_span. If the end of t_span is reached, the
    % first element of the t_span will be later in this loop redefined so
    % that it equals user defined tvect's last element. This will terminate
    % the while loop.
    %
    % Otherwise ode45 returns te, ye and ie:
    % te is the time when the event occured,
    % ye is the y-vector at time te
    % ie is the index of diameter(s) that grew over limit

    % Get the time vector's size
    nt = length(t);
    
    % Define the initial conditions as y(nt,:), so ode45 will continue from the
    % same point where the event occured.
    y0 = y(nt,:);
    
    % If there are more than one index, several Dp:s have grown over their
    % limits. In that case, the program will take care of the Dp with
    % lowest index. The other Dp:s will be handled in the next loop if
    % needed.
    if(length(ie)>1)
        i = 1;
        while( i < length(ye(:,1)) && isequal(ye(i,:),ye(i+1,:)))
            i = i+1;
        end
        ie=ie(1:i);
        clear i;

        % If there are two events that occur in sections that are next to
        % each other, discard the other one of them, so that only the first
        % one of consecutive sections is handled.
        difference = diff(ie);
        difference = [2;difference];
        ie=ie(difference ~= 1)
        
        y0=ye(1,:);
        te=te(1);
        
        % Redefine t so that it ends to the time point where the first
        % event occurs:
        for i=1:length(t)
            if(t(i) <= te)
                t_temp(i) = t(i);
            else
                break;
            end
        end
        t = t_temp';
        clear t_temp;
        nt = length(t);
    end

    % Set the initial step of ode to one second, so the first step of
    % solver after the event will not be too big.
    options = odeset(options, 'InitialStep', 1);
    
    % If there are particles in the section that has grown over limit, move
    % the particles to next section and calculate the new diameter inside
    % the next section.
    for i=1:length(ie)
        if(ie(i) == nSec)
            % If the ie is nSec, it means that particles are to
            % nucleate into an empty section. If this happens, the diameter
            % of the section must be moved to correspond the diameter of
            % nucleating particles.
            for j=1:length(params.part_source(1,1,:))                
                % Index of the section where particles nucleate is stored in
                % params.part_source, as well as the diameter of nucleating
                % particles.
                index = params.part_source(1,3,j);
                diam = params.part_source(1,4,j);
                if(y0(1+index) < eps)
                    y0(nSec+1+index) = diam;
                end
            end
            
            % Set the value of y(3*nSec+6) back to eps to indicate that the
            % diameter of the section is now all right.
            y0(3*nSec+6) = eps;
            clear index diam;
       
        elseif(ie(i) == 1)
            % Retracking.
%             % Omatekoinen "interpolaatio" alkaa:
%             y0_temp=zeros(1,length(y0));
%             for j=1:nSec-1
%                 if(y0(1+j) > 0)
%                     new_ind = j;
%                     for k=j:nSec-1
%                         if(y0(nSec+1+j) > obj.Dplims(k))
%                             new_ind = k+1;
%                         end
%                     end
%                     Dp0=y0(nSec+1+j);
%                     Dp1=obj.center_diameters(new_ind);
%                     Dp2=obj.center_diameters(new_ind+1);
%                     N0=y0(1+j);
%                     N1 = (N0*(Dp0^3-Dp2^3))/(Dp1^3-Dp2^3);
%                     N2 = N0-N1;
%                     y0_temp(1+new_ind) = y0_temp(1+new_ind)+N1;
%                     y0_temp(1+1+new_ind) = y0_temp(1+1+new_ind)+N2;
%                 end
%             end
%             % Omatekoinen "interpolaatio p‰‰ttyy

            retrack_moment = retrack_interval + te
            
            % Jakauman interpolaatio alkaa:
%             y0(2:nSec+1) = interp1(y0(nSec+2:2*nSec+1),y0(2:nSec+1),obj.center_diameters,'pchip');
            % Jakauman interpolaatio p‰‰ttyy
            
%             y0(2:nSec+1)=y0_temp;
%             y0(2:nSec+1) = y0_temp(2:nSec+1);
%             clear y0_temp;

            % DndlogDp interpolaatio alkaa
            dn = obj.N_to_dlog(y0(nSec+2:2*nSec+1),y0(2:nSec+1));
            dn_new = interp1(log10(y0(nSec+2:2*nSec+1)),dn,log10(obj.center_diameters),'pchip');
            if(~all(dn_new>=0))
                pause;
            end
            [~, y0(2:nSec+1)] = obj.Dlog_to_N_vect(obj.center_diameters,dn_new);
            % DndlogDp interpolaatio p‰‰ttyy
            
            y0(nSec+2:2*nSec+1) = obj.center_diameters;
        end
    end
    % The t and y vectors from ode will be saved to cumulative output
    % vectors tout and yout. The first row of y and t is the same as the
    % last row of previous run, so it will not be saved twice.
    %
    % If mod(t(end),delta_t) ~= 0, it means that the event has occured at a
    % time that does not equal any element of the time vector. In that case
    % this last time event shall not be saved into output, because we want
    % to get tout that equals the user input time vector.
    %
    % Example: Time vector is [0 60 120 180 ...] which means that
    % delta_t = 60. Event occurs at time 125, so ode returns 
    % t=[0 60 120 125]. 
    % Then mod(t(end), delta_t) = 5 ~= 0.
    % Only values 60 and 120 of t and respective values of y are saved to
    % tout and yout.
    % Ode will continue calculating from time 125, and the first time step
    % is delta_t-mod(t(nt),delta_t) = 60-5 = 55. The next time steps are 
    % delta_t. 
    % So the time span vector that is input to ode is [125 180 240 ...] and
    % ode will eventually return t = [125 180 ...] and in this way only the
    % values of t that equal the values of user input time vector are saved
    % to the output data.
    modulo_time = mod(t(end),delta_t);
    
    % If t consists of only one element, it has the same value as the last
    % element of t during previous round.
    if(length(t) < 2)
        t_span = [te, te+(delta_t-modulo_time):delta_t:t_end];
        continue;
    end
    
%     if(modulo_time ~= 0)
%         % If the last element of t does not equal any of the elements in
%         % the time vector, save all other elements but the first and last,
%         % because the first one is already saved in previous round and the
%         % elements between the first and last do equal time vector's
%         % elements.
%         length_addition = length(t(2:nt-1));
%         tout(run_ind:run_ind+length_addition-1) = t(2:nt-1);
%         yout(run_ind:run_ind+length_addition-1,:) = y(2:nt-1,:);
%         run_ind = run_ind + length_addition;
%         
%         % TODO: It's possible that event happens between y(nt-2) and y(nt-1),
%         % so the diameters in y(nt-1) may be too big, the correct values are
%         % in y0, but represent the y at time t_span(1).
%         
%         
%         % Redefine the t_span for ode so that the first step is
%         % delta_t-modulo_time. In this way the second element of t_span
%         % (and thus t) will equal some of the elements of time vector. The
%         % next elements will have spacing of delta_t.
%         t_span = [t(nt), t(nt)+(delta_t-modulo_time):delta_t:t_end];
%         
%     elseif(length(t) > length(t_span))
%         % If t is longer than t_span, it means that for the previous round
%         % t_span has had only two elements and in that case ode will return
%         % t with denser spacing than delta_t. This can happen only at the
%         % end of time vector, so this means we have reached the end. That
%         % is why only the last element of t will be saved; other elements
%         % do not equal any of the elements in the time vector.
%         
%         length_addition = 1;
%         tout(run_ind:run_ind+length_addition-1) = t(nt);
%         yout(run_ind:run_ind+length_addition-1,:) = y(nt,:);
%         run_ind = run_ind + length_addition;
%         
%         % Let t_span equal t_end, so the while loop will be terminated.
%         t_span = t_end;
        
    if(te ~= t_end)
        % Now the event has happened at such time that equals some
        % element of the time vector OR ode has reached the end of the time
        % vector. In that case, save all the values of t and y except the
        % first row.
        length_addition = length(t(2:nt-1));
        tout(run_ind:run_ind+length_addition-1) = t(2:nt-1);
        yout(run_ind:run_ind+length_addition-1,:) = y(2:nt-1,:);
        run_ind = run_ind + length_addition;

        % And redefine the t_span for ode beginning from current time point
        % to the end of the time vector with spacing of delta_t.
        t_span = [t(nt):delta_t:t_end];
        if(length(t_span) < 2)
            t_span = [t(nt), t(nt)+delta_t];
        end
    else
        
        length_addition = length(t(2:nt));
        tout(run_ind:run_ind+length_addition-1) = t(2:nt);
        yout(run_ind:run_ind+length_addition-1,:) = y(2:nt,:);
        run_ind = run_ind + length_addition;

        % And redefine the t_span for ode beginning from current time point
        % to the end of the time vector with spacing of delta_t.
        t_span = [t(nt):delta_t:t_end];
    end
end

% If the length of tout is longer than the length of tvect, it means that
% the spacing of tout is different to tvect, so we need to interpolate the
% results to fit tvect. This should happen rarely or not at all.
if(length(tout) > length(tvect))
    % Y2 will be similar to Y, but only the variable diameters will be
    % saved to get as much information as possible. The fixed diameters
    % will be skipped.
    Y2 = interp1(tout,yout,tvect);
    tout = tvect;
else
    % Now tout equals tvect, so there is no need for interpolation.
    Y2 = yout;
end

out_Y = Y2;
out_t = tout;

% Close waitbar:
close(h);


function dy = chamberODE(t,y)
    dy = zeros(size(y));
    
    part_source = params.part_source;
    Source = params.gas_source ; % the condensing vapor source rate (1/cm^3/s) (scalar or 2-column array)
    Dilu   = params.dilu_coeff ; % the dilution coefficient (scalar or 2-column array)
    Vap_wallsink = params.vap_wallsink;
    Csat   = params.satu_conc ; % The condensing vapor saturation concentration
    lambda = params.lambda ; % the condensing vapor mean free path
    diffu  = params.diff_coeff ; % the condensing vapor diffusion coefficient
    mv     = params.vap_molmass ; % the condensing vapor molecular weight
    rool   = params.particle_dens ; % the particle density
    alfa   = params.stick_coeff ; % the sticking coefficient
    T      = params.T;               % Temperature
    
    Df     = params.Df;  % Fractal dimension of agglomerates. Used only if agglomeration is on.
    r0     = params.r0;  % Radius of agglomerate primary particles. Used only if agglomeration is on.


    CX     = params.coag_on; % Coagulation switch.
    coag_possible = CX;
    
    % Get the coagulation mode (coagulation or agglomeration).
    % 1 = coagulation
    % 0 = agglomeration
    coagmode = params.coag_num;
    
    NA = 6.022e23; % Avogadro constant
    
    % If dilu_coeff is defined as vector, interpolate it to find the value of
    % dilu_coeff for the current t.
    if params.dilu_vect_on, 
        Dilu = interp1(Dilu(:,1),Dilu(:,2),t,'linear',0);
    end
    
    % If gas_source is defined as vector, interpolate it to find the value of
    % gas_source for the current t.
    if params.gas_source_is_vect,
        Source = interp1(Source(:,1),Source(:,2),t,'linear',0);
    end
    
    % Nucleation:
    dy = obj.add_nucleation(dy, y,t, part_source);
    
    % Calculation of coagulation kernels:
    if(all(diff(y(nSec+2:2*nSec+1))>0))
        % Calculate the coagulation kernels only if the diameter vector is
        % increasing. In this case it is possible to calculate the
        % coagulation, so set coag_possible to 1.
        coag_possible = 1;
        % Make coagulation kernel. Different functions for coagulation and
        % agglomeration.
        kk=zeros(nSec,length(y((nSec+2):(2*nSec+1)))); % Preallocate
        if(coagmode == 1) % coagmode == 1 => particles coagulate.
            for i = 1:nSec,
                kk(i,:) = obj.koag_kernel(y(nSec+1+i),y((nSec+2):(2*nSec+1)),rool,T).*1e6;
            end
        else    % Else coagmode == 0 => particles agglomerate.
            for i = 1:nSec,
                kk(i,:) = obj.aggl_kernel(y(nSec+1+i),y((nSec+2):(2*nSec+1)),rool,T,Df,r0).*1e6;
            end
        end
    else
        % If diameter vector is not increasing, coagulation is not possible
        % (because coagulationMatrix() will return NaN and Inf).
        coag_possible = 0;
    end
        
    % Show the time evolution in Matlab command window:
%     time = t
    % Or update the visual waitbar:
    perc=round(t/params.tvect(end)*100);
    waitbar(perc/100,h,sprintf('%d %%',perc))
    
    % Increase condensing vapor concentration by vapor coming from Source
    dy(1) = dy(1) + Source;

    % Dilution of vapor
    if params.dilu_on,
        % Dilute the vapor:
        dy(1) = dy(1)-Dilu.*y(1);
        % And save the information about lost vapor molecules:
        dy(2*nSec+4)  = dy(2*nSec+4)+Dilu.*y(1);
    end
    
    % Deposition of vapor on walls
    if params.vap_wallsink_on,
        % Dilute the vapor by wall sink coefficient:
        dy(1) = dy(1)-Vap_wallsink.*y(1);
        % And save the information about lost vapor molecules:
        dy(2*nSec+5) = dy(2*nSec+5)+Vap_wallsink.*y(1);
    end
        
    
    % Go through all particle diameters and calculate the effect of dilution,
    % coagulation, condensation and sedimentation on particle
    % concentrations.
    for i = 1:nSec
        % Dilution of particles
        if params.dilu_on
            dy(i+1) = dy(i+1)-Dilu.*y(i+1);   % Decrease particle concentration

            % Calculate the molecules lost to dilution in aerosol phase and
            % save information to dy(2*nSec+3).
            dy(2*nSec+3) = dy(2*nSec+3)+ Dilu.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(nSec+1+i).^3))./(6.*mv);
        end
        
        % coagulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if (CX && coag_possible)    
            % calculate a coagulation matrix
            % this tells how to partition the particles 
            cM = obj.coagulationMatrix(y((nSec+2):(2*nSec+1)),i);
            for j = 1:i,
                if i == j
                    dy(j+1) = dy(j+1)-y(i+1).*y(j+1).*kk(i,j); % loss             
                    dy(2:(nSec+1)) = dy(2:(nSec+1))+cM(j,:)'.*0.5.*y(i+1).*y(j+1).*kk(i,j); % gain
                else
                    dy(i+1) = dy(i+1)-y(i+1).*y(j+1).*kk(i,j); %loss
                    dy(j+1) = dy(j+1)-y(i+1).*y(j+1).*kk(i,j); %loss           
                    dy(2:(nSec+1)) = dy(2:(nSec+1))+cM(j,:)'.*y(i+1).*y(j+1).*kk(i,j); % gain
                end
            end
        end %if CX
        % end coagulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % condensation %%%%
        if(y(i+1) ~= 0)
            dy = obj.add_condensation(dy, y, params, i);
        end
        
        % wall losses and sedimentation according to T. Anttila model...fitted
        % by M. Dal Maso; only usable for SAPPHIR chamber!!
        if (params.sedi_on && coag_possible)
            beta = obj.sapphir_beta2(y(nSec+1+i),T);
            %beta = 3.5e-5; % 0th order approx
            dy(i+1) = dy(i+1)-beta.*y(i+1);

            %molecules lost to the wall = volume lost to wall (in aerosol phase)
            dy(2*nSec+2) = dy(2*nSec+2)+ beta.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(nSec+1+i).^3))./(6.*mv);               
        end
    end %% End for loop.
    
    % If vapor concentration is kept constant, reset the value of dy(1) to
    % initial.Cvap0.
    if(params.Cvap_const == 1)
        dy(1) = 0;
    end
end  % End function chamberODE


function[value,isterminal,direction] = events(t,y)   
%     Dps = y(nSec+2:2*nSec+1);
% 
%     % limits(i) is the upper limit of Dp(i) and the lower limit of Dp(i+1)
%     limits = obj.Dplims'; 
%     
%     val_1 = limits-Dps(1:end-1); % Goes below zero if Dp(i) > limits(i)
%     val_2 = Dps(2:end)-limits;   % Goes below zero if Dp(i+1) < limits(i)
% 
%     indices_1 = find(val_1 < 0); % Find the indices of negative elements.
%     indices_2 = find(val_2 < 0); % Find the indices of negative elements.
%     
%     value=ones(length(val_1),1); % Set the event's values initially as ones.
%     
%     % Add the negative values of val_1 and val_2 to value vector.
%     value(indices_1) = val_1(indices_1); 
%     value(indices_2) = val_2(indices_2);
    
    retracking = retrack_moment - t;
    
    % Add y(3*nSec+6) to the end of value vector. This value indicates if
    % particles are to nucleate into an empty section. If this happens, ode
    % must be stopped and the diameter of the section must be moved to
    % correspond the diameter of nucleating particles.
%     value = [value; 2*nSec+6;retracking];
    value = retracking;
    isterminal = ones(length(value),1);
    direction = zeros(length(value),1);
end

end

