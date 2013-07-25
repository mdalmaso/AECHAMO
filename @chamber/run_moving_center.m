function run_moving_center(obj)
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
% -Make obj.sections stand only for number of sections. Create a new
%  variable for user-defined section spacing.
% -Constant particle sources.
% -Remove nested functions.

initials=obj.initials;

% Make a logarithmically spaced vector between 10^(Dp_min) and 10^(Dpmax).
% Number of cells is initials.sections.
% If initials.sections is not a scalar, the user has already defined the
% sections. Otherwise initials.sections tells the number of sections.
if(isscalar(initials.sections))
    Dp=logspace(initials.Dp_min,initials.Dp_max,initials.sections);
else
    Dp = unique(initials.sections);
    initials.sections = length(Dp);
end


% Dp_variable keeps track of the diameter inside a section. In the
% beginning it is the same as Dp.
Dp_variable = Dp;

% All but the last section have a limit.
Dplims=zeros(1,length(Dp)-1); 

for i = 1:length(Dplims)
    exp1 = log10(Dp(i))/log10(10);  % Find the exponents of current and
    exp2 = log10(Dp(i+1))/log10(10);% next Dp.
    % Then create a logarithmically spaced 3-element vector between these Dp:s.
    temp = logspace(exp1,exp2,3);
    
    % Now the second value of the vector is the logarithmic center between
    % Dp:s and will be the upper limit of section i.
    Dplims(i)=temp(2);
end

% Form a lognormal distribution dN/dlogDp where the minimum diameter is
% vector Dp's first value and maximum is Dp's last value. Number of cells
% is the same as Dp's length (i.e. initials.sections).
dNdlogDp = zeros(size(Dp));
for i=1:length(initials.mu)
    dNdlogDp=dNdlogDp + obj.log_normal(Dp,initials.mu(i),initials.sigma(i),initials.N(i));
end

% If the particle source is defined, search the index that corresponds the
% diameter defined in part_source. Replace the diameter in part_source with
% the index.
if(initials.part_source_is_vect)
    % If part_source is a 3d-matrix, length(part_source(1,1,:)) > 1, and
    % the loop will go through all particle sources.
    for i=1:length(initials.part_source(1,1,:))
        ind = 1;
        while(Dplims(ind) < initials.part_source(1,3,i))
            ind = ind+1;
        end
        initials.part_source(:,4,i) = initials.part_source(1,3,i);
        initials.part_source(:,3,i) = ind;        
    end
    clear ind;
end


% Get the concentration of particles (N) in each cell (=section) of the
% distribution (Dp,dNdlogDp). N(i) is the concentration of particles in
% section i.
[~, N]  = obj.Dlog_to_N_vect(Dp,dNdlogDp);
N = abs(N); % Make sure that concentrations are positive 
            % (Dlog_to_N_vect may make near-zero concentrations negative).


AE_Wall  = 0; % aerosol lost to wall (molecules)
AE_dilu  = 0; % aerosol lost to dilution (molecules)
Vap_dilu = 0; % vapor lost to dilution (molecules)
Vap_Wall = 0; % vapor lost to walls (molecules)

% Initial conditions:
y0 = [initials.Cvap0 N Dp AE_Wall AE_dilu Vap_dilu Vap_Wall Dp_variable eps];

% [y0, C3] = add_vapor(y0, C1, C2, 1e5, 'constant');
% C3.constant = 0;
% C3.source = [...];

% Error tolerance options for ode45:
absTol = ones(size(y0)).*1e-6; % Preallocate the tolerance vector.
absTol(1) = initials.Cvap_tol; % Vapor concentration tolerance
absTol(2:initials.sections+1) = initials.N_tol; % Particle concentration tolerance
absTol(initials.sections+2:(2*initials.sections+1)) = initials.Dp_tol; % Particle diameter tolerance
absTol(2*initials.sections+6:(3*initials.sections+5)) = initials.Dp_tol; % Particle diameter tolerance
% absTol(3*initials.sections+6) = 20*eps;

% Apply ode45 options: the tolerance settings and the function 'events' as
% Event-function.
options = odeset('absTol',absTol,'Events',@events, 'NonNegative', 1:initials.sections+1,'NonNegative', 2*initials.sections+6:3*initials.sections+5); 

% Summary:
% concentrations: Y(2:nSec+1)
% diameters:      Y((nSec+2):(2*nSec+1))
% lost:           Y(2*nSec+2 2*nSec+3 2*nSec+4 2*nSec+5)
% variable diams: Y((2*nSec+6 : 3*nSec+5))
yout = [];
tout = [];
ye = 0;

tvect = initials.tvect;

tout=tvect(1);
yout=y0;



t_start=tvect(1);
delta_t = tvect(2)-tvect(1);
t_end = tvect(end);

options = odeset(options, 'MaxStep', delta_t);

% Define the time span for ode so that it equals the user defined tvect.
t_span = t_start:delta_t:t_end;

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
        
        % Discard events that occur in sections that are next to each
        % other:
        difference = diff(ie);
        difference = [2;difference];
        ie=ie(difference ~= 1)
        
        y0=ye(1,:); % And the first row of ye as well.
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
        if(ie(i) == initials.sections)
            % If the ie is nSec, it means that particles are to
            % nucleate into an empty section. If this happens, the diameter
            % of the section must be moved to correspond the diameter of
            % nucleating particles.
            for j=1:length(initials.part_source(1,1,:))                
                % Index of the section where particles nucleate is stored in
                % initials.part_source, as well as the diameter of nucleating
                % particles.
                index = initials.part_source(1,3,j);
                diam = initials.part_source(1,4,j);
                if(y0(1+index) < eps)
                    y0(2*nSec+5+index) = diam;
                end
            end
            
            % Set the value of y(3*nSec+6) back to eps to indicate that the
            % diameter of the section is now all right.
            y0(3*nSec+6) = eps;
            clear index;
            clear diam;
        
        % Otherwise particles just need to be moved to another section.
        elseif(y0(1+ie(i)) > 0)    % Is there particles in section?
            Ni1=y0(1+ie(i));   % Number of particles in section ie
            Ni2=y0(1+ie(i)+1); % Number of particles in section ie+1
            Ntot = Ni1+Ni2;
            v1 = pi/6*y0(2*nSec+5+ie(i))^3*Ni1;    % Total vol of particles in section ie
            v2 = pi/6*y0(2*nSec+1+5+ie(i))^3*Ni2;  % Total vol of particles in section ie+1
            vtot = (v1+v2)/Ntot;           % Average vol of particles in section ie+1 when the particles from ie are moved there.
            if(vtot > 0)
                y0(2*nSec+5+1+ie(i)) = (6/pi*vtot)^(1/3);  % New average diameter inside section ie+1
            end
            y0(1+ie(i)+1)=Ntot; % Add particles from section ie to ie+1
            y0(1+ie(i)) = 0;       % Delete particles from section ie.
            
            % Reset the diameter of the section that has now zero
            % particles. If the diameter is not reset, it will grow
            % immediatly over its limit if particles coagulate into it and
            % grow by condensation.
            y0(2*nSec+5+ie(i)) = y0(nSec+1+ie(i));
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
    
    if(modulo_time ~= 0)
        % If the last element of t does not equal any of the elements in
        % the time vector, save all other elements but the first and last,
        % because the first one is already saved in previous round and the
        % elements between the first and last do equal time vector's
        % elements.
        tout = [tout; t(2:nt-1)];
        yout = [yout; y(2:nt-1,:)];
        
        %TODO: It's possible that event happens between y(nt-2) and y(nt-1),
        %so the diameters in y(nt-1) may be too big, the correct values are
        %in y0, but represent the y at time t_span(1).
        
        
        % Redefine the t_span for ode so that the first step is
        % delta_t-modulo_time. In this way the second element of t_span
        % (and thus t) will equal some of the elements of time vector. The
        % next elements will have spacing of delta_t.
        t_span = [t(nt), t(nt)+(delta_t-modulo_time):delta_t:t_end];
        
    elseif(length(t) > length(t_span))
        % If t is longer than t_span, it means that for the previous round
        % t_span has had only two elements and in that case ode will return
        % t with denser spacing than delta_t. This can happen only at the
        % end of time vector, so this means we have reached the end. That
        % is why only the last element of t will be saved; other elements
        % do not equal any of the elements in the time vector.
        tout = [tout; t(nt)];
        yout = [yout; y(nt,:)];
        
        % Let t_span equal t_end, so the while loop will be terminated.
        t_span = t_end;
        
    else
        % Now the event has happened at such time that equals some
        % element of the time vector OR ode has reached the end of the time
        % vector. In that case, save all the values of t and y except the
        % first row.
        yout = [yout; y(2:nt,:)];
        tout = [tout; t(2:nt)];
        
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
    temp = [yout(:,1:nSec+1),yout(:,2*nSec+6:3*nSec+5),yout(:,2*nSec+2:2*nSec+5)];        
    Y2=interp1(tout,temp,tvect);
    tout = tvect;
else
    % Now tout equals tvect, so there is no need for interpolation.
    
    % Y2 will be similar to Y, but only the variable diameters will be
    % saved to get as much information as possible. The fixed diameters
    % will be skipped.
    Y2 = [yout(:,1:nSec+1),yout(:,2*nSec+6:3*nSec+5),yout(:,2*nSec+2:2*nSec+5)];
end

display('Ode45 finished, processing data...');

% This makes a handy structure of the results:
obj.output_data = obj.model_convert(tout,Y2);


function dy = chamberODE(t,y)
    dy = zeros(size(y));
    
    part_source = initials.part_source;
    Source = initials.gas_source ; % the condensing vapor source rate (1/cm^3/s) (scalar or 2-column array)
    Dilu   = initials.dilu_coeff ; % the dilution coefficient (scalar or 2-column array)
    Vap_wallsink = initials.vap_wallsink;
    Csat   = initials.satu_conc ; % The condensing vapor saturation concentration
    lambda = initials.lambda ; % the condensing vapor mean free path
    diffu  = initials.diff_coeff ; % the condensing vapor diffusion coefficient
    mv     = initials.vap_molmass ; % the condensing vapor molecular weight
    rool   = initials.particle_dens ; % the particle density
    alfa   = initials.stick_coeff ; % the sticking coefficient
    nSec   = initials.sections; % the number of size sections
    T      = initials.T;               % Temperature
    
    Df     = initials.Df;  % Fractal dimension of agglomerates. Used only if agglomeration is on.
    r0     = initials.r0;  % Radius of agglomerate primary particles. Used only if agglomeration is on.


    CX     = initials.coag_on; % Coagulation switch.
    coag_possible = CX;
    
    % Get the coagulation mode (coagulation or agglomeration).
    % 1 = coagulation
    % 0 = agglomeration
    coagmode = initials.coag_num;
    
    NA = 6.022e23; % Avogadro constant
    
    % If dilu_coeff is defined as vector, interpolate it to find the value of
    % dilu_coeff for the current t.
    if initials.dilu_vect_on, 
        Dilu = interp1(Dilu(:,1),Dilu(:,2),t,'linear',0);
    end
    
    % If gas_source is defined as vector, interpolate it to find the value of
    % gas_source for the current t.
    if initials.gas_source_is_vect,
        Source = interp1(Source(:,1),Source(:,2),t,'linear',0);
    end
    
    % Nucleation:
    dy = obj.add_nucleation(dy, y,t, part_source);
    
    % Calculation of coagulation kernels:
    if(all(diff(y(2*nSec+6:3*nSec+5))>0))
        % Calculate the coagulation kernels only if the diameter vector is
        % increasing. In this case it is possible to calculate the
        % coagulation, so set coag_possible to 1.
        coag_possible = 1;
        % Make coagulation kernel. Different functions for coagulation and
        % agglomeration.
        kk=zeros(nSec,length(y((2*nSec+6):(3*nSec+5)))); % Preallocate
        if(coagmode == 1) % coagmode == 1 => particles coagulate.
            for i = 1:nSec,
                kk(i,:) = obj.koag_kernel(y(2*nSec+5+i),y((2*nSec+6):(3*nSec+5)),rool,T).*1e6;
            end
        else    % Else coagmode == 0 => particles agglomerate.
            for i = 1:nSec,
                kk(i,:) = obj.aggl_kernel(y(2*nSec+5+i),y((2*nSec+6):(3*nSec+5)),rool,T,Df,r0).*1e6;
            end
        end
    else
        % If diameter vector is not increasing, coagulation is not possible
        % (because coagulationMatrix() will return NaN and Inf).
        coag_possible = 0;
    end
        
    % Show the time evolution in Matlab command window:
    time = t
    
    % Increase condensing vapor concentration by vapor coming from Source
    dy(1) = dy(1) + Source;

    % Dilution of vapor
    if initials.dilu_on,
        % Dilute the vapor:
        dy(1) = dy(1)-Dilu.*y(1);
        % And save the information about lost vapor molecules:
        dy(2*nSec+4)  = dy(2*nSec+4)+Dilu.*y(1);
    end
    
    % Deposition of vapor on walls
    if initials.vap_wallsink_on,
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
        if initials.dilu_on
            dy(i+1) = dy(i+1)-Dilu.*y(i+1);   % Decrease particle concentration

            % Calculate the molecules lost to dilution in aerosol phase and
            % save information to dy(2*nSec+3).
            dy(2*nSec+3) = dy(2*nSec+3)+ Dilu.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(2*nSec+5+i).^3))./(6.*mv);
        end
        
        % coagulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if (CX && coag_possible)    
            % calculate a coagulation matrix
            % this tells how to partition the particles 
            cM = obj.coagulationMatrix(y((2*nSec+6):(3*nSec+5)),i);
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
            dy = obj.add_condensation(dy, y, initials, i);
        end
        
        % wall losses and sedimentation according to T. Anttila model...fitted
        % by M. Dal Maso; only usable for SAPPHIR chamber!!
        if (initials.sedi_on && coag_possible)
            beta = obj.sapphir_beta2(y(2*nSec+5+i),T);
            %beta = 3.5e-5; % 0th order approx
            dy(i+1) = dy(i+1)-beta.*y(i+1);

            %molecules lost to the wall = volume lost to wall (in aerosol phase)
            dy(2*nSec+2) = dy(2*nSec+2)+ beta.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(2*nSec+5+i).^3))./(6.*mv);               
        end
    end %% End for loop.
    
    % If vapor concentration is kept constant, reset the value of dy(1) to
    % initial.Cvap0.
    if(initials.Cvap_const == 1)
        dy(1) = 0;
    end
end  % End function chamberODE


function[value,isterminal,direction] = events(t,y)
    nSec = initials.sections;
    
    Dps = y(2*nSec+6:3*nSec+4);
    limits = Dplims';
    
    value = limits-Dps; % Run ode45 as long as Dp(i) < Dplim(i). When value <= 0, ode45 is terminated.
    
    % Add y(3*nSec+6) to the end of value vector. This value indicates if
    % particles are to nucleate into an empty section. If this happens, ode
    % must be stopped and the diameter of the section must be moved to
    % correspond the diameter of nucleating particles.
    value = [value; y(3*nSec+6)];
    
    isterminal = ones(length(value),1);
    direction = ones(length(value),1).*(-1);
end




end