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

% TODO:
% -Make obj.sections stand only for number of sections. Create a new
%  variable for user-defined section spacing.
% -Constant particle sources.

function run_moving_center(obj)
% RUN_MOVING_CENTER Runs the simulation with fixed sections using ode45.
% Saves the results to chamber.output_data.
%  
% A private function of class chamber, used by public method chamber.run.

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
        while(Dp(ind) < initials.part_source(1,3,i))
            ind = ind+1;
        end
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

% Initial conditions:
y0 = [initials.Cvap0 N Dp AE_Wall AE_dilu Vap_dilu Dp_variable];


% Error tolerance options for ode45:
absTol = ones(size(y0)).*1e-6; % Preallocate the tolerance vector.
absTol(1) = initials.Cvap_tol; % Vapor concentration tolerance
absTol(2:initials.sections+1) = initials.N_tol; % Particle concentration tolerance
absTol(initials.sections+1:(2*initials.sections+1)) = initials.Dp_tol; % Particle diameter tolerance

% Apply ode45 options: the tolerance settings and the function 'events' as
% Event-function.
options = odeset('absTol',absTol,'Events',@events); 


% Summary:
% concentrations: Y(2:nSec+1)
% diameters:      Y((nSec+2):(2*nSec+1))
% lost:           Y(2*nSec+2 2*nSec+3 2*nSec+4)
% variable diams: Y((2*nSec+5 : 3*nSec+4))
yout = [];
tout = [];
ye = 0;

tvect = initials.tvect;

tout=tvect(1);
yout=y0;



t_start=tvect(1);
delta_t = tvect(2)-tvect(1);
t_end = tvect(end);

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
        display('ie:ssa useampi alkio');
        pause;
        ie=ie(1);   % Take only the first index.
        y0=ye(1,:); % And the first row of ye as well.
    end
    
    % If there are particles in the section that has grown over limit, move
    % the particles to next section and calculate the new diameter inside
    % the next section.
    if(y0(1+ie) > 0)    % Is there particles in section?
        Ni1=y0(1+ie);   % Number of particles in section ie
        Ni2=y0(1+ie+1); % Number of particles in section ie+1
        v1 = pi/6*y0(2*nSec+4+ie)^3*Ni1;    % Total vol of particles in section ie
        v2 = pi/6*y0(2*nSec+1+4+ie)^3*Ni2;  % Total vol of particles in section ie+1
        vtot = (v1+v2)/(Ni1+Ni2);           % Average vol of particles in section ie+1 when the particles from ie are moved there.
        y0(2*nSec+4+1+ie) = (6/pi*vtot)^(1/3);  % New average diameter inside section ie+1
        y0(2*nSec+4+ie) = y0(nSec+1+ie);        % Reset the diameter of section ie.
        y0(1+ie+1)=Ni1+Ni2; % Add particles from section ie to ie+1
        y0(1+ie) = 0;       % Delete particles from section ie.
    else
        y0(2*nSec+4+ie) = y0(nSec+1+ie); % Reset the diameter of section ie.
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
    % t=[0 60 72 120 125]. 
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
    if(modulo_time ~= 0)
        % If the last element of t does not equal any of the elements in
        % the time vector, save all other elements but the first and last,
        % because the first one is already saved in previous round and the
        % elements between the first and last do equal time vector's
        % elements.
        tout = [tout; t(2:nt-1)];
        yout = [yout; y(2:nt-1,:)];
        
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
    temp = [yout(:,1:nSec+1),yout(:,2*nSec+5:3*nSec+4),yout(:,2*nSec+2:2*nSec+4)];        
    Y2=interp1(tout,temp,tvect);
    tout = tvect;
else
    % Now tout equals tvect, so there is no need for interpolation.
    
    % Y2 will be similar to Y, but only the variable diameters will be
    % saved to get as much information as possible. The fixed diameters
    % will be skipped.
    Y2 = [yout(:,1:nSec+1),yout(:,2*nSec+5:3*nSec+4),yout(:,2*nSec+2:2*nSec+4)];
end

display('Ode45 finished, processing data...');

% This makes a handy structure of the results:
obj.output_data = obj.model_convert(tout,Y2);


function dy = chamberODE(t,y)
    dy = zeros(size(y));
    
    part_source = initials.part_source;
    Source = initials.gas_source ; % the condensing vapor source rate (1/cm^3/s) (scalar or 2-column array)
    Dilu   = initials.dilu_coeff ; % the dilution coefficient (scalar or 2-column array)
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
    
    % Get the coagulation mode (coagulation or agglomeration).
    % 1 = coagulation
    % 0 = agglomeration
    coagmode = initials.coag_num;
    
    NA = 6.022e23; % Avogadro constant
    
    if initials.dilu_vect_on, 
        Dilu = interp1(Dilu(:,1),Dilu(:,2),t,'linear',0);
    end
    
    % If gas_source is defined as vector, interpolate it to find the value of
    % gas_source for the current t.
    if initials.gas_source_is_vect,
        Source = interp1(Source(:,1),Source(:,2),t,'linear',0);
    end
    
    if(initials.part_source_is_vect)
        for i=1:length(part_source(1,1,:))
            part_source_temp = interp1(part_source(:,1,i), part_source(:,2:3,i), t,'linear',0);
            index = part_source_temp(2);
            dy(1+index) = dy(1+index)+part_source_temp(1);
        end
%     else % Korjaa: index yms, part_source on aina vektori tai nolla
%         dy(1+index) = dy(1+index)+part_source(1);
    end
    
    
    
    % Make coagulation kernel. Different functions for coagulation and
    % agglomeration.
    kk=zeros(nSec,length(y((2*nSec+5):(3*nSec+4)))); % Preallocate
    if(coagmode == 1) % coagmode == 1 => particles coagulate.
        for i = 1:nSec,
            kk(i,:) = obj.koag_kernel(y(2*nSec+4+i),y((2*nSec+5):(3*nSec+4)),rool,T).*1e6;
        end
    else    % Else coagmode == 0 => particles agglomerate.
        for i = 1:nSec,
            kk(i,:) = obj.aggl_kernel(y(2*nSec+4+i),y((2*nSec+5):(3*nSec+4)),rool,T,Df,r0).*1e6;
        end
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
    
    % Go through all particle diameters and calculate the effect of dilution,
    % coagulation, condensation and sedimentation on particle
    % concentrations.
    for i = 1:nSec
        % Dilution of particles
        if initials.dilu_on
            dy(i+1) = dy(i+1)-Dilu.*y(i+1);   % Decrease particle concentration

            % Calculate the molecules lost to dilution in aerosol phase and
            % save information to dy(2*nSec+3).
            dy(2*nSec+3) = dy(2*nSec+3)+ Dilu.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(2*nSec+4+i).^3))./(6.*mv);
        end
        
        % coagulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if CX,    
        % calculate a coagulation matrix
        % this tells how to partition the particles 
        cM = obj.coagulationMatrix(y((2*nSec+5):(3*nSec+4)),i);
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
        Kn = (2.*lambda)./y(2*nSec+4+i); % Knudsen number
        betam = (Kn+1)./((0.377.*Kn)+1+(4/(3.*alfa)).*(Kn.^2)+(4/(3.*alfa)).*Kn);
    
        % I is the flux of molecules to the particle phase    
        I = 2.*pi.*max([y(2*nSec+4+i) 0]).*1e2.*diffu.*(y(1)-Csat).*betam; %1/s
        
        % Move the variable diameters by condensation.
        dy(2*nSec+4+i) = dy(2*nSec+4+i)+(2.*mv.*I)./(pi.*rool.*y(2*nSec+4+i).^2.*NA.*1e6); % particle diameter (m/s)
        dy(1) = dy(1) - y(i+1).*I;
        % end condensation %%%%
        
        % wall losses and sedimentation according to T. Anttila model...fitted
        % by M. Dal Maso; only usable for SAPPHIR chamber!!
        if initials.sedi_on,
            beta = obj.sapphir_beta2(y(2*nSec+4+i),T);
            %beta = 3.5e-5; % 0th order approx
            dy(i+1) = dy(i+1)-beta.*y(i+1);

            %molecules lost to the wall = volume lost to wall (in aerosol phase)
            dy(2*nSec+2) = dy(2*nSec+2)+ beta.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(2*nSec+4+i).^3))./(6.*mv);               
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
    
    Dps = y(2*nSec+5:3*nSec+3);
    limits = Dplims';
    
%     value = limits-Dps+eps(Dps); % Run ode45 as long as Dp(i) <= Dplim(i)
    value = limits-Dps; % Run ode45 as long as Dp(i) < Dplim(i). When value == 0, ode45 is terminated.
    isterminal = ones(length(value),1);
    direction = ones(length(value),1).*(-1);
end




end