function run_movsec(obj)
% RUN_MOVSEC runs the simulation with moving sections.
% 
% Saves the results to chamber.output_data.
% 
% A private function of class chamber, used by public method chamber.run.

% (c) Miikka Dal Maso & Pauli Simonen 2013
%
% Version history:
% 2013-05-24    0.1.0
% 
% TODO:
% -User defined initial sections in the same way as in run_moving_center.
%  Maybe initials.sections should not be used in this meaning for clarity.


initials=obj.initials;
% Make a logarithmically spaced vector between 10^(Dp_min) and 10^(Dpmax).
% Number of cells is initials.sections.
Dp=logspace(initials.Dp_min,initials.Dp_max,initials.sections);

% Form a lognormal distribution dN/dlogDp where the minimum diameter is
% vector Dp's first value and maximum is Dp's last value. Number of cells
% is the same as Dp's length (i.e. initials.sections).
dNdlogDp = zeros(size(Dp));
for i=1:length(initials.mu)
    dNdlogDp=dNdlogDp + obj.log_normal(Dp,initials.mu(i),initials.sigma(i),initials.N(i));
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
y0 = [initials.Cvap0 N Dp AE_Wall AE_dilu Vap_dilu];

% The basic structure of the aerosol dynamics model is the following:
% there are 2*nSec+4 variables, where nSec is the number of size sections.
% The model stores all data to vector Y which length is 2*nSec+4.
% Basically, the model keeps track of each size section's size (central
% diameter) and concentration. The concentrations of sections are stored
% in Y(2:nSec+1) and the diamters in Y((nSec+2):(2*nSec+1)).
% 
% Additionally, the model keeps track of the condensing vapor 
% concentration (in Y(1)) and the molecules lost to the wall as aerosols
% (Y(2*nSec+2)), diluted as aerosols Y(2*nSec+3) and diluted in the gas
% phase y(2*nSec+4).

% Summary:
% concentrations: Y(2:nSec+1)
% diameters:      Y((nSec+2):(2*nSec+1))
% lost:           Y(2*nSec+2 2*nSec+3 2*nSec+4)

% The model forms an array [t Y] using ode45. The first row of array is
% [t0 y0], that is, the initial values of Y. Ode45 solves the other values
% as a function of time using the function dy = chamberODE that defines
% the change of y in a time step (dy/dt = chamberODE).
% So finally we get [t Y] = |t0  Y0| = Y as a function of time.
%                           |.   . |
%                           |.   . |
%                           |.   . |
%                           |tf  Yf|



% Error tolerance options for ode45:
absTol = ones(size(y0)).*1e-6; % Preallocate the tolerance vector.
absTol(1) = initials.Cvap_tol; % Vapor concentration tolerance
absTol(2:initials.sections+1) = initials.N_tol; % Particle concentration tolerance
absTol(initials.sections+1:(2*initials.sections+1)) = initials.Dp_tol; % Particle diameter tolerance
opts = odeset('absTol',absTol); % Assign tolerances




% Create a visual waitbar (slows the program down a little)
% h = waitbar(0,'0 %');

% Solve Y:
[t Y]= ode45(@chamberODE,initials.tvect,y0,opts,obj);

display('Ode45 finished, processing data...');

% Close waitbar
% close(h);

% This makes a handy structure of the results:
obj.output_data = obj.model_convert(t,Y);


% Now
% output_data.Y           % the Y array
% output_data.distr =     % the put as a "sum-file" (Distribution of particles for
%                                             all time points)
% output_data.vap =       % the vapor concentration 
% output_data.tim =       % the time vector
% output_data.Ntot =      % the total number of particles
% output_data.Vtot =      % the total volume of particles
% output_data.Dpmean=     % the mean diameter of particles
% output_data.Mtot =      % the total mass of particles
% output_data.Mwall=      % mass lost to wall
% output_data.Mdilu=      % mass diluted as aerosol
% output_data.Mvdilu=     % mass diluted as vapor

% toc

function dy = chamberODE(t,y,chamb)
    
dy = zeros(size(y));

initial = chamb.initials;

% load the parameters:
Source = initial.gas_source ; % the condensing vapor source rate (1/cm^3/s) (scalar or 2-column array)
Dilu   = initial.dilu_coeff ; % the dilution coefficient (scalar or 2-column array)
Csat   = initial.satu_conc ; % The condensing vapor saturation concentration
lambda = initial.lambda ; % the condensing vapor mean free path
diffu  = initial.diff_coeff ; % the condensing vapor diffusion coefficient
mv     = initial.vap_molmass ; % the condensing vapor molecular weight
rool   = initial.particle_dens ; % the particle density
alfa   = initial.stick_coeff ; % the sticking coefficient
nSec   = initial.sections; % the number of size sections
T      = initial.T;               % Temperature

Df     = initial.Df;  % Fractal dimension of agglomerates. Used only if agglomeration is on.
r0     = initial.r0;  % Radius of agglomerate primary particles. Used only if agglomeration is on.


CX     = initial.coag_on; % Coagulation switch.

% Get the coagulation mode (coagulation or agglomeration).
% 1 = coagulation
% 0 = agglomeration
coagmode = initial.coag_num;


NA = 6.022e23; % Avogadro constant

% If dilution coefficient is defined as vector, interpolate it to find the
% value of coefficient for the current t.
if initial.dilu_vect_on, 
    Dilu = interp1(Dilu(:,1),Dilu(:,2),t);
end

% If gas_source is defined as vector, interpolate it to find the value of
% gas_source for the current t.
if initial.gas_source_is_vect,
    Source = interp1(Source(:,1),Source(:,2),t);
end


% Make coagulation kernel. Different functions for coagulation and
% agglomeration.
kk=zeros(nSec,length(y((nSec+2):(2*nSec+1)))); % Preallocate
if(coagmode == 1) % coagmode == 1 => particles coagulate.
    for i = 1:nSec,
        kk(i,:) = chamb.koag_kernel(y(nSec+1+i),y((nSec+2):(2*nSec+1)),rool,T).*1e6;
    end
else    % Else coagmode == 0 => particles agglomerate.
    for i = 1:nSec,
        kk(i,:) = chamb.aggl_kernel(y(nSec+1+i),y((nSec+2):(2*nSec+1)),rool,T,Df,r0).*1e6;
    end
end
 
% Show the time evolution in Matlab command window:
time = t

% Or update the visual waitbar:
% perc=round(t/chamb.tvect(end)*100);
% waitbar(perc/100,h,sprintf('%d %%',perc))



% Locations of different values in dy-vector:
% vapor concentration:      dy(1)
% particle concentrations:  dy(2:Nsec+1)
% particle diameters:       dy((NSec+2):(2*nSec+1))
% lost molecules:           dy([2*nSec+2         2*nSec+3             2*nSec+4])
%                               lost to wall     diluted as aerosols  diluted as vapor molecules   


% Increase condensing vapor concentration by vapor coming from Source
dy(1) = dy(1) + Source;

% Dilution of vapor
if initial.dilu_on,
    % Dilute the vapor:
    dy(1) = dy(1)-Dilu.*y(1);
    % And save the information about lost vapor molecules:
    dy(2*nSec+4)  = dy(2*nSec+4)+Dilu.*y(1);
end


% Go through all particle diameters and calculate the effect of dilution,
% coagulation, condensation and sedimentation on particle
% concentrations.
for i = 1:nSec,
    Kn = (2.*lambda)./y(nSec+1+i); % Knudsen number
%     Kn = lambda./y(nSec+1+i);
    betam = (Kn+1)./((0.377.*Kn)+1+(4/(3.*alfa)).*(Kn.^2)+(4/(3.*alfa)).*Kn);
    
    % the flux of molecules to the particle phase    
    I = 2.*pi.*max([y(nSec+1+i) 0]).*1e2.*diffu.*(y(1)-Csat).*betam; %1/s

    % Dilution of particles
    if initial.dilu_on
        dy(i+1) = dy(i+1)-Dilu.*y(i+1);   % Decrease particle concentration
        
        % Calculate the molecules lost to dilution in aerosol phase and
        % save information to dy(2*nSec+3).
        dy(2*nSec+3) = dy(2*nSec+3)+ Dilu.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(nSec+1+i).^3))./(6.*mv);
    end



    % coagulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if CX,    
    % calculate a coagulation matrix
    % this tells how to partition the particles 
    cM = chamb.coagulationMatrix(y((nSec+2):(2*nSec+1)),i);
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
    
    % condensation
    dy(nSec+1+i) = dy(nSec+1+i)+(2.*mv.*I)./(pi.*rool.*y(nSec+1+i).^2.*NA.*1e6); % particle diameter (m/s)
    dy(1) = dy(1) - y(i+1).*I;
    % end condensation
        
    % wall losses and sedimentation according to T. Anttila model...fitted
    % by M. Dal Maso; only usable for SAPPHIR chamber!!
    if initial.sedi_on,
        beta = chamb.sapphir_beta2(y(nSec+1+i),T);
        %beta = 3.5e-5; % 0th order approx
        dy(i+1) = dy(i+1)-beta.*y(i+1);
        
        %molecules lost to the wall = volume lost to wall (in aerosol phase)
        dy(2*nSec+2) = dy(2*nSec+2)+ beta.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(nSec+1+i).^3))./(6.*mv);               
    end
    
    % If vapor concentration is kept constant, reset the value of dy(1) to
    % initial.Cvap0.
    if(initial.Cvap_const == 1)
        dy(1) = 0;
    end
    
    
end % for i

end % function ChamberODE

end