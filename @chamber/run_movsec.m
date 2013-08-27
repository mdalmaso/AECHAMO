function [t,Y] = run_movsec(obj)
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

% Get the initial values:
initials=obj.initials;
initials.sections = obj.sections;

% And diameter vector:
Dp = obj.Dps;

% And the particle number distribution. N(i) represents the particle
% concentration in section i, and the initial diameter of section i is
% Dp(i).
N = obj.number_distribution;

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
absTol(2:obj.sections+1) = initials.N_tol; % Particle concentration tolerance
absTol(obj.sections+1:(2*obj.sections+1)) = initials.Dp_tol; % Particle diameter tolerance
opts = odeset('absTol',absTol); % Assign tolerances

opts = odeset(opts, 'NonNegative', 1:obj.sections+1);

if(initials.max_timestep)
    opts = odeset(opts, 'MaxStep', initials.max_timestep);
end


% Create a visual waitbar (slows the program down a little, a few seconds)
h = waitbar(0,'0 %','Name','Running simulation...');

% Solve and return Y:
% [t, Y]= ode45(@chamberODE,initials.tvect,y0,opts,obj);
[t, Y]= ode45(@chamberODE,initials.tvect,y0,opts);

% Close waitbar:
close(h);

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

% function dy = chamberODE(t,y,chamb)
function dy = chamberODE(t,y)
        
dy = zeros(size(y));

% initial = chamb.initials;

% load the parameters:
Source = initials.gas_source ; % the condensing vapor source rate (1/cm^3/s) (scalar or 2-column array)
Dilu   = initials.dilu_coeff ; % the dilution coefficient (scalar or 2-column array)
Csat   = initials.satu_conc ; % The condensing vapor saturation concentration
lambda = initials.lambda ; % the condensing vapor mean free path
diffu  = initials.diff_coeff ; % the condensing vapor diffusion coefficient
mv     = initials.vap_molmass ; % the condensing vapor molecular weight
rool   = initials.particle_dens ; % the particle density
alfa   = initials.stick_coeff ; % the sticking coefficient
nSec   = obj.sections; % the number of size sections
T      = initials.T;               % Temperature

Df     = initials.Df;  % Fractal dimension of agglomerates. Used only if agglomeration is on.
r0     = initials.r0;  % Radius of agglomerate primary particles. Used only if agglomeration is on.


CX     = initials.coag_on; % Coagulation switch.

% Get the coagulation mode (coagulation or agglomeration).
% 1 = coagulation
% 0 = agglomeration
coagmode = initials.coag_num;


NA = 6.022e23; % Avogadro constant

% If dilution coefficient is defined as vector, interpolate it to find the
% value of coefficient for the current t.
if initials.dilu_vect_on, 
    Dilu = interp1(Dilu(:,1),Dilu(:,2),t);
end

% If gas_source is defined as vector, interpolate it to find the value of
% gas_source for the current t.
if initials.gas_source_is_vect,
    Source = interp1(Source(:,1),Source(:,2),t);
end


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
 
% Show the time evolution in Matlab command window:
% time = t

% Or update the visual waitbar:
perc=round(t/initials.tvect(end)*100);
waitbar(perc/100,h,sprintf('%d %%',perc))



% Locations of different values in dy-vector:
% vapor concentration:      dy(1)
% particle concentrations:  dy(2:Nsec+1)
% particle diameters:       dy((NSec+2):(2*nSec+1))
% lost molecules:           dy([2*nSec+2         2*nSec+3             2*nSec+4])
%                               lost to wall     diluted as aerosols  diluted as vapor molecules   


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
for i = 1:nSec,
    Kn = (2.*lambda)./y(nSec+1+i); % Knudsen number
%     Kn = lambda./y(nSec+1+i);
    betam = (Kn+1)./((0.377.*Kn)+1+(4/(3.*alfa)).*(Kn.^2)+(4/(3.*alfa)).*Kn);
    
    % the flux of molecules to the particle phase    
    I = 2.*pi.*max([y(nSec+1+i) 0]).*1e2.*diffu.*(y(1)-Csat).*betam; %1/s

    % Dilution of particles
    if initials.dilu_on
        dy(i+1) = dy(i+1)-Dilu.*y(i+1);   % Decrease particle concentration
        
        % Calculate the molecules lost to dilution in aerosol phase and
        % save information to dy(2*nSec+3).
        dy(2*nSec+3) = dy(2*nSec+3)+ Dilu.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(nSec+1+i).^3))./(6.*mv);
    end



    % coagulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if CX,    
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
    
    % condensation
    dy(nSec+1+i) = dy(nSec+1+i)+(2.*mv.*I)./(pi.*rool.*y(nSec+1+i).^2.*NA.*1e6); % particle diameter (m/s)
    dy(1) = dy(1) - y(i+1).*I;
    % end condensation
        
    % wall losses and sedimentation according to T. Anttila model...fitted
    % by M. Dal Maso; only usable for SAPPHIR chamber!!
    if initials.sedi_on,
        beta = obj.sapphir_beta2(y(nSec+1+i),T);
        %beta = 3.5e-5; % 0th order approx
        dy(i+1) = dy(i+1)-beta.*y(i+1);
        
        %molecules lost to the wall = volume lost to wall (in aerosol phase)
        dy(2*nSec+2) = dy(2*nSec+2)+ beta.*y(i+1).*(NA.*1e6.*rool.*pi.*(y(nSec+1+i).^3))./(6.*mv);               
    end
    
    % If vapor concentration is kept constant, reset the value of dy(1) to
    % initials.Cvap0.
    if(initials.Cvap_const == 1)
        dy(1) = 0;
    end
    
    
end % for i

end % function ChamberODE

end