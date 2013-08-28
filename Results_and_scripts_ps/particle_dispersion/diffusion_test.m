function [ out ] = diffusion_test(N0,num_of_circles, R, sections, time, interval)
%DIFFUSION_TEST Simulates particle dispersion and simultaneous processes
% Simulaatio tekee num_of_circles määrän sisäkkäisiä palloja. Uloimman
% pallon säde on R (m), ja muiden pallojen säde määräytyy siten, että
% jokaisen pallokuoren sisältämä tilavuus on sama.
% 
% Alussa kaikki hiukkaset ovat sisimmässä pallossa, josta ne leviävät
% Fickin lain mukaisesti ulompiin. N0 tarkoittaa tässä
% hiukkaskonsentraatiota (1/cm^3).
% 
% Jokaisessa pallossa on 'sections' määrä sektioita, joilla simuloidaan
% dynaamisia prosesseja. Simulaatio ajetaan 'interval'-muuttujan määrämissä
% jaksoissa. Jos interval on 2, ensin ajetaan 2 sekuntia kammiomallia,
% jonka jälkeen 2 sekuntia diffuusiota. Tätä jatketaan loppuun asti.
%
% Diffuusio vaikuttaa liian hitaalta, kun käytetään realistisia arvoja
% diffuusiokertoimelle (5.4e-8 10nm hiukkaselle). Nyt käytössä oleva
% kerroin on 0.9 m^2/s.
% 
% Testaa esim. out = diffusion_test(5e7,12,18,15,0:2:160,2);
% ja sen jälkeen diffusion_script.m:llä voi katsoa tuloksia. Tämä
% simulaatio kestää n. 900 sekuntia.

tic
kam(num_of_circles) = chamber;
kam(1).initialize('N',[N0],...
                  'mu', [5e-9],...
                  'sigma', [1.2]);
kam(1).initialize('fixed_sections', 1,...
                  'dilu_on', 0,...
                  'sections', sections,...
                  'coag_mode', 'aggl');

for i=2:length(kam)
    kam(i).initialize('N',0);
    kam(i).initialize('fixed_sections', 1,...
                  'dilu_on', 0,...
                  'sections', sections,...
                  'coag_mode', 'aggl');
end

r = zeros(length(num_of_circles),1);
V = zeros(length(num_of_circles),1);
r(1) = R/sqrt(num_of_circles);
V(1) = 4/3.*pi.*r(1)^3;
for i=2:num_of_circles
    r(i) = i^(1/3)*r(1);
    V(i) = 4/3*pi*r(i)^3-sum(V(1:i-1));
end


y0 = [N0*V(1), zeros(1,num_of_circles-1), R];


Y = zeros(length(time),2*sections+5,num_of_circles);

Ns = zeros(length(kam), sections);
Dps = zeros(length(kam),sections);
Ntot = zeros(1,length(kam));
Ns_deleted = zeros(length(kam),sections);
h2 = waitbar(0,'0 %','Name','Total progress','Units', 'normalized', 'Position', [0.5 0.4 0.25 0.2]);

for k=time(1):interval:time(end)
    
    perc=round(k/time(end)*100);
    waitbar(perc/100,h2,sprintf('%d %%',perc))

    for i=1:length(kam)
        kam(i).initialize('tvect',k:k+interval);
        [~,Ytemp] = kam(i).run_moving_center;
        Y(k+1:k+interval,:,i) = Ytemp(1:end-1,:);
        Ns(i,:) = Ytemp(end,2:sections+1);
        Dps(i,:) = Ytemp(end,sections+2:2*sections+1);
        Ntot(i) = sum(Ns(i,:));
        y0(i) = Ntot(i)*V(i)*1e6;
    end
    
    [~, N] = ode45(@diffusion, k:k+interval, y0);
    
    difference = Ntot-N(end,1:end-1)./(V.*1e6);
        
    for i=1:length(kam)
        Ns_deleted(i,:) = Ns(i,:)./Ntot(i).*difference(i);
        if(difference(i) > 0)
            Ns(i,:) = Ns(i,:)-Ns_deleted(i,:);
        end
        
        if((i > 1) && (difference(i-1) > 0))
            [Dps(i,:), Ns(i,:)] = add_particles(Dps(i-1,:), Ns_deleted(i-1,:), Dps(i,:), Ns(i,:));
        end
        kam(i).initialize('number_distr',Ns(i,:),'center_diameters',Dps(i,:));
    end
end

close(h2);

for i=1:length(kam)
    kam(i).output_data = kam(i).model_convert(time,Y(:,:,i));
end

out.kam = kam;
out.R = R;
out.r = r;
out.Ntot0 = N0;
toc
end

function [dy] = diffusion(t, y)
% Fick's law.

num_of_spheres = length(y)-1;
dy=zeros(length(y),1);

diff_coeff = 9e-1;

R = y(num_of_spheres+1);

% The radii of spheres is such that every sphere has equal area.
r(1) = R/((num_of_spheres)^(1/3));
V(1) = 4/3.*pi.*r(1)^3;

for i=2:num_of_spheres
    r(i) = i^(1/3)*r(1);
   
    V(i) = 4/3*pi*r(i)^3-sum(V(1:i-1));
end

% A is the area of sphere shell and the area through which the particle
% flux goes.
A = 4.*pi.*r.^2;

% Fick: J = -D*dfii/dx.
% dfii = particles/m^3, particle concentration difference between two
% sphere shells. dx is the average distance between two sphere shells.
% dy = J*A, where A is the area through which the particle flux, J, goes.
J = zeros(1,num_of_spheres-1);
for i=1:num_of_spheres-1
    dfii = y(i+1)/V(i+1)-y(i)/V(i);

    if(i == 1)
        distance = r(1)+(r(2)-r(1))/2;
    else
        distance = .5.*(r(i+1)-r(i-1));
    end
    dfiidx = dfii/distance;
    J(i) = -diff_coeff*dfiidx;
    
    if(i>1)
        dy(i) = -J(i)*A(i) + J(i-1)*A(i-1);
    else
        dy(i) = -J(i)*A(i);
    end
end

% At the boundary, dfii = -(concentration of last circle). This is obtained
% by assuming that the concentration outside the outer circle is zero.
% dfii_bound = -y(num_of_spheres)/A(num_of_spheres);
dfii_bound = -y(num_of_spheres)/V(num_of_spheres);
distance_bound = .5.*(r(num_of_spheres)-r(num_of_spheres-2));
J_bound = -diff_coeff*dfii_bound/distance_bound;
dy(num_of_spheres) = -J_bound*A(end) + J(end)*A(end-1);
end

function [Dp_new, N_new] = add_particles(Dp1, N1, Dp2, N2)
% Calculates the new diameter and N vector when N1 particles with
% corresponding diameters Dp1 are moved to sections of diameter Dp2 and
% particle concentration N2.
% It's assumed that the diameter vectors are initially identical, so that
% they are fixed and their limits are identical. However, the center
% diameters of sections can differ.

for i=1:length(Dp1)
    if(N1(i) > 0)
        if(N2(i) <= eps)
            Dp_new(i) = Dp1(i);
        else
            Ntot = N1(i)+N2(i);
            v1 = pi/6*Dp1(i)^3*N1(i);
            v2 = pi/6*Dp2(i)^3*N2(i);
            vtot = (v1+v2)/Ntot;
            Dp_new(i) = (6/pi*vtot)^(1/3);
        end
    else
        Dp_new(i) = Dp2(i);
    end
    N_new(i) = N1(i)+N2(i);
end  
    
    
end

