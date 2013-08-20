function [ out ] = diffusion_test(N0,num_of_circles, R, sections, time, interval)
%DIFFUSION_TEST Summary of this function goes here
%   Detailed explanation goes here
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

for k=time(1):interval:time(end)
    
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
        kam(i).initialize('distr',Ns(i,:),'center_diameters',Dps(i,:));
    end
end

for i=1:length(kam)
    kam(i).output_data = kam(i).model_convert(time,Y(:,:,i));
end

out.kam = kam;
out.R = R;
out.r = r;
out.Ntot0 = N0;
end

function [dy] = diffusion(t, y)

num_of_sects = length(y)-1;
dy=zeros(length(y),1);

diff_coeff = 8e-2;

R = y(num_of_sects+1);

% The radii of circles is such that every circle has equal area.
% r(1) = R/sqrt(num_of_sects);
% A(1) = pi*r(1)^2;

r(1) = R/((num_of_sects)^(1/3));
V(1) = 4/3.*pi.*r(1)^3;

% for i=2:num_of_sects
%     r(i) = sqrt(i)*r(1);
%    
%     % Area of circle is: pi*radius^2 - (area of inner circles).
%     A(i) = pi*r(i)^2-sum(A(1:i-1));
% end

for i=2:num_of_sects
    r(i) = i^(1/3)*r(1);
   
    V(i) = 4/3*pi*r(i)^3-sum(V(1:i-1));
end

% % l is the circumference of circles.
% l=2.*pi.*r;

A = 4.*pi.*r.^2;


for i=1:num_of_sects-1
%     dN = y(i+1)/A(i+1)-y(i)/A(i);
    dN = y(i+1)/V(i+1)-y(i)/V(i);

    if(i == 1)
        distance = r(1)+(r(2)-r(1))/2;
    else
        distance = .5.*(r(i+1)-r(i-1));
    end
    
    J(i) = -diff_coeff*dN/distance;
    
    if(i>1)
%         dy(i) = -J(i)*l(i)+J(i-1)*l(i-1);
        dy(i) = -J(i)*A(i) + J(i-1)*A(i-1);
    else
%         dy(i) = -J(i)*l(i);
        dy(i) = -J(i)*A(i);
    end
end

% At the boundary, dN = -(concentration of last circle). This is obtained
% by assuming that the concentration outside the outer circle is zero.
% dN_bound = -y(num_of_sects)/A(num_of_sects);
dN_bound = -y(num_of_sects)/V(num_of_sects);
distance_bound = .5.*(r(num_of_sects)-r(num_of_sects-2));
J_bound = -diff_coeff*dN_bound/distance_bound;
% dy(num_of_sects) = -J_bound*l(end) + J(end)*l(end-1);
dy(num_of_sects) = -J_bound*A(end) + J(end)*A(end-1);
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
% 
% 

