function [C] = gauss_puff(Q,r,x,y,z,time,H,u)
%GAUSS_PUFF Summary of this function goes here
%   Detailed explanation goes here

% L‰hde: Stockie: The Mathematics of Atmospheric Dispersion Modelling

%  [Q] = kg, hiukkasm‰‰r‰ alkuhetkell‰ l‰hdepisteess‰.

% Testiarvoja
% x=0:10
% y=-5:5
% z=0:10
% u=1
% K=1
% Q=10
% time=0:20
% H=3
% r(i)=K/u*x(i) %Vakio K


% C = [C, time, r, y, z]
C = zeros(length(time),length(r),length(y),length(z));

for t=1:length(time)
    for i=1:length(r)
        for j=1:length(y)
            for k=1:length(z)
                C(t,i,j,k) = gauss(Q,r(i),x(i),y(j),z(k),time(t),H,u);
            end
        end
    end
end

% Plot the view from top

% timeindex = floor(length(time)/2);
timeindex = 1;
z_index = 3;
[X, Y] = meshgrid(x,y);
intensity = zeros(size(X));
intensity(:,:) = squeeze(C(timeindex,:,:,z_index))';
figure;
pcolor(X,Y,intensity);

timeindex = 11;
z_index = 3;
[X, Y] = meshgrid(x,y);
intensity = zeros(size(X));
intensity(:,:) = squeeze(C(timeindex,:,:,z_index))';
figure;
pcolor(X,Y,intensity);

% timeindex = 20;
% z_index = 4;
% [X, Y] = meshgrid(x,y);
% intensity = zeros(size(X));
% intensity(:,:) = squeeze(C(timeindex,:,:,z_index))';
% figure;
% pcolor(X,Y,intensity);
end

function [C] = gauss(Q,r,x,y,z,time,H,u)

C = Q/(8*(pi*r)^(3/2))*exp(-((x-u*time)^2+y^2)/(4*r))*(exp(-(z-H)^2/(4*r))+exp(-(z+H)^2/(4*r)));

end


