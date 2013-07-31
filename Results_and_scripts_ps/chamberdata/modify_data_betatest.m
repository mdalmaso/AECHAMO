day=19;

% time_start = v(2,1)+(day-18)+0.5832;
time_start = datenum(2009,9,day,0,0,0);

% time_end = v(2,1)+(day-18+1)+0.5832;
time_end = datenum(2009,9,day+1,0,0,0);

time = v(2:end,1);
A=find(time <= time_start);
begin = A(end);
A=find(time >= time_end);
last = A(1);

% begin = begin + (day-18)*width1;

% time = v(begin:begin+width2,1);
time = v(begin:last,1);
% time = (time-v(2,1)-0.5832-(day-18))*24; % Time in hours

% dist = [v(1,:);v(begin:begin+width2,:)];
dist = [v(1,:);v(begin:last,:)];
% dist(2:end,1) = time;

time_a = a(:,1);
% time_a = (time_a-v(2,1)-0.5832-(day-18))*24;
A=find(time_a <= time(1));
ind_first = A(end);
A=find(time_a >= time(end));
ind_last = A(1);
time_a = time_a(ind_first:ind_last);
matrixdata = a(ind_first:ind_last,:);

data.time = time_a;
data.cpc = matrixdata(:,6);
data.CS = matrixdata(:,8);
data.MT_plant = matrixdata(:,12); % ppb
data.inflow = matrixdata(:,14).*(1000/60); % cm^3/s
data.dilu = matrixdata(:,16).*(1000/60); % cm^3/s
data.MT_rc = matrixdata(:,10); % ppb
data.isoprene_plant = matrixdata(:,13); % ppb
min_isoprene = min(matrixdata(:,13));
if(min_isoprene < 0)
    data.isoprene_plant = data.isoprene_plant + abs(min_isoprene);
end


V_rc=1.45e6; % cm^3
NA = 6.022e23; % 1/mol
mol_in_cm3 = 41.6*1e-6; % mol/cm^3
inflow_mols = data.inflow.*mol_in_cm3; % cm^3/s*mol/cm^3 = mol/s
N_tot = inflow_mols.*NA; % mol/s*1/mol = 1/s Includes all molecules
N_mt = N_tot.*data.MT_plant./1e9; % 1/s
N_mt_rc = data.MT_rc.*2.505152e10; % ppb*1/cm3 = 1/cm3
N_isoprene = N_tot.*data.isoprene_plant./1e9;

% alfa = 0.83;
alfa = 0.35;
beta = 0.45;

% alfa=0.5;
cond_vap_rc = alfa.*N_mt_rc; % Mitattu MT-konsentr. kertaa alfa.

Q1=alfa.*N_mt./V_rc; % 1/s*1/cm^3 = 1/(cm^3s)
Q2=beta.*N_mt./V_rc;

% test:
% Q=alfa.*(N_mt+N_isoprene)./V_pc;
% Q=(alfa.*N_mt + beta.*N_isoprene)./V_pc;

dilu = data.dilu./V_rc;

Dps=dist(1,3:end);
% 
% dist_fixed=dist;
% dist_fixed(:,1)=dist_fixed(:,1).*3600;

for i=begin:last
    Vtot(i-begin+1)=trapz(log10(Dps),v(i,3:end).*(pi/6.*v(1,3:end).^3));
end
Ntot = v(begin:last,2);


clear time_start;
clear time_end;
clear begin;
clear last;
clear A;
clear ind_first;
clear ind_last;




