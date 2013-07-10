
day=18;

time_start = v(2,1)+(day-18)+0.5832;
time_end = v(2,1)+(day-18+1)+0.5832;
time = v(2:end,1);
A=find(time <= time_start);
begin = A(end);
A=find(time >= time_end);
last = A(1);

% begin = begin + (day-18)*width1;

% time = v(begin:begin+width2,1);
time = v(begin:last,1);
time = (time-v(2,1)-0.5832-(day-18))*24; % Time in hours

% dist = [v(1,:);v(begin:begin+width2,:)];
dist = [v(1,:);v(begin:last,:)];
dist(2:end,1) = time;

time_a = a(:,1);
time_a = (time_a-v(2,1)-0.5832-(day-18))*24;
A=find(time_a <= time(1));
ind_first = A(end);
A=find(time_a >= time(end));
ind_last = A(1);
time_a = time_a(ind_first:ind_last);
matrixdata = a(ind_first:ind_last,:);

data.time = time_a;
data.cpc = matrixdata(:,6);
data.CS = matrixdata(:,8);
data.MT_plant = matrixdata(:,12); % dimensionless
data.inflow = matrixdata(:,14).*(1000/60); % cm^3/s
data.dilu = matrixdata(:,16).*(1000/60);

V_pc=1.46e6; % cm^3
NA = 6.022e23; % 1/mol
mol_in_cm3 = 41.6*1e-6; % mol/cm^3
inflow_mols = data.inflow.*mol_in_cm3; % cm^3/s*mol/cm^3 = mol/s
N_tot = inflow_mols.*NA; % mol/s*1/mol = 1/s Includes all molecules
N_mt = N_tot./(1./data.MT_plant + 1); % 1/s

alfa = 0.3;
Q=alfa.*N_mt./V_pc; % 1/s*1/cm^3 = 1/(cm^3s)
dilu = data.dilu./V_pc;

clear time_start;
clear time_end;
clear begin;
clear last;
clear A;
clear ind_first;
clear ind_last;




