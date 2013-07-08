clear all;
load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\Batch\10nm\run_20130705T204109.mat')

for i = 1:32
N(i) = chamb(i).output_data.Ntot(end);
end
for i = 1:32
N(i) = chamb(i).output_data.Ntot(189);
end
load('K:\603_L\60304\Users\Poikkimäki\GitHub\AECHAMO\Results and scripts_mp\SOA formation\Batch\80nm_isoM\run_20130708T141228.mat')
for i = 1:32
N(i+32) = chamb(i).output_data.Ntot(117);
end
sum(N)/64
ans =

   1.4599e+04

median(N)

ans =

   1.0737e+04
