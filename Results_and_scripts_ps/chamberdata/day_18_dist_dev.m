sections = kam.initials.sections;
d_orig=kam.output_data.distr_original;
% plot(dist(1,3:end),dist(109,3:end),'bx-');
% set(gca,'xscale','log');
% hold on;
% plot(d_orig(45,sections+3:end),d_orig(45,3:sections+2),'rx-');
% 
% plot(dist(1,3:end),dist(117,3:end),'bx-');

plot(d_orig(45,sections+3:end),d_orig(45,3:sections+2),'rx-');
hold on;
plot(dist(1,3:end),dist(109,3:end),'bx-');
set(gca,'xscale','log');
figure('Color',[1 1 1]);
plot(dist(1,3:end),dist(117,3:end),'bx-');
set(gca,'xscale','log');
hold on;
plot(d_orig(45,sections+3:end),d_orig(130,3:sections+2),'rx-');
figure('Color',[1 1 1]);
plot(dist(1,3:end),dist(148,3:end),'bx-');
set(gca,'xscale','log');
hold on;
plot(d_orig(45,sections+3:end),d_orig(336,3:sections+2),'rx-');
