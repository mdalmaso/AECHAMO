
dist1=kam.output_data.distr_original;
dist2=kam2.output_data.distr_original;
dist3=kam3.output_data.distr_original;

%% Puolivälin jakauma
figure;
plot(dist2(floor(end/2),33:end),dist2(floor(end/2),3:32),'rx-');
set(gca,'xscale','log');
hold on;
plot(dist1(floor(end/2),33:end),dist1(floor(end/2),3:32),'bx-');

plot(dist3(floor(end/2),33:end),dist3(floor(end/2),3:32),'mx-');

%% Lopuujakauma
figure;
plot(dist1(end,33:end),dist1(end,3:32),'bx-');
set(gca,'xscale','log');
hold on;
plot(dist2(end,33:end),dist2(end,3:32),'rx-');

plot(dist3(end,33:end),dist3(end,3:32),'mx-');