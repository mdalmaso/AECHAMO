

d = dir('*.mat');
cols = 'rgbcmykrgbcmyk'
for i = 1:length(d),
    clear a N50
    load(d(i).name);
    l{i} = d(i).name;
    for t = 1:length(a.output_data.tim),
        N50(t) = integrate_distribution(a.output_data.distr_original(t+1,(end-2)/2+3:end),a.output_data.distr_original(t+1,3:(end-2)/2+2),50e-9,1e-6);
    end
    figure(1)
    hold on
    plot(a.output_data.tim,N50,[cols(i)])
end
% legend(l)