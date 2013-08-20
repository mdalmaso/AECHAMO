%%
kam=out.kam;
Ntot=zeros(length(kam(1).output_data.tim),1);
figure;
for i=1:length(kam)
    plot(kam(i).output_data.tim, kam(i).output_data.Ntot, 'Linewidth', 2);
    hold on;
    Ntot = Ntot + kam(i).output_data.Ntot;
end

plot(kam(1).output_data.tim, Ntot,'--','LineWidth', 2);

%%
plot(kam(1).output_data.tim,kam(1).output_data.Ntot,'LineWidth',2);
hold on;
plot(kam(2).output_data.tim,kam(2).output_data.Ntot,'r','LineWidth',2);
hold on;
plot(kam(3).output_data.tim,kam(3).output_data.Ntot,'m','LineWidth',2);
Nsum = kam(1).output_data.Ntot + kam(2).output_data.Ntot + kam(3).output_data.Ntot;
plot(kam(1).output_data.tim, Nsum,'--','LineWidth',2);
leg=legend('Section 1', 'Section 2', 'Section 3', 'Total');
kam(1).plot('dist','original');
kam(2).plot('dist','original');
kam(3).plot('dist','original');


%%

for i=1:length(kam)
    kam(i).plot('dist','original');
    caxis([1 7]);
    gcc=colorbar('horiz');
    set(gcc,'xlim',[1 7],'xtick',[1 2 3 4 5 6 7],'fontsize',10)
    set(gcc,'xticklabel',[10,100,1000,10000 100000, 1000000, 10000000]')
end


%%
kam=out.kam;
r = out.r;
Ntot0 = out.Ntot0;
time = kam(1).output_data.tim;
for t=1:50:length(time)
    figure;
    for i=length(r):-1:1
        intensity = kam(i).output_data.Ntot(t)/Ntot0;
        filledCircle([0,0],r(i),1000,[1-intensity 1-intensity 1-intensity]);
        hold on;
    end
end
