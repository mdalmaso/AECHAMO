%%
% Plottaa hiukkasten kokonaism‰‰r‰n ajan funktiona (katkoviiva) sek‰
% jokaisen pallon sis‰lt‰m‰n hiukkasm‰‰‰r‰n (yhtein‰iset viivat).
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
% Plottaa jokaisen pallon hiukkasjakauman ajan funktiona eri kuviin
kam = out.kam;
for i=1:length(kam)
    kam(i).plot('dist','original');
    caxis([1 7]);
    gcc=colorbar('horiz');
    set(gcc,'xlim',[1 7],'xtick',[1 2 3 4 5 6 7],'fontsize',10)
    set(gcc,'xticklabel',[10,100,1000,10000 100000, 1000000, 10000000]')
end


%%
% Plottaa pallot ylh‰‰lt‰p‰in katsottuna siten, ett‰ pallon (tai 
% pallokuoren) sis‰lt‰m‰ hiukkasm‰‰r‰ ilmaistaan harmaan s‰vyn‰. Musta =
% paljon hiukkasia, valkoinen = ei yht‰‰n hiukkasia. Tilanne plotataan
% time_interval v‰lein eri kuviin. Vaatii filledCircle-funktion 
% (Draw a filled circle, file exchange).
time_interval = 50;
kam=out.kam;
r = out.r;
Ntot0 = out.Ntot0;
time = kam(1).output_data.tim;
for t=1:time_interval:length(time)
    figure;
    for i=length(r):-1:1
        intensity = kam(i).output_data.Ntot(t)/Ntot0;
        filledCircle([0,0],r(i),1000,[1-intensity 1-intensity 1-intensity]);
        hold on;
    end
end
