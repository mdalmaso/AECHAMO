ajoja = 3;

for i = 1:ajoja
    chamb(i).plot
end

CS = zeros(2881,ajoja);
for i = 1:ajoja
    CS(1:end,i) = CS_tot(chamb(i).output_data.distr);
end
