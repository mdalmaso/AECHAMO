alfa = 0.38;
gamma = kam.initials.vap_wallsink;
sections = kam.initials.sections;

delta_M=kam.output_data.Mtot'+kam.output_data.Mdilu+kam.output_data.Mwall;
for i=2:length(source_time)
    delta_P(i) = trapz(source_time(1:i),Q_MTOH(1:i).*kam.initials.vap_molmass./6.022e23);
end
delta_P_interp = interp1(source_time, delta_P, kam.output_data.tim);
% delta_P = trapz(source_time,Q_MTOH.*kam.initials.vap_molmass./6.022e23);
% delta_P = trapz(timehours(begin:last).*3600,Q(begin:last)./alfa.*kam.initials.vap_molmass./6.022e23);
Yield1 = delta_M./delta_P_interp;
for i=1:length(kam.output_data.Y(:,1))
    CS(i)=CS_general(kam.output_data.Y(i,2+sections:1+2*sections),kam.output_data.Y(i,2:1+sections),273.15+16,1.0);
%     CS(i)=CS(i)+kam.initials.dilu_coeff(i,2);
end
gamma_2=kam.initials.dilu_coeff(:,2)+gamma;
Yield2=alfa./(1+gamma_2'./CS);