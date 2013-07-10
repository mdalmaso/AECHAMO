clear Dps;
clear Ns;
clear gsd;
clear v_avg;

secs = chamb.initials.sections;

Dps=chamb.output_data.Y(:,secs+2:2*secs+1);
Ns=chamb.output_data.Y(:,2:secs+1);
for i=1:length(Dps(:,1))
gsd(i)=geostd(Dps(i,:),Ns(i,:));
end

 v_avg=chamb.output_data.Vtot./chamb.output_data.Ntot;
 