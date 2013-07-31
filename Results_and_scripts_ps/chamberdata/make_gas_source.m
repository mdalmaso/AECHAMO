kO3=15e-17;
kOH=6e-11;
Q_MTOH = kOH.*res.model.RC_MT.*res.meas.OH;
Q_MTO3 = kO3.*res.model.RC_MT.*res.meas.O3;
Q_MTOH(isnan(Q_MTOH))=0;
Q_MTO3(isnan(Q_MTO3))=0;
min(Q_MTO3)
source_tvect=datevec(res.meas.time);
source_time=3600.*source_tvect(:,4)+60.*source_tvect(:,5)+source_tvect(:,6);