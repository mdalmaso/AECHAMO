function [ out ] = gsd( chamb )
%GSD Summary of this function goes here
%   Detailed explanation goes here

secs = chamb.initials.sections;

Dps=chamb.output_data.Y(:,secs+2:2*secs+1);
Ns=chamb.output_data.Y(:,2:secs+1);
Ns(find(Ns < 0)) = 0;
for i=1:length(Dps(:,1))
out(i)=geostd(Dps(i,:),Ns(i,:));
end



end

