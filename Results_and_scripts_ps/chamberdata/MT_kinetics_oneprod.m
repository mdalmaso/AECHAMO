function[out] = MT_kinetics_oneprod(day)
% test new theory
close all
% dilu_gamma = dilution flow

load smps.mat % into v
load allmatrixdata_new_gc.mat
% into aa

st_time = datenum(2009,9,day);
en_time = datenum(2009,9,day+1);


ppb = M_dens(1013,288).*1e-9;

ix = find(aa(:,1)>st_time & aa(:,1)<en_time);
time = aa(ix,1);
RC = aa(ix,10);
PC = aa(ix,12);
CS = aa(ix,8);
MCPC = aa(ix,7);
outflow = aa(ix,16);
PC_flow = aa(ix,14);


figure(1)
subplot(311)
plot(time,PC,'g.')
hold on
plot(time,RC,'b.')
datetick('x','keeplimits');
hold off
legend('Plant','Recation')
subplot(312)
plot(time,CS,'r.')
datetick('x','keeplimits');
subplot(313)
plot(time,MCPC,'k.')
datetick('x','keeplimits');

clear ix time RC PC CS MCPC outflow PC_flow
% refine the time selection
[x y] = ginput(2);

ix = find(aa(:,1)>x(1) & aa(:,1)<x(2));
time = aa(ix,1);
RC = aa(ix,10);
PC = aa(ix,12);
CS = aa(ix,8);
MCPC = aa(ix,7);
outflow = aa(ix,16);
PC_flow = aa(ix,14);
O3_flow = aa(ix,15);
O3 = aa(ix,3);
OH = aa(ix,9);

figure(1)
subplot(311)
plot(time,PC,'g.')
hold on
plot(time,RC,'b.')
datetick('x','keeplimits');
hold off
legend('Plant','Recation')
subplot(312)
plot(time,O3,'r.')
datetick('x','keeplimits');
subplot(313)
plot(time,MCPC,'k.')
datetick('x','keeplimits');

[UVx yy] = ginput(2);

UV = zeros(size(time));
UVix = find(time>UVx(1) & time<UVx(2));
UV(UVix) = 1;
subplot(311)
hold on
plot(time,UV,'y*-')
waitforbuttonpress


odepara.tim = (time-time(1)).*(24.*3600);
odepara.O3 = O3;
odepara.O3_flow = O3_flow;
odepara.PC_flow = PC_flow;
odepara.PC = PC;
odepara.outflow = outflow;
odepara.UV = UV;
odepara.ppb = ppb
odepara.OH = OH;



tspan = odepara.tim;
y0 = RC(1).*ppb; 

[T Y] = ode45(@(t,y) kinefun(t,y,odepara),tspan,y0);




figure(2)
plot(odepara.tim,RC.*ppb,'b.')
hold on
plot(T,Y,'r.')
plot(odepara.tim,PC.*ppb,'g.')
plot(odepara.tim,UV.*5e10,'y-')
plot(odepara.tim,O3.*ppb./10,'k.')


hold off


kO3 = 15e-17;
kOH = 6e-11;

figure(3)
plot(T,Y(:).*OH(:).*kOH,'r*')
hold on
plot(T,Y(:).*O3(:).*ppb.*kO3,'k*')
plot(T,UV.*2e7+1e4,'y.-')
hold off

legend('OH reaction','O3 reaction','UV')
imgnam = sprintf('SOURCES_%02i.png',day)
title(imgnam)
set(gca,'yscale','log')
export_fig(3, imgnam,'-transparent')



out.meas.time = time;
out.meas.PC = PC.*ppb;
out.meas.RC = RC.*ppb;
out.meas.O3 = O3.*ppb;
out.ppb = ppb;
out.meas.UV = UV;
out.meas.OH = OH; 
out.model.RC_MT = Y(:);
out.model.time= T./(24.*3600)+time(1);

end



















function dydt = kinefun(t,y,odepara)

t
iO3 = interp1(odepara.tim,odepara.O3,t)
iO3_flow = interp1(odepara.tim,odepara.O3_flow,t);
iPC_flow = interp1(odepara.tim,odepara.PC_flow,t);
iPC = interp1(odepara.tim,odepara.PC,t);
ioutflow = interp1(odepara.tim,odepara.outflow,t);
iUV = interp1(odepara.tim,odepara.UV,t);

iOH = interp1(odepara.tim,odepara.OH,t);


% inflow of MT
Q = iPC_flow./60./1460.*iPC.*odepara.ppb % molecules/second;

kO3 = 15e-17;
% MT + O3 reaction rate
rO3 = iO3.*odepara.ppb.*kO3;





% MT + OH reaction rate
if iUV>0,
    if isfinite(iOH),
        cOH = iOH;
    else
        cOH = max(odepara.OH);
        if ~isfinite(cOH)
            cOH = 5e7;
        end
    end
else
    cOH = iO3.*1e-6.*odepara.ppb
    cOH = 0;
end

kOH = 6e-11;
rOH = cOH.*kOH;

dilu = ioutflow./60./1460;

dydt = Q - rOH.*y - rO3.*y - dilu.*y;

end














