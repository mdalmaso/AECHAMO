P0 = 100e9.*2.6908e19;
t = 0:600:86400;
P = P0.*exp(-2e-16*60e9*2.6908e19.*t);