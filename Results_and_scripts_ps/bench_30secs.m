% Benchmark using 30 sections

%% mu 1:
# 1
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [1e3 1e2];
mu=[50e-9, 140e-9];
dilu_on = 0;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 1.4e7;

gas_source = 0;

# 2
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [1e3 1e2];
mu=[50e-9, 140e-9];
dilu_on = 1;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 1.4e7;

gas_source = 0;


# 3
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [1e3 1e2];
mu=[50e-9, 140e-9];
dilu_on = 1;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e7;

gas_source = 0;

# 4
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [1e3 1e2];
mu=[50e-9, 140e-9];
dilu_on = 0;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e7;

gas_source = 0;

%% Change the N values

# 5
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [5e4 3e3];
mu=[50e-9, 140e-9];
dilu_on = 0;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 1.4e7;

gas_source = 0;

# 6
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [5e4 3e3];
mu=[50e-9, 140e-9];
dilu_on = 1;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 1.4e7;

gas_source = 0;

# 7
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [5e4 3e3];
mu=[50e-9, 140e-9];
dilu_on = 1;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e7;

gas_source = 0;

# 8
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.3];

N = [5e4 3e3];
mu=[50e-9, 140e-9];
dilu_on = 0;
sections = 30;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e7;

gas_source = 0;

#