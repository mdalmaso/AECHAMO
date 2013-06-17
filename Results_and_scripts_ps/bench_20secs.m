% Benchmark using 20 sections

%% Part 1: N values are constant
%% Part 1.1: Cvap0, sigma and mu are constants
# 1.11
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
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 2e7;

gas_source = 0;

# 1.12
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
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 2e7;

gas_source = 0;

%% Part 1.2: change the Cvap0, sigma and mu
# 1.21
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.4];

N = [1e3 1e2];
mu=[40e-9, 200e-9];
dilu_on = 0;
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e8;

gas_source = 0;

# 1.22
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.4];

N = [1e3 1e2];
mu=[40e-9, 200e-9];
dilu_on = 1;
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e8;

gas_source = 0;

%% Part 2: change the N values
%% Part 2.1: Cvap0, sigma and mu are constants
# 2.11
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
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 2e7;

gas_source = 0;

# 2.12
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
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 2e7;

gas_source = 0;

%% Part 2.2: change the Cvap0, sigma and mu
# 2.21
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.4];

N = [5e4 3e3];
mu=[40e-9, 200e-9];
dilu_on = 0;
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e8;

gas_source = 0;

# 2.22
fixed_sections = 1;
sedi_on = 0;
coag_on = 1;
Dp_min = -9;
Dp_max = -6;
tvect = 0:60:32400;
sigma = [1.25, 1.4];

N = [5e4 3e3];
mu=[40e-9, 200e-9];
dilu_on = 1;
sections = 20;
output_sections = 10*sections;
Cvap_const = 1;
Cvap0 = 5e8;

gas_source = 0;
#
