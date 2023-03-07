% Test system data for the Strategic Offering for a Wind Power Producer in EN and GC

%% system MVA base
baseMVA = 100;

%% branch data
% The index and meaning of each column of the branch matrix is given below:
%   1   "from" bus number
%   2   "to" bus number
%   3   susceptance (p.u.)
%   4   transmission capacity (MW)

%   fbus    tbus    b   rateA
branch = [
    1   2   5   100;
    1   3   5   100;
    2   3   5   100;
];

% change unit of branch susceptance to siemens (S)
branch(:, 3) = baseMVA * branch(:, 3);

%% generator data
% The index and meaning of each column of the gen matrix is givne below:
%   1   bus number
%   2   maximum size of the 1st block (MW)
%   3   maximum size of the 2nd block (MW)
%   4   maximum size of the 3rd block (MW)
%   5   maximum size of the 4th block (MW)
%   6   marginal cost of the 1st block ($/MWh)
%   7   marginal cost of the 2nd block ($/MWh)
%   8   marginal cost of the 3rd block ($/MWh)
%   9   marginal cost of the 4th block ($/MWh)

%   bus   PG1max   PG2max   PG3max   PG4max   lambdaG1   lambdaG2   lambdaG3   lambdaG4
gen = [
    1   60  60  60  60  49  55  67  79;
    2   60  60  60  60  50  54  65  70;
    3   70  60  60  60  44  60  66  80;
];

%% wind unit data
% The index and meaning of each column of the wind matrix is given below:
%   1   bus number
%   2   wind unit capacity (MW)
%   3   marginal cost of wind unit ($/MWh)
%   4   marginal cost of green certificate of wind unit ($/MWh or $/p)

%   bus   PWmax   lambdaW   lambdaWGC
wind = [
    2   300   0   22.18;
];

%% demand data
% The index and meaning of each column of the load matrix is given below:
%   1   bus number
%   2   maximum size of the 1st block (MW)
%   3   maximum size of the 2nd block (MW)
%   4   maximum size of the 3rd block (MW)
%   5   marginal utility of the 1st block ($/MWh)
%   6   marginal utility of the 2nd block ($/MWh)
%   7   marginal utility of the 3rd block ($/MWh)
%   8   renewable energy consumption weight (%)
%   9   marginal utility of green certificate of the 1st block ($/MWh or $/p)
%   10  margianl utility of green certificate of the 2nd block ($/MWh or $/p)
%   11  marginal utility of green certificate of the 3rd block ($/MWh or $/p)

%   bus   PD1max   PD2max   PD3max   lambdaD1   lambdaD2   lambdaD3   ratio   lambdaDGC1   lambdaDGC2   lambdaDGC3
load = [
    1   130   50   50   120   71   59   0.1944   44.49   22.24   11.12;
    2   110   40   40   125   69   50   0.1842   44.63   22.31   11.15;
    3   120   30   30   121   60   40   0.1793   44.19   22.09   11.04;
];
