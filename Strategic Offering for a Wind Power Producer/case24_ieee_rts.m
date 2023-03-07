% CASE24_IEEE_RTS  Power flow data for the IEEE RELIABILITY TEST SYSTEM.
% 
% This system data is from the IEEE RELIABILITY TEST SYSTEM, see
% 
% IEEE Reliability Test System Task Force of Applications of
% Probability Methods Subcommittee, "IEEE reliability test system-96,"
% IEEE Transactions on Power Systems, Vol. 14, No. 3, Aug. 1999,
% pp. 1010-1020.

%% system MVA base
baseMVA = 100;

%% branch data
% The index and meaning of each column of the branch matrix is given below:
%   1   "from" bus number
%   2   "to" bus number
%   3   reactance (p.u.)
%   4   transmission capacity (MW)

%   fbus   tbus   x   rateA
branch = [
    1   2   0.014   175;
    1   3   0.211   175;
    1   5   0.085   175;
    2   4   0.127   175;
    2   6   0.192   175;
    3   9   0.119   175;
    3   24  0.084   400;
    4   9   0.104   175;
    5   10  0.088   175;
    6   10  0.061   175;
    7   8   0.061   175;
    8   9   0.161   175;
    8   10  0.165   175;
    9   11  0.084   400;
    9   12  0.084   400;
    10  11  0.084   400;
    10  12  0.084   400;
    11  13  0.048   500;
    11  14  0.042   500;
    12  13  0.048   500;
    12  23  0.087   500;
    13  23  0.075   500;
    14  16  0.059   500;
    15  16  0.017   500;
    15  21  0.049   500;
    15  21  0.049   500;
    15  24  0.052   500;
    16  17  0.026   500;
    16  19  0.023   500;
    17  18  0.014   500;
    17  22  0.105   500;
    18  21  0.026   500;
    18  21  0.026   500;
    19  20  0.040   500;
    19  20  0.040   500;
    20  23  0.022   500;
    20  23  0.022   500;
    21  22  0.068   500;
];
% transform the branch reactance to susceptance (S)
branch(:, 3) = baseMVA * (1 ./ branch(:, 3));

%% generator data
% The index and meaning of each column of the gen matrix is given below:
%   1   bus number
%   2   maximum size of the 1st block (MW)
%   3   maximum size of the 2nd block (MW)
%   4   maximum size of the 3rd block (MW)
%   5   marginal cost of the 1st block ($/MWh)
%   6   marginal cost of the 2nd block ($/MWh)
%   7   marginal cost of the 3rd block ($/MWh)

%   bus   PG1max   PG2max   PG3max   lambdaG1   lambdaG2   lambdaG3
gen = [
    1   172   90   82   75   83   91;
    2   172   90   82   77   84   93;
    7   240   140  100  75   85   92;
    13  285   150  135  70   79   82;
    14  200   100  100  72   81   86;
    15  215   130  85   67   76   80;
    16  155   100  55   69   77   87;
    18  400   250  150  71   75   84;
    21  400   250  150  68   79   89;
    22  300   200  100  70   80   75;
    23  660   400  260  65   81   90;
];
% reduce the offer price of the blocks of the gen units by 25%
gen(:, [5 6 7]) = 0.75 * gen(:, [5 6 7]);

%% wind unit data
% The index and meaning of each column of the wind matrix is given below:
%   1   bus number
%   3   wind unit capacity (MW)
%   3   marginal cost of wind unit ($/MWh)
%   4   marginal cost of green certificate of wind unit ($/MWh or $/p)

%   bus   PWmax   lambdaW   lambdaWGC
wind = [
    7   800   0   22.18;
    18  400   0   22.30;
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
    1   100   50   50   124   77   53   0.1944   44.60   22.30   11.15;
    2   100   40   40   129   76   54   0.1842   44.19   22.09   11.04;
    3   150   100  100  135   72   54   0.1793   44.69   22.34   11.17;
    4   50    30   30   144   78   63   0.2141   44.66   22.33   11.16;
    5   50    30   30   139   82   46   0.1444   44.53   22.26   11.13;
    6   150   50   50   127   77   45   0.2187   44.69   22.34   11.17;
    7   150   50   50   135   71   57   0.1690   44.67   22.33   11.16;
    8   200   100  100  137   79   47   0.2929   44.64   22.32   11.16;
    9   200   100  100  130   70   56   0.2340   44.70   22.35   11.17;
    10  200   100  100  127   77   59   0.3245   44.68   22.34   11.17;
    13  300   100  100  133   78   48   0.1771   44.64   22.32   11.16;
    14  200   100  100  140   72   52   0.1946   44.70   22.35   11.17;
    15  400   100  100  125   67   58   0.1734   44.72   22.36   11.18;
    16  100   50   50   132   83   60   0.1996   45.05   22.52   11.26;
    18  400   100  100  131   84   53   0.3239   44.88   22.44   11.22;
    19  200   100  100  138   65   61   0.2277   44.71   22.35   11.17;
    20  150   50   50   137   73   52   0.3750   44.72   22.36   11.18;
];

load(:, 8) = load(:, 8) * 0.4;
