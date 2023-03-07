function mpc = case24_ieee_rts_update
% An updated version of the IEEE RTS 24-bus system for electricity market and power system operation studies
%
% For more details, refer to the technique note
% url:https://orbit.dtu.dk/en/publications/an-updated-version-of-the-ieee-rts-24-bus-system-for-electricity-
%
%% system MVA base
mpc.baseMVA = 100;

%% branch data
% The index, name and meaning of each column of the branch matrix is given below:
%   1   fbus   from bus number
%   2   tbus   to bus number
%   3   x      reactance (p.u.)
%   4   rateA  transmission capacity (MW)

%   fbus   tbus   x   rateA
mpc.branch = [
    1   2   0.0146  175;
    1   3   0.2253  175;
    1   5   0.0907  350;
    2   4   0.1356  175;
    2   6   0.205   175;
    3   9   0.1271  175;
    3   24  0.084   400;
    4   9   0.111   175;
    5   10  0.094   350;
    6   10  0.0642  175;
    7   8   0.0652  350;
    8   9   0.1762  175;
    8   10  0.1762  175;
    9   11  0.084   400;
    9   12  0.084   400;
    10  11  0.084   400;
    10  12  0.084   400;
    11  13  0.0488  500;
    11  14  0.0426  500;
    12  13  0.0488  500;
    12  23  0.0985  500;
    13  23  0.0884  500;
    14  16  0.0594  500;
    15  16  0.0172  500;
    15  21  0.0249  1000;
    15  24  0.0529  500;
    16  17  0.0263  500;
    16  19  0.0234  500;
    17  18  0.0143  500;
    17  22  0.1069  500;
    18  21  0.0132  1000;
    19  20  0.0203  1000;
    20  23  0.0112  1000;
    21  22  0.0692  500;
];

%% gen data
% The index, name and meaning of each column of the gen matrix is given below:
%   1   bus     bus number
%   2   PGmax   maximum power output of generating unit i (MW)
%   3   PGmin   minumum power output of generating unit i (MW)
%   4   R+      maximum up reserve capacity of generating unit i (MW)
%   5   R-      maximum down reserve capacity of generating unit i (MW)
%   6   RU      ramp up rate of generating unit i (MW/h)
%   7   RD      ramp down rate of generating unit i (MW/h)
%   8   UT      minimum up time of generating unit i (h)
%   9   UD      minimum down time of generating unit i (h)
%   10  CG      day-ahead offer price of generating unit i ($/MWh)
%   11  CU      upward reserve capacity cost of generating unit i ($/MWh)
%   12  CD      downward reserve capacity cost of generating unit i ($/MWh)
%   13  C+      up regulation offer price of generating unit i ($/MWh)
%   14  C-      down regulation offer price of generating unit i ($/MWh)
%   15  CSU     start-upcost of generating unit i ($)
%   16  PINI    initial power output of generating unit i when t=0 (MW)
%   17  UINI    stating whether generating unit i is online/pffline when t=0 (0/1)
%   18  TINI    number of hours of which the generating unit i was in/out at the beginning of scheduling horizon (h)

%   bus   PGmax   PGmin   R+   R-   RU   RD   UT   UD   CG   CU   CD   C+   C-   CSU   PINI   UINI   TINI
mpc.gen = [
    1   152   30.4    40   40   120   120   8   4   13.32   15   14   15   11   1430.4   76   1   22;
    2   152   30.4    40   40   120   120   8   4   13.32   15   14   15   11   1430.4   76   1   22;
    7   350   75      70   70   350   350   8   8   20.7    10   9    24   16   1725     0    0   -2;
    13  591   206.85  180  180  240   240   12  10  20.93   8    7    25   17   3056.7   0    0   -1;
    15  60    12      60   60   60    60    4   2   26.11   7    5    28   23   437      0    0   -1;
    15  155   54.25   30   30   155   155   8   8   10.52   16   14   16   7    312      0    0   -2;
    16  155   54.25   30   30   155   155   8   8   10.52   16   14   16   7    312      124  1   10;
    18  400   100     0    0    280   280   1   1   6.02    0    0    0    0    0        240  1   50;
    21  400   100     0    0    280   280   1   1   5.47    0    0    0    0    0        240  1   16;
    22  300   300     0    0    300   300   0   0   0       0    0    0    0    0        240  1   24;
    23  310   108.5   60   60   180   180   8   8   10.52   17   16   14   8    624      248  1   10;
    23  350   140     40   40   240   240   8   8   10.89   16   14   16   8    2298     280  1   50;
];

%% load data
sys = [1775.835, 1669.815, 1590.3, 1563.795, 1563.795, 1590.3, 1961.37, 2279.43, 2517.975, 2544.48, 2544.48, 2517.975, ...
    2517.975, 2517.975, 2464.965, 2464.965, 2623.995, 2650.5 2650.5 2544.48, 2411.955, 2199.915, 1934.865, 1669.815];

% The index, name and meaning of each column of the load matrix is given below:
%   1   bus     bus number
%   2   ratio   percentage of the total system load (%)

%   bus   ratio
mpc.load  = [
    1   3.8
    2   3.4
    3   6.3
    4   2.6
    5   2.5
    6   4.8
    7   4.4
    8   6
    9   6.1
    10  6.8
    13  9.3
    14  6.8
    15  11.1
    16  3.5
    18  11.7
    19  6.4
    20  4.5
];

mpc.load(:, 2:25) = repmat(mpc.load(:,2), 1, 24) .* repmat(sys, 17, 1) / 100; % bus load profile from t=1 to t=24 (MW)

end