% Strategic Bidding Under Uncertainty: A Binary Expansion Approach
% 
% For more details, refer to the paper IEEE TRANS ON POWER SYSTEMS
% url: https://ieeexplore.ieee.org/abstract/document/1388508
% doi: 10.1109/TPWRS.2004.840397
%
clc; clear; close all;

%% Define Constant
ns = 5;             % number of strategic    generator
no = 20;            % number of nonstrategic generator
M  = 128;           % number of segment         for discreting price bids of strategic generator
K  = log2(M) + 1;   % number of binary variable for discreting price bids of strategic generator

%% Load Data
Gmaxs   = [76;   30;   1149; 1480; 4186];                             % generator capacities of strategic    generator (MW)
Gmaxo   = [331;  1955; 28;   130;  378; 1022; 1347; 124; 924; 363; ...
           1205; 171;  200;  78;   66;  598;  369;  436; 94;  567];   % generator capacities of nonstrategic generator (MW)
Es      = Gmaxs;                                                      % quantity bids of        strategic    generator (MW)
Eo      = Gmaxo;                                                      % quantity bids of        nonstrategic generator (MW)
Cs      = [23;  163; 177; 224; 279];                                  % operation cost of       strategic    generator ($/MWh)
Co      = [23;  40;  89;  95;  95;  139; 146; 172; 180; 236; ...
           238; 339; 352; 282; 245; 361; 93;  109; 278; 151];         % operation cost of       nonstrategic generator ($/MWh)
lambdao = Co;                                                         % price bids of           nonstrategic generator ($/MWh)
d       = 0.6 * (sum(Gmaxs) + sum(Gmaxo));                            % system load                                    (MWh)

%% Binary Expansion Scheme
lambdasLB  = Cs;                                                      % lower bound of price bids of strategic generator       ($/MWh)
lambdasUB  = 3 * Cs;                                                  % upper bound of price bids of strategic generator       ($/MWh)
lambdasInt = (lambdasUB - lambdasLB) / M;                             % interval between each neighboring discreted price bids ($/MWh)
coeff      = 2 .^ [ 0: log2(M)];                                      % coefficients of each auxiliary binary variable
G          = Gmaxs + 1;                                               % upper bound for strategic generator production gs      (MW)

%% Define variables
% upper level variable
lambdas = sdpvar(ns, 1, 'full'); % price bids of strategic generator                                 ($/MWh)

% lower level primal variable
gs      = sdpvar(ns, 1, 'full'); % energy produced by strategic generator                            (MW)
go      = sdpvar(no, 1, 'full'); % energy produceed by nonstrategic generaator                       (MW)

% lower level dual variable
pid     = sdpvar(1,  1, 'full'); % spot price                                                        ($/MWh)
pigs    = sdpvar(ns, 1, 'full'); % marginal benefit of increaing strategic generator quantity bid    ($/MWh);
pigo    = sdpvar(no, 1, 'full'); % marginal benefit of increaing nonstrategic generator quantity bid ($/MWh)

% auxiliary variable for linearzing lambdas * gs
x_sk    = binvar(ns, K, 'full'); % binary variable for approximating lambdas
z_sk    = sdpvar(ns, K, 'full'); % continuous variable for defining x_ks * gs                        (MW)

%% Constraints and Model Definition of Strategic Bidding using BE
sum_x_sk    = sum(repmat(coeff, ns, 1) .* x_sk, 2); % sum{2^k * x_sk}, k = 0, ... ,K for each strategic generator (ns x 1)
sum_z_sk    = sum(repmat(coeff, ns, 1) .* z_sk, 2); % sum{2^k * z_sk}, k = 0, ... ,K for each strategic generator (ns x 1)

obj         = sum(lambdasLB .* gs + lambdasInt .* sum_z_sk - pigs .* Es - Cs .* gs); % maximize strategic agent's total net revenue ($)
const       = [];
const       = [const, ( lambdasLB + lambdasInt .* sum_x_sk <= lambdasUB ):               'upper bound for the discreted bid price'  ];
const       = [const, ( sum(gs) + sum(go) == d ):                                        'primal feasibility: supply demand balance'];
const       = [const, ( 0 <= gs <= Es ):                                                 'primal feasibility: energy bounds'        ];
const       = [const, ( 0 <= go <= Eo ):                                                 'primal feasibility: energy bounds'        ];
const       = [const, ( pid + pigs - lambdasLB - lambdasInt .* sum_x_sk <= 0 ):          'dual feasibility: gs >= 0'                ];
const       = [const, ( pid + pigo <= lambdao ):                                         'dual feasibility: go >= 0'                ];
const       = [const, ( pigs <= 0 ):                                                     'dual feasibility: gs <= Es'               ];
const       = [const, ( pigo <= 0 ):                                                     'dual feasibility: go <= Eo'               ];
const       = [const, ( sum(lambdasLB .* gs + lambdasInt .* sum_z_sk - pigs .* Es) ...
                        + sum(lambdao .* go - pigo .* Eo) - pid * d == 0 ):              'strong duality'                           ];
const       = [const, ( 0 <= repmat(gs, 1, K) - z_sk <= repmat(G, 1, K) .* (1 - x_sk) ): 'IF-THEN relation of z_sk = x_sk * gs'     ];
const       = [const, ( 0 <= z_sk <= repmat(G, 1, K) .* x_sk ):                          'IF-THEN relation of z_sk = x_sk * gs'     ];
opts        = sdpsettings('solver', 'gurobi', 'verbose', 2);
diag        = optimize(const, -obj, opts);

G_cort = all(value(gs) < G, 'all');
if ~G_cort
    disp('rechoose a valid bound for gs');
end

result.Prod = value(gs);
result.Prof = value(obj);
