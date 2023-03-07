% Pool Strategy of a Producer With Endogenous Formulation of Locational Marginal Prices
% 
% For more details, refer to the paper IEEE TRANS ON POWER SYSTEMS
% url: https://ieeexplore.ieee.org/abstract/document/5272241
% doi: 10.1109/TPWRS.2009.2030378
% 
clc; clear; close all;

%% Define Constant
nt = 24;                    % number of time periods
ns = 4;                     % number of strategic    gens
no = 4;                     % number of nonstrategic gens
nd = 4;                     % number of demands
ng = 4;                     % number of generation blocks
nk = 5;                     % number of demand     blocks;
nb = 6;                     % number of buses
nl = 16;                    % number of lines

%% Network Topology
% system MVA base
baseMVA = 100;

%   fbus    tbus    b (p.u.)    RATE (MW)
branch  = [
    1   2   9.412   500;
    1   3   9.412   500;
    2   3   9.412   500;
    2   4   9.412   500;
    2   1   9.412   500;
    3   6   9.412   500;
    3   1   9.412   500;
    3   2   9.412   500;
    4   5   9.412   500;
    4   6   9.412   500;
    4   2   9.412   500;
    5   6   9.412   500;
    5   4   9.412   500;
    6   3   9.412   500;
    6   4   9.412   500;
    6   5   9.412   500;
];

gsbus = [1 2 3 6];                                              % buses where strategic gens located at
Cgs   = sparse(gsbus, 1 : ns, 1, nb, ns);                       % connection matrix element (i, j) is 1 if strategic gen j is at bus i
gobus = [1 2 3 5];                                              % buses where nonstrategic gens located at
Cgo   = sparse(gobus, 1 : no, 1, nb, no);                       % connection matrix element (i, j) is 1 if nonstrategic gen j is at bus i
dbus  = [3 4 5 6];                                              % buses where demand located at
Cd    = sparse(dbus,  1 : nd, 1, nb, nd);                       % connection matrix element (i, j) is 1 if demand j is at bus i
f     = branch(:, 1);                                           % list of "from" buses
t     = branch(:, 2);                                           % list of "to" buses
i     = [(1 : nl)'; (1:nl)'];                                   % double set of row indices
Cft   = sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb); % connection matrix for line and from-to buses
Bf    = spdiags(branch(:, 3), 0, nl, nl) * Cft;                 % build Bf such that Bf * Va is the vector of branch flow injected at each branch's "from" bus (p.u.)
Bbus  = 1/2 * (Cft' * Bf);                                      % build Bbus such that Bbus * Va is the vector of real power injections at each bus (p.u.)

%% Load Data
% data for the strategic gen unit
lambdaS = [9.92 10.25 10.68 11.26;
          18.60 20.03 21.67 22.72;
          9.92 10.25 10.68 11.26;
          10.08 10.66 11.09 11.72];  % (ns x ng)
lambdaS = repmat(lambdaS, [1 1 nt]); % (ns x ng x nt)
lambdaS = permute(lambdaS, [3 1 2]); % marginal cost of block b of unit i of the strategic producer in period t ($/MWh) (nt x ns x ng)

RUP     = [90; 210; 90; 90];         % ramp-up   limit for generation unit i of the strategic producer (MW)
RLO     = [90; 210; 90; 90];         % ramp-down limit for generation unit i of the strategic producer (MW)

PS_MAX  = [54.25 38.75 31.00 31.00;
           25.00 25.00 20.00 20.00;
           54.25 38.75 31.00 31.00;
           68.95 49.25 39.40 39.40]; % (ns x ng)
PS_MAX  = repmat(PS_MAX, [1 1 nt]);  % (ns x ng x nt)
PS_MAX  = permute(PS_MAX, [3 1 2]);  % upper limit of block b of generation unit i of the strategic producer in period t (MWh) (nt x ns x ng)

% data for the nonstrategic gen unit
lambdaO = [19.20 20.32 21.22 22.13;
           10.08 10.66 11.09 11.72;
           10.08 10.66 11.09 11.72;
           9.92  10.25 10.68 11.26]; % (no x ng)
lambdaO = repmat(lambdaO, [1 1 nt]); % (no x ng x nt)
lambdaO = permute(lambdaO, [3 1 2]); % marginal cost of block b of unit j of the nonstrategic producer in period t ($/MWh) (nt x no x ng)

PO_MAX  = [140.0 97.50 52.50 70.00;
          68.95 49.25 39.40 39.40;
          68.95 49.25 39.40 39.40;
          54.25 38.75 31.00 31.00];  % (no x ng)
PO_MAX  = repmat(PO_MAX, [1 1 nt]);  % (no x ng x nt)
PO_MAX  = permute(PO_MAX, [3 1 2]);  % upper limit of block b of generation unit j of the nonstrategic producer in period t (MWh) (nt x no x ng)

% data for the demand
lambdaD = [17.430 17.250 17.216 16.886 16.790;
           17.250 17.216 16.886 16.790 16.380;
           17.216 16.886 16.790 16.380 16.320;
           17.216 16.886 16.790 16.380 16.320;
           16.886 16.790 16.380 16.320 16.130;
           16.886 16.790 16.380 16.320 16.130;
           17.250 17.216 16.886 16.790 16.380;
           17.940 17.612 17.430 17.250 17.216;
           19.232 18.932 18.806 18.344 18.152;
           20.378 19.992 19.532 19.232 18.932;
           24.968 22.628 20.876 20.606 20.378;
           25.000 24.968 22.628 20.876 20.606;
           25.000 24.968 22.628 20.876 20.606;
           24.968 22.628 20.876 20.606 20.378;
           20.378 19.922 19.532 19.232 18.932;
           20.378 19.922 19.532 19.232 18.932;
           20.876 20.606 20.378 19.992 19.532;
           25.000 24.968 22.628 20.876 20.606;
           25.000 24.968 22.628 20.876 20.606;
           25.000 24.968 22.628 20.876 20.606;
           25.000 24.968 22.628 20.876 20.606;
           24.968 22.628 20.876 20.606 20.378;
           19.532 19.232 18.932 18.806 18.344;
           17.940 17.612 17.430 17.250 17.216];  % (nt x nk)
lambdaD = repmat(lambdaD, [1 1 nd]);             % (nt x nk x nd)
lambdaD = permute(lambdaD, [1 3 2]);             % marginal utility of block k of demand d in period t ($/MWh) (nt x nd x nk)

PD_MAX  = reshape(cat(2, 0.19 * ones(nt, nk), 0.27 * ones(nt, nk), 0.27 * ones(nt, nk), 0.27 * ones(nt, nk)) .* ...
                 repmat(repmat([900 25 25 25 25], nt, 1), 1, nd), [nt nk nd]); % (nt x nk x nd)
PD_MAX  = permute(PD_MAX, [1 3 2]);              % upper limit of block k of demand d in period t (MWh) (nt x nd x nk)

% data for the network
C_MAX   = repmat(branch(:, 4), 1, nt);           % transmission capacity of line n-m in period t (MW) (nl x nt)

%% Define Variables
% upper level primal variable
alphas_tib      = sdpvar(nt, ns, ng, 'full'); % price offer of block b of unit i of the strategic producer in period t       ($/MWh)

% lower level primal variable
ps_tib          = sdpvar(nt, ns, ng, 'full'); % power produced by block b of unit i of the strategic producer in period t    (MW)
po_tjb          = sdpvar(nt, no, ng, 'full'); % power produced by block b of unit j of the nonstrategic producer in period t (MW)
pd_tdk          = sdpvar(nt, nd, nk, 'full'); % power consumer by block k of demand d in period t                            (MW)
delta_nt        = sdpvar(nb, nt,     'full'); % voltage angle of bus n in period t                                           (rad)

% lower level dual variable
lambda          = sdpvar(nb, nt,     'full'); % dual of generation-demand equilibrium at bus n    in period t                ($/MWh)
musUB           = sdpvar(nt, ns, ng, 'full'); % dual of capacity of           block b of unit i   in period t                ($/MWh)
musLB           = sdpvar(nt, ns, ng, 'full'); % dual of minimum production of block b of unit i   in period t                ($/MWh)
muoUB           = sdpvar(nt, no, ng, 'full'); % dual of capacity of           block b of unit j   in period t                ($/MWh)
muoLB           = sdpvar(nt, no, ng, 'full'); % dual of minimum production of block b of unit j   in period t                ($/MWh)
mudUB           = sdpvar(nt, nd, nk, 'full'); % dual of capacity of           block k of demand d in period t                ($/MWh)
mudLB           = sdpvar(nt, nd, nk, 'full'); % dual of minimum power of      block k of demand d in period t                ($/MWh)
nuUB            = sdpvar(nl, nt,     'full'); % dual of transmission capacity of line n-m in period t and direction n-m      ($/MWh)
nuLB            = sdpvar(nl, nt,     'full'); % dual of transmission capacity of line n-m in period t and direction n-m      ($/MWh)
xiUB            = sdpvar(nb, nt,     'full'); % dual of upper bound of the voltage angle at bus n in period t                ($/h)
xiLB            = sdpvar(nb, nt,     'full'); % dual of lower bound of the voltage angle at bus n in period t                ($/h)
xi1             = sdpvar(nb, nt,     'full'); % dual of voltage angle at bus n = 1 in period t                               ($/h)
xi1(2 : end, :) = 0;

% auxiliary binary variable
omegasLB        = binvar(nt, ns, ng, 'full'); % binary variables for the complementarity condition           ps_tib                      * musLB = 0
omegasUB        = binvar(nt, ns, ng, 'full'); % binary variables for the complementarity condition (PS_MAX - ps_tib)                     * musUB = 0
omegaoLB        = binvar(nt, no, ng, 'full'); % binary variables for the complementarity condition           po_tjb                      * muoLB = 0
omegaoUB        = binvar(nt, no, ng, 'full'); % binary variables for the complementarity condition (PO_MAX - po_tjb)                     * muoUB = 0
omegadLB        = binvar(nt, nd, nk, 'full'); % binary variables for the complementarity condition           pd_tdk                      * mudLB = 0
omegadUB        = binvar(nt, nd, nk, 'full'); % binary variables for the complementarity condition (PD_MAX - pd_tdk)                     * mudUB = 0
psiLB           = binvar(nl, nt,     'full'); % binary variables for the complementarity condition (C_MAX + Bnm * (delta_nt - delta_nm)) * nuLB  = 0
psiUB           = binvar(nl, nt,     'full'); % binary variables for the complementarity condition (C_MAX - Bnm * (delta_nt - delta_nm)) * nuUB  = 0
phiLB           = binvar(nb, nt,     'full'); % binary variables for the complementarity condition (pi + delta_nt)                       * xiLB  = 0
phiUB           = binvar(nb, nt,     'full'); % binary variables for the complementarity condition (pi - delta_nt)                       * xiUB  = 0

%% Constraints and Model Definition of Single-level Model
objS         = sum(lambdaS .* ps_tib, [1 2 3]) + sum(lambdaO .* po_tjb, [1 2 3]) - sum(lambdaD .* pd_tdk, [1 2 3]); % minimize minus social welfare ($)

constS       = [];

re_ps_tib    = permute(ps_tib, [2 3 1]);          % (ns x ng x nt)
re_ps_tib    = reshape(re_ps_tib, [ns, ng * nt]); % (ns x (ng * nt))
sum_ps_tib   = Cgs * re_ps_tib;                   % (nb x (ng * nt))
sum_ps_tib   = reshape(sum_ps_tib, [nb ng nt]);   % (nb x ng x nt)
sum_ps_tib   = squeeze(sum(sum_ps_tib, 2));       % sum of block b of unit i of the strategic    producer at bus n in period t (MW) (nb x nt)

re_po_tjb    = permute(po_tjb, [2 3 1]);          % (no x ng x nt)
re_po_tjb    = reshape(re_po_tjb, [no, ng * nt]); % (no x (ng * nt))
sum_po_tjb   = Cgo * re_po_tjb;                   % (nb x (ng * nt))
sum_po_tjb   = reshape(sum_po_tjb, [nb ng nt]);   % (nb x ng x nt)
sum_po_tjb   = squeeze(sum(sum_po_tjb, 2));       % sum of block b of unit j of the nonstrategic producer at bus n in period t (MW) (nb x nt)

re_pd_tdk    = permute(pd_tdk, [2 3 1]);          % (nd x nk x nt)
re_pd_tdk    = reshape(re_pd_tdk, [nd, nk * nt]); % (nd x (nk * nt))
sum_pd_tdk   = Cd * re_pd_tdk;                    % (nb x (nk * nt))
sum_pd_tdk   = reshape(sum_pd_tdk, [nb nk nt]);   % (nb x nk x nt)
sum_pd_tdk   = squeeze(sum(sum_pd_tdk, 2));       % sum of block k of dmeand d at bus n in period t (MW) (nb x nt)

sum_bus_inj  = Bbus * delta_nt * baseMVA;         % sum of branch flow at branch's "from" bus n in period t (MW) (nb x nt)

constS       = [constS, ( sum_ps_tib + sum_po_tjb - sum_pd_tdk == sum_bus_inj ): 'power balance'        ];
constS       = [constS, ( 0 <= ps_tib <= PS_MAX ):                               'power bounds'         ];
constS       = [constS, ( 0 <= po_tjb <= PO_MAX ):                               'power bounds'         ];
constS       = [constS, ( 0 <= pd_tdk <= PD_MAX ):                               'power bounds'         ];
constS       = [constS, ( -C_MAX <= Bf * delta_nt * baseMVA <= C_MAX ):          'transmission capacity'];
constS       = [constS, ( -pi <= delta_nt <= pi ):                               'voltage angle bounds' ];
constS       = [constS, ( delta_nt(1, :) == 0 ):                                 'slack bus'            ];
opts         = sdpsettings('solver', 'Gurobi', 'verbose', 0);
diag         = optimize(constS, objS, opts);

% production of the strategic producer when not exercising market power (MWh)
resultS.Prod = sum(value(ps_tib), [1 3]);
% profit of the strategic producer when not exercising market power ($)
resultS.Prof = sum(-dual(constS(1)) .* value(sum_ps_tib), [1 2]) - sum(lambdaS .* value(ps_tib), [1 2 3]);

%% Big-M Value
% constants related to dual variables
MUS           = dual(constS(2));
MUS_MIN       = MUS(1 : nt * ns * ng);
MUS_MIN       = reshape(MUS_MIN, [nt ns ng]);
MUS_MIN       = (MUS_MIN + 1) * 100;            % linearization constants for the complementarity condition ps_tib * musLB = 0                               ($/MWh)
MUS_MAX       = MUS(nt * ns * ng + 1 : end);
MUS_MAX       = reshape(MUS_MAX, [nt ns ng]);
MUS_MAX       = (MUS_MAX + 1) * 100;            % linearization constants for the complementarity condition (PS_MAX - ps_tib) * musUB = 0                    ($/MWh)

MUO           = dual(constS(3));
MUO_MIN       = MUO(1 : nt * no * ng);
MUO_MIN       = reshape(MUO_MIN, [nt no ng]);
MUO_MIN       = (MUO_MIN + 1) * 100;            % linearization constants for the complementarity condition po_tjb * muoLB = 0                               ($/MWh)
MUO_MAX       = MUO(nt * no * ng + 1 : end);
MUO_MAX       = reshape(MUO_MAX, [nt no ng]);
MUO_MAX       = (MUO_MAX + 1) * 100;            % linearization constants for the complementarity condition (PO_MAX - po_tjb) * muoUB = 0                    ($/MWh)

MUD           = dual(constS(4));
MUD_MIN       = MUD(1 : nt * nd * nk);
MUD_MIN       = reshape(MUD_MIN, [nt nd nk]);
MUD_MIN       = (MUD_MIN + 1) * 100;            % linearization constants for the complementarity condition pd_tdk * mudLB = 0                               ($/MWh)
MUD_MAX       = MUD(nt * nd * nk + 1 : end);
MUD_MAX       = reshape(MUD_MAX, [nt nd nk]);
MUD_MAX       = (MUD_MAX + 1) * 100;            % linearization constants for the complementarity condition (PD_MAX - pd_tdk) * mudUB = 0                    ($/MWh)

NU            = dual(constS(5));
NU_MIN        = NU(1 : nl * nt);
NU_MIN        = reshape(NU_MIN, [nl nt]);
NU_MIN        = (NU_MIN + 1) * 100;             % linearization constants for the complementarity condition (C_MAX + Bnm * (delta_nt - delta_nm)) * nuLB = 0 ($/MWh)
NU_MAX        = NU(nl * nt + 1 : end);
NU_MAX        = reshape(NU_MAX, [nl nt]);
NU_MAX        = (NU_MAX + 1) * 100;             % linearization constants for the complementarity condition (C_MAX - Bnm * (delta_nt - delta_nm)) * nuUB = 0 ($/MWh)

XI            = dual(constS(6));
XI_MIN        = XI(1 : nb * nt);
XI_MIN        = reshape(XI_MIN, [nb nt]);
XI_MIN        = (XI_MIN + 1) * 100;             % linearization constants for the complementarity condition (pi + delta_nt) * xiLB = 0 ($/h)
XI_MAX        = XI(nb * nt + 1 : end);
XI_MAX        = reshape(XI_MAX, [nb nt]);
XI_MAX        = (XI_MAX + 1) * 100;             % linearization constants for the complementarity condition (pi - delta_nt) * xiUB = 0 ($/h)

% constants related to primal variables
PS_LIN_MIN    = PS_MAX + 1;                     % linearization constants for the complementarity condition           ps_tib                      * musLB = 0 (MWh)
PS_LIN_MAX    = PS_MAX + 1;                     % linearization constants for the complementarity condition (PS_MAX - ps_tib)                     * musUB = 0 (MWh)
PO_LIN_MIN    = PO_MAX + 1;                     % linearization constants for the complementarity condition           po_tjb                      * muoLB = 0 (MWh)
PO_LIN_MAX    = PO_MAX + 1;                     % linearization constants for the complementarity condition (PO_MAX - po_tjb)                     * muoUB = 0 (MWh)
PD_LIN_MIN    = PD_MAX + 1;                     % linearization constants for the complementarity condition           pd_tdk                      * mudLB = 0 (MWh)
PD_LIN_MAX    = PD_MAX + 1;                     % linearization constants for the complementarity condition (PD_MAX - pd_tdk)                     * mudUB = 0 (MWh)
F_LIN_MIN     = 2 * C_MAX;                      % linearization constants for the complementarity condition (C_MAX + Bnm * (delta_nt + delta_mt)) * nuLB  = 0 (MW)
F_LIN_MAX     = 2 * C_MAX;                      % linearization constants for the complementarity condition (C_MAX - Bnm * (delta_nt - delta_mt)) * nuUB  = 0 (MW)
DELTA_LIN_MIN = 2 * pi;                         % linearization constants for the complementarity condition (pi + delta_nt)                       * xiLB  = 0 (rad)
DELTA_LIN_MAX = 2 * pi;                         % linearization constants for the complementarity condition (pi - delta_nt)                       * xiUB  = 0 (rad)

%% Constraints and Model Definition of Bilevel Model
% Uncongested Network
objB         = sum(lambdaS .* ps_tib, [1 2 3]) + sum(lambdaO .* po_tjb, [1 2 3]) - sum(lambdaD .* pd_tdk, [1 2 3]) + sum(muoUB .* PO_MAX, [1 2 3]) ...
               + sum(mudUB .* PD_MAX, [1 2 3]) + sum(nuLB .* C_MAX, [1 2]) + sum(nuUB .* C_MAX, [1 2]) + pi * sum(xiUB, [1 2]) + pi * sum(xiLB, [1 2]);

constB       = [];

sum_ps_b     = squeeze(sum(ps_tib, 3));
sum_ps_b     = transpose(sum_ps_b);
constB       = [constB, ( -RLO <= sum_ps_b(:, 1) <= RUP ): 'ramp limit'                                                                                           ];
constB       = [constB, ( -repmat(RLO, 1, nt - 1) <= sum_ps_b(:, 2 : end) - sum_ps_b(:, 1 : end - 1) <= repmat(RUP, 1, nt - 1) ): 'ramp limit'                    ];

cov_lambdas  = Cgs' * lambda;                 % (ns x nt)
cov_lambdas  = transpose(cov_lambdas);        % (nt x ns)
cov_lambdas  = repmat(cov_lambdas, [1 1 ng]); % (nt x ns x ng)
constB       = [constB, ( alphas_tib - cov_lambdas + musUB - musLB == 0 ): 'stationary condition'                                                                 ];

cov_lambdao  = Cgo' * lambda;                 % (no x nt)
cov_lambdao  = transpose(cov_lambdao);        % (nt x no)
cov_lambdao  = repmat(cov_lambdao, [1 1 ng]); % (nt x no x ng)
constB       = [constB, ( lambdaO - cov_lambdao + muoUB - muoLB == 0 ):    'stationary condition'                                                                 ];

cov_lambdad = Cd' * lambda;                   % (nd x nt)
cov_lambdad = transpose(cov_lambdad);         % (nt x nd)
cov_lambdad = repmat(cov_lambdad, [1 1 nk]);  % (nt x nd x nk)
constB      = [constB, ( -lambdaD + cov_lambdad + mudUB - mudLB == 0 ):    'stationary condition'                                                                 ];

constB       = [constB, ( Bbus * lambda * baseMVA + Bf' * nuUB * baseMVA - Bf' * nuLB * baseMVA + xiUB - xiLB + xi1 == 0 ): 'stationary condition'                ];
constB       = [constB, ( sum_ps_tib + sum_po_tjb - sum_pd_tdk == sum_bus_inj ):        'power balance'];
constB       = [constB, ( delta_nt(1, :) == 0 ):                                        'slack bus'];
constB       = [constB, ( 0 <= ps_tib ):                                                'complementarity condition of lower bound ps_tib'                         ];
constB       = [constB, ( ps_tib <= (1 - omegasLB) .* PS_LIN_MIN ):                     'complementarity condition of upper bound ps_tib'                         ];
constB       = [constB, ( 0 <= musLB ):                                                 'complementarity condition of lower bound musLB'                          ];
constB       = [constB, ( musLB <= omegasLB .* MUS_MIN ):                               'complementarity condition of upper bound musLB'                          ];
constB       = [constB, ( 0 <= po_tjb ):                                                'complementarity condition of lower bound po_tjb'                         ];
constB       = [constB, ( po_tjb <= (1 - omegaoLB) .* PO_LIN_MIN ):                     'complementarity condition of upper bound po_tjb'                         ];
constB       = [constB, ( 0 <= muoLB ):                                                 'complementarity condition of lower bound muoLB'                          ];
constB       = [constB, ( muoLB <= omegaoLB .* MUO_MIN ):                               'complementarity condition of upper bound muoLB'                          ];
constB       = [constB, (0 <= pd_tdk ):                                                 'complementarity condition of lower bound pd_tdk'                         ];
constB       = [constB, ( pd_tdk <= (1 - omegadLB) .* PD_LIN_MIN ):                     'complementarity condition of upper bound pd_tdk'                         ];
constB       = [constB, ( 0 <= mudLB ):                                                 'complementarity condition of lower bound mudLB'                          ];
constB       = [constB, ( mudLB <= omegadLB .* MUD_MIN ):                               'complementarity condition of upper bound mudLB'                          ];
constB       = [constB, ( 0 <= PS_MAX - ps_tib ):                                       'complementarity condition of lower bound PS_MAX - ps_tib'                ];
constB       = [constB, ( PS_MAX - ps_tib <= (1 - omegasUB) .* PS_LIN_MAX ):            'complementarity condition of upper bound PS_MAX - ps_tib'                ];
constB       = [constB, ( 0 <= musUB ):                                                 'complementarity condition of lower bound musUB'                          ];
constB       = [constB, ( musUB <= omegasUB .* MUS_MAX ):                               'complementarity condition of upper bound musUB'                          ];
constB       = [constB, ( 0 <= PO_MAX - po_tjb ):                                       'complementarity condition of lower bound PO_MAX - po_tjb'                ];
constB       = [constB, ( PO_MAX - po_tjb <= (1 - omegaoUB) .* PO_LIN_MAX ):            'complementarity condition of upper bound PO_MAX - po_tjb'                ];
constB       = [constB, ( 0 <= muoUB ):                                                 'complementarity condition of lower bound muoUB'                          ];
constB       = [constB, ( muoUB <= omegaoUB .* MUO_MAX ):                               'complementarity condition of upper bound muoUB'                          ];
constB       = [constB, ( 0 <= PD_MAX - pd_tdk ):                                       'complementarity condition of lower bound PD_MAX - pd_tdk'                ];
constB       = [constB, ( PD_MAX - pd_tdk <= (1 - omegadUB) .* PD_LIN_MAX ):            'complementarity condition of upper bound PD_MAX - pd_tdk'                ];
constB       = [constB, ( 0 <= mudUB ):                                                 'complementarity condition of lower bound mudUB'                          ];
constB       = [constB, ( mudUB <= omegadUB .* MUD_MAX ):                               'complementarity condition of upper bound mudUB'                          ];
constB       = [constB, ( 0 <= C_MAX + Bf * delta_nt * baseMVA):                        'transmission capacity of lower bound C_MAX + Bnm * (delta_nt - delta_mt)'];
constB       = [constB, ( C_MAX + Bf * delta_nt * baseMVA <= (1 - psiLB) .* F_LIN_MIN): 'transmission capacity of upper bound C_MAX + Bnm * (delta_nt - delta_mt)'];
constB       = [constB, ( 0 <= nuLB ):                                                  'complementarity condition of lower bound nuLB'                           ];
constB       = [constB, ( nuLB <= psiLB .* NU_MIN ):                                    'complementarity condition of upper bound nuLB'                           ];
constB       = [constB, ( 0 <= C_MAX - Bf * delta_nt * baseMVA):                        'transmission capacity of lower bound C_MAX - Bnm * (delta_nt - delta_mt)'];
constB       = [constB, ( C_MAX - Bf * delta_nt * baseMVA <= (1 - psiUB) .* F_LIN_MAX): 'transmission capacity of upper bound C_MAX - Bnm * (delta_nt - delta_mt)'];
constB       = [constB, ( 0 <= nuUB ):                                                  'complementarity condition of lower bound nuUB'                           ];
constB       = [constB, ( nuUB <= psiUB .* NU_MAX ):                                    'complementarity condition of upper bound nuUB'                           ];
constB       = [constB, ( 0 <= pi + delta_nt ):                                         'complementarity condition of lower bound pi + delta_nt'                  ];
constB       = [constB, ( pi + delta_nt <= (1 - phiLB) .* DELTA_LIN_MIN ):              'complementarity condition of upper bound delta_nt'                       ];
constB       = [constB, ( 0 <= xiLB ):                                                  'complementarity condition of lower bound xiLB'                           ];
constB       = [constB, ( xiLB <= phiLB .* XI_MIN ):                                    'complementarity condition of upper bound xiLB'                           ];
constB       = [constB, ( 0 <= pi - delta_nt ):                                         'complementarity condition of lower bound pi - delta_nt'                  ];
constB       = [constB, ( pi - delta_nt <= (1 - phiUB) .* DELTA_LIN_MAX ):              'complementarity condition of upper bound pi - delta_nt'                  ];
constB       = [constB, ( 0 <= xiUB ):                                                  'complementarity condition of lower bound xiUB'                           ];
constB       = [constB, ( xiUB <= phiUB .* XI_MAX ):                                    'complementarity condition of upper bound xiUB'                           ];
diag         = optimize(constB, objB, opts);

musLB_cort   = all(value(musLB) < MUS_MIN, 'all');
musUB_cort   = all(value(musUB) < MUS_MAX, 'all');
muoLB_cort   = all(value(muoLB) < MUO_MIN, 'all');
muoUB_cort   = all(value(muoUB) < MUO_MAX, 'all');
mudLB_cort   = all(value(mudLB) < MUD_MIN, 'all');
mudUB_cort   = all(value(mudUB) < MUD_MAX, 'all');
nuLB_cort    = all(value(nuLB)  < NU_MIN,  'all');
nuUB_cort    = all(value(nuUB)  < NU_MAX,  'all');
xiLB_cort    = all(value(xiLB)  < XI_MIN,  'all');
xiUB_cort    = all(value(xiUB)  < XI_MAX,  'all');
if ~all([musLB_cort musUB_cort muoLB_cort muoUB_cort mudLB_cort mudUB_cort nuLB_cort nuUB_cort xiLB_cort xiUB_cort])
    disp('The dual variables are bounded by dual big-Ms');
end

resultB.Prod = sum(value(ps_tib), [1 3]); % production of the strategic producer when exercising market power (MWh)
resultB.Prof = -value(objB);              % profit of the strategic producer when exercising market power     ($)

figure
plot([1 : nt], value(lambda(1, :)), '-s');
axis([1 24 16 19.5]);
xlabel('t (h)');
ylabel('\lambda ($/MWh');
title('Fig.2. Clearing prices in the uncongested network case');

FLOW = Bf * value(delta_nt) * baseMVA;
f_24 = max(FLOW(4, :));
f_36 = max(FLOW(6, :));
f_64 = max(FLOW(15, :));

% Congested Network
CONG               = C_MAX;
CONG(6,  :)        = 230;
CONG(14, :)        = 230;

objCONG            = sum(lambdaS .* ps_tib, [1 2 3]) + sum(lambdaO .* po_tjb, [1 2 3]) - sum(lambdaD .* pd_tdk, [1 2 3]) + sum(muoUB .* PO_MAX, [1 2 3]) ...
                     + sum(mudUB .* PD_MAX, [1 2 3]) + sum(nuLB .* CONG, [1 2]) + sum(nuUB .* CONG, [1 2]) + pi * sum(xiUB, [1 2]) + pi * sum(xiLB, [1 2]);
constCONG          = constB(setdiff([1 : 48], [33 34 37 38]));
constCONG          = [constCONG, ( 0 <= CONG + Bf * delta_nt * baseMVA ):                        'transmission capacity of lower bound C_MAX + Bnm * (delta_nt - delta_mt)'];
constCONG          = [constCONG, ( CONG + Bf * delta_nt * baseMVA <= (1 - psiLB) .* F_LIN_MIN ): 'transmission capacity of upper bound C_MAX + Bnm * (delta_nt - delta_mt)'];
constCONG          = [constCONG, ( 0 <= CONG - Bf * delta_nt * baseMVA ):                        'transmission capacity of lower bound C_MAX - Bnm * (delta_nt - delta_mt)'];
constCONG          = [constCONG, ( CONG - Bf * delta_nt * baseMVA <= (1 - psiUB) .* F_LIN_MAX ): 'transmission capacity of upper bound C_MAX - Bnm * (delta_nt - delta_mt)'];
diag               = optimize(constCONG, objCONG, opts);

resultCONG.Prod_36 = sum(value(ps_tib), [1 3]); % production of the strategic producer when capacity limit of line 3-6 to 230MW (MWh)
resultCONG.Prof_36 = -value(objCONG);           % profit     of the strategic producer when capacity limit of line 3-6 to 230MW ($)

figure
plot([1 : nt], value(lambda(1, :)), '-*', ...
     [1 : nt], value(lambda(2, :)), '-d', ...
     [1 : nt], value(lambda(3, :)), '-h', ...
     [1 : nt], value(lambda(4, :)), '-o', ...
     [1 : nt], value(lambda(5, :)), '-s', ...
     [1 : nt], value(lambda(6, :)), '-x');
axis([1 nt 16 21]);
legend({'N1', 'N2', 'N3', 'N4', 'N5', 'N6'}, 'Location', 'southeast');
xlabel('t (h)');
ylabel('\lambda ($/MWh)');
title('Fig.3. LMPs with line 3-6 limited');

CONG               = C_MAX;
CONG(10, :)        = 39;
CONG(15, :)        = 39;

objCONG            = sum(lambdaS .* ps_tib, [1 2 3]) + sum(lambdaO .* po_tjb, [1 2 3]) - sum(lambdaD .* pd_tdk, [1 2 3]) + sum(muoUB .* PO_MAX, [1 2 3]) ...
                     + sum(mudUB .* PD_MAX, [1 2 3]) + sum(nuLB .* CONG, [1 2]) + sum(nuUB .* CONG, [1 2]) + pi * sum(xiUB, [1 2]) + pi * sum(xiLB, [1 2]);
constCONG          = constB(setdiff([1 : 48], [33 34 37 38]));
constCONG          = [constCONG, ( 0 <= CONG + Bf * delta_nt * baseMVA ):                        'transmission capacity of lower bound C_MAX + Bnm * (delta_nt - delta_mt)'];
constCONG          = [constCONG, ( CONG + Bf * delta_nt * baseMVA <= (1 - psiLB) .* F_LIN_MIN ): 'transmission capacity of upper bound C_MAX + Bnm * (delta_nt - delta_mt)'];
constCONG          = [constCONG, ( 0 <= CONG - Bf * delta_nt * baseMVA ):                        'transmission capacity of lower bound C_MAX - Bnm * (delta_nt - delta_mt)'];
constCONG          = [constCONG, ( CONG - Bf * delta_nt * baseMVA <= (1 - psiUB) .* F_LIN_MAX ): 'transmission capacity of upper bound C_MAX - Bnm * (delta_nt - delta_mt)'];
diag               = optimize(constCONG, objCONG, opts);

resultCONG.Prod_46 = sum(value(ps_tib), [1 3]); % production of the strategic producer when capacity limit of line 4-6 to 39MW (MWh)
resultCONG.Prof_46 = -value(objCONG);           % profit     of the strategic producer when capacity limit of line 4-6 to 39MW ($)

figure
plot([1 : nt], value(lambda(1, :)), '-d');
axis([1 nt 16 19.5]);
xlabel('t (h)');
ylabel('\lambda ($/MWh)');
title('Fig.4. LMPs with line 4-6 limited');
