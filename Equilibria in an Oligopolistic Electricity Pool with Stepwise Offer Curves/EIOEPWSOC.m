% Equilibria in an Oligopolistic Electricity Pool with Stepwise Offer Curves
%
% For more details, refer to the paper IEEE TRANS ON POWER SYSTEMS
% url: https://ieeexplore.ieee.org/abstract/document/6064916
% doi: 10.1109/TPWRS.2011.2170439
%
clc; clear; close all;

%% Define Constant
np   = 2;     % number of producers
ng   = 2;     % number of generating units
ng1  = 1;     % number of generating units of producer 1
ng2  = 1;     % number of generating units of producer 2
p1   = [1];   % index of generating units of producer 1
p2   = [2];   % index of generating units of producer 2
no   = 2;     % number of generation blocks
nd   = 2;     % number of demands
nk   = 2;     % number of demand blocks
nb   = 2;     % number of buses
nl   = 2;     % number of lines
CONG = false; % trnasmission congestion

%% Network Topology
% load system data
case2;

% congested network
if CONG
    branch(:, 4) = 10;
    disp('limit transmission capacity to 10MW');
end

% connection matrix of bus-gen
gbus = gen(:, 1);
Cg   = sparse(gbus, 1 : ng, 1, nb, ng);
Cg1  = sparse([1 0; 0 0]);
Cg2  = sparse([0 0; 0 1]);

% connection matrix of bus-demand
dbus = load(:, 1);
Cd   = sparse(dbus, 1 : nd, 1, nb, nd);

% connection matrix of line-bus
f    = branch(:, 1);
t    = branch(:, 2);
i    = [(1 : nl)'; (1 : nl)'];
Cft  = sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb);

%% Load Data
% data for the gen unit
lambdaG = gen(:, [4 5]);                          % marginal cost of block b of unit i ($/MWh) (ng x no)
PG_MAX  = gen(:, [2 3]);                          % capacity of block b of unit i (MW) (ng x no)

% data for the demand
lambdaD = load(:, [4 5]);                         % marginal utility of block k of demand d ($/MWh) (nd x nk)
PD_MAX  = load(:, [2 3]);                         % capacity of block k of demand d (MW) (nd x nk)

% data for the network
Bf      = spdiags(branch(:, 3), 0, nl, nl) * Cft; % Bf * Va is the vector of branch flow injected ...
                                                  % at each branch's "from" bus (S) (nl x nb)
Bbus    = 1/2 * Cft' * Bf;                        % Bbus * Va is the vector of power injection at each bus (S) (nb x nb)
F_MAX   = branch(:, 4);                           % transmission capacity of line n-m (MW) (nl x 1)

%% Define Variable
% upper-level variable of monopoly producer
alpha  = sdpvar(ng,  no, 'full'); % price offer for block b of unit i ($/MWh)

% upper-level variable of duopoly producer 1
alpha1 = sdpvar(ng1, no, 'full'); % price offer for block b of unit i ($/MWh)

% upper-level variable of duopoly producer 2
alpha2 = sdpvar(ng2, no, 'full'); % price offer for block b of unit i ($/MWh)

% lower-level variable
pg     = sdpvar(ng,  no, 'full'); % power produced through block b of unit i (MW)
pd     = sdpvar(nd,  nk, 'full'); % power consumed through block k of demand d (MW)
delta  = sdpvar(nb,  1,  'full'); % voltage angle of bus n (rad)

lambda = sdpvar(nb,  1,  'full'); % dual of power balance at bus n ($/MWh)
mugLB  = sdpvar(ng,  no, 'full'); % dual of lower bound of block b of unit i ($/MWh)
mugUB  = sdpvar(ng,  no, 'full'); % dual of upper bound of block b of unit i ($/MWh)
mudLB  = sdpvar(nd,  nk, 'full'); % dual of lower bound of block k of demand d ($/MWh)
mudUB  = sdpvar(nd,  nk, 'full'); % dual of upper bound of block k of demand d ($/MWh)
nuUB   = sdpvar(nl,  1,  'full'); % dual of capacity limit of transmission line ($/MWh)
xiLB   = sdpvar(nb,  1,  'full'); % dual of lower bound of voltage angle at bus n ($/h)
xiUB   = sdpvar(nb,  1,  'full'); % dual of upper bound of voltage angle at bus n ($/h)
xi1    = sdpvar(nb,  1,  'full'); % dual of voltage angle at bus n = 1 ($/h)
xi1(2:end, :) = 0;

% dual variable of producer 1
beta_alpha11  = sdpvar(ng1, 1,    'full'); % dual of price offer for block 1 of unit i (MW)
beta_alphab1  = sdpvar(ng1, no-1, 'full'); % dual of price offer for block b of unit i (MW)
gamma_bal1    = sdpvar(nb,  1,    'full'); % dual of power balance at bus n ($/MWh)
beta_gLB1     = sdpvar(ng,  no,   'full'); % dual of lower bound of block b of unit i ($/MWh)
beta_gUB1     = sdpvar(ng,  no,   'full'); % dual of upper bound of block b of unit i ($/MWh)
beta_dLB1     = sdpvar(nd,  nk,   'full'); % dual of lower bound of block k of demand d ($/MWh)
beta_dUB1     = sdpvar(nd,  nk,   'full'); % dual of upper bound of block k of demand d ($/MWh)
beta_fUB1     = sdpvar(nl,  1,    'full'); % dual of capacity limit of transmission line ($/MWh)
beta_deltaLB1 = sdpvar(nb,  1,    'full'); % dual of lower bound of voltage angle at bus n ($/h)
beta_deltaUB1 = sdpvar(nb,  1,    'full'); % dual of upper bound of voltage angle at bus n ($/h)
gamma_delta11 = sdpvar(nb,  1,    'full'); % dual of voltage angle at bus n = 1 ($/h)
gamma_delta11(2:end, :) = 0;
% gamma_DT1     = sdpvar(1,   1,    'full'); % dual of strong duality equality (p.u.)
gamma_DT1     = 0.9;
gamma_g1      = sdpvar(ng,  no,   'full'); % dual of stationary condition w.r.t. pg (MW)
gamma_d1      = sdpvar(nd,  nk,   'full'); % dual of stationary condition w.r.t. pd (MW)
gamma_delta1  = sdpvar(nb,  1,    'full'); % dual of stationary condition w.r.t. delta (rad)
eta_gLB1      = sdpvar(ng,  no,   'full'); % dual of nonnegative of mugLB (MW)
eta_gUB1      = sdpvar(ng,  no,   'full'); % dual of nonnegative of mugUB (MW)
eta_dLB1      = sdpvar(nd,  nk,   'full'); % dual of nonnegative of mudLB (MW)
eta_dUB1      = sdpvar(nd,  nk,   'full'); % dual of nonnegative of mudUB (MW)
eta_nuUB1     = sdpvar(nl,  1,    'full'); % dual of nonnegative of nuUB (MW)
eta_xiLB1     = sdpvar(nb,  1,    'full'); % dual of nonnegative of xiLB (rad)
eta_xiUB1     = sdpvar(nb,  1,    'full'); % dual of nonnegative of xiUB (rad)

% dual variable of producer 2
beta_alpha12  = sdpvar(ng2, 1,    'full'); % dual of price offer for block 1 of unit i (MW)
beta_alphab2  = sdpvar(ng2, no-1, 'full'); % dual of price offer for block b of unit i (MW)
gamma_bal2    = sdpvar(nb,  1,    'full'); % dual of power balance at bus n ($/MWh)
beta_gLB2     = sdpvar(ng,  no,   'full'); % dual of lower bound of block b of unit i ($/MWh)
beta_gUB2     = sdpvar(ng,  no,   'full'); % dual of upper bound of block b of unit i ($/MWh)
beta_dLB2     = sdpvar(nd,  nk,   'full'); % dual of lower bound of block k of demand d ($/MWh)
beta_dUB2     = sdpvar(nd,  nk,   'full'); % dual of upper bound of block k of demand d ($/MWh)
beta_fUB2     = sdpvar(nl,  1,    'full'); % dual of capacity limit of transmission line ($/MWh)
beta_deltaLB2 = sdpvar(nb,  1,    'full'); % dual of lower bound of voltage angle at bus n ($/h)
beta_deltaUB2 = sdpvar(nb,  1,    'full'); % dual of upper bound of voltage angle at bus n ($/h)
gamma_delta12 = sdpvar(nb,  1,    'full'); % dual of voltage angle at bus n = 1 ($/h)
gamma_delta12(2:end, :) = 0;
% gamma_DT2     = sdpvar(1,   1,    'full'); % dual of strong duality equality (p.u.)
gamma_DT2     = 1.4;
gamma_g2      = sdpvar(ng,  no,   'full'); % dual of stationary condition w.r.t. pg (MW)
gamma_d2      = sdpvar(nd,  nk,   'full'); % dual of stationary condition w.r.t. pd (MW)
gamma_delta2  = sdpvar(nb,  1,    'full'); % dual of stationary condition w.r.t. delta (rad)
eta_gLB2      = sdpvar(ng,  no,   'full'); % dual of nonnegative of mugLB (MW)
eta_gUB2      = sdpvar(ng,  no,   'full'); % dual of nonnegative of mugUB (MW)
eta_dLB2      = sdpvar(nd,  nk,   'full'); % dual of nonnegative of mudLB (MW)
eta_dUB2      = sdpvar(nd,  nk,   'full'); % dual of nonnegative of mudUB (MW)
eta_nuUB2     = sdpvar(nl,  1,    'full'); % dual of nonnegative of nuUB (MW)
eta_xiLB2     = sdpvar(nb,  1,    'full'); % dual of nonnegative of xiLB (rad)
eta_xiUB2     = sdpvar(nb,  1,    'full'); % dual of nonnegative of xiUB (rad)

% auxiliary binary variable
auxpgLB    = binvar(ng, no, 'full'); % binary variables for the CP pg * mugLB = 0
auxpgUB    = binvar(ng, no, 'full'); % binary variables for the CP (PG_MAX - pg) * mugUB = 0
auxpdLB    = binvar(nd, nk, 'full'); % binary variables for the CP pd * mudLB = 0
auxpdUB    = binvar(nd, nk, 'full'); % binary variables for the CP (PD_MAX - pd) * mudUB = 0
auxfUB     = binvar(nl, 1,  'full'); % binary variables for the CP (F_MAX - Bf * delta) * nuUB = 0
auxdeltaLB = binvar(nb, 1,  'full'); % binary variables for the CP (pi + delta) * xiLB = 0
auxdeltaUB = binvar(nb, 1,  'full'); % binary variables for the CP (pi - delta) * xiUB = 0

auxalpha11  = binvar(ng1, 1,    'full'); % binary variables for the CP alpha1(:, 1) * beta_alpha11 = 0
auxalphab1  = binvar(ng1, no-1, 'full'); % binary variables for the CP (alpha1(:, 2:end) - alpha1(:, 1:end-1)) * beta_alphab1 = 0
auxpgLB1    = binvar(ng,  no,   'full'); % binary variables for the CP pg * beta_gLB1 = 0
auxpgUB1    = binvar(ng,  no,   'full'); % binary variables for the CP (PG_MAX - pg) * beta_gUB1 = 0
auxpdLB1    = binvar(nd,  nk,   'full'); % binary variables for the CP pd * beta_dLB1 = 0
auxpdUB1    = binvar(nd,  nk,   'full'); % binary variables for the CP (PD_MAX - pd) * beta_dUB1 = 0
auxfUB1     = binvar(nl,  1,    'full'); % binary variables for the CP (F_MAX - Bf * delta) * beta_fUB1 = 0
auxdeltaLB1 = binvar(nb,  1,    'full'); % binary variables for the CP (pi + delta) * beta_deltaLB1 = 0
auxdeltaUB1 = binvar(nb,  1,    'full'); % binary variables for the CP (pi - delta) * beta_deltaUB1 = 0
auxmugLB1   = binvar(ng,  no,   'full'); % binary variables for the CP mugLB * eta_gLB1 = 0
auxmugUB1   = binvar(ng,  no,   'full'); % binary variables for the CP mugUB * eta_gUB1 = 0
auxmudLB1   = binvar(nd,  nk,   'full'); % binary variables for the CP mudLB * eta_dLB1 = 0
auxmudUB1   = binvar(nd,  nk,   'full'); % binary variables for the CP mudUB * eta_dUB1 = 0
auxnuUB1    = binvar(nl,  1,    'full'); % binary variables for the CP nuUB * eta_nuUB1 = 0
auxxiLB1    = binvar(nb,  1,    'full'); % binary variables for the CP xiLB * eta_xiLB1 = 0
auxxiUB1    = binvar(nb,  1,    'full'); % binary variables for the CP xiUB * eta_xiUB1 = 0

auxalpha12  = binvar(ng2, 1,    'full'); % binary variables for the CP alpha2(:, 1) * beta_alpha12 = 0
auxalphab2  = binvar(ng2, no-1, 'full'); % binary variables for the CP (alpha2(:, 2:end) - alpha2(:, 1:end-1)) * beta_alphab2 = 0
auxpgLB2    = binvar(ng,  no,   'full'); % binary variables for the CP pg * beta_gLB2 = 0
auxpgUB2    = binvar(ng,  no,   'full'); % binary variables for the CP (PG_MAX - pg) * beta_gUB2 = 0
auxpdLB2    = binvar(nd,  nk,   'full'); % binary variables for the CP pd * beta_dLB2 = 0
auxpdUB2    = binvar(nd,  nk,   'full'); % binary variables for the CP (PD_MAX - pd) * beta_dUB2 = 0
auxfUB2     = binvar(nl,  1,    'full'); % binary variables for the CP (F_MAX - Bf * delta) * beta_fUB2 = 0
auxdeltaLB2 = binvar(nb,  1,    'full'); % binary variables for the CP (pi + delta) * beta_deltaLB2 = 0
auxdeltaUB2 = binvar(nb,  1,    'full'); % binary variables for the CP (pi - delta) * beta_deltaUB2 = 0
auxmugLB2   = binvar(ng,  no,   'full'); % binary variables for the CP mugLB * eta_gLB2 = 0
auxmugUB2   = binvar(ng,  no,   'full'); % binary variables for the CP mugUB * eta_gUB2 = 0
auxmudLB2   = binvar(nd,  nk,   'full'); % binary variables for the CP mudLB * eta_dLB2 = 0
auxmudUB2   = binvar(nd,  nk,   'full'); % binary variables for the CP mudUB * eta_dUB2 = 0
auxnuUB2    = binvar(nl,  1,    'full'); % binary variables for the CP nuUB * eta_nuUB2 = 0
auxxiLB2    = binvar(nb,  1,    'full'); % binary variables for the CP xiLB * eta_xiLB2 = 0
auxxiUB2    = binvar(nb,  1,    'full'); % binary variables for the CP xiUB * eta_xiUB2 = 0

%% Solve the Competitive Market
objC   = sum(lambdaG .* pg, 'all') - sum(lambdaD .* pd, 'all'); % minus social welfare ($)

constC = [];

re_pg  = Cg * pg;       % (nb x no)
sum_pg = sum(re_pg, 2); % sum of block b of unit i at bus n (MW) (nb x 1)
re_pd  = Cd * pd;       % (nb x nk)
sum_pd = sum(re_pd, 2); % sum of block k of demand d at bus n (MW) (nb x 1)
sum_f  = Bbus * delta;  % sum of power flow leaving at bus n (MW) (nb x 1)

constC = [constC, (sum_pd - sum_pg + sum_f == 0): 'power balance'];
constC = [constC, (0 <= pg <= PG_MAX): 'power bounds on pg'];
constC = [constC, (0 <= pd <= PD_MAX): 'power bounds on pd'];
constC = [constC, (Bf * delta <= F_MAX): 'transmission capacity limit'];
constC = [constC, (-pi <= delta <= pi): 'bounds on delta'];
constC = [constC, (delta(1, 1) == 0): 'reference bus'];

optsC  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagC  = optimize(constC, objC, optsC);

re_lambda      = Cg' * dual(constC(1));     % LMP at bus n in which unit i located at ($/MWh) (ng x 1)
re_lambda      = repmat(re_lambda, [1 no]); % (ng x no)
resultC.profit = sum(re_lambda .* value(pg), 'all') - sum(lambdaG .* value(pg), 'all');
                                            % total profit of all producers ($)
resultC.social = -value(objC);              % social welfare of the market ($)

%% Big-M Value
% constant related to primal constraints
pgLB_LIN    = PG_MAX + 1; % upper bounds of the primal constraint 0 <= pg                 (MW)  (ng x no)
pgUB_LIN    = PG_MAX + 1; % upper bounds of the primal constraint 0 <= PG_MAX - pg        (MW)  (ng x no)
pdLB_LIN    = PD_MAX + 1; % upper bounds of the primal constraint 0 <= pd                 (MW)  (nd x nk)
pdUB_LIN    = PD_MAX + 1; % upper bounds of the primal constraint 0 <= PD_MAX - pd        (MW)  (nd x nk)
fUB_LIN     = 2 * F_MAX;  % upper bounds of the primal constraint 0 <= F_MAX - Bf * delta (MW)  (nl x 1)
deltaLB_LIN = 2 * pi;     % upper bounds of the primal constraint 0 <= pi + delta         (rad) (nb x 1)
deltaUB_LIN = 2 * pi;     % upper bounds of the primal constraint 0 <= pi - delta         (rad) (nb x 1)

% constant related dual variables
mugDL       = dual(constC(2));
mugLB_LIN   = mugDL(1 : ng * no, 1);
mugLB_LIN   = reshape(mugLB_LIN, [ng no]);
mugLB_LIN   = (mugLB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= mugLB ($/MWh) (ng x no)
mugUB_LIN   = mugDL(ng * no + 1 : end, 1);
mugUB_LIN   = reshape(mugUB_LIN, [ng no]);
mugUB_LIN   = (mugUB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= mugUB ($/MWh) (ng x no)

mudDL       = dual(constC(3));
mudLB_LIN   = mudDL(1 : nd * nk, 1);
mudLB_LIN   = reshape(mudLB_LIN, [nd nk]);
mudLB_LIN   = (mudLB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= mudLB ($/MWh) (nd x nk)
mudUB_LIN   = mudDL(nd * nk + 1 : end);
mudUB_LIN   = reshape(mudUB_LIN, [nd nk]);
mudUB_LIN   = (mudUB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= mudUB ($/MWh) (nd x nk)

nuUB_LIN    = dual(constC(4));
nuUB_LIN    = (nuUB_LIN + 1) * 100;         % upper bounds of the dual variable 0 <= nuUB  ($/MWh) (nl x 1)

xiDL        = dual(constC(5));
xiLB_LIN    = xiDL(1 : nb, 1);
xiLB_LIN    = (xiLB_LIN + 1) * 100;         % upper bounds of the dual variable 0 <= xiLB  ($/h) (nb x 1)
xiUB_LIN    = xiDL(nb + 1 : end);
xiUB_LIN    = (xiUB_LIN + 1) * 100;         % upper bounds of the dual variable 0 <= xiUB  ($/h) (nb x 1)

%% Solve the Monopoly Market
% minimize minus total profit of monopoly producer ($)
objM   = -sum(lambdaD .* pd, 'all') + sum(mudUB .* PD_MAX, 'all') + sum(nuUB .* F_MAX, 'all') + ...
         pi * sum(xiLB, 'all') + pi * sum(xiUB, 'all') + sum(lambdaG .* pg, 'all');

constM = [];
% upper-level constraints
constM = [constM, (alpha(:, 1) >= 0): 'nonnegative of offer price'];
constM = [constM, (alpha(:, 2:end) >= alpha(:, 1:end-1)): 'stepwise offer price'];

% primal constraints of lower-level problem
constM = [constM, (sum_pd - sum_pg + sum_f == 0): 'power balance'];
constM = [constM, (0 <= pg <= PG_MAX): 'power bounds on pg'];
constM = [constM, (0 <= pd <= PD_MAX): 'power bounds on pd'];
constM = [constM, (Bf * delta <= F_MAX): 'transmission capacity limit'];
constM = [constM, (-pi <= delta <= pi): 'bounds on delta'];
constM = [constM, (delta(1, 1) == 0): 'reference bus'];

% dual constraints of lower-level problem
re_lambdag = Cg' * lambda;               % LMP at bus n in which unit i located at ($/MWh) (ng x 1)
re_lambdag = repmat(re_lambdag, [1 no]); % (ng x no)
constM     = [constM, (alpha - re_lambdag + mugUB - mugLB == 0): 'stationary condition w.r.t. pg'];

re_lambdad = Cd' * lambda;               % LMP at bus n in which demand d located at ($/MWh) (nd x 1)
re_lambdad = repmat(re_lambdad, [1 nk]); % (nd x nk)
constM     = [constM, (re_lambdad - lambdaD + mudUB - mudLB ==0): 'stationary condition w.r.t. pd'];

constM = [constM, (Bbus * lambda + Bf' * nuUB + xiUB - xiLB + xi1 == 0): 'stationary condition w.r.t. delta'];

constM = [
    constM, ...
    (0 <= mugLB): 'nonnegative of mugLB', (0 <= mugUB): 'nonnegative of mugUB', ...
    (0 <= mudLB): 'nonnegative of mudLB', (0 <= mudUB): 'nonnegative of mudUB', ...
    (0 <= nuUB):  'nonnegative of nuUB', ...
    (0 <= xiLB):  'nonnegative of xiLB',  (0 <= xiUB):  'nonnegative of xiUB'
];

% complementarity conditions of lower-level problem
constM = [constM, (pg                 <= pgLB_LIN    .* (1 - auxpgLB)   ): 'upper bounds of the primal constraint 0 <= pg'];
constM = [constM, (PG_MAX - pg        <= pgUB_LIN    .* (1 - auxpgUB)   ): 'upper bounds of the primal constraint 0 <= PG_MAX - pg'];
constM = [constM, (pd                 <= pdLB_LIN    .* (1 - auxpdLB)   ): 'upper bounds of the primal constraint 0 <= pd'];
constM = [constM, (PD_MAX - pd        <= pdUB_LIN    .* (1 - auxpdUB)   ): 'upper bounds of the primal constraint 0 <= PD_MAX - pd'];
constM = [constM, (F_MAX - Bf * delta <= fUB_LIN     .* (1 - auxfUB)    ): 'upper bounds of the primal constraint 0 <= F_MAX - Bf * delta'];
constM = [constM, (pi + delta         <= deltaLB_LIN .* (1 - auxdeltaLB)): 'upper bounds of the primal constraint 0 <= pi + delta'];
constM = [constM, (pi - delta         <= deltaUB_LIN .* (1 - auxdeltaUB)): 'upper bounds of the primal constraint 0 <= pi - delta'];
constM = [constM, (mugLB              <= mugLB_LIN   .* auxpgLB         ): 'upper bounds of the dual variable 0 <= mugLB'];
constM = [constM, (mugUB              <= mugUB_LIN   .* auxpgUB         ): 'upper bounds of the dual variable 0 <= mugUB'];
constM = [constM, (mudLB              <= mudLB_LIN   .* auxpdLB         ): 'upper bounds of the dual variable 0 <= mudLB'];
constM = [constM, (mudUB              <= mudUB_LIN   .* auxpdUB         ): 'upper bounds of the dual variable 0 <= mudUB'];
constM = [constM, (nuUB               <= nuUB_LIN    .* auxfUB          ): 'upper bounds of the dual variable 0 <= nuUB'];
constM = [constM, (xiLB               <= xiLB_LIN    .* auxdeltaLB      ): 'upper bounds of the dual variable 0 <= xiLB'];
constM = [constM, (xiUB               <= xiUB_LIN    .* auxdeltaUB      ): 'upper bounds of the dual variable 0 <= xiUB'];

optsM  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagM  = optimize(constM, objM, optsM);

% check the correctness of big-M of the primal constraints
valpgLB    = all(value(pg)                 < pgLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pg
valpgUB    = all(PG_MAX - value(pg)        < pgUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= PG_MAX - pg
valpdLB    = all(value(pd)                 < pdLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pd
valpdUB    = all(PD_MAX - value(pd)        < pdUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= PD_MAX - pd
valfUB     = all(F_MAX - Bf * value(delta) < fUB_LIN,     'all'); % valid upper bounds of the primal constraint 0 <= F_MAX - Bf * delta
valdeltaLB = all(pi + value(delta)         < deltaLB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pi + delta
valdeltaUB = all(pi - value(delta)         < deltaUB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pi - delta

% check the correctness of big-M of the dual variables
valmugLB = all(value(mugLB) < mugLB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= mugLB
valmugUB = all(value(mugUB) < mugUB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= mugUB
valmudLB = all(value(mudLB) < mudLB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= mudLB
valmudUB = all(value(mudUB) < mudUB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= mudUB
valnuUB  = all(value(nuUB)  < nuUB_LIN,  'all'); % valid upper bounds of the dual variable 0 <= nuUB
valxiLB  = all(value(xiLB)  < xiLB_LIN,  'all'); % valid upper bounds of the dual variable 0 <= xiLB
valxiUB  = all(value(xiUB)  < xiUB_LIN,  'all'); % valid upper bounds of the dual variable 0 <= xiUB

if ~all([valpgLB valpgUB valpdLB valpdUB valfUB valdeltaLB valdeltaUB]) || ...
        ~all([valmugLB valmugUB valmudLB valmudUB valnuUB valxiLB valxiUB])
    disp('INCORRECT BIG-M');
end

resultM.profit = -value(objM); % totla profit of monopoly producer ($)
resultM.social = sum(lambdaD .* value(pd), 'all') - sum(lambdaG .* value(pg), 'all'); % social welfare of the market ($)

%% Big-M Value
% constant related to primal constraints
alpha11_LIN = 1000 * ones(ng1, 1);    % upper bounds of the primal constraint 0 <= alpha1(:, 1) ($/MWh) (ng1 x 1)
alphab1_LIN = 1000 * ones(ng1, no-1); % upper bounds of the primal constraint 0 <= alpha1(:, 2:end) - alpha1(:, 1:end-1) ($/MWh) (ng1 x no-1)
alpha12_LIN = 1000 * ones(ng2, 1);    % upper bounds of the primal constraint 0 <= alpha2(:, 1) ($/MWh) (ng2 x 1)
alphab2_LIN = 1000 * ones(ng2, no-1); % upper bounds of the primal constraint 0 <= alpha2(:, 2:end) - alpha2(:, 1:end-1) ($/MWh) (ng2 x no-1)

% constant related to dual variables
beta_alpha11_LIN = 1000 * ones(ng1, 1);    % upper bounds of the dual variable 0 <= beta_alpha11 (MW) (ng1 x 1)
beta_alphab1_LIN = 1000 * ones(ng1, no-1); % upper bounds of the dual variable 0 <= beta_alphab1 (MW) (ng1 x no-1)
beta_alpha12_LIN = 1000 * ones(ng2, 1);    % upper bounds of the dual variable 0 <= beta_alpha12 (MW) (ng2 x 1)
beta_alphab2_LIN = 1000 * ones(ng2, no-1); % upper bounds of the dual variable 0 <= beta_alphab2 (MW) (ng2 x no-1)
beta_gLB_LIN     = 10 * mugLB_LIN;         % upper bounds of the dual variable 0 <= beta_gLB1 & beta_gLB2 ($/MWh) (ng x no)
beta_gUB_LIN     = 10 * mugUB_LIN;         % upper bounds of the dual variable 0 <= beta_gUB1 & beta_gUB2 ($/MWh) (ng x no)
beta_dLB_LIN     = 10 * mudLB_LIN;         % upper bounds of the dual variable 0 <= beta_dLB1 & beta_dLB2 ($/MWh) (nd x nk)
beta_dUB_LIN     = 10 * mudUB_LIN;         % upper bounds of the dual variable 0 <= beta_dUB1 & beta_dUB2 ($/MWh) (nd x nk)
beta_fUB_LIN     = 10 * nuUB_LIN;          % upper bounds of the dual variable 0 <= beta_fUB1 & beta_fUB2 ($/MWh) (nl x 1)
beta_deltaLB_LIN = 10 * xiLB_LIN;          % upper bounds of the dual variable 0 <= beta_deltaLB1 & beta_deltaLB2 ($/h) (nb x 1)
beta_deltaUB_LIN = 10 * xiUB_LIN;          % upper bounds of the dual variable 0 <= beta_deltaUB1 & beta_deltaUB2 ($/h) (nb x 1)
eta_gLB_LIN      = 1000 * ones(ng, no);    % upper bounds of the dual variable 0 <= eta_gLB1 & eta_gLB2 (MW) (ng x no)
eta_gUB_LIN      = 1000 * ones(ng, no);    % upper bounds of the dual variable 0 <= eta_gUB1 & eta_gUB2 (MW) (ng x no)
eta_dLB_LIN      = 1000 * ones(nd, nk);    % upper bounds of the dual variable 0 <= eta_dLB1 & eta_dLB2 (MW) (nd x nk)
eta_dUB_LIN      = 1000 * ones(nd, nk);    % upper bounds of the dual variable 0 <= eta_dUB1 & eta_dUB2 (MW) (nd x nk)
eta_nuUB_LIN     = 1000 * ones(nl, 1);     % upper bounds of the dual variable 0 <= eta_nuUB1 & eta_nuUB2 (MW) (nl x 1)
eta_xiLB_LIN     = 1000 * ones(nb, 1);     % upper bounds of the dual variable 0 <= eta_xiLB1 & eta_xiLB2 (rad) (nb x 1)
eta_xiUB_LIN     = 1000 * ones(nb, 1);     % upper bounds of the dual variable 0 <= eta_xiUB1 & eta_xiUB2 (rad) (nb x 1)

%% Solve the Duopoly Market
% maximize total profit of all producers ($)
objTP = sum(lambdaD .* pd, 'all') - sum(mudUB .* PD_MAX, 'all') - sum(nuUB .* F_MAX, 'all') - ...
        pi * sum(xiLB, 'all') - pi * sum(xiUB, 'all') - sum(lambdaG .* pg, 'all');

% maximize social welfare of the market ($)
objSW = sum(lambdaD .* pd, 'all') - sum(lambdaG .* pg, 'all');

% primal constraints of MPECs
constP1 = [];
constP1 = [constP1, (alpha1(:, 1) >= 0): 'nonnegative of offer price'];
constP1 = [constP1, (alpha1(:, 2:end) >= alpha1(:, 1:end-1)): 'stepwise offer price'];
constP1 = [constP1, (alpha1 == alpha(1, :)): 'offer price of producer 1'];
constP1 = [constP1, (sum_pd - sum_pg + sum_f == 0): 'power balance'];
constP1 = [constP1, (0 <= pg <= PG_MAX): 'power bounds on pg'];
constP1 = [constP1, (0 <= pd <= PD_MAX): 'power bounds on pd'];
constP1 = [constP1, (Bf * delta <= F_MAX): 'transmission capacity limit'];
constP1 = [constP1, (-pi <= delta <= pi): 'bounds on delta'];
constP1 = [constP1, (delta(1, 1) == 0): 'reference bus'];
constP1 = [constP1, (alpha - re_lambdag + mugUB -mugLB == 0): 'stationary condition w.r.t. pg'];
constP1 = [constP1, (re_lambdad - lambdaD + mudUB - mudLB == 0): 'stationary condition w.r.t. pd'];
constP1 = [constP1, (Bbus * lambda + Bf' * nuUB + xiUB - xiLB + xi1 == 0): 'stationary condition w.r.t. delta'];
constP1 = [
    constP1, ...
    (0 <= mugLB): 'nonnegative of mugLB', (0 <= mugUB): 'nonnegative of mugUB', ...
    (0 <= mudLB): 'nonnegative of mudLB', (0 <= mudUB): 'nonnegative of mudUB', ...
    (0 <= nuUB ): 'nonnegative of nuUB', ...
    (0 <= xiLB ): 'nonnegative of xiLB',  (0 <= xiUB ): 'nonnegative of xiUB'
];
constP1 = [constP1, (pg                 <= pgLB_LIN    .* (1 - auxpgLB)   ): 'upper bounds of the primal constraint 0 <= pg'];
constP1 = [constP1, (PG_MAX - pg        <= pgUB_LIN    .* (1 - auxpgUB)   ): 'upper bounds of the primal constraint 0 <= PG_MAX - pg'];
constP1 = [constP1, (pd                 <= pdLB_LIN    .* (1 - auxpdLB)   ): 'upper bounds of the primal constraint 0 <= pd'];
constP1 = [constP1, (PD_MAX - pd        <= pdUB_LIN    .* (1 - auxpdUB)   ): 'upper bounds of the primal constraint 0 <= PD_MAX - pd'];
constP1 = [constP1, (F_MAX - Bf * delta <= fUB_LIN     .* (1 - auxfUB)    ): 'upper bounds of the primal constraint 0 <= F_MAX - Bf * delta'];
constP1 = [constP1, (pi + delta         <= deltaLB_LIN .* (1 - auxdeltaLB)): 'upper bounds of the primal constraint 0 <= pi + delta'];
constP1 = [constP1, (pi - delta         <= deltaUB_LIN .* (1 - auxdeltaUB)): 'upper bounds of the primal constraint 0 <= pi - delta'];
constP1 = [constP1, (mugLB              <= mugLB_LIN   .* auxpgLB         ): 'upper bounds of the dual variable 0 <= mugLB'];
constP1 = [constP1, (mugUB              <= mugUB_LIN   .* auxpgUB         ): 'upper bounds of the dual variable 0 <= mugUB'];
constP1 = [constP1, (mudLB              <= mudLB_LIN   .* auxpdLB         ): 'upper bounds of the dual variable 0 <= mudLB'];
constP1 = [constP1, (mudUB              <= mudUB_LIN   .* auxpdUB         ): 'upper bounds of the dual variable 0 <= mudUB'];
constP1 = [constP1, (nuUB               <= nuUB_LIN    .* auxfUB          ): 'upper bounds of the dual variable 0 <= nuUB'];
constP1 = [constP1, (xiLB               <= xiLB_LIN    .* auxdeltaLB      ): 'upper bounds of the dual variable 0 <= xiLB'];
constP1 = [constP1, (xiUB               <= xiUB_LIN    .* auxdeltaUB      ): 'upper bounds of the dual variable 0 <= xiUB'];

constP2 = [];
constP2 = [constP2, (alpha2(:, 1) >= 0): 'nonnegative of offer price'];
constP2 = [constP2, (alpha2(:, 2:end) >= alpha2(:, 1:end-1)): 'stepwise offer price'];
constP2 = [constP2, (alpha2 == alpha(2, :)): 'offer price of producer 2'];

constP  = [constP1, constP2];

% dual constraints of MPECs
constD1 = [];
constD1 = [constD1, (gamma_g1(p1, 1) + gamma_DT1 * pg(p1, 1) - beta_alpha11 + beta_alphab1(:, 1) == 0): ...
           'differentiating L1 w.r.t. alpha1'];
constD1 = [constD1, (gamma_g1(p1, end) + gamma_DT1 * pg(p1, end) - beta_alphab1(:, end) == 0): ...
           'differentiating L1 w.r.t. alpha1'];
re_gamma_bal1g = Cg' * gamma_bal1;               % dual of power balance at bus n in which unit i located at ($/MWh) (ng x 1)
re_gamma_bal1g = repmat(re_gamma_bal1g, [1 no]); % (ng x no)
constD1 = [constD1, (-re_lambdag(p1, :) + lambdaG(p1, :) - re_gamma_bal1g(p1, :) - beta_gLB1(p1, :) + beta_gUB1(p1, :) + ...
           gamma_DT1 * alpha(p1, :) == 0): 'differentiating L1 w.r.t. pg'];
constD1 = [constD1, (-re_gamma_bal1g(p2, :) - beta_gLB1(p2, :) + beta_gUB1(p2, :) + gamma_DT1 * alpha(p2, :) == 0): ...
           'differentiating L1 w.r.t. pg'];
re_gamma_bal1d = Cd' * gamma_bal1;               % dual of power balance at bus n in which demand d located at ($/MWh) (nd x1)
re_gamma_bal1d = repmat(re_gamma_bal1d, [1 nk]); % (nd x nk)
constD1 = [constD1, (re_gamma_bal1d - beta_dLB1 + beta_dUB1 - gamma_DT1 * lambdaD == 0): 'differentiating L1 w.r.t. pd'];
constD1 = [constD1, (Bbus * gamma_bal1 + Bf' * beta_fUB1 + beta_deltaUB1 - beta_deltaLB1 + gamma_delta11 == 0): ...
           'differentiating L1 w.r.t. delta'];
re_pg            = Cg1 * pg;            % (nb x no)
sum_pg           = sum(re_pg, 2);       % (nb x 1)
re_gamma_g1      = Cg * gamma_g1;       % (nb x no)
sum_gamma_g1     = sum(re_gamma_g1, 2); % (nb x 1)
re_gamma_d1      = Cd * gamma_d1;       % (nb x nk)
sum_gamma_d1     = sum(re_gamma_d1, 2); % (nb x 1)
sum_gamma_delta1 = Bbus * gamma_delta1; % (nb x 1)
constD1 = [constD1, (-sum_pg - sum_gamma_g1 + sum_gamma_d1 + sum_gamma_delta1 == 0): 'differentiating L1 w.r.t. lambda'];
constD1 = [constD1, (-gamma_g1 - eta_gLB1 == 0): 'differentiating L1 w.r.t. mugLB'];
constD1 = [constD1, (gamma_g1 + gamma_DT1 * PG_MAX - eta_gUB1 ==0): 'differentiating L1 w.r.t. mugUB'];
constD1 = [constD1, (-gamma_d1 - eta_dLB1 == 0): 'differentiating L1 w.r.t. mudLB'];
constD1 = [constD1, (gamma_d1 + gamma_DT1 * PD_MAX - eta_dUB1 == 0): 'differentiating L1 w.r.t. mudUB'];
constD1 = [constD1, (Bf * gamma_delta1 + gamma_DT1 * F_MAX - eta_nuUB1 == 0): 'differentiating L1 w.r.t. nuUB'];
constD1 = [constD1, (-gamma_delta1 + gamma_DT1 * pi - eta_xiLB1 ==0): 'differentiating L1 w.r.t. xiLB'];
constD1 = [constD1, (gamma_delta1 + gamma_DT1 * pi - eta_xiUB1 ==0):' differentiating L1 w.r.t. xiUB'];
constD1 = [constD1, (gamma_delta1(1, 1) == 0): 'differentiating L1 w.r.t. xi1'];
constD1 = [
    constD1, ...
    (0 <= beta_alpha11 ): 'nonnegative of beta_alpha11',  (0 <= beta_alphab1 ): 'nonnegative of beta_alphab1', ...
    (0 <= beta_gLB1    ): 'nonnegative of beta_gLB1',     (0 <= beta_gUB1    ): 'nonnegative of beta_gUB1', ....
    (0 <= beta_dLB1    ): 'nonnegative of beta_dLB1',     (0 <= beta_dUB1    ): 'nonnegative of beta_dUB1', ...
    (0 <= beta_fUB1    ): 'nonnegative of beta_fUB1', ...
    (0 <= beta_deltaLB1): 'nonnegative of beta_deltaLB1', (0 <= beta_deltaUB1): 'nonnegative of beta_deltaUb1', ...
    (0 <= eta_gLB1     ): 'nonnegative of eta_gLB1',      (0 <= eta_gUB1     ): 'nonnegative of eta_gUB1', ...
    (0 <= eta_dLB1     ): 'nonnegative of eta_dLB1',      (0 <= eta_dUB1     ): 'nonnegative of eta_dUB1', ...
    (0 <= eta_nuUB1    ): 'nonnegative of eta_nuUB1', ...
    (0 <= eta_xiLB1    ): 'nonnegative of eta_xiLB1',     (0 <= eta_xiUB1    ): 'nonnegative of eta_xiUB1'
];

constD2 = [];
constD2 = [constD2, (gamma_g2(p2, 1) + gamma_DT2 * pg(p2, 1) - beta_alpha12 + beta_alphab2(:, 1) == 0): ...
           'differentiating L2 w.r.t. alpha2'];
constD2 = [constD2, (gamma_g2(p2, end) + gamma_DT2 * pg(p2, end) - beta_alphab2(:, end) == 0): ...
           'differentiating L2 w.r.t. alpha2'];
re_gamma_bal2g = Cg' * gamma_bal2;               % dual of power balance at bus n in which unit i located at ($/MWh) (ng x 1)
re_gamma_bal2g = repmat(re_gamma_bal2g, [1 no]); % (ng x no)
constD2 = [constD2, (-re_lambdag(p2, :) + lambdaG(p2, :) - re_gamma_bal2g(p2, :) - beta_gLB2(p2, :) + beta_gUB2(p2, :) + ...
           gamma_DT2 * alpha(p2, :) == 0): 'differentiating L2 w.r.t. pg'];
constD2 = [constD2, (-re_gamma_bal2g(p1, :) - beta_gLB2(p1, :) + beta_gUB2(p1, :) + gamma_DT2 * alpha(p1, :) == 0): ...
           'differentiating L2 w.r.t. pg'];
re_gamma_bal2d = Cd' * gamma_bal2;               % dual of power balance at bus n in which demand d located at ($/MWh) (nd x 1)
re_gamma_bal2d = repmat(re_gamma_bal2d, [1 nk]); % (nd x nk)
constD2 = [constD2, (re_gamma_bal2d - beta_dLB2 + beta_dUB2 - gamma_DT2 * lambdaD == 0): 'differentiating L2 w.r.t. pd'];
constD2 = [constD2, (Bbus * gamma_bal2 + Bf' * beta_fUB2 + beta_deltaUB2 - beta_deltaLB2 + gamma_delta12 == 0): ...
           'differentiating L2 w.r.t. delta'];
re_pg            = Cg2 * pg;            % (nb x no)
sum_pg           = sum(re_pg, 2);       % (nb x 1)
re_gamma_g2      = Cg * gamma_g2;       % (nb x no)
sum_gamma_g2     = sum(re_gamma_g2, 2); % (nb x 1)
re_gamma_d2      = Cd * gamma_d2;       % (nb x nk)
sum_gamma_d2     = sum(re_gamma_d2, 2); % (nb x 1)
sum_gamma_delta2 = Bbus * gamma_delta2; % (nb x 1)
constD2 = [constD2, (-sum_pg - sum_gamma_g2 + sum_gamma_d2 + sum_gamma_delta2 == 0): 'differentiating L2 w.r.t. lambda'];
constD2 = [constD2, (-gamma_g2 - eta_gLB2 == 0): 'differentiating L2 w.r.t. mugLB'];
constD2 = [constD2, (gamma_g2 + gamma_DT2 * PG_MAX - eta_gUB2 ==0): 'differentiating L2 w.r.t. mugUB'];
constD2 = [constD2, (-gamma_d2 - eta_dLB2 == 0): 'differentiating L2 w.r.t. mudLB'];
constD2 = [constD2, (gamma_d2 + gamma_DT2 * PD_MAX - eta_dUB2 == 0): 'differentiating L2 w.r.t. mudUB'];
constD2 = [constD2, (Bf * gamma_delta2 + gamma_DT2 * F_MAX - eta_nuUB2 == 0): 'differentiating L2 w.r.t. nuUB'];
constD2 = [constD2, (-gamma_delta2 + gamma_DT2 * pi - eta_xiLB2 == 0): 'differentiating L2 w.r.t. xiLB'];
constD2 = [constD2, (gamma_delta2 + gamma_DT2 * pi - eta_xiUB2 == 0): 'differentiating L2 w.r.t. xiUB'];
constD2 = [constD2, (gamma_delta2(1, 1) == 0): 'differentiating L2 w.r.t. xi1'];
constD2 = [
    constD2, ...
    (0 <= beta_alpha12 ): 'nonnegative of beta_alpha12',  (0 <= beta_alphab2 ): 'nonnegative of beta_alphab2', ...
    (0 <= beta_gLB2    ): 'nonnegative of beta_gLB2',     (0 <= beta_gUB2    ): 'nonnegative of beta_gUB2', ....
    (0 <= beta_dLB2    ): 'nonnegative of beta_dLB2',     (0 <= beta_dUB2    ): 'nonnegative of beta_dUB2', ...
    (0 <= beta_fUB2    ): 'nonnegative of beta_fUB2', ...
    (0 <= beta_deltaLB2): 'nonnegative of beta_deltaLB2', (0 <= beta_deltaUB2): 'nonnegative of beta_deltaUb2', ...
    (0 <= eta_gLB2     ): 'nonnegative of eta_gLB2',      (0 <= eta_gUB2     ): 'nonnegative of eta_gUB2', ...
    (0 <= eta_dLB2     ): 'nonnegative of eta_dLB2',      (0 <= eta_dUB2     ): 'nonnegative of eta_dUB2', ...
    (0 <= eta_nuUB2    ): 'nonnegative of eta_nuUB2', ...
    (0 <= eta_xiLB2    ): 'nonnegative of eta_xiLB2',     (0 <= eta_xiUB2    ): 'nonnegative of eta_xiUB2'
];

constD  = [constD1, constD2];

% complementarity conditions of MPECs
constCP1 = [];
constCP1 = [constCP1, (alpha1(:, 1)       <= alpha11_LIN .* (1 - auxalpha11) ): 'upper bounds of the primal constraint 0 <= alpha1(:, 1)'];
constCP1 = [constCP1, (alpha1(:, 2:end) - alpha1(:, 1:end-1) <= alphab1_LIN .* (1 - auxalphab1)): ...
            'upper bounds of the primal constraint 0 <= alpha1(:, 2:end) - alpha1(:, 1:end-1)'];
constCP1 = [constCP1, (pg                 <= pgLB_LIN    .* (1 - auxpgLB1)   ): 'upper bounds of the primal constraint 0 <= pg'];
constCP1 = [constCP1, (PG_MAX - pg        <= pgUB_LIN    .* (1 - auxpgUB1)   ): 'upper bounds of the primal constraint 0 <= PG_MAX - pg'];
constCP1 = [constCP1, (pd                 <= pdLB_LIN    .* (1 - auxpdLB1)   ): 'upper bounds of the primal constraint 0 <= pd'];
constCP1 = [constCP1, (PD_MAX - pd        <= pdUB_LIN    .* (1 - auxpdUB1)   ): 'upper bounds of the primal constraint 0 <= PD_MAX - pd'];
constCP1 = [constCP1, (F_MAX - Bf * delta <= fUB_LIN     .* (1 - auxfUB1)    ): 'upper bounds of the primal constraint 0 <= F_MAX - Bf * delta'];
constCP1 = [constCP1, (pi + delta         <= deltaLB_LIN .* (1 - auxdeltaLB1)): 'upper bounds of the primal constraint 0 <= pi + delta'];
constCP1 = [constCP1, (pi - delta         <= deltaUB_LIN .* (1 - auxdeltaUB1)): 'upper bounds of the primal constraint 0 <= pi - delta'];
constCP1 = [constCP1, (mugLB              <= mugLB_LIN   .* (1 - auxmugLB1)  ): 'upper bounds of the primal constraint 0 <= mugLB'];
constCP1 = [constCP1, (mugUB              <= mugUB_LIN   .* (1 - auxmugUB1)  ): 'upper bounds of the primal constraint 0 <= mugUB'];
constCP1 = [constCP1, (mudLB              <= mudLB_LIN   .* (1 - auxmudLB1)  ): 'upper bounds of the primal constraint 0 <= mudLB'];
constCP1 = [constCP1, (mudUB              <= mudUB_LIN   .* (1 - auxmudUB1)  ): 'upper bounds of the primal constraint 0 <= mudUB'];
constCP1 = [constCP1, (nuUB               <= nuUB_LIN    .* (1 - auxnuUB1)   ): 'upper bounds of the primal constraint 0 <= nuUB'];
constCP1 = [constCP1, (xiLB               <= xiLB_LIN    .* (1 - auxxiLB1)   ): 'upper bounds of the primal constraint 0 <= xiLB'];
constCP1 = [constCP1, (xiUB               <= xiUB_LIN    .* (1 - auxxiUB1)   ): 'upper bounds of the primal constraint 0 <= xiUB'];
constCP1 = [constCP1, (beta_alpha11       <= beta_alpha11_LIN .* auxalpha11  ): 'upper bounds of the dual variable 0 <= beta_alpha11'];
constCP1 = [constCP1, (beta_alphab1       <= beta_alphab1_LIN .* auxalphab1  ): 'upper bounds of the dual variable 0 <= beta_alphab1'];
constCP1 = [constCP1, (beta_gLB1          <= beta_gLB_LIN     .* auxpgLB1    ): 'upper bounds of the dual variable 0 <= beta_gLB1'];
constCP1 = [constCP1, (beta_gUB1          <= beta_gUB_LIN     .* auxpgUB1    ): 'upper bounds of the dual variable 0 <= beta_gUB1'];
constCP1 = [constCP1, (beta_dLB1          <= beta_dLB_LIN     .* auxpdLB1    ): 'upper bounds of the dual variable 0 <= beta_dLB1'];
constCP1 = [constCP1, (beta_dUB1          <= beta_dUB_LIN     .* auxpdUB1    ): 'upper bounds of the dual variable 0 <= beta_dUB1'];
constCP1 = [constCP1, (beta_fUB1          <= beta_fUB_LIN     .* auxfUB1     ): 'upper bounds of the dual variable 0 <= beta_fUB1'];
constCP1 = [constCP1, (beta_deltaLB1      <= beta_deltaLB_LIN .* auxdeltaLB1 ): 'upper bounds of the dual variable 0 <= beta_deltaLB1'];
constCP1 = [constCP1, (beta_deltaUB1      <= beta_deltaUB_LIN .* auxdeltaUB1 ): 'upper bounds of the dual variable 0 <= beta_deltaUB1'];
constCP1 = [constCP1, (eta_gLB1           <= eta_gLB_LIN      .* auxmugLB1   ): 'upper bounds of the dual variable 0 <= eta_gLB1'];
constCP1 = [constCP1, (eta_gUB1           <= eta_gUB_LIN      .* auxmugUB1   ): 'upper bounds of the dual variable 0 <= eta_gUB1'];
constCP1 = [constCP1, (eta_dLB1           <= eta_dLB_LIN      .* auxmudLB1   ): 'upper bounds of the dual variable 0 <= eta_dLB1'];
constCP1 = [constCP1, (eta_dUB1           <= eta_dUB_LIN      .* auxmudUB1   ): 'upper bounds of the dual variable 0 <= eta_dUB1'];
constCP1 = [constCP1, (eta_nuUB1          <= eta_nuUB_LIN     .* auxnuUB1    ): 'upper bounds of the dual variable 0 <= eta_nuUB1'];
constCP1 = [constCP1, (eta_xiLB1          <= eta_xiLB_LIN     .* auxxiLB1    ): 'upper bounds of the dual variable 0 <= eta_xiLB1'];
constCP1 = [constCP1, (eta_xiUB1          <= eta_xiUB_LIN     .* auxxiUB1    ): 'upper bounds of the dual variable 0 <= eta_xiUB1'];

constCP2 = [];
constCP2 = [constCP2, (alpha2(:, 1)       <= alpha12_LIN .* (1 - auxalpha12) ): 'upper bounds of the primal constraint 0 <= alpha2(:, 1)'];
constCP2 = [constCP2, (alpha2(:, 2:end) - alpha2(:, 1:end-1) <= alphab2_LIN .* (1 - auxalphab2)): ...
            'upper bounds of the primal constraint 0 <= alpha2(:, 2:end) - alpha2(:, 1:end-1)'];
constCP2 = [constCP2, (pg                 <= pgLB_LIN    .* (1 - auxpgLB2)   ): 'upper bounds of the primal constraint 0 <= pg'];
constCP2 = [constCP2, (PG_MAX - pg        <= pgUB_LIN    .* (1 - auxpgUB2)   ): 'upper bounds of the primal constraint 0 <= PG_MAX - pg'];
constCP2 = [constCP2, (pd                 <= pdLB_LIN    .* (1 - auxpdLB2)   ): 'upper bounds of the primal constraint 0 <= pd'];
constCP2 = [constCP2, (PD_MAX - pd        <= pdUB_LIN    .* (1 - auxpdUB2)   ): 'upper bounds of the primal constraint 0 <= PD_MAX - pd'];
constCP2 = [constCP2, (F_MAX - Bf * delta <= fUB_LIN     .* (1 - auxfUB2)    ): 'upper bounds of the primal constraint 0 <= F_MAX - Bf * delta'];
constCP2 = [constCP2, (pi + delta         <= deltaLB_LIN .* (1 - auxdeltaLB2)): 'upper bounds of the primal constraint 0 <= pi + delta'];
constCP2 = [constCP2, (pi - delta         <= deltaUB_LIN .* (1 - auxdeltaUB2)):' upper bounds of the primal constraint 0 <= pi - delta'];
constCP2 = [constCP2, (mugLB              <= mugLB_LIN   .* (1 - auxmugLB2)  ): 'upper bounds of the primal constraint 0 <= mugLB'];
constCP2 = [constCP2, (mugUB              <= mugUB_LIN   .* (1 - auxmugUB2)  ): 'upper bounds of the primal constraint 0 <= mugUB'];
constCP2 = [constCP2, (mudLB              <= mudLB_LIN   .* (1 - auxmudLB2)  ): 'upper bounds of the primal constraint 0 <= mudLB'];
constCP2 = [constCP2, (mudUB              <= mudUB_LIN   .* (1 - auxmudUB2)  ): 'upper bounds of the primal constraint 0 <= mudUB'];
constCP2 = [constCP2, (nuUB               <= nuUB_LIN    .* (1 - auxnuUB2)   ): 'upper bounds of the primal constraint 0 <= nuUB'];
constCP2 = [constCP2, (xiLB               <= xiLB_LIN    .* (1 - auxxiLB2)   ): 'upper bounds of the primal constraint 0 <= xiLB'];
constCP2 = [constCP2, (xiUB               <= xiUB_LIN    .* (1 - auxxiUB2)   ): 'upper bounds of the primal constraint 0 <= xiUB'];
constCP2 = [constCP2, (beta_alpha12       <= beta_alpha12_LIN .* auxalpha12  ): 'upper bounds of the dual variable 0 <= beta_alpha12'];
constCP2 = [constCP2, (beta_alphab2       <= beta_alphab2_LIN .* auxalphab2  ): 'upper bounds of the dual variable 0 <= beta_alphab2'];
constCP2 = [constCP2, (beta_gLB2          <= beta_gLB_LIN     .* auxpgLB2    ): 'upper bounds of the dual variable 0 <= beta_gLB2'];
constCP2 = [constCP2, (beta_gUB2          <= beta_gUB_LIN     .* auxpgUB2    ): 'upper bounds of the dual variable 0 <= beta_gUB2'];
constCP2 = [constCP2, (beta_dLB2          <= beta_dLB_LIN     .* auxpdLB2    ): 'upper bounds of the dual variable 0 <= beta_dLB2'];
constCP2 = [constCP2, (beta_dUB2          <= beta_dUB_LIN     .* auxpdUB2    ): 'upper bounds of the dual variable 0 <= beta_dUB2'];
constCP2 = [constCP2, (beta_fUB2          <= beta_fUB_LIN     .* auxfUB2     ): 'upper bounds of the dual variable 0 <= beta_fUB2'];
constCP2 = [constCP2, (beta_deltaLB2      <= beta_deltaLB_LIN .* auxdeltaLB2 ): 'upper bounds of the dual variable 0 <= beta_deltaLB2'];
constCP2 = [constCP2, (beta_deltaUB2      <= beta_deltaUB_LIN .* auxdeltaUB2 ): 'upper bounds of the dual variable 0 <= beta_deltaUB2'];
constCP2 = [constCP2, (eta_gLB2           <= eta_gLB_LIN      .* auxmugLB2   ): 'upper bounds of the dual variable 0 <= eta_gLB2'];
constCP2 = [constCP2, (eta_gUB2           <= eta_gUB_LIN      .* auxmugUB2   ): 'upper bounds of the dual variable 0 <= eta_gUB2'];
constCP2 = [constCP2, (eta_dLB2           <= eta_dLB_LIN      .* auxmudLB2   ): 'upper bounds of the dual variable 0 <= eta_dLB2'];
constCP2 = [constCP2, (eta_dUB2           <= eta_dUB_LIN      .* auxmudUB2   ): 'upper bounds of the dual variable 0 <= eta_dUB2'];
constCP2 = [constCP2, (eta_nuUB2          <= eta_nuUB_LIN     .* auxnuUB2    ): 'upper bounds of the dual variable 0 <= eta_nuUB2'];
constCP2 = [constCP2, (eta_xiLB2          <= eta_xiLB_LIN     .* auxxiLB2    ): 'upper bounds of the dual variable 0 <= eta_xiLB2'];
constCP2 = [constCP2, (eta_xiUB2          <= eta_xiUB_LIN     .* auxxiUB2    ): 'upper bounds of the dual variable 0 <= eta_xiUB2'];

constCP  = [constCP1, constCP2];

% strong stationarity conditions of all MPECs
constD = [constP, constD, constCP];
optsD  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagD  = optimize(constD, -objTP, optsD);

% check the correctness of big-M of the primal constraints
valalpha11 = all(value(alpha1(:, 1))       < alpha11_LIN, 'all'); % 0<= alpha1(:, 1)
valalphab1 = all(value(alpha1(:, 2:end)) - value(alpha1(:, 1:end-1)) < alphab1_LIN, 'all'); % 0 <= alpha1(:, 2:end) - alpha1(:, 1:end-1)
valalpha12 = all(value(alpha2(:, 1))       < alpha12_LIN, 'all'); % 0 <= alpha1(:, 1)
valalphab2 = all(value(alpha2(:, 2:end)) - value(alpha2(:, 1:end-1)) < alphab2_LIN, 'all'); % 0 <= alpha2(:, 2:end) - alpha2(:, 1:end-1)
valpgLB    = all(value(pg)                 < pgLB_LIN,    'all'); % 0 <= pg
valpgUB    = all(PG_MAX - value(pg)        < pgUB_LIN,    'all'); % 0 <= PG_MAX - pg
valpdLB    = all(value(pd)                 < pdLB_LIN,    'all'); % 0 <= pd
valpdUB    = all(PD_MAX - value(pd)        < pdUB_LIN,    'all'); % 0 <= PD_MAX - pd
valfUB     = all(F_MAX - Bf * value(delta) < fUB_LIN,     'all'); % 0 <= F_MAX - Bf * delta
valdeltaLB = all(pi + value(delta)         < deltaLB_LIN, 'all'); % 0 <= pi + delta
valdeltaUB = all(pi - value(delta)         < deltaUB_LIN, 'all'); % 0 <= pi - delta
valmugLB   = all(value(mugLB)              < mugLB_LIN,   'all'); % 0 <= mugLB
valmugUB   = all(value(mugUB)              < mugUB_LIN,   'all'); % 0 <= mugUB
valmudLB   = all(value(mudLB)              < mudLB_LIN,   'all'); % 0 <= mudLB
valmudUB   = all(value(mudUB)              < mudUB_LIN,   'all'); % 0 <= mudUB
valnuUB    = all(value(nuUB)               < nuUB_LIN,    'all'); % 0 <= nuUB
valxiLB    = all(value(xiLB)               < xiLB_LIN,    'all'); % 0 <= xiLB
valxiUB    = all(value(xiUB)               < xiUB_LIN,    'all'); % 0 <= xiUB

% check the correctness of big-M of the dual variables
valbeta_alpha11 = all(value(beta_alpha11) < beta_alpha11_LIN, 'all');                             % 0 <= beta_alpha11
valbeta_alphab1 = all(value(beta_alphab1) < beta_alphab1_LIN, 'all');                             % 0 <= beta_alphab1
valbeta_alpha12 = all(value(beta_alpha12) < beta_alpha12_LIN, 'all');                             % 0 <= beta_alpha12
valbeta_alphab2 = all(value(beta_alphab2) < beta_alphab2_LIN, 'all');                             % 0 <= beta_alphab2
valbeta_gLB     = all([value(beta_gLB1) < beta_gLB_LIN, value(beta_gLB2) < beta_gLB_LIN], 'all'); % 0 <= beta_gLB1 & beta_gLB2
valbeta_gUB     = all([value(beta_gUB1) < beta_gUB_LIN, value(beta_gUB2) < beta_gUB_LIN], 'all'); % 0 <= beta_gUB1 & beta_gUB2
valbeta_dLB     = all([value(beta_dLB1) < beta_dLB_LIN, value(beta_dLB2) < beta_dLB_LIN], 'all'); % ) <= beta_dLB1 & beta_dLB2
valbeta_dUB     = all([value(beta_dUB1) < beta_dUB_LIN, value(beta_dUB2) < beta_dUB_LIN], 'all'); % 0 <= beta_dUB1 & beta_dUB2
valbeta_fUB     = all([value(beta_fUB1) < beta_fUB_LIN, value(beta_fUB2) < beta_fUB_LIN], 'all'); % 0 <= beta_fUB1 & beta_fUB2
valbeta_deltaLB = all([value(beta_deltaLB1) < beta_deltaLB_LIN, value(beta_deltaLB2) < beta_deltaLB_LIN], 'all');
                                                                                                  % 0 <= beta_deltaLB1 & beta_deltaLB2
valbeta_deltaUB = all([value(beta_deltaUB1) < beta_deltaUB_LIN, value(beta_deltaUB2) < beta_deltaUB_LIN], 'all');
                                                                                                  % 0 <= beta_deltaUB1 & beta_deltaUB2
valeta_gLB      = all([value(eta_gLB1)  < eta_gLB_LIN, value(eta_gLB2)   < eta_gLB_LIN],  'all'); % 0 <= eta_gLB1 & eta_gLB2
valeta_gUB      = all([value(eta_gUB1)  < eta_gUB_LIN, value(eta_gUB2)   < eta_gUB_LIN],  'all'); % 0 <= eta_gUB1 & eta_gUB2
valeta_dLB      = all([value(eta_dLB1)  < eta_dLB_LIN, value(eta_dLB2)   < eta_dLB_LIN],  'all'); % 0 <= eta_dLB1 & eta_dLB2
valeta_dUB      = all([value(eta_dUB1)  < eta_dUB_LIN, value(eta_dUB2)   < eta_dUB_LIN],  'all'); % 0 <= eta_dUB1 & eta_dUB2
valeta_nuUB     = all([value(eta_nuUB1) < eta_nuUB_LIN, value(eta_nuUB2) < eta_nuUB_LIN], 'all'); % 0 <= eta_nuUB1 & eta_nuUB2
valeta_xiLB     = all([value(eta_xiLB1) < eta_xiLB_LIN, value(eta_xiLB2) < eta_xiLB_LIN], 'all'); % 0 <= eta_xiLB1 & eta_xiLB2
valeta_xiUB     = all([value(eta_xiUB1) < eta_xiUB_LIN, value(eta_xiUB2) < eta_xiUB_LIN], 'all'); % 0 <= eta_xiUB1 & eta_xiUB2

if ~all([valalpha11 valalphab1 valalpha12 valalphab2 valpgLB valpgUB valpdLB valpdUB valfUB valdeltaLB valdeltaUB ...
        valmugLB valmugUB valmudLB valmudUB valnuUB valxiLB valxiUB ...
        valbeta_alpha11 valbeta_alphab1 valbeta_alpha12 valbeta_alphab2 valbeta_gLB valbeta_gUB valbeta_dLB valbeta_dUB ...
        valbeta_fUB valbeta_deltaLB valbeta_deltaUB valeta_gLB valeta_gUB valeta_dLB valeta_dUB valeta_nuUB valeta_xiLB valeta_xiUB])
    disp('INCORRECT BIG -M!');
end

re_lambda = Cg' * value(lambda);       % LMP at bus n in which unit i located at ($/MWh) (ng x 1)
re_lambda = repmat(re_lambda, [1 no]); % (ng x no)
resultD.profit = sum(re_lambda .* value(pg), 'all') - sum(lambdaG .* value(pg), 'all'); % total profit of all producers ($)
resultD.social = sum(lambdaD .* value(pd), 'all') - sum(lambdaG .* value(pg), 'all');   % social welfare of the market ($)

Parameterized = optimizer(constD, -objTP, optsD, {gamma_DT1, gamma_DT2}, {gamma_bal1, gamma_bal2, pg, lambda});
[X, Y] = meshgrid(0.5:0.05:2, 0.5:0.05:3);
for i = 1:size(X,1)
    for j = 1:size(X,2)
        sol = Parameterized({X(i,j), Y(i,j)});
        vau_gamma_bal1 = sol{1,1};
        vau_gamma_bal11(i,j) = vau_gamma_bal1(1,1);
        vau_gamma_bal2 = sol{1,2};
        vau_gamma_bal21(i,j) = vau_gamma_bal2(1,1);
        vau_pg = sol{1,3};
        vau_lambda = sol{1,4};
        re_lambda = Cg' * vau_lambda;
        re_lambda = repmat(re_lambda, [1 no]);
        vau_profit(i,j) = sum(re_lambda .* vau_pg, 'all') - sum(lambdaG .* vau_pg, 'all');
    end
end
figure(1);
mesh(X, Y, vau_gamma_bal11);
xlim([0.5 2]); ylim([0.5 3]);
xlabel('\gamma^{DT}_{1}'); ylabel('\gamma^{DT}_{2}'); zlabel('\gamma^{bal}_{11}');
print('gamma_bal11', '-dsvg');

figure(2);
mesh(X, Y, vau_gamma_bal21);
xlim([0.5 2]); ylim([0.5 3]);
xlabel('\gamma^{DT}_{1}'); ylabel('\gamma^{DT}_{2}'); zlabel('\gamma^{bal}_{21}');
print('gamma_bal21', '-dsvg');

figure(3);
mesh(X, Y, vau_profit);
xlim([0.5 2]); ylim([0.5 3]);
xlabel('\gamma^{DT}_{1}'); ylabel('\gamma^{DT}_{2}'); zlabel('total profit');
print('total profit', '-dsvg');

disp('Done!');
