% Pricing Non-Convexities in an Electricity Pool
%
% For more details, refer to the paper IEEE TRANS ON POWER SYSTEMS
% url: https://ieeexplore.ieee.org/abstract/document/6153412
% doi: 10.1109/TPWRS.2012.2184562
%
clc; clear; close all;

%% Define Constant
ng = 18;    % number of generating units
no = 4;     % number of generation blocks
nd = 17;    % number of demands
nk = 5;     % number of demand blocks
nb = 24;    % number of buses
nl = 76;    % number of lines
nq = 30;    % number of discretized generation blocks

%% Network Topology
% load system data
case24_ieee_rts;

% connection matrix of bus-gen
gbus = gen(:, 1);
Cg   = sparse(gbus, 1:ng, 1, nb, ng);

% connection matrix of bus-demand
dbus = load(:, 1);
Cd   = sparse(dbus, 1:nd, 1, nb, nd);

% connection matrix of line-bus
f    = branch(:, 1);
t    = branch(:, 2);
i    = [(1:nl)'; (1:nl)'];
Cft  = sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb);

%% Load Data
% data for the gen unit
lambdaG = gen(:, [8 9 10 11]);          % price offer for power block b of unit i ($/MWh) (ng x no)
PGO_MAX = gen(:, [3 4 5 6]);            % upper limit of power block b of unit i (MW) (ng x no)

PG_MIN  = gen(:, 7);                    % minimum power output of unit i (MW) (ng x 1)
PG_MAX  = gen(:, 2);                    % maximum power output of unit i (MW) (ng x 1)

deltaG  = PG_MAX / (nq - 1);            % parameter used to discretize the profit constraint of unit i (MW) (ng x 1)
PG_BAR  = repmat([0:nq-1], [ng 1]) .* repmat(deltaG, [1 nq]); 
                                        % discrete power produced by block q of unit i (MW) (ng x nq)

KU      = gen(:, 12);                   % startup cost of unit i ($) (ng x 1)
KD      = gen(:, 13);                   % shutdown cost of unit i ($) (ng x 1)

v0      = gen(:, 14);                   % initial on/off status of unit i (ng x 1)

% data for the demand
lambdaD = load(:, [7 8 9 10 11]);       % price bid for power block k of demand j ($/MWh) (nd x nk)
PD_MAX  = load(:, [2 3 4 5 6]);         % upper limit of power block k of demand j (MW) (nd x nk)

% data for the network
Bf      = spdiags(branch(:, 3), 0, nl, nl) * Cft;
                                        % Bf * Va is the vector of branch flow injected at each branch's "from" bus (S) (nl x nb)
Bbus    = 1/2 * (Cft' * Bf);            % Bbus * Va is the vector of real power injections at each bus (S) (nb x nb)
F_MAX   = branch(:, 4);                 % transmission capacity of line n-m (MW) (nl x 1)

% data for linearization
G       = 10 * max(max(lambdaG, [], 'all'), max(lambdaD, [], 'all'));

%% Define Variable
% primal variable
pg     = sdpvar(ng, no, 'full'); % power produced by block b of unit i (MW)
pd     = sdpvar(nd, nk, 'full'); % power consumed by block k of demand j (MW)
delta  = sdpvar(nb, 1,  'full'); % voltage angle of bus n (rad)
v      = binvar(ng, 1,  'full'); % on/off status of unit i
cu     = sdpvar(ng, 1,  'full'); % startup cost of unit i ($)
cd     = sdpvar(ng, 1,  'full'); % shutdown cost of unit i ($)

% dual variable
lambda = sdpvar(nb, 1,  'full'); % dual of power balance at bus n ($/MWh)
mugLB  = sdpvar(ng, no, 'full'); % dual of lower bound of block b of unit i ($/MWh)
mugUB  = sdpvar(ng, no, 'full'); % dual of upper bound of block b of unit i ($/MWh)
mudLB  = sdpvar(nd, nk, 'full'); % dual of lower bound of block k of demand j ($/MWh)
mudUB  = sdpvar(nd, nk, 'full'); % dual of upper bound of block k of demand j ($/MWh)
muLB   = sdpvar(ng, 1,  'full'); % dual of lower bound of unit i ($/MWh)
nuUB   = sdpvar(nl, 1,  'full'); % dual of transmission capacity limit of line n-m ($/MWh)
xiLB   = sdpvar(nb, 1,  'full'); % dual of lower bound of voltage angle at bus n ($/h)
xiUB   = sdpvar(nb, 1,  'full'); % dual of upper bound of voltage angle at bus n ($/h)
xi1    = sdpvar(nb, 1,  'full'); % dual of voltage angle at reference bus n = 1 ($/h)
xi1(2:end, :) = 0;
etaU   = sdpvar(ng, 1,  'full'); % dual of startup cost of unit i (p.u.)
etaU0  = sdpvar(ng, 1,  'full'); % dual of nonnegative of startup cost of unit i (p.u.)
etaD   = sdpvar(ng, 1,  'full'); % dual of shutdown cost of unit i (p.u.)
etaD0  = sdpvar(ng, 1,  'full'); % dual of nonnegative of shutdown cost of unit i (p.u.)
betaLB = sdpvar(ng, 1,  'full'); % dual of lower bound of relaxed on/off status of unit i ($/h)
betaUB = sdpvar(ng, 1,  'full'); % dual of upper bound of relaxed on/off status of unit i ($/h)

% auxiliary variable
x      = binvar(ng, nq, 'full'); % binary variable used to discretize the power produced by block q of unit i
z      = sdpvar(ng, nq, 'full'); % supplementary variable used to discretize the power produced by block q of unit i

%% Conventional Market-clearing
objC = sum(lambdaG .* pg, 'all') + sum(cu, 'all') + sum(cd, 'all') - sum(lambdaD .* pd, 'all');

constC = [];
constC = [constC, (sum(Cd * pd, 2) -sum(Cg * pg, 2) + Bbus * delta == 0): 'power balacne'];
constC = [constC, (0 <= pg <= repmat(v, [1 no]) .* PGO_MAX): 'bound on pg'];
constC = [constC, (0 <= pd <= PD_MAX): 'bound on pd'];
constC = [constC, (v .* PG_MIN <= sum(pg, 2)): 'bound on pg'];
constC = [constC, (Bf * delta <= F_MAX): 'transmission capacity limit'];
constC = [constC, (-pi <= delta <= pi): 'bound on delta'];
constC = [constC, (delta(1, 1) == 0): 'reference bus'];
constC = [constC, (cu >= (v - v0) .* KU): 'startup cost'];
constC = [constC, (cu >= 0): 'startup cost'];
constC = [constC, (cd >= (v0 - v) .* KD): 'shutdown cost'];
constC = [constC, (cd >= 0): 'shutdown cost'];

optsC = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagC = optimize(constC, objC, optsC);

fixv = value(v);
resultC.v = fixv;
resultC.socwel = -value(objC);

constF = [];
constF = [constF, constC];
constF = [constF, v == fixv];

optsF = sdpsettings('solver', 'gurobi', 'verbose', 2, 'relax', 2);
diagF = optimize(constF, objC, optsF);

resultC.profit = sum(repmat(Cg' * dual(constF(1)), [1 no]) .* value(pg), 2) - sum(lambdaG .* value(pg), 2) - value(cu) - value(cd);
resultC.prod = sum(value(pg), 2);

%% Relaxed Duality
objD = -sum(PD_MAX .* mudUB, 'all') - sum(F_MAX .* nuUB, 'all') - pi * sum(xiLB, 'all') - pi * sum(xiUB, 'all') ...
    - sum(betaUB, 'all') - sum(v0 .* KU .* etaU, 'all') + sum(v0 .* KD .* etaD, 'all');

constD = [];
constD = [constD, (lambdaG - repmat(Cg' * lambda, [1 no]) + mugUB - mugLB - repmat(muLB, [1 no]) == 0): 'diff L w.r.t. pg'];
constD = [constD, (-lambdaD + repmat(Cd' * lambda, [1 nk]) + mudUB - mudLB == 0): 'diff L w.r.t. pd'];
constD = [constD, (Bbus * lambda + Bf' * nuUB + xiUB - xiLB + xi1 == 0): 'diff L w.r.t. delta'];
constD = [constD, (1 - etaU - etaU0 == 0): 'diff L w.r.t. cu'];
constD = [constD, (1 - etaD - etaD0 == 0): 'diff L w.r.t. cd'];
constD = [constD, (muLB .* PG_MIN - sum(mugUB .* PGO_MAX, 2) + KU .* etaU - KD .* etaD + betaUB - betaLB == 0): 'diff L w.r.t. v'];
constD = [constD, ...
    (mugUB >= 0): 'mugUB >= 0', (mugLB >= 0): 'mugLB >= 0', ...
    (mudUB >= 0): 'mudUB >= 0', (mudLB >= 0): 'mudLB >= 0', ...
    (muLB >= 0): 'muLB >= 0', ...
    (nuUB >= 0): 'nuUB >= 0', ...
    (xiUB >= 0): 'xiUB >= 0', (xiLB >= 0): 'xiLB >= 0', ...
    (etaU >= 0): 'etaU >= 0', (etaU0 >= 0): 'etaU0 >= 0', ...
    (etaD >= 0): 'etaD >= 0', (etaD0 >= 0): 'etaD0 >= 0', ...
    (betaUB >= 0): 'betaUB >= 0', (betaLB >= 0): 'betaLB >= 0'];

optsD = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagD = optimize(constD, -objD, optsD);

constP = [];
constP = [constP, (sum(pg, 2) - deltaG <= sum(PG_BAR .* x, 2) <= sum(pg, 2)): '(11a)'];
constP = [constP, (sum(x, 2) == 1): '(11b)'];
constP = [constP, (0 <= repmat(Cg' * lambda, [1 nq]) - z <= G * (1 - x)): '(11c)'];
constP = [constP, (0 <= z <= G * x): '(11d)'];
constP = [constP, (sum(z .* PG_BAR, 2) - sum(lambdaG .* pg, 2) - cu - cd >= 0): '(11e)'];
% constP = [sum(repmat(Cg' * lambda, [1 no]) .* pg, 2) - sum(lambdaG .* pg, 2) - cu - cd >= 0]; % bilinear

objR = objC - objD;
constR = [constC, constD, constP];
optsR = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagR = optimize(constR, objR, optsR);

valG = all(repmat(Cg' * value(lambda), [1 nq]) - value(z) < G, 'all') & all(value(z) < G, 'all');

resultR.v = value(v);
resultR.socwel = -value(objC);
resultR.prod = sum(value(pg), 2);
resultR.gap = value(objR);
resultR.profit = sum(repmat(Cg' * value(lambda), [1 no]) .* value(pg), 2) - sum(lambdaG .* value(pg), 2) - value(cu) - value(cd);


disp('DONE!');