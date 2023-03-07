% Strategic Offering for a Wind Power Producer
%
% For more details, refer to the paper IEEE TRANS ON POWER SYSTEMS
% url: https://ieeexplore.ieee.org/abstract/document/6571288
% doi: 10.1109/TPWRS.2013.2273276
%
clc; clear; close all;

%% Define Constant
nt   = 24;              % number of time periods
ns   = 10 ;             % number of scenarios
nw   = 1;               % number of wind units
ng   = 3;               % number of gen units
nd   = 3;               % number of demands
no   = 4;               % number of gen blocks
nj   = 3;               % number of demand blocks
nb   = 3;               % number of buses
nl   = 3;               % number of lines
prob = 1 / ns;          % weight of scenario
CONG = false;           % transmission congestion

%% Network Topology
% load system data
case3;

% congested network
if CONG
    branch(:, 4) = 30;
    disp('limit transmission capacity to cause congestion');
end

% connection matrix of bus-wind
wbus    = wind(:, 1);
Cw      = sparse(wbus, 1 : nw, 1, nb, nw);

% connection matrix of bus-gen
gbus    = gen(:, 1);
Cg      = sparse(gbus, 1 : ng, 1, nb, ng);

% connection matrix of bus-demand
dbus    = load(:, 1);
Cd      = sparse(dbus, 1 : nd, 1, nb, nd);

% connection matrix of line-bus
f       = branch(:, 1);
t       = branch(:, 2);
i       = [(1 : nl)'; (1 : nl)'];
Cft     = sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb);

%% Load Data
% data for the wind units
lambdaW = zeros(nt, nw);                 % marginal cost of the lth wind unit in period t ($/MWh) (nt x nw)

PW_MAX  = 300 * ones(nt, nw);            % maximum wind power output of the lth wind unit in period t (MW) (nt x nw)

factW(:, 1, :) = xlsread('SourceData', 'Wind Power Capacity Factors', 'B2:K25');
                                         % wind power capacity factor of the lth wind unit ...
                                         % in scenario s in period t (p.u.) (nt x nw x ns)

avefactW = sum(factW, 3) ./ ns;          % averaged wind power capacity factor over scenarios ...
                                         % of the lth wind unit in period t (p.u.) (nt x nw)

PW_EXP   = avefactW .* PW_MAX;           % expected wind power output of the lth wind unit in period t (MW) (nt x nw)

PW_SCEN  = factW .* repmat(PW_MAX, [1 1 ns]);
                                         % wind power production by lth wind unit ...
                                         % in scenario s in period t (MW) (nt x nw x ns)

lambdaWGC = wind(:, 4);                  % (nw x 1)
lambdaWGC = repmat(lambdaWGC, [1 nt]);   % (nw x nt)
lambdaWGC = transpose(lambdaWGC);        % marginal cost of green certificate of the lth wind unit in period t ($/MWh or $/p) (nt x nw)

% data for the gen units
lambdaG = gen(:, [6 7 8 9]);             % (ng x no)
lambdaG = repmat(lambdaG, [1 1 nt]);     % (ng x no x nt)
lambdaG = permute(lambdaG, [3 1 2]);     % marginal cost of the bth block of the ith gen unit in period t ($/MWh) (nt x ng x no)

PG_MAX  = gen(:, [2 3 4 5]);             % (ng x no)
PG_MAX  = repmat(PG_MAX, [1 1 nt]);      % (ng x no x nt)
PG_MAX  = permute(PG_MAX, [3 1 2]);      % upper limit of the bth block of the ith gen unit in period t (MW) (nt x ng x no)

% data for the demands
lambdaD = load(:, [5 6 7]);              % (nd x nj)
lambdaD = repmat(lambdaD, [1 1 nt]);     % (nd x nj x nt)
lambdaD = permute(lambdaD, [3 1 2]);     % marginal utility of the jth block of the dth demand in period t ($/MWh) (nt x nd x nj)

PD_MAX  = load(:, [2 3 4]);              % (nd x nj)
PD_MAX  = repmat(PD_MAX, [1 1 nt]);      % (nd x nj x nt)
PD_MAX  = permute(PD_MAX, [3 1 2]);      % (nt x nd x nj)
factD   = [
    0.829 0.777 0.745 0.736 0.728 0.732 0.736 0.829 0.892 0.908 0.932 0.935 ...
    0.917 0.905 0.870 0.869 0.864 0.870 0.881 0.974 1.000 0.981 0.926 0.915
];                                       % hourly demand factros (p.u.) (1 x nt)
factD   = repelem(factD, nd, nj);        % (nd x (nj * nt))
factD   = reshape(factD, [nd nj nt]);    % (nd x nj x nt)
factD   = permute(factD, [3 1 2]);       % (nt x nd x nj)
PD_MAX  = factD .* PD_MAX;               % upper limit of the jth block of the dth demand in period t (MW) (nt x nd x nj)

ratio = load(:, 8);                      % (nd x 1)
ratio = repmat(ratio, [1 nt nj]);        % (nd x nt x nj)
ratio = permute(ratio, [2 1 3]);         % renewable energy consumption weight (%) (nt x nd x nj)

lambdaDGC = load(:, [9 10 11]);          % (nd x nj)
lambdaDGC = repmat(lambdaDGC, [1 1 nt]); % (nd x nj x nt)
lambdaDGC = permute(lambdaDGC, [3 1 2]); % marginal utility of green certificate of 
                                         % the jth block of the dth demand in period t ($/MWh or $/p) (nt x nd x nj)
% data for the network
B       = branch(:, 3);                  % susceptance of link (S) (nl x 1)

F_MAX   = branch(:, 4);                  % (nl x 1)
F_MAX   = repmat(F_MAX, [1 nt]);         % transmission capacity of line k in period t (MW) (nl x nt)

sla     = 1;
numsla  = length(sla);

nsla    = setdiff([1 : nb]', sla);
numnsla = length(nsla);

% data for the balancing market
lambdaB = xlsread('SourceData', 'Balancing Market Prices', 'B2:K25');
                                         % balancing market price in scenario s in period t ($/MWh) (nt x ns)
lambdaB = repmat(lambdaB, [1 1 nw]);     % (nt x ns x nw)
lambdaB = permute(lambdaB, [1 3 2]);     % (nt x nw x ns)

%% Define Variable
% upper-level variables
alphaw     = sdpvar(nt, nw,      'full'); % offer price of the lth wind unit in period t ($/MWh)
alphawGC   = sdpvar(nt, nw,  ns, 'full'); % offer price of green certificate of the lth wind unit in period t ($/MWh or $/p)
pwof       = sdpvar(nt, nw,      'full'); % offer quantity of the lth wind unit in period t (MW)
pwp        = sdpvar(nt, nw,  ns, 'full'); % wind power produced by the lth wind power unit in scenario s in period t (MW)
pwbal      = sdpvar(nt, nw,  ns, 'full'); % power brought from/sold to the balancing market by the lth wind power unit ...
                                          % in scenario s in period t (MW)
% lower-level variables of energy market
pw         = sdpvar(nt, nw,      'full'); % wind power cleared in the DA market for the lth wind unit in period t (MW)
pg         = sdpvar(nt, ng,  no, 'full'); % power scheduled to be produced by the bth block of the ith gen unit in period t (MW)
pd         = sdpvar(nt, nd,  nj, 'full'); % power scheduled to be consumed by the jth block of the dth demand in period t (MW)
f          = sdpvar(nl, nt,      'full'); % power flow through line k in period t (MW)
delta      = sdpvar(nb, nt,      'full'); % voltage angle at bus n in period t (rad)

lambda     = sdpvar(nb, nt,      'full'); % dual of power balance at bus n in period t ($/MWh)
phi        = sdpvar(nl, nt,      'full'); % dual of power flow at line k in period t ($/MWh)
phiLB      = sdpvar(nl, nt,      'full'); % dual of transmission capacity of line k in period t ($/MWh)
phiUB      = sdpvar(nl, nt,      'full'); % dual of transmission capacity of line k in period t ($/MWh)
etaLB      = sdpvar(nt, nd,  nj, 'full'); % dual of lower bound of the jth block of the dth demand in period t ($/MWh)
etaUB      = sdpvar(nt, nd,  nj, 'full'); % dual of upper bound of the jth block of the dth demand in period t ($/MWh)
psiLB      = sdpvar(nt, ng,  no, 'full'); % dual of lower bound of the bth block of the lth gen unit in period t ($/MWh)
psiUB      = sdpvar(nt, ng,  no, 'full'); % dual of upper bound of the bth block of the lth gen unit in period t ($/MWh)
sigmaLB    = sdpvar(nt, nw,      'full'); % dual of lower bound of the lth wind unit in period t ($/MWh)
sigmaUB    = sdpvar(nt, nw,      'full'); % dual of upper bound of the lth wind unit in period t ($/MWh)
xiLB       = sdpvar(numnsla, nt, 'full'); % dual of lower bound of the voltage angle at bus n in period t ($/h)
xiUB       = sdpvar(numnsla, nt, 'full'); % dual of upper bound of the voltage angle at bus n in period t ($/h)
xi         = sdpvar(numsla,  nt, 'full'); % dual of voltage angle at slack bus in period t ($/h)

% lower-level variables of green certificate market
tw         = sdpvar(nt, nw,      'full'); % green certificate sold by the lth wind unit in period t (MW or p)
td         = sdpvar(nt, nd,  nj, 'full'); % green certificate brought by the jth block of the dth demand in period t (MW or p)

lambdaGC   = sdpvar(nt, 1,       'full'); % dual of green certificate balance in period t ($/MWh or $/p)
muLB       = sdpvar(nt, nw,      'full'); % dual of lower bound of green certificate of the lth wind unit in period t ($/MWh or $/p)
muUB       = sdpvar(nt, nw,      'full'); % dual of upper bound of green certificate of the lth wind unit in period t ($/MWh or $/p)
rhoLB      = sdpvar(nt ,nd,  nj, 'full'); % dual of lower bound of green certificate of the jth block of the dth demand in period t ($/MWh ot $/p)
rhoUB      = sdpvar(nt, nd,  nj, 'full'); % dual of upper bound of green certificate of the jth block of the dth demand in period t ($/MWh or $/p)

% auxiliart binary variables of energy market
auxfLB     = binvar(nl, nt,      'full'); % binary variables for the CP (F_MAX + f) * phiLB = 0
auxfUB     = binvar(nl, nt,      'full'); % binary variables for the CP (F_MAX - f) * phiUB = 0
auxpdLB    = binvar(nt, nd,  nj, 'full'); % binary variables for the CP pd * etaLB = 0
auxpdUB    = binvar(nt, nd,  nj, 'full'); % binary variables for the CP (PD_MAX - pd) * etaUB = 0
auxpgLB    = binvar(nt, ng,  no, 'full'); % binary variables for the CP pg * psiLB = 0
auxpgUB    = binvar(nt, ng,  no, 'full'); % binary variables for the CP (PG_MAX - pg) * psiUB = 0
auxpwLB    = binvar(nt, nw,      'full'); % binary variables fot the CP pw * sigmaLB = 0
auxpwUB    = binvar(nt, nw,      'full'); % binary variables for the CP (pwof - pw) * sigmaUB = 0
auxdeltaLB = binvar(numnsla, nt, 'full'); % binary variables for the CP (pi + delta(nsla, :)) * xiLB = 0
auxdeltaUB = binvar(numnsla, nt, 'full'); % binary variables for the CP (pi - delta(nsla, :)) * xiUB = 0

% auxiliary binary variables of green certificate market
auxtwLB    = binvar(nt, nw,      'full'); % binary variables for the CP tw * muLB = 0
auxtwUB    = binvar(nt, nw,      'full'); % binary variables for the CP (pwp - tw) * muUB = 0
auxtdLB    = binvar(nt, nd,  nj, 'full'); % binary variables for the CP td * rhoLB = 0
auxtdUB    = binvar(nt, nd,  nj, 'full'); % binary variables for the CP (ratio * PD_MAX - td) * rhoUB = 0

%% Solve the Primal Lower-level Problem
% 1): energy market
% minus social welfare ($)
objP   = sum(lambdaG .* pg, 'all') + sum(lambdaW .* pw, 'all') - sum(lambdaD .* pd, 'all');

constP = [];

re_pg  = permute(pg, [2 3 1]);          % (ng x no x nt)
re_pg  = reshape(re_pg, [ng, no * nt]); % (ng x (no * nt))
sum_pg = Cg * re_pg;                    % (nb x (no * nt))
sum_pg = reshape(sum_pg, [nb no nt]);   % (nb x no x nt)
sum_pg = squeeze(sum(sum_pg, 2));       % sum of the bth block of the ith gen unit at bus n in period t (MW) (nb x nt)

sum_pw = Cw * pw';                      % sum of the lth wind production at bus n in period t (MW) (nb x nt)

re_pd  = permute(pd, [2 3 1]);          % (nd x nj x nt)
re_pd  = reshape(re_pd, [nd, nj * nt]); % (nd x (nj * nt))
sum_pd = Cd * re_pd;                    % (nb x (nj * nt))
sum_pd = reshape(sum_pd, [nb nj nt]);   % (nb x nj x nt)
sum_pd = squeeze(sum(sum_pd, 2));       % sum of the jth block of the dth demand at bus n in period t (MW) (nb x nt)

sum_f  = Cft' * f;                      % sum of power flow leaving bus n in period t (MW) (nb x nt)

constP = [constP, (sum_pg + sum_pw - sum_f == sum_pd): 'power balance at bus n in period t in DA market'];
constP = [constP, (f == diag(B) * Cft * delta): 'def of power flow'];
constP = [constP, (-F_MAX <= f <= F_MAX): 'transmission capacity'];
constP = [constP, (0 <= pd <= PD_MAX): 'power bounds on pd'];
constP = [constP, (0 <= pg <= PG_MAX): 'power bounds on pg'];
constP = [constP, (0 <= pw <= PW_EXP): 'power bounds on pw'];
constP = [constP, (-pi <= delta(nsla, :) <= pi): 'bounds on delta'];
constP = [constP, (delta(sla, :) == 0): 'def of slack bus'];

optsP  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagP  = optimize(constP, objP, optsP);

% wind power cleared in DA market (MW) (nt x nw)
resultP.prod = value(pw);
% difference of wind power level in the DA market and real wind power output of the lth wind in scenario s in period t (MW) (nt x nw x ns)
diffbal = repmat(value(pw), [1 1 ns]) - PW_SCEN;
% wind power profit when submited expected wind power output to the DA market and anticipating balancing market ($)
resultP.profit = sum(-dual(constP(1)) .* value(sum_pw), 'all') - prob * sum(lambdaB .* diffbal, 'all');

fprintf('total production of wind units behave unstrategically solved using primal form is %d MW\n', sum(resultP.prod));
fprintf('total profit     of wind units behave unstrategically solved using primal form is %d $\n\n', resultP.profit);

% 2): green certificate market
% minus social welfare ($)
objPGC   = sum(lambdaWGC .* tw, 'all') - sum(lambdaDGC .* td, 'all');

constPGC = [];
constPGC = [constPGC, (sum(tw, 2) == sum(td, [2 3])): 'green certificate balance in period t'];
constPGC = [constPGC, (0 <= tw <= PW_EXP): 'green certificate bounds on tw'];
constPGC = [constPGC, (0 <= td <= ratio .* PD_MAX): 'green certificate bounds on td'];

optsPGC  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagPGC  = optimize(constPGC, objPGC, optsPGC);

% green certificate traded between wind unit and demand (MW or p)
resultPGC.prod = value(tw);
% wind power profit of trading green certificate ($)
resultPGC.profit = sum(repmat(-dual(constPGC(1)), [1 nw]) .* value(tw), 'all') - sum(lambdaWGC .* value(tw), 'all');

fprintf('total traded green certificate of wind units behave unstrategically solved using primal form is %d MW or p\n', sum(resultPGC.prod));
fprintf('total profit of green certificate of wind units behave unstragically solved using primal form is %d $\n\n', resultPGC.profit);

%% Solve the Dual Lower-level Problem
% 1): energy market
% positive social welfare ($)
objD       = - sum((phiLB + phiUB) .* F_MAX, 'all') - sum(etaUB .* PD_MAX, 'all') - sum(psiUB .* PG_MAX, 'all') ...
             - sum(sigmaUB .* PW_EXP, 'all') - pi * sum(xiLB + xiUB, 'all');

constD     = [];

re_lambdag = Cg' * lambda;                 % n(i): node at which gen i located at ($/MWh) (ng x nt)
re_lambdag = repmat(re_lambdag, [1 1 no]); % (ng x nt x no)
re_lambdag = permute(re_lambdag, [2 1 3]); % (nt x ng x no)
constD     = [constD, (lambdaG - re_lambdag + psiUB - psiLB == 0): 'stationary condtion w.r.t. pg'];

re_lambdaw = Cw' * lambda;                 % n(l): node at which wind unit l located at ($/MWh) (nw x nt)
re_lambdaw = transpose(re_lambdaw);        % (nt x nw)
constD     = [constD, (lambdaW - re_lambdaw + sigmaUB - sigmaLB == 0): 'stationary condition w.r.t. pw'];

re_lambdad = Cd' * lambda;                 % n(d): node at which demand d located at ($/MWh) (nd x nt)
re_lambdad = repmat(re_lambdad, [1 1 nj]); % (nd x nt x nj)
re_lambdad = permute(re_lambdad, [2 1 3]); % (nt x nd x nj)
constD     = [constD, (-lambdaD + re_lambdad + etaUB - etaLB == 0):'stationary condition w.r.t. pd'];

constD     = [constD, (Cft * lambda - phi + phiUB - phiLB == 0): 'stationary condition w.r.t. f'];

re_phi     = Cft' * diag(B) * phi;         % sum_k|s(k)=n{B_k * phi_kt} - sum_k|r(k)=n{b_k*phi_kt} ($/h) (nb x nt)
constD     = [constD, (re_phi(nsla, :) + xiUB - xiLB == 0): 'stationary condition w.r.t. delta'];
constD     = [constD, (re_phi(sla, :) - xi == 0): 'stationary condition w.r.t. delta'];

% nonnegativity of dual variables associated with inequality constraints
constD     = [
    constD, ...
    (0 <= phiLB):   'nonnegative of phiLB',   (0 <= phiUB):   'nonnegative of phiUB', ...
    (0 <= etaLB):   'nonnegative of etaLB',   (0 <= etaUB):   'nonnegative of etaUB', ...
    (0 <= psiLB):   'nonnegative of psiLB',   (0 <= psiUB):   'nonnegative of psiUB', ...
    (0 <= sigmaLB): 'nonnegative of sigmaLB', (0 <= sigmaUB): 'nonnegative of sigmaUB', ...
    (0 <= xiLB):    'nonnegative of xiLB',    (0 <= xiUB):    'nonnegative of xiUB'
];

optsD      = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagD      = optimize(constD, -objD, optsD);

% wind power cleared in DA market (MW) (nt x nw)
resultD.prod = -dual(constD(2));
% difference of wind power level in the DA market and real wind power output of the lth wind in scenario s in period t (MW) (nt x nw x ns)
diffbal = repmat(-dual(constD(2)), [1 1 ns]) - PW_SCEN;
% wind power proft when submited expected wind power output to the DA market and anticipating balancing market ($)
resultD.profit = sum(value(re_lambdaw) .* -dual(constD(2)), 'all') - prob * sum(lambdaB .* diffbal, 'all');

fprintf('total production of wind units behave unstrategically solved using dual form is %d MW\n', sum(resultD.prod));
fprintf('total profit     of wind units behave unstrategically solved using dual form is %d $\n\n', resultD.profit);

% 2): green certificate market
% positive social welfare ($)
objDGC   = -sum(muUB .* PW_EXP, 'all') - sum(rhoUB .* ratio .* PD_MAX, 'all');

constDGC = [];
constDGC = [constDGC, (lambdaWGC - repmat(lambdaGC, [1 nw]) + muUB - muLB == 0): 'stationary condition w.r.t. tw'];
constDGC = [constDGC, (-lambdaDGC + repmat(lambdaGC, [1 nd nj]) + rhoUB - rhoLB == 0): 'stationary condition w.r.t. td'];
constDGC = [
    constDGC, ...
    (0 <= muLB): 'nonnegative of muLB', (0 <= muUB): 'nonnegative of muUB', ...
    (0 <= rhoLB): 'nonnegative of rhoLB', (0 <= rhoUB): 'nonnegative of rhoUB'
];

optsDGC = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagDGC = optimize(constDGC, -objDGC, optsDGC);

% green certificate traded between wind unit and demand (MW or p)
resultDGC.prod = -dual(constDGC(1));
% wind power profit of trading green certificate ($)
resultDGC.profit = sum(repmat(value(lambdaGC), [1 nw]) .* -dual(constDGC(1)), 'all') - sum(lambdaWGC .* -dual(constDGC(1)), 'all');

fprintf('total traded green certificate of wind units behave unstrategically solved using dual form is %d MW or p\n', sum(resultDGC.prod));
fprintf('total profit of green certificate of wind units behave unstragically solved using dual form is %d $\n\n', resultDGC.profit);

%% Big-M Value
% constant related to primal constraints of energy market
fLB_LIN     = 2 * F_MAX + 1; % upper bounds of the primal constraint 0 <= F_MAX + f           (MW)  (nl x nt)
fUB_LIN     = 2 * F_MAX + 1; % upper bounds of the primal constraint 0 <= F_MAX - f           (MW)  (nl x nt)
pdLB_LIN    = PD_MAX    + 1; % upper bounds of the primal constraint 0 <= pd                  (MW)  (nt x nd x nj)
pdUB_LIN    = PD_MAX    + 1; % upper bounds of the primal constraint 0 <= PD_MAX - pd         (MW)  (nt x nd x nj)
pgLB_LIN    = PG_MAX    + 1; % upper bounds of the primal constraint 0 <= pg                  (MW)  (nt x ng x no)
pgUB_LIN    = PG_MAX    + 1; % upper bounds of the primal constraint 0 <= PG_MAX - pg         (MW)  (nt x ng x no)
pwLB_LIN    = PW_MAX    + 1; % upper bounds of the primal constraint 0 <= pw                  (MW)  (nt x nw)
pwUB_LIN    = PW_MAX    + 1; % upper bounds of the primal constraint 0 <= pwof - pw           (MW)  (nt x nw)
deltaLB_LIN = 2 * pi;        % upper bounds of the primal constraint 0 <= pi + delta(nsla, :) (rad) (numnsla x nt)
deltaUB_LIN = 2 * pi;        % upper bounds of the primal constraint 0 <= pi - delta(nsla, :) (rad) (numnsla x nt)

% constants related to primal constraints of green certificate market
twLB_LIN    = PW_MAX + 1;    % upper bounds of the primal constraint 0 <= tw                  (MW or p) (nt x nw)
twUB_LIN    = PW_MAX + 1;    % upper bounds of the primal constraint 0 <= pwp - tw            (MW or p) (nt x nw)
tdLB_LIN    = PD_MAX + 1;    % upper bounds of the primal constraint 0 <= td                  (MW or p) (nt x nd x nj)
tdUB_LIN    = PD_MAX + 1;    % upper bounds of the primal constraint 0 <= ratio * PD_MAX - td (MW or p) (nt x nd x nj)

% constant related to dual variables of energy market
phiLB_LIN   = (value(phiLB) + 1)   * 100; % upper bounds of the dual variable 0 <= phiLB   ($/MWh) (nl x nt)
phiUB_LIN   = (value(phiUB) + 1)   * 100; % upper bounds of the dual variable 0 <= phiUB   ($/MWh) (nl x nt)
etaLB_LIN   = (value(etaLB) + 1)   * 100; % upper bounds of the dual variable 0 <= etaLB   ($/MWh) (nt x nd x nj)
etaUB_LIN   = (value(etaUB) + 1)   * 100; % upper bounds of the dual variable 0 <= etaUB   ($/MWh) (nt x nd x nj)
psiLB_LIN   = (value(psiLB) + 1)   * 100; % upper bounds of the dual variable 0 <= psiLB   ($/MWh) (nt x ng x no)
psiUB_LIN   = (value(psiUB) + 1)   * 100; % upper bounds of the dual variable 0 <= psiUB   ($/MWh) (nt x ng x no)
sigmaLB_LIN = (value(sigmaLB) + 1) * 100; % upper bounds of the dual variable 0 <= sigmaLB ($/MWh) (nt x nw)
sigmaUB_LIN = (value(sigmaUB) + 1) * 100; % upper bounds of the dual variable 0 <= sigmaUB ($/MWh) (nt x nw)
xiLB_LIN    = (value(xiLB) + 1)    * 100; % upper bounds of the dual variable 0 <= xiLB    ($/h)   (numnsla x nt)
xiUB_LIN    = (value(xiUB) + 1)    * 100; % upper bounds of the dual variable 0 <= xiUB    ($/h)   (numnsla x nt)

% constant related to dual variables of green certificate market
muLB_LIN    = (value(muLB) + 1)    * 100; % upper bounds of the dual variable 0 <= muLB  ($/MWh or $/p) (nt x nw)
muUB_LIN    = (value(muUB) + 1)    * 100; % upper bounds of the dual variable 0 <= muUB  ($/MWh or $/p) (nt x nw)
rhoLB_LIN   = (value(rhoLB) + 1)   * 100; % upper bounds of the dual variable 0 <= rhoLB ($/MWh or $/p) (nt x nd x nj)
rhoUB_LIN   = (value(rhoUB) + 1)   * 100; % upper bounds of the dual variable 0 <= rhoUB ($/MWh or $/p) (nt x nd x nj)

%% Solve KKT Systems of the Lower-level Problem
% 1): energy market
% strong duality
constSD = [];
constSD = [constSD, (constP): 'primal feasibility'];
constSD = [constSD, (constD): 'dual feasibility'];
constSD = [constSD, (objP == objD): 'strong duality'];
optsSD  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagSD  = optimize(constSD, [], optsSD);

% wind power cleared in DA market (MW) (nt x nw)
resultSD.prod = value(pw);
% locational marginal prices of bus n at which wind unit l located at ($/MWh)
LMPw = Cw' * value(lambda); % (nw x nt)
LMPw = transpose(LMPw);     % (nt x nw)
% difference of wind power level in the DA market and real wind power output of the lth wind in scenario s in period t (MW) (nt x nw x ns)
diffbal = repmat(value(pw), [1 1 ns]) - PW_SCEN;
% wind power proft when submited expected wind power output to the DA market and anticipating balancing market ($)
resultSD.profit = sum(LMPw .* value(pw), 'all') - prob * sum(lambdaB .* diffbal, 'all');

fprintf('total production of wind units behave unstrategically solved using strong duality is %d MW\n', sum(resultSD.prod));
fprintf('total profit     of wind units behave unstrategically solved using strong duality is %d $\n\n', resultSD.profit);

% 2): green certificate market
% strong duality
constSDGC = [];
constSDGC = [constSDGC, (constPGC): 'primal feasibility'];
constSDGC = [constSDGC, (constDGC): 'dual feasibility'];
constSDGC = [constSDGC, (objPGC == objDGC): 'strong duality'];
optsSDGC  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagSDGC  = optimize(constSDGC, [], optsSDGC);

% green certificate traded between wind unit and demand (MW or p)
resultSDGC.prod = value(tw);
% wind power profit of trading green certificate ($)
resultSDGC.profit = sum(repmat(value(lambdaGC), [1 nw]) .* value(tw), 'all') - sum(lambdaWGC .* value(tw), 'all');

fprintf('total traded green certificate of wind units behave unstrategically solved using strong duality is %d MW or p\n', sum(resultSDGC.prod));
fprintf('total profit of green certificate of wind units behave unstragically solved using strong duality is %d $\n\n', resultSDGC.profit);

% 1): energy market
% complementarity condition
constCP    = [];
constCP    = [constCP, ((F_MAX + f)           <= fLB_LIN     .* (1 - auxfLB)    ): 'upper bounds of the primal constraint 0 <= F_MAX + f'];
constCP    = [constCP, ((F_MAX - f)           <= fUB_LIN     .* (1 - auxfUB)    ): 'upper bounds of the primal constraint 0 <= F MAX - f'];
constCP    = [constCP, (pd                    <= pdLB_LIN    .* (1 - auxpdLB)   ): 'upper bounds of the primal constraint 0 <= pd'];
constCP    = [constCP, ((PD_MAX - pd)         <= pdUB_LIN    .* (1 - auxpdUB)   ): 'upper bounds of the primal constraint 0 <= PD_MAX - pd'];
constCP    = [constCP, (pg                    <= pgLB_LIN    .* (1 - auxpgLB)   ): 'upper bounds of the primal constraint 0 <= pg'];
constCP    = [constCP, ((PG_MAX - pg)         <= pgUB_LIN    .* (1 - auxpgUB)   ): 'upper bounds of the primal constraint 0 <= PG_MAX - pg'];
constCP    = [constCP, (pw                    <= pwLB_LIN    .* (1 - auxpwLB)   ): 'upper bounds of the primal constraint 0 <= pw'];
constCP    = [constCP, ((PW_EXP - pw)         <= pwUB_LIN    .* (1 - auxpwUB)   ): 'upper bounds of the primal constraint 0 <= pwof - pw'];
constCP    = [constCP, ((pi + delta(nsla, :)) <= deltaLB_LIN .* (1 - auxdeltaLB)): 'upper bounds of the primal constraint 0 <= pi + delta'];
constCP    = [constCP, ((pi - delta(nsla, :)) <= deltaUB_LIN .* (1 - auxdeltaUB)): 'upper bounds of the primal constraint 0 <= pi - delta'];
constCP    = [constCP, (phiLB                 <= phiLB_LIN   .* auxfLB          ): 'upper bounds of the dual variable 0 <= phiLB'];
constCP    = [constCP, (phiUB                 <= phiUB_LIN   .* auxfUB          ): 'upper bounds of the dual variable 0 <= phiUB'];
constCP    = [constCP, (etaLB                 <= etaLB_LIN   .* auxpdLB         ): 'upper bounds of the dual variable 0 <= etaLB'];
constCP    = [constCP, (etaUB                 <= etaUB_LIN   .* auxpdUB         ): 'upper bounds of the dual variable 0 <= etaUB'];
constCP    = [constCP, (psiLB                 <= psiLB_LIN   .* auxpgLB         ): 'upper bounds of the dual variable 0 <= psiLB'];
constCP    = [constCP, (psiUB                 <= psiUB_LIN   .* auxpgUB         ): 'upper bounds of the dual variable 0 <= psiUB'];
constCP    = [constCP, (sigmaLB               <= sigmaLB_LIN .* auxpwLB         ): 'upper bounds of the dual variable 0 <= sigmaLB'];
constCP    = [constCP, (sigmaUB               <= sigmaUB_LIN .* auxpwUB         ): 'upper bounds of the dual variable 0 <= sigmaUB'];
constCP    = [constCP, (xiLB                  <= xiLB_LIN    .* auxdeltaLB      ): 'upper bounds of the dual variable 0 <= xiLB'];
constCP    = [constCP, (xiUB                  <= xiUB_LIN    .* auxdeltaUB      ): 'upper bounds of the dual variable 0 <= xiUB'];

constKKT   = [constP, constD, constCP];
optsKKT    = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagKKT    = optimize(constKKT, [], optsKKT);

% check the correctness of big-M of the primal inequalities
valfLB     = all(F_MAX + value(f)           < fLB_LIN,     'all'); % valid upper bounds of the primal constraint 0 <= F_MAX + f
valfUB     = all(F_MAX - value(f)           < fUB_LIN,     'all'); % valid upper bounds of the primal constraint 0 <= F_MAX - f
valpdLB    = all(value(pd)                  < pdLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pd
valpdUB    = all(PD_MAX - value(pd)         < pdUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= PD_MAX - pd
valpgLB    = all(value(pg)                  < pgLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pg
valpgUB    = all(PG_MAX - value(pg)         < pgUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= PG_MAX - pg
valpwLB    = all(value(pw)                  < pwLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pw
valpwUB    = all(PW_EXP - value(pw)         < pwUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pwof - pw
valdeltaLB = all(pi + value(delta(nsla, :)) < deltaLB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pi + delta(nsla, :)
valdeltaUB = all(pi - value(delta(nsla, :)) < deltaUB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pi - delta(nsla, :)

% check the correctness of big-M of the dual variables
valphiLB   = all(value(phiLB)   < phiLB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= phiLB
valphiUB   = all(value(phiUB)   < phiUB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= phiUB
valetaLB   = all(value(etaLB)   < etaLB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= etaLB
valetaUB   = all(value(etaUB)   < etaUB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= etaUB
valpsiLB   = all(value(psiLB)   < psiLB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= psiLB
valpsiUB   = all(value(psiUB)   < psiUB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= psiUB
valsigmaLB = all(value(sigmaLB) < sigmaLB_LIN, 'all'); % valid upper bounds of the dual variable 0 < sigmaLB
valsigmaUB = all(value(sigmaUB) < sigmaUB_LIN, 'all'); % valid upper bounds of the dual variable 0 < sigmaUB
valxiLB    = all(value(xiLB)    < xiLB_LIN,    'all'); % valid upper bounds of the dual variable 0 < xiLB
valxiUB    = all(value(xiUB)    < xiUB_LIN,    'all'); % valid upper bounds of the dual variable 0 < xiUB

if ~all([valfLB valfUB valpdLB valpdUB valpgLB valpgUB valpwLB valpwUB valdeltaLB valdeltaUB]) || ...
        ~all([valphiLB valphiUB valetaLB valetaUB valpsiLB valpsiUB valsigmaLB valsigmaUB valxiLB valxiUB])
    disp('INCORRECT BIG-M');
    return;
end

% wind power cleared in DA market (MW) (nt x nw)
resultKKT.prod = value(pw);
% locational marginal prices of bus n at which wind unit l located at ($/MWh)
LMPw = Cw' * value(lambda); % (nw x nt)
LMPw = transpose(LMPw);     % (nt x nw)
% difference of wind power level in the DA market and real wind power output of the lth wind in scenario s in period t (MW) (nt x nw x ns)
diffbal = repmat(value(pw), [1 1 ns]) - PW_SCEN;
% wind power proft when submited expected wind power output to the DA market and anticipating balancing market ($)
resultKKT.profit = sum(LMPw .* value(pw), 'all') - prob * sum(lambdaB .* diffbal, 'all');

fprintf('total production of wind units behave unstrategically solved using KKT conditions is %d MW\n', sum(resultKKT.prod));
fprintf('total profit     of wind units behave unstrategically solved using KKT conditions is %d $\n\n', resultKKT.profit);

% 2): green certificate market
% complementarity condition
constCPGC  = [];
constCPGC  = [constCPGC, (tw                   <= twLB_LIN  .* (1 - auxtwLB)): 'upper bounds of the primal constraint 0 <= tw'];
constCPGC  = [constCPGC, (PW_EXP - tw          <= twUB_LIN  .* (1 - auxtwUB)): 'upper bounds of the primal constraint 0 <= pwp - tw'];
constCPGC  = [constCPGC, (td                   <= tdLB_LIN  .* (1 - auxtdLB)): 'upper bounds of the primal constraint 0 <= td'];
constCPGC  = [constCPGC, (ratio .* PD_MAX - td <= tdUB_LIN  .* (1 - auxtdUB)): 'upper bounds of the primal constraint 0 <= ratio * PD_MAX - td'];
constCPGC  = [constCPGC, (muLB                 <= muLB_LIN  .* auxtwLB      ): 'upper bounds of the dual variable 0 <= muLB'];
constCPGC  = [constCPGC, (muUB                 <= muUB_LIN  .* auxtwUB      ): 'upper bounds of the dual variable 0 <= muUB'];
constCPGC  = [constCPGC, (rhoLB                <= rhoLB_LIN .* auxtdLB      ): 'upper bounds of the dual variable 0 <= rhoLB'];
constCPGC  = [constCPGC, (rhoUB                <= rhoUB_LIN .* auxtdUB      ): 'upper bounds of the dual variable 0 <= rhoUB'];

constKKTGC = [constPGC, constDGC, constCPGC];
optsKKTGC  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagKKTGC  = optimize(constKKTGC, [], optsKKTGC);

% check the correctness of big-M of the primal inequalities
valtwLB    = all(value(tw)                   < twLB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= tw
valtwUB    = all(PW_EXP - value(tw)          < twUB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pwp - tw
valtdLB    = all(value(td)                   < tdLB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= td
valtdUB    = all(ratio .* PD_MAX - value(td) < tdUB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= ratio * PD_MAX - td

% check the correctness of big-M of the dual variables
valmuLB    = all(value(muLB)  < muLB_LIN,  'all'); % valid upper bounds of the dual variable 0 <= muLB
valmuUB    = all(value(muUB)  < muUB_LIN,  'all'); % valid upper bounds of the dual variable 0 <= muUB
valrhoLB   = all(value(rhoLB) < rhoLB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= rhoLB
valrhoUB   = all(value(rhoUB) < rhoUB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= rhoUB

if ~all([valtwLB valtwUB valtdLB valtdUB]) || ...
        ~all([valmuLB valmuUB valrhoLB valrhoUB])
    disp('INCORRECT BIG-M');
end

% green certificate traded between wind unit and demand (MW or p)
resultKKTGC.prod = value(tw);
% wind power profit of trading green certificate ($)
resultKKTGC.profit = sum(repmat(value(lambdaGC), [1 nw]) .* value(tw), 'all') - sum(lambdaWGC .* value(tw), 'all');

fprintf('total traded green certificate of wind units behave unstrategically solved using KKT conditions is %d MW or p\n', sum(resultKKTGC.prod));
fprintf('total profit of green certificate of wind units behave unstragically solved using KKT conditions is %d $\n\n', resultKKTGC.profit);

fprintf('PRIMAL & DUAL & SD & KKT INDEED SAME!\n\n');

%% Solve the Bi-level Problem
% upper-level constraints
constUL     = [];
constUL     = [constUL, (pwp + pwbal == repmat(pw, [1 1 ns])): 'power balance of the lth wind unit in scenario s in period t'];
constUL     = [constUL, (0 <= pwp <= PW_SCEN): 'bounds on wind power production pwp'];
constUL     = [constUL, (0 <= pwof <= PW_MAX): 'bounds on wind power offered pwof'];

% 1): energy market
% change lambdaW to accommodate scenarios ($/MWh) (nt x nw) -> (nt x nw x ns)
lambdaW     = repmat(lambdaW, [1 1 ns]);

% maximize wind unit expected profit in energy market ($)
objBEN      = sum(lambdaD .* pd, 'all') - sum(lambdaG .* pg, 'all') ...
            - sum((phiLB + phiUB) .* F_MAX, 'all') - sum(etaUB .* PD_MAX, 'all') - sum(psiUB .* PG_MAX, 'all') - pi * sum(xiLB + xiUB, 'all') ...
            - prob * sum(lambdaW .* pwp, 'all') - prob * sum(lambdaB .* pwbal, 'all');

diffpwof    = setdiff([1:8], 6);
constP      = constP(diffpwof);
constP      = [constP, (0 <= pw <= pwof): 'power bounds on pw'];

diffalphaw  = setdiff([1:16], 2);
constD      = constD(diffalphaw);
constD      = [constD, (alphaw - re_lambdaw + sigmaUB - sigmaLB == 0): 'stationary condition w.r.t. pw'];

diffpwof    = setdiff([1:20], 8);
constCP     = constCP(diffpwof);
constCP     = [constCP, ((pwof - pw) <= pwUB_LIN .* (1 - auxpwUB)): 'upper bounds of the primal constraint 0 <= pwof - pw'];

% 2): green certificate market
% redefine lower-level variables of green certificate market
tw          = repmat(tw, [1 1 ns]);   % (nt x nw x ns)
td          = repmat(td, [1 1 1 ns]); % (nt x nd x nj x ns)

lambdaGC    = repmat(lambdaGC, [1 ns]);  % (nt x ns)
muLB        = repmat(muLB, [1 1 ns]);    % (nt x nw x ns)
muUB        = repmat(muUB, [1 1 ns]);    % (nt x nw x ns)
rhoLB       = repmat(rhoLB, [1 1 1 ns]); % (nt x nd x nj x ns)
rhoUB       = repmat(rhoUB, [1 1 1 ns]); % (nt x nd x nj x ns)

% redefine auxiliary binary variables of green certificate market
auxtwLB     = repmat(auxtwLB, [1 1 ns]);   % (nt x nw x ns)
auxtwUB     = repmat(auxtwUB, [1 1 ns]);   % (nt x nw x ns)
auxtdLB     = repmat(auxtdLB, [1 1 1 ns]); % (nt x nd x nj x ns)
auxtdUB     = repmat(auxtdUB, [1 1 1 ns]); % (nt x nd x nj x ns)

% redefine constants related to primal constraints of green certificate market
twLB_LIN    = repmat(twLB_LIN, [1 1 ns]);   % (nt x nw x ns)
twUB_LIN    = repmat(twUB_LIN, [1 1 ns]);   % (nt x nw x ns)
tdLB_LIN    = repmat(tdLB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)
tdUB_LIN    = repmat(tdUB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)

% redefine constants related to dual variables of green certificate market
muLB_LIN    = repmat(muLB_LIN, [1 1 ns]);    % (nt x nw x ns)
muUB_LIN    = repmat(muUB_LIN, [1 1 ns]);    % (nt x nw x ns)
rhoLB_LIN   = repmat(rhoLB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)
rhoUB_LIN   = repmat(rhoUB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)

% maximize wind unit expected profit in green certificate market ($)
objBGC      = prob * (-sum(rhoUB .* repmat(ratio .* PD_MAX,[1 1 1 ns]), 'all') + sum(repmat(lambdaDGC, [1 1 1 ns]) .* td, 'all') ...
            - sum(repmat(lambdaWGC, [1 1 ns]) .* tw, 'all'));

constPGC    = [];
constPGC    = [constPGC, (sum(tw, 2) == sum(td, [2 3])): 'green certificate balance in period t in scenario s'];
constPGC    = [constPGC, (0 <= tw): 'green certificate lower bounds on tw'];
constPGC    = [constPGC, (tw <= pwp): 'green certificate upper bounds on tw'];
constPGC    = [constPGC, (0 <= td <= repmat(ratio .* PD_MAX, [1 1 1 ns])): 'green certificate bounds on td'];

constDGC    = [];

re_lambdaGC = repmat(lambdaGC, [1 1 nw]);       % (nt x ns x nw)
re_lambdaGC = permute(re_lambdaGC, [1 3 2]);    % (nt x nw x ns)
constDGC    = [constDGC, (alphawGC - re_lambdaGC + muUB - muLB == 0): 'stationary condition w.r.t. tw'];

re_lambdaGC = repmat(lambdaGC, [1 1 nd nj]);    % (nt x ns x nd x nj)
re_lambdaGC = permute(re_lambdaGC, [1 3 4 2]);  % (nt x nd x nj x ns)
constDGC    = [constDGC, (-repmat(lambdaDGC, [1 1 1 ns]) + re_lambdaGC + rhoUB - rhoLB == 0): 'stationary condition w.r.t. td'];

constDGC    = [
               constDGC, ...
               (0 <= muLB):  'nonnegative of muLB',  (0 <= muUB):  'nonnegative of muUB', ...
               (0 <= rhoLB): 'nonnegative of rhoLB', (0 <= rhoUB): 'nonnegative of rhoUB'
];

constCPGC   = [];
constCPGC   = [constCPGC, (tw                   <= twLB_LIN  .* (1 - auxtwLB)): 'upper bounds of the primal constraint 0 <= tw'];
constCPGC   = [constCPGC, (pwp - tw             <= twUB_LIN  .* (1 - auxtwUB)): 'upper bounds of the primal constraint 0 <= pwp - tw'];
constCPGC   = [constCPGC, (td                   <= tdLB_LIN  .* (1 - auxtdLB)): 'upper bounds of the primal constraint 0 <= td'];
constCPGC   = [constCPGC, (repmat(ratio .* PD_MAX, [1 1 1 ns]) ...
                                           - td <= tdUB_LIN  .* (1 - auxtdUB)): 'upper bounds of the primal constraint 0 <= ratio * PD_MAX - td'];
constCPGC   = [constCPGC, (muLB                 <= muLB_LIN  .* auxtwLB      ): 'upper bounds of the dual variable 0 <= muLB'];
constCPGC   = [constCPGC, (muUB                 <= muUB_LIN  .* auxtwUB      ): 'upper bounds of the dual variable 0 <= muUB'];
constCPGC   = [constCPGC, (rhoLB                <= rhoLB_LIN .* auxtdLB      ): 'upper bounds of the dual variable 0 <= rhoLB'];
constCPGC   = [constCPGC, (rhoUB                <= rhoUB_LIN .* auxtdUB      ): 'upper bounds of the dual variable 0 <= rhoUB'];

constB      = [constUL, constP, constD, constCP, constPGC, constDGC, constCPGC];
objB        = objBEN + objBGC;
optsB       = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagB       = optimize(constB, -objB, optsB);

% check the correctness of big-M of the primal inequalities
valfLB      = all(F_MAX + value(f)           < fLB_LIN,     'all'); % valid upper bounds of the primal constraint 0 <= F_MAX + f
valfUB      = all(F_MAX - value(f)           < fUB_LIN,     'all'); % valid upper bounds of the primal constraint 0 <= F_MAX - f
valpdLB     = all(value(pd)                  < pdLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pd
valpdUB     = all(PD_MAX - value(pd)         < pdUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= PD_MAX - pd
valpgLB     = all(value(pg)                  < pgLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pg
valpgUB     = all(PG_MAX - value(pg)         < pgUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= PG_MAX - pg
valpwLB     = all(value(pw)                  < pwLB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pw
valpwUB     = all(value(pwof) - value(pw)    < pwUB_LIN,    'all'); % valid upper bounds of the primal constraint 0 <= pwof - pw
valdeltaLB  = all(pi + value(delta(nsla, :)) < deltaLB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pi + delta(nsla, :)
valdeltaUB  = all(pi - value(delta(nsla, :)) < deltaUB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pi - delta(nsla, :)

valtwLB     = all(value(tw)                   < twLB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= tw
valtwUB     = all(value(pwp) - value(tw)      < twUB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= pwp - tw
valtdLB     = all(value(td)                   < tdLB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= td
valtdUB     = all(repmat(ratio .* PD_MAX, [1 1 1 ns]) ...
                                  - value(td) < tdUB_LIN, 'all'); % valid upper bounds of the primal constraint 0 <= ratio * PD_MAX - td

% check the correctness of big-M of the dual variables
valphiLB    = all(value(phiLB)   < phiLB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= phiLB
valphiUB    = all(value(phiUB)   < phiUB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= phiUB
valetaLB    = all(value(etaLB)   < etaLB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= etaLB
valetaUB    = all(value(etaUB)   < etaUB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= etaUB
valpsiLB    = all(value(psiLB)   < psiLB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= psiLB
valpsiUB    = all(value(psiUB)   < psiUB_LIN,   'all'); % valid upper bounds of the dual variable 0 <= psiUB
valsigmaLB  = all(value(sigmaLB) < sigmaLB_LIN, 'all'); % valid upper bounds of the dual variable 0 < sigmaLB
valsigmaUB  = all(value(sigmaUB) < sigmaUB_LIN, 'all'); % valid upper bounds of the dual variable 0 < sigmaUB
valxiLB     = all(value(xiLB)    < xiLB_LIN,    'all'); % valid upper bounds of the dual variable 0 < xiLB
valxiUB     = all(value(xiUB)    < xiUB_LIN,    'all'); % valid upper bounds of the dual variable 0 < xiUB

valmuLB     = all(value(muLB)  < muLB_LIN,  'all'); % valid upper bounds of the dual variable 0 <= muLB
valmuUB     = all(value(muUB)  < muUB_LIN,  'all'); % valid upper bounds of the dual variable 0 <= muUB
valrhoLB    = all(value(rhoLB) < rhoLB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= rhoLB
valrhoUB    = all(value(rhoUB) < rhoUB_LIN, 'all'); % valid upper bounds of the dual variable 0 <= rhoUB

if ~all([valfLB valfUB valpdLB valpdUB valpgLB valpgUB valpwLB valpwUB valdeltaLB valdeltaUB valtwLB valtwUB valtdLB valtdUB]) || ...
        ~all([valphiLB valphiUB valetaLB valetaUB valpsiLB valpsiUB valsigmaLB valsigmaUB valxiLB valxiUB valmuLB valmuUB valrhoLB valrhoUB])
    disp('INCORRECT BIG-M');
    return;
end

%% Post Precessing
% 1): energy market
% wind power cleared in DA market (MW) (nt x nw)
resultB.prodEN = value(pw);
% wind units profit when strategically behave in DA market and anticipating balancing market($)
resultB.profitEN = value(objBEN);

fprintf('total production of wind units behave strategically solved using mip is %d MW\n', sum(resultB.prodEN));
fprintf('total profit     of wind units behave strategically solved using mip is %d $\n\n', resultB.profitEN);

figure
plot(value(pw), '-o');
xlabel('time (h)'); ylabel('wind power cleared in DA (MW)');

figure
plot(value(alphaw), '-o');
xlabel('time (h)'); ylabel('wind power offer price in DA ($/MWh)');

figure
vau_pw    = value(pw);                   % wind power cleared in day-ahead market (MW) (nt x nw)
vau_pw    = vau_pw(:, 1);                % (nt x 1)
vau_pw    = repmat(vau_pw, [1 ns]);      % (nt x ns)
vau_pwp   = value(pwp);                  % wind power produced in balancing market (MW) (nt x nw x ns)
vau_pwp   = squeeze(vau_pwp(:, 1, :));   % (nt x ns)
vau_pwbal = value(pwbal);                % power brought/sold to the balancing market (MW) (nt x nw x ns)
vau_pwbal = squeeze(vau_pwbal(:, 1, :)); % (nt x ns)
plot([1 : ns], vau_pw(9, :), '-*', [1 : ns], vau_pwbal(9, :), '-d', [1 : ns], vau_pwp(9, :), '-o');
xlabel('scenario'); ylabel('wind power (MW)'); ylim([-100 100]);
legend('cleared', 'balacing', 'produced');

figure
vau_lambda  = value(lambda);                 % LMP (nb x nt)
vau_lambda  = Cw' * vau_lambda;              % offer price of wind unit (nw x nt)
vau_lambda	= vau_lambda(1, :)';             % (nt x 1)
vau_lambda  = repmat(vau_lambda, [1 ns]);    % (nt x ns)
vau_lambdaB = lambdaB;                       % balancing market price (nt x nw x ns)
vau_lambdaB = squeeze(vau_lambdaB(:, 1, :)); % (nt x ns)
plot([1 : ns], vau_lambda(9, :), '-*', [1 : ns], vau_lambdaB(9, :), '-d');
xlabel('scenario'); ylabel('price ($/MWh)');
legend('offer', 'balacing');

% 2): green certificate market
% green certificate traded between wind unit and demand (MW or p)
resultB.prodGC = value(tw);
% wind power profit of trading green certificate ($)
resultB.profitGC = value(objBGC);

fprintf('total traded green certificate of wind units behave strategically solved using mip is %d MW or p\n', sum(resultB.prodGC));
fprintf('total profit of green certificate of wind units behave stragically solved using mip is %d $\n\n', resultB.profitGC);

figure
vau_tw = value(tw);
vau_tw = squeeze(vau_tw(:, 1, 1));
plot(vau_tw, '-o');
xlabel('time (h)'); ylabel('green certificate cleared (MW or p)');

figure
plot(value(alphawGC(:, :, 1)), '-o');
xlabel('time (h)'); ylabel('wind power offer price in GC maarket ($/MWh or $/p)');

% to get a picture of hourly profit day-ahead, balacing, green certificate and total profit


disp('DONE!');