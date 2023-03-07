% Strategic Offering of a Wind Producer in Energy Market and Green Certificate Market

clc; clear; close all;

%% Define Constant
nt   = 24;              % number of time periods
ns   = 10;              % number of scenarios
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
wbus = wind(:, 1);
Cw   = sparse(wbus, 1 : nw, 1, nb, nw);

% connection matrix of bus-gen
gbus = gen(:, 1);
Cg   = sparse(gbus, 1 : ng, 1, nb, ng);

% connection matrix of bus-demand
dbus = load(:, 1);
Cd   = sparse(dbus, 1 : nd, 1, nb, nd);

% connection matrix of line-bus
f    = branch(:, 1);
t    = branch(:, 2);
i    = [(1 : nl)'; (1 : nl)'];
Cft  = sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb);

%% Load Data
% data for the wind units
lambdaW   = zeros(nt, nw);               % marginal cost of the lth wind unit in period t ($/MWh) (nt x nw)

PW_MAX    = 300 * ones(nt, nw);          % maximum wind power output of the lth wind unit in period t (MW) (nt x nw)

factW(:, 1, :) = xlsread('SourceData', 'Wind Power Capacity Factors', 'B2:K25');
                                         % wind power capacity factor of the lth wind unit ...
                                         % in scenario s in period t (p.u.) (nt x nw x ns)

avefactW  = sum(factW, 3) ./ ns;         % averaged wind power capacity factor over scenarios ...
                                         % of the lth wind unit in period t (p.u.) (nt x nw)

PW_EXP    = avefactW .* PW_MAX;          % expected wind power output of the lth wind unit in period t (MW) (nt x nw)

PW_SCEN   = factW .* repmat(PW_MAX, [1 1 ns]);
                                         % wind power production by lth wind unit ...
                                         % in scenario s in period t (MW) (nt x nw x ns)

lambdaWGC = wind(:, 4);                  % (nw x 1)
lambdaWGC = repmat(lambdaWGC, [1 nt]);   % (nw x nt)
lambdaWGC = transpose(lambdaWGC);        % marginal cost of green certificate of the lth wind unit in period t ($/MWh or $/p) (nt x nw)

% data for the gen units
lambdaG   = gen(:, [6 7 8 9]);           % (ng x no)
lambdaG   = repmat(lambdaG, [1 1 nt]);   % (ng x no x nt)
lambdaG   = permute(lambdaG, [3 1 2]);   % marginal cost of the bth block of the ith gen unit in period t ($/MWh) (nt x ng x no)

PG_MAX    = gen(:, [2 3 4 5]);           % (ng x no)
PG_MAX    = repmat(PG_MAX, [1 1 nt]);    % (ng x no x nt)
PG_MAX    = permute(PG_MAX, [3 1 2]);    % upper limit of the bth block of the ith gen unit in period t (MW) (nt x ng x no)

% data for the demands
lambdaD   = load(:, [5 6 7]);            % (nd x nj)
lambdaD   = repmat(lambdaD, [1 1 nt]);   % (nd x nj x nt)
lambdaD   = permute(lambdaD, [3 1 2]);   % marginal utility of the jth block of the dth demand in period t ($/MWh) (nt x nd x nj)

PD_MAX    = load(:, [2 3 4]);            % (nd x nj)
PD_MAX    = repmat(PD_MAX, [1 1 nt]);    % (nd x nj x nt)
PD_MAX    = permute(PD_MAX, [3 1 2]);    % (nt x nd x nj)
factD     = [
    0.829 0.777 0.745 0.736 0.728 0.732 0.736 0.829 0.892 0.908 0.932 0.935 ...
    0.917 0.905 0.870 0.869 0.864 0.870 0.881 0.974 1.000 0.981 0.926 0.915
];                                       % hourly demand factros (p.u.) (1 x nt)
factD     = repelem(factD, nd, nj);      % (nd x (nj * nt))
factD     = reshape(factD, [nd nj nt]);  % (nd x nj x nt)
factD     = permute(factD, [3 1 2]);     % (nt x nd x nj)
PD_MAX    = factD .* PD_MAX;             % upper limit of the jth block of the dth demand in period t (MW) (nt x nd x nj)

ratio     = load(:, 8);                  % (nd x 1)
ratio     = repmat(ratio, [1 nt nj]);    % (nd x nt x nj)
ratio     = permute(ratio, [2 1 3]);     % renewable energy consumption weight (%) (nt x nd x nj)

lambdaDGC = load(:, [9 10 11]);          % (nd x nj)
lambdaDGC = repmat(lambdaDGC, [1 1 nt]); % (nd x nj x nt)
lambdaDGC = permute(lambdaDGC, [3 1 2]); % marginal utility of green certificate of ...
                                         % the jth block of the dth demand in period t ($/MWh or $/p) (nt x nd x nj)
% data for the network
B         = branch(:, 3);                % susceptance of line (S) (nl x 1)

F_MAX     = branch(:, 4);                % (nl x 1)
F_MAX     = repmat(F_MAX, [1 nt]);       % transmission capacity of line k in period t (MW) (nl x nt)

sla       = 1;
numsla    = length(sla);

nsla      = setdiff([1 : nb]', sla);
numnsla   = length(nsla);

% data for the balancing market
lambdaB   = xlsread('SourceData', 'Balancing Market Prices', 'B2:K25');
                                         % balancing market price in scenario s in period t ($/MWh) (nt x ns)
lambdaB   = repmat(lambdaB, [1 1 nw]);   % (nt x ns x nw)
lambdaB   = permute(lambdaB, [1 3 2]);   % (nt x nw x ns)

%% Define Variable
% upper-level variables
alphaw     = sdpvar(nt, nw,      'full'); % offer price of wind power of the lth wind unit in period t ($/MWh)
alphawGC   = sdpvar(nt, nw,  ns, 'full'); % offer price of green certificate of the lth wind unit in period t in scenario s ($/MWh or $/p)
pwof       = sdpvar(nt, nw,      'full'); % offer quantity of the lth wind unit in period t (MW)
pwp        = sdpvar(nt, nw,  ns, 'full'); % wind power produced by the lth wind power unit in period t in scenario s (MW)
pwbal      = sdpvar(nt, nw,  ns, 'full'); % power brought from/sold to the balancing market by the lth wind power unit ...
                                          % in period t in scenario s (MW)
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
psiLB      = sdpvar(nt, ng,  no, 'full'); % dual of lower bound of the bth block of the ith gen unit in period t ($/MWh)
psiUB      = sdpvar(nt, ng,  no, 'full'); % dual of upper bound of the bth block of the ith gen unit in period t ($/MWh)
sigmaLB    = sdpvar(nt, nw,      'full'); % dual of lower bound of the lth wind unit in period t ($/MWh)
sigmaUB    = sdpvar(nt, nw,      'full'); % dual of upper bound of the lth wind unit in period t ($/MWh)
xiLB       = sdpvar(numnsla, nt, 'full'); % dual of lower bound of the voltage angle at nonslack bus n in period t ($/h)
xiUB       = sdpvar(numnsla, nt, 'full'); % dual of upper bound of the voltage angle at nonslack bus n in period t ($/h)
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

%% Solve the Single-level Energy Market
% minus social welfare ($)
objSEN   = sum(lambdaG .* pg, 'all') + sum(lambdaW .* pw, 'all') - sum(lambdaD .* pd, 'all');

constSEN = [];

re_pg    = permute(pg, [2 3 1]);          % (ng x no x nt)
re_pg    = reshape(re_pg, [ng, no * nt]); % (ng x (no * nt))
sum_pg   = Cg * re_pg;                    % (nb x (no * nt))
sum_pg   = reshape(sum_pg, [nb no nt]);   % (nb x no x nt)
sum_pg   = squeeze(sum(sum_pg, 2));       % sum of the bth block of the ith gen unit at bus n in period t (MW) (nb x nt)

sum_pw   = Cw * pw';                      % sum of the lth wind production at bus n in period t (MW) (nb x nt)

re_pd    = permute(pd, [2 3 1]);          % (nd x nj x nt)
re_pd    = reshape(re_pd, [nd, nj * nt]); % (nd x (nj * nt))
sum_pd   = Cd * re_pd;                    % (nb x (nj * nt))
sum_pd   = reshape(sum_pd, [nb nj nt]);   % (nb x nj x nt)
sum_pd   = squeeze(sum(sum_pd, 2));       % sum of the jth block of the dth demand at bus n in period t (MW) (nb x nt)

sum_f    = Cft' * f;                      % sum of power flow leaving bus n in period t (MW) (nb x nt)

constSEN = [constSEN, (sum_pg + sum_pw - sum_f == sum_pd): 'power balance in DA market'];
constSEN = [constSEN, (f == diag(B) * Cft * delta): 'def of power flow'];
constSEN = [constSEN, (-F_MAX <= f <= F_MAX): 'transmission capacity'];
constSEN = [constSEN, (0 <= pd <= PD_MAX): 'power bounds on pd'];
constSEN = [constSEN, (0 <= pg <= PG_MAX): 'power bounds on pg'];
constSEN = [constSEN, (0 <= pw <= PW_EXP): 'power bounds on pw'];
constSEN = [constSEN, (-pi <= delta(nsla, :) <= pi): 'bounds on delta'];
constSEN = [constSEN, (delta(sla, :) == 0): 'def of slack bus'];

optsSEN  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagSEN  = optimize(constSEN, objSEN, optsSEN);

% wind power cleared in DA market (MW)
resultSEN.prod = value(pw);
% difference of wind power level in the DA market and real wind power output of the lth wind in scenario s in period t (MW)
diffbal = repmat(value(pw), [1 1 ns]) - PW_SCEN;
% wind power profit when submited expected wind power output to the DA market and anticipating balancing market ($)
resultSEN.profit = sum(-dual(constSEN(1)) .* value(sum_pw), 'all') - prob * sum(lambdaB .* diffbal, 'all');

fprintf('total production of wind units behave unstrategically in energy market is %d MW\n', sum(resultSEN.prod));
fprintf('total profit of wind units behave unstrategically in energy market is %d $\n\n', resultSEN.profit);

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

% constant related to dual variables of energy market
phiDL       = dual(constSEN(3));
phiLB_LIN   = phiDL(1 : nl * nt);
phiLB_LIN   = reshape(phiLB_LIN, [nl nt]);
phiLB_LIN   = (phiLB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= phiLB   ($/MWh) (nl x nt)
phiUB_LIN   = phiDL(nl * nt + 1 : end);
phiUB_LIN   = reshape(phiUB_LIN, [nl nt]);
phiUB_LIN   = (phiUB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= phiUB   ($/MWh) (nl x nt)

etaDL       = dual(constSEN(4));
etaLB_LIN   = etaDL(1 : nt * nd * nj);
etaLB_LIN   = reshape(etaLB_LIN, [nt nd nj]);
etaLB_LIN   = (etaLB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= etaLB   ($/MWh) (nt x nd x nj)
etaUB_LIN   = etaDL(nt * nd * nj + 1 : end);
etaUB_LIN   = reshape(etaUB_LIN, [nt nd nj]);
etaUB_LIN   = (etaUB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= etaUB   ($/MWh) (nt x nd x nj)

psiDL       = dual(constSEN(5));
psiLB_LIN   = psiDL(1 : nt * ng * no);
psiLB_LIN   = reshape(psiLB_LIN, [nt ng no]);
psiLB_LIN   = (psiLB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= psiLB   ($/MWh) (nt x ng x no)
psiUB_LIN   = psiDL(nt * ng * no + 1 : end);
psiUB_LIN   = reshape(psiUB_LIN, [nt ng no]);
psiUB_LIN   = (psiUB_LIN + 1) * 100;        % upper bounds of the dual variable 0 <= psiUB   ($/MWh) (nt x ng x no)

sigmaDL = dual(constSEN(6));
sigmaLB_LIN = sigmaDL(1 : nt * nw);
sigmaLB_LIN = reshape(sigmaLB_LIN, [nt nw]);
sigmaLB_LIN = (sigmaLB_LIN + 1) * 100;      % upper bounds of the dual variable 0 <= sigmaLB ($/MWh) (nt x nw)
sigmaUB_LIN = sigmaDL(nt * nw + 1 : end);
sigmaUB_LIN = reshape(sigmaUB_LIN, [nt nw]);
sigmaUB_LIN = (sigmaUB_LIN + 1) * 100;      % upper bounds of the dual variable 0 <= sigmaUB ($/MWh) (nt x nw)

xiDL        = dual(constSEN(7));
xiLB_LIN    = xiDL(1 : numnsla * nt);
xiLB_LIN    = reshape(xiLB_LIN, [numnsla nt]);
xiLB_LIN    = (xiLB_LIN + 1) * 100;         % upper bounds of the dual variable 0 <= xiLB    ($/h)   (numnsla x nt)
xiUB_LIN    = xiDL(numnsla * nt + 1 : end);
xiUB_LIN    = reshape(xiUB_LIN, [numnsla nt]);
xiUB_LIN    = (xiUB_LIN + 1) * 100;         % upper bounds of the dual variable 0 <= xiUB    ($/h)   (numnsla x nt)

%% Solve the Single-level Green Certificate Market
% minus social welfare ($)
objSGC   = sum(lambdaWGC .* tw, 'all') - sum(lambdaDGC .* td, 'all');

constSGC = [];
constSGC = [constSGC, (sum(tw, 2) == sum(td, [2 3])): 'green certificate balance in GC market'];
constSGC = [constSGC, (0 <= tw <= PW_EXP): 'green certificate bounds on tw'];
constSGC = [constSGC, (0 <= td <= ratio .* PD_MAX): 'green certificate bounds on td'];

optsSGC  = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagSGC  = optimize(constSGC, objSGC, optsSGC);

% green certificate traded between wind unit and demand (MW or p)
resultSGC.prod = value(tw);
% wind power profit of trading green certificate ($)
resultSGC.profit = sum(repmat(-dual(constSGC(1)), [1 nw]) .* value(tw), 'all') - sum(lambdaWGC .* value(tw), 'all');

fprintf('total traded green certificate of wind units behave unstrategically is %d p\n', sum(resultSGC.prod));
fprintf('total profit of green certificate of wind units behave unstragically is %d $\n\n', resultSGC.profit);

% constant related to primal constraints of green certificate market
twLB_LIN  = PW_MAX + 1;    % upper bounds of the primal constraint 0 <= tw                  (MW or p) (nt x nw)
twUB_LIN  = PW_MAX + 1;    % upper bounds of the primal constraint 0 <= pwp - tw            (MW or p) (nt x nw)
tdLB_LIN  = PD_MAX + 1;    % upper bounds of the primal constraint 0 <= td                  (MW or p) (nt x nd x nj)
tdUB_LIN  = PD_MAX + 1;    % upper bounds of the primal constraint 0 <= ratio * PD_MAX - td (MW or p) (nt x nd x nj)

% constant related to dual variables of green certificate market
muDL      = dual(constSGC(2));
muLB_LIN  = muDL(1 : nt * nw);
muLB_LIN  = reshape(muLB_LIN, [nt nw]);
muLB_LIN  = (muLB_LIN + 1) * 100;           % upper bounds of the dual variable 0 <= muLB  ($/MWh or $/p) (nt x nw)
muUB_LIN  = muDL(nt * nw + 1 : end);
muUB_LIN  = reshape(muUB_LIN, [nt nw]);
muUB_LIN  = (muUB_LIN + 1) * 100;           % upper bounds of the dual variable 0 <= muUB  ($/MWh or $/p) (nt x nw)

rhoDL     = dual(constSGC(3));
rhoLB_LIN = rhoDL(1 : nt * nd * nj);
rhoLB_LIN = reshape(rhoLB_LIN, [nt nd nj]); % upper bounds of the dual variable 0 <= rhoLB ($/MWh or $/p) (nt x nd x nj)
rhoLB_LIN = (rhoLB_LIN + 1) * 100;
rhoUB_LIN = rhoDL(nt * nd * nj + 1 : end);
rhoUB_LIN = reshape(rhoUB_LIN, [nt nd nj]);
rhoUB_LIN = (rhoUB_LIN + 1) * 100;          % upper bounds of the dual variable 0 <= rhoUB ($/MWh or $/p) (nt x nd x nj)

% redefine lower-level variables of green certificate market
tw        = repmat(tw, [1 1 ns]);   % (nt x nw x ns)
td        = repmat(td, [1 1 1 ns]); % (nt x nd x nj x ns)

lambdaGC  = repmat(lambdaGC, [1 ns]);  % (nt x ns)
muLB      = repmat(muLB, [1 1 ns]);    % (nt x nw x ns)
muUB      = repmat(muUB, [1 1 ns]);    % (nt x nw x ns)
rhoLB     = repmat(rhoLB, [1 1 1 ns]); % (nt x nd x nj x ns)
rhoUB     = repmat(rhoUB, [1 1 1 ns]); % (nt x nd x nj x ns)

% redefine auxiliary binary variables of green certificate market
auxtwLB   = repmat(auxtwLB, [1 1 ns]);   % (nt x nw x ns)
auxtwUB   = repmat(auxtwUB, [1 1 ns]);   % (nt x nw x ns)
auxtdLB   = repmat(auxtdLB, [1 1 1 ns]); % (nt x nd x nj x ns)
auxtdUB   = repmat(auxtdUB, [1 1 1 ns]); % (nt x nd x nj x ns)

% redefine constants related to primal constraints of green certificate market
twLB_LIN  = repmat(twLB_LIN, [1 1 ns]);   % (nt x nw x ns)
twUB_LIN  = repmat(twUB_LIN, [1 1 ns]);   % (nt x nw x ns)
tdLB_LIN  = repmat(tdLB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)
tdUB_LIN  = repmat(tdUB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)

% redefine constants related to dual variables of green certificate market
muLB_LIN  = repmat(muLB_LIN, [1 1 ns]);    % (nt x nw x ns)
muUB_LIN  = repmat(muUB_LIN, [1 1 ns]);    % (nt x nw x ns)
rhoLB_LIN = repmat(rhoLB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)
rhoUB_LIN = repmat(rhoUB_LIN, [1 1 1 ns]); % (nt x nd x nj x ns)

%% Solve the Bi-level Problem
% change lambdaW to accommodate scenarios ($/MWh) (nt x nw) -> (nt x nw x ns)
lambdaW     = repmat(lambdaW, [1 1 ns]);

% maximize wind unit expected profit in energy market ($)
objBEN      = sum(lambdaD .* pd, 'all') - sum(lambdaG .* pg, 'all') ...
            - sum((phiLB + phiUB) .* F_MAX, 'all') - sum(etaUB .* PD_MAX, 'all') - sum(psiUB .* PG_MAX, 'all') - pi * sum(xiLB + xiUB, 'all') ...
            - prob * sum(lambdaW .* pwp, 'all') - prob * sum(lambdaB .* pwbal, 'all');

% maximize wind unit expected profit in green certificate market ($)
objBGC      = prob * (-sum(rhoUB .* repmat(ratio .* PD_MAX,[1 1 1 ns]), 'all') + sum(repmat(lambdaDGC, [1 1 1 ns]) .* td, 'all') ...
            - sum(repmat(lambdaWGC, [1 1 ns]) .* tw, 'all'));

% upper-level constraints
constUL     = [];
constUL     = [constUL, (pwp + pwbal == repmat(pw, [1 1 ns])): 'power balance in balancing market'];
constUL     = [constUL, (0 <= pwp <= PW_SCEN): 'bounds on wind power production pwp'];
constUL     = [constUL, (0 <= pwof <= PW_MAX): 'bounds on wind power offered pwof'];

% primal feasibility of energy market
constPEN    = [];
constPEN    = [constPEN, (sum_pg + sum_pw - sum_f == sum_pd): 'power balance in DA market'];
constPEN    = [constPEN, (f == diag(B) * Cft * delta): 'def of power flow'];
constPEN    = [constPEN, (-F_MAX <= f <= F_MAX): 'transmission capacity'];
constPEN    = [constPEN, (0 <= pd <= PD_MAX): 'power bounds on pd'];
constPEN    = [constPEN, (0 <= pg <= PG_MAX): 'power bounds on pg'];
constPEN    = [constPEN, (0 <= pw <= pwof): 'power bounds on pw'];
constPEN    = [constPEN, (-pi <= delta(nsla, :) <= pi): 'bounds on delta'];
constPEN    = [constPEN, (delta(sla, :) == 0): 'def of slack bus'];

% dual feasibility of energy market
constDEN    = [];

re_lambdag  = Cg' * lambda;                 % n(i): node at which gen i located at ($/MWh) (ng x nt)
re_lambdag  = repmat(re_lambdag, [1 1 no]); % (ng x nt x no)
re_lambdag  = permute(re_lambdag, [2 1 3]); % (nt x ng x no)
constDEN    = [constDEN, (lambdaG - re_lambdag + psiUB - psiLB == 0): 'stationary condtion w.r.t. pg'];

re_lambdaw  = Cw' * lambda;                 % n(l): node at which wind unit l located at ($/MWh) (nw x nt)
re_lambdaw  = transpose(re_lambdaw);        % (nt x nw)
constDEN    = [constDEN, (alphaw - re_lambdaw + sigmaUB - sigmaLB == 0): 'stationary condition w.r.t. pw'];

re_lambdad  = Cd' * lambda;                 % n(d): node at which demand d located at ($/MWh) (nd x nt)
re_lambdad  = repmat(re_lambdad, [1 1 nj]); % (nd x nt x nj)
re_lambdad  = permute(re_lambdad, [2 1 3]); % (nt x nd x nj)
constDEN    = [constDEN, (-lambdaD + re_lambdad + etaUB - etaLB == 0):'stationary condition w.r.t. pd'];

constDEN    = [constDEN, (Cft * lambda - phi + phiUB - phiLB == 0): 'stationary condition w.r.t. f'];

re_phi      = Cft' * diag(B) * phi;         % sum_k|s(k)=n{B_k * phi_kt} - sum_k|r(k)=n{b_k*phi_kt} ($/h) (nb x nt)
constDEN    = [constDEN, (re_phi(nsla, :) + xiUB - xiLB == 0): 'stationary condition w.r.t. delta'];
constDEN    = [constDEN, (re_phi(sla, :) - xi == 0): 'stationary condition w.r.t. delta'];

constDEN    = [
    constDEN, ...
    (0 <= phiLB):   'nonnegative of phiLB',   (0 <= phiUB):   'nonnegative of phiUB', ...
    (0 <= etaLB):   'nonnegative of etaLB',   (0 <= etaUB):   'nonnegative of etaUB', ...
    (0 <= psiLB):   'nonnegative of psiLB',   (0 <= psiUB):   'nonnegative of psiUB', ...
    (0 <= sigmaLB): 'nonnegative of sigmaLB', (0 <= sigmaUB): 'nonnegative of sigmaUB', ...
    (0 <= xiLB):    'nonnegative of xiLB',    (0 <= xiUB):    'nonnegative of xiUB'
];

% complementarity condition of energy market
constCPEN   = [];
constCPEN   = [constCPEN, ((F_MAX + f)           <= fLB_LIN     .* (1 - auxfLB)    ): 'upper bounds of the primal constraint 0 <= F_MAX + f'];
constCPEN   = [constCPEN, ((F_MAX - f)           <= fUB_LIN     .* (1 - auxfUB)    ): 'upper bounds of the primal constraint 0 <= F MAX - f'];
constCPEN   = [constCPEN, (pd                    <= pdLB_LIN    .* (1 - auxpdLB)   ): 'upper bounds of the primal constraint 0 <= pd'];
constCPEN   = [constCPEN, ((PD_MAX - pd)         <= pdUB_LIN    .* (1 - auxpdUB)   ): 'upper bounds of the primal constraint 0 <= PD_MAX - pd'];
constCPEN   = [constCPEN, (pg                    <= pgLB_LIN    .* (1 - auxpgLB)   ): 'upper bounds of the primal constraint 0 <= pg'];
constCPEN   = [constCPEN, ((PG_MAX - pg)         <= pgUB_LIN    .* (1 - auxpgUB)   ): 'upper bounds of the primal constraint 0 <= PG_MAX - pg'];
constCPEN   = [constCPEN, (pw                    <= pwLB_LIN    .* (1 - auxpwLB)   ): 'upper bounds of the primal constraint 0 <= pw'];
constCPEN   = [constCPEN, ((pwof - pw)           <= pwUB_LIN    .* (1 - auxpwUB)   ): 'upper bounds of the primal constraint 0 <= pwof - pw'];
constCPEN   = [constCPEN, ((pi + delta(nsla, :)) <= deltaLB_LIN .* (1 - auxdeltaLB)): 'upper bounds of the primal constraint 0 <= pi + delta'];
constCPEN   = [constCPEN, ((pi - delta(nsla, :)) <= deltaUB_LIN .* (1 - auxdeltaUB)): 'upper bounds of the primal constraint 0 <= pi - delta'];
constCPEN   = [constCPEN, (phiLB                 <= phiLB_LIN   .* auxfLB          ): 'upper bounds of the dual variable 0 <= phiLB'];
constCPEN   = [constCPEN, (phiUB                 <= phiUB_LIN   .* auxfUB          ): 'upper bounds of the dual variable 0 <= phiUB'];
constCPEN   = [constCPEN, (etaLB                 <= etaLB_LIN   .* auxpdLB         ): 'upper bounds of the dual variable 0 <= etaLB'];
constCPEN   = [constCPEN, (etaUB                 <= etaUB_LIN   .* auxpdUB         ): 'upper bounds of the dual variable 0 <= etaUB'];
constCPEN   = [constCPEN, (psiLB                 <= psiLB_LIN   .* auxpgLB         ): 'upper bounds of the dual variable 0 <= psiLB'];
constCPEN   = [constCPEN, (psiUB                 <= psiUB_LIN   .* auxpgUB         ): 'upper bounds of the dual variable 0 <= psiUB'];
constCPEN   = [constCPEN, (sigmaLB               <= sigmaLB_LIN .* auxpwLB         ): 'upper bounds of the dual variable 0 <= sigmaLB'];
constCPEN   = [constCPEN, (sigmaUB               <= sigmaUB_LIN .* auxpwUB         ): 'upper bounds of the dual variable 0 <= sigmaUB'];
constCPEN   = [constCPEN, (xiLB                  <= xiLB_LIN    .* auxdeltaLB      ): 'upper bounds of the dual variable 0 <= xiLB'];
constCPEN   = [constCPEN, (xiUB                  <= xiUB_LIN    .* auxdeltaUB      ): 'upper bounds of the dual variable 0 <= xiUB'];

% KKT optimality condition of energy market
constKKTEN  = [constPEN, constDEN, constCPEN];

% primal feasibility of green certificate market
constPGC    = [];
constPGC    = [constPGC, (sum(tw, 2) == sum(td, [2 3])): 'green certificate balance in period t in scenario s'];
constPGC    = [constPGC, (0 <= tw): 'green certificate lower bounds on tw'];
constPGC    = [constPGC, (tw <= pwp): 'green certificate upper bounds on tw'];
constPGC    = [constPGC, (0 <= td <= repmat(ratio .* PD_MAX, [1 1 1 ns])): 'green certificate bounds on td'];

% dual feasibility of green certificate market
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

% complementarity condition of green certificate market
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

% KKT optimality condition of green certificate market
constKKTGC  = [constPGC, constDGC, constCPGC];

% solve the bi-level problem containing energy and green certificate market
constB      = [constUL, constKKTEN, constKKTGC];
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

%% Post Processing
% wind power cleared in DA market (MW) (nt x nw)
resultB.prodEN = value(pw);
% wind units profit when strategically behave in DA market and anticipating balancing market($)
resultB.profitEN = value(objBEN);

fprintf('total production of wind units behave strategically is %d MW\n', sum(resultB.prodEN));
fprintf('total profit     of wind units behave strategically is %d $\n\n', resultB.profitEN);

% green certificate traded between wind unit and demand (MW or p)
resultB.prodGC = value(tw);
% wind power profit of trading green certificate ($)
resultB.profitGC = value(objBGC);

fprintf('total traded green certificate of wind units behave strategically is %d p\n', sum(resultB.prodGC));
fprintf('total profit of green certificate of wind units behave stragically is %d $\n\n', resultB.profitGC);

figure(1);
plot(value(pw), '-o');
xlabel('time (h)'); ylabel('wind power cleared in DA market (MW)');
print('WP cleared', '-dsvg');

figure(2)
plot(value(alphaw), '-o');
xlabel('time (h)'); ylabel('wind unit offer price in DA market ($/MWh)');
print('offer price in DA', '-dsvg');

figure(3)
vau_tw = value(tw);
vau_tw = squeeze(vau_tw(:, 1, 1));
plot(vau_tw, '-o');
xlabel('time (h)'); ylabel('green certificate cleared in GC market (p)');
print('GC cleared', '-dsvg');

figure(4)
plot(value(alphawGC(:, :, 1)), '-o');
xlabel('time (h)'); ylabel('wind unit offer price in GC maarket ($/p)');
print('offer price in GC', '-dsvg');

figure(5)
profitDA  = sum(lambdaD .* pd, [2 3]) - sum(lambdaG .* pg, [2 3]) ...
          - sum((phiLB + phiUB) .* F_MAX, 1)' - sum(etaUB .* PD_MAX, [2 3]) - sum(psiUB .* PG_MAX, [2 3]) - pi * sum(xiLB + xiUB, 1)' ...
          - prob * sum(lambdaW .* pwp, [2 3]);
profitBL  = - prob * sum(lambdaB .* pwbal, [2 3]);
profitGC  = prob * (-sum(rhoUB .* repmat(ratio .* PD_MAX,[1 1 1 ns]), [2 3 4]) + sum(repmat(lambdaDGC, [1 1 1 ns]) .* td, [2 3 4]) ...
          - sum(repmat(lambdaWGC, [1 1 ns]) .* tw, [2 3]));
profitTol = profitDA + profitBL + profitGC;
plot([1 : nt], value(profitDA),  '-o', ...
     [1 : nt], value(profitBL),  '-*', ...
     [1 : nt], value(profitGC),  '-d', ...
     [1 : nt], value(profitTol), '-s');
xlabel('time (h)'); ylabel('profit ($)');
legend('day-ahead', 'balacing', 'green certificate', 'total', 'Location', 'north');
print('wind unit profit', '-dsvg');

disp('DONE!');