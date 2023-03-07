% Strategic Investment in Generation and Storage Units
% 
% For more details please refer to Chapter 5 of the following book:
% Bilevel Optimization: Advances and Next Challenges
% URL: https://link.springer.com/chapter/10.1007/978-3-030-52119-6_5
% 
clc; clear; close all;

% Define Constant
ng = 20;              % number of generating units to be built
th = 10;              % number of thermal units within the planned generating units
so = 10;              % number of solar   units within the planned generating units
ns = 10;              % number of storage units to be built
nt = 24;              % number of hours of a representative day

% Load Data
ETA                   = 4 * ones(ns, 1);                             % energy capacity of storage unit s                 (h)
RHO                   = zeros(ng, nt);                               % capacity factor of generating unit g at time t    (p.u.)
RHO(1 : th, :)        = 1;                                           % capacity factor of thermal    unit g at time t    (p.u.)
soRHO                 = [0.00 0.00 0.00 0.00 0.03 0.35 0.51 0.59 ... % solar capacity factors for the representative day (p.u.)
                         0.58 0.51 0.23 0.54 0.28 0.34 0.45 0.69 ...
                         0.70 0.61 0.32 0.02 0.00 0.00 0.00 0.00];
RHO(th + (1 : so), :) = repmat(soRHO, so, 1);                        % capacity factor of solar unit g at time t  (p.u.)
CG                    = zeros(ng, 1);                                % linear cost parameter of generating unit g ($/MWh)
CG(1 : th, :)         = 60;                                          % linear cost parameter of thermal    unit g ($/MWh)
CG(th + (1 : so), :)  = 0;                                           % linear cost parameter of solar      unit g ($/MWh)
D                     = [0.65 0.60 0.50 0.28 0.31 0.46 0.65 0.74 ... % demand level at time period t              (p.u.)
                         0.79 0.86 0.88 0.82 0.69 0.59 0.56 0.66 ...
                         0.79 0.94 1.00 0.98 0.88 0.75 0.69 0.65]';
IG                    = zeros(ng, 1);                                % annualized investment cost of generating unit g ($/kW/year)
IG(1 : th, :)         = 42;                                          % annualized investment cost of thermal    unit g ($/kW/year)
IG(th + (1 : so), :)  = 85;                                          % annualized investment cost of solar      unit g ($/kW/year)
IS                    = 4 * ones(ns, 1);                             % annualized investment cost of storage    unit s ($/kW/year)
PG                    = zeros(ng, 1);                                % capacity of generating unit g                   (MW)
PG(1 : th, :)         = 100;                                         % capacity of thermal    unit g                   (MW)
PG(th + (1 : so), :)  = 100;                                         % capacity of solar      unit g                   (MW)
PS                    = 100 * ones(ns, 1);                           % capacity of storage    unit s                   (MW)
CS                    = 300;                                         % load shedding cost                              ($/MWh)

% Unit Conversion
CS   = CS * 1e-3;       % M$/GWh
DMAX = 1.0;             % GW
CG   = CG .* 1e-3;      % M$/GWh
PG   = PG .* 1e-3;      % GW
IG   = IG .* PG ./ 365; % M$ per unit of generating units per hour
PS   = PS .* 1e-3;      % GW
IS   = IS .* PS ./ 365; % M$ per unit of storage    units per hour
D    = D .* DMAX;       % GW

% Big-M Values
BETAUB_MIN  =      zeros(ng, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables betaUB  * u_g
BETAUB_MAX  = 1e2 * ones(ng, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables betaUB  * u_g
GAMMALB_MIN =      zeros(ns, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables gammaLB * v_s
GAMMALB_MAX = 1e2 * ones(ns, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables gammaLB * v_s
GAMMAUB_MIN =      zeros(ns, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables gammaUB * v_s
GAMMAUB_MAX = 1e2 * ones(ns, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables gammaUB * v_s
MUUB_MIN    =      zeros(ns, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables muUB    * v_s
MUUB_MAX    = 1e2 * ones(ns, nt); % auxiliary large constants used for linearization of the product of continuous and binary variables muUB    * v_s

% Define Variables
u_g         = binvar(ng, 1);  % binary variable equal to 1 if generating unit g is built in the current planning period, and 0 otherwise
v_s         = binvar(ns, 1);  % binvar variable equal to 1 if storage    unit s is built in the current planning period, and 0 otherwise
p_gt        = sdpvar(ng, nt); % output of generating    unit g at time t (GW)
p_st        = sdpvar(ns, nt); % output of storage       unit s at time t (GW)
e_st        = sdpvar(ns, nt); % energy level of storage unit s at time t (GWh)
d_t         = sdpvar(nt, 1);  % satisfied demand at time t               (GW)
lambda      = sdpvar(nt, 1);  % electricity price at time t              (M$/GWh)
alphaLB     = sdpvar(nt, 1);  % dual of lower bound on demand
alphaUB     = sdpvar(nt, 1);  % dual of upper bound on demand
betaLB      = sdpvar(ng, nt); % dual of lower bound on output of      generating unit g at time t
betaUB      = sdpvar(ng, nt); % dual of upper bound on output of      generating unit g at time t
gammaLB     = sdpvar(ns, nt); % dual of lower bound on output of      storage    unit s at time t
gammaUB     = sdpvar(ns, nt); % dual of upper bound on output of      storage    unit s at time t
kappa       = sdpvar(ns, nt); % dual of definition of storage balance of storage unit s at time t (M$/GWh)
muLB        = sdpvar(ns, nt); % dual of lower bound on energy level of storage   unit s at time t
muUB        = sdpvar(ns, nt); % dual of upper bound on energy level of storage   unit s at time t
betaUB_aux  = sdpvar(ng, nt); % dual of upper bound on output of      generating unit g at time t
gammaLB_aux = sdpvar(ns, nt); % dual of lower bound on output of      storage    unit s at time t
gammaUB_aux = sdpvar(ns, nt); % dual of upper bound on output of      storage    unit s at time t
muUB_aux    = sdpvar(ns, nt); % dual of upper bound on energy level of storage   unit s at time t

% Constraints and Model definition of Centralized approach
objC   = sum(sum(repmat(CG, 1, nt) .* p_gt)) + sum(IG .* u_g) + sum(IS .* v_s) + CS * sum(D - d_t); % minimize total operating and investment costs (M$ per day)
constC = [];
constC = [constC, ( sum(p_gt, 1)' + sum(p_st, 1)' == d_t ):                          'Demand balance (1f)'];
constC = [constC, ( 0 <= d_t <= D ):                                                 'Lower and upper bound on satisfied demand (1g)'];
constC = [constC, ( 0 <= p_gt <= repmat(u_g .* PG, 1, nt) .* RHO ):                  'Lower and upper bound on generating unit production (1h)'];
constC = [constC, ( -repmat(v_s .* PS, 1, nt) <= p_st <= repmat(v_s .* PS, 1, nt) ): 'Lower and upper bound on storage unit generation (1i)'];
constC = [constC, ( e_st(:, 1) == 0 - p_st(:, 1) ):                                  'Storage balance constraint at time t = 1 (1j)'];
constC = [constC, ( e_st(:, 2 : end) == e_st(:, 1 : end - 1) - p_st(:, 2 : end) ):   'Storage balance constraint at time t = 2 ~ 24 (1j)'];
constC = [constC, ( 0 <= e_st <= repmat(ETA .* v_s .* PS, 1, nt) ):                  'Lower and upper bound on storage unit level (1k)'];
opts   = sdpsettings('solver', 'Gurobi', 'verbose', 0);
diag   = optimize(constC, objC, opts);

u_g_vau  = value(u_g);
v_s_vau  = value(v_s);
p_gt_vau = value(p_gt);
p_st_vau = value(p_st);
d_t_vau  = value(d_t);

resC.TolCost  = value(objC) * 365;                                                         % total costs per year (M$)
resC.InveCost = ( sum(IG .* u_g_vau) + sum(IS .* v_s_vau) ) * 365;                         % investment costs per unit of generating units per year (M$)
resC.OperCost = ( sum(sum(repmat(CG, 1, nt) .* p_gt_vau)) + CS * sum(D - d_t_vau) ) * 365; % operating costs per year (M$)
resC.LoadShed = sum(D - d_t_vau) / sum(D) * 1e2;                                           % load shedding (%)
resC.ThInv    = sum(u_g_vau(1 : th) .* PG(1 : th)) * 1e3;                                  % thermal capacity (MW)
resC.RenInv   = sum(u_g_vau(th + (1 : so)) .* PG(th + (1 : so))) * 1e3;                    % solar investment (MW)
resC.SoInv    = sum(v_s_vau .* PS) * 1e3;                                                  % storage investment (MW)

% Fix the Investment decisions to Determine Marginal prices
objC   = sum(sum(repmat(CG, 1, nt) .* p_gt)) + sum(IG .* u_g_vau) + sum(IS .* v_s_vau) + CS * sum(D - d_t);
constC = [];
constC = [constC, ( sum(p_gt, 1)' + sum(p_st, 1)' == d_t ):                                  'Demand balance (1f)'];
constC = [constC, ( 0 <= d_t <= D ):                                                         'Lower and upper bound on satisfied demand (1g)'];
constC = [constC, ( 0 <= p_gt <= repmat(u_g_vau .* PG, 1, nt) .* RHO ):                      'Lower and upper bound on generating unit production (1h)'];
constC = [constC, ( -repmat(v_s_vau .* PS, 1, nt) <= p_st <= repmat(v_s_vau .* PS, 1, nt) ): 'Lower and upper bound on storage unit generation (1i)'];
constC = [constC, ( e_st(:, 1) == 0 - p_st(:, 1) ):                                          'Storage balance constraint at time t = 1 (1j)'];
constC = [constC, ( e_st(:, 2 : end) == e_st(:, 1 : end - 1) - p_st(:, 2 : end) ):           'Storage balance constraint at time t = 2 ~ 24 (1j)'];
constC = [constC, ( 0 <= e_st <= repmat(ETA .* v_s_vau .* PS, 1, nt) ):                      'Lower and upper bound on storage unit level (1k)'];
opts   = sdpsettings('solver', 'Gurobi', 'verbose', 0);
diag   = optimize(constC, objC, opts);

resC.Profit = ( sum(sum((repmat(-dual(constC('Demand balance (1f)'))', ng, 1) - repmat(CG, 1, nt)) .* p_gt_vau)) + ...
                sum(sum(repmat(-dual(constC('Demand balance (1f)'))', ns, 1) .* p_st_vau)) - sum(IG .* u_g_vau) - sum(IS .* v_s_vau) ) * 365; % power producer profit per year (M$)
resC.AvePrc = ( sum(-dual(constC('Demand balance (1f)')) .* d_t_vau) / sum(d_t_vau) ) * 1e3;                                                  % average price ($/MWh)

% Constraints and Model definition of Bilevel approach
objB   = sum(sum((betaUB - betaUB_aux) .* RHO .* repmat(PG, 1, nt))) + sum(sum((gammaLB - gammaLB_aux) .* repmat(PS, 1, nt))) + ... % maximize the revenue of merchant investors (M$ per day)
         sum(sum((gammaUB - gammaUB_aux) .* repmat(PS, 1, nt))) + sum(sum((muUB - muUB_aux) .* repmat(PS, 1, nt))) - sum(IG .* u_g) - sum(IS .* v_s);
constB = [];
constB = [constB, ( sum(p_gt, 1)' + sum(p_st, 1)' == d_t ):                                                              'Deamnd balance (3a)'];
constB = [constB, ( 0 <= d_t <= D ):                                                                                     'Lower and upper bound on satisfied demand (3b)'];
constB = [constB, ( 0 <= p_gt <= repmat(u_g .* PG, 1, nt) .* RHO ):                                                      'Lower and upper bound on generating unit production (3c)'];
constB = [constB, ( -repmat(v_s .* PS, 1, nt) <= p_st <= repmat(v_s .* PS, 1, nt) ):                                     'Lower and upper bound on storage unit production (3d)'];
constB = [constB, ( e_st(:, 1) == 0 - p_st(:, 1) ):                                                                      'Storage balance constraint at time t = 1 (3e)'];
constB = [constB, ( e_st(:, 2 : end) == e_st(:, 1 : end - 1) - p_st(:, 2 : end) ):                                       'Storage balance constraint at time t = 2 ~ 24 (3e)'];
constB = [constB, ( 0 <= e_st <= repmat(ETA .* v_s .* PS, 1, nt) ):                                                      'Lower and upper bound on storage unit level (3f)'];
constB = [constB, ( -CS + lambda - alphaLB + alphaUB == 0 ):                                                             'Derivative of Lagrangian with respect to d_t (3g)'];
constB = [constB, ( repmat(CG, 1, nt) - repmat(lambda', ng, 1) -betaLB + betaUB == 0 ):                                  'Derivative of Lagrangian with respect to p_gt (3h)'];
constB = [constB, ( -repmat(lambda', ns, 1) - gammaLB + gammaUB + kappa == 0 ):                                          'Derivative of Lagraigian with respect to p_st (3i)'];
constB = [constB, ( kappa(:, 1 : end - 1) - kappa(:, 2 : end) - muLB(:, 1 : end - 1) + muUB(:, 1 : end - 1) == 0 ):      'Derivative of Lagrangian with respect to e_st (3j)'];
constB = [constB, ( kappa(:, end) - muLB(:, end) + muUB(:, end) == 0 ):                                                  'Derivative of Lagraigian with respect to e_st (3k)'];
constB = [constB, ( alphaLB >= 0 ):                                                                                      'Nonnegativity of dual varable (3l)'];
constB = [constB, ( alphaUB >= 0 ):                                                                                      'Nonnegativity of dual varable (3l)'];
constB = [constB, ( betaLB >= 0 ):                                                                                       'Nonnegativity of dual varable (3l)'];
constB = [constB, ( betaUB >= 0 ):                                                                                       'Nonnegativity of dual varable (3l)'];
constB = [constB, ( gammaLB >= 0 ):                                                                                      'Nonnegativity of dual varable (3l)'];
constB = [constB, ( gammaUB >= 0 ):                                                                                      'Nonnegativity of dual varable (3l)'];
constB = [constB, ( muLB >= 0 ):                                                                                         'Nonnegativity of dual varable (3l)'];
constB = [constB, ( muUB >= 0 ):                                                                                         'Nonnegativity of dual varable (3l)'];
constB = [constB, ( sum(sum(repmat(CG, 1, nt) .* p_gt)) + CS * sum(D - d_t) == -sum((alphaUB - CS) .* D) - sum(sum((betaUB - betaUB_aux) .* RHO .* repmat(PG, 1, nt))) - ...
    sum(sum((gammaLB - gammaLB_aux) .* repmat(PS, 1, nt))) - sum(sum((gammaUB - gammaUB_aux) .* repmat(PS, 1, nt))) - sum(sum((muUB - muUB_aux) .* repmat(ETA .* PS, 1, nt))) ): ...
                                                                                                                         'Linearized complementarity equality constraint (5a)'];
constB = [constB, ( repmat(u_g, 1, nt) .* BETAUB_MIN <= betaUB - betaUB_aux <= repmat(u_g, 1, nt) .* BETAUB_MAX ):       'Linearzied complementarity lower and upper bound (5b)'];
constB = [constB, ( (1 - repmat(u_g, 1, nt)) .* BETAUB_MIN <= betaUB_aux <= (1 - repmat(u_g, 1, nt)) .* BETAUB_MAX ):    'Linearzied complementarity lower and upper bound (5c)'];
constB = [constB, ( repmat(v_s, 1, nt) .* GAMMALB_MIN <= gammaLB - gammaLB_aux <= repmat(v_s, 1, nt) .* GAMMALB_MAX ):   'Linearized complementarity lower and upper bound (5d)'];
constB = [constB, ( (1 - repmat(v_s, 1, nt)) .* GAMMALB_MIN <= gammaLB_aux <= (1 - repmat(v_s, 1, nt)) .* GAMMALB_MAX ): 'Linearzied complementarity lower and upper bound (5e)'];
constB = [constB, ( repmat(v_s, 1, nt) .* GAMMAUB_MIN <= gammaUB - gammaUB_aux <= repmat(v_s, 1, nt) .* GAMMAUB_MAX ):   'Linearized complementarity lower and upper bound (5f)'];
constB = [constB, ( (1 - repmat(v_s, 1, nt)) .* GAMMAUB_MIN <= gammaUB_aux <= (1 - repmat(v_s, 1, nt)) .* GAMMAUB_MAX ): 'Linearized complementarity lower and upper bound (5g)'];
constB = [constB, ( repmat(v_s, 1, nt) .* MUUB_MIN <= muUB - muUB_aux <= repmat(v_s, 1, nt) .* MUUB_MAX ):               'Linearized complementarity lower and upper bound (5h)'];
constB = [constB, ( (1 - repmat(v_s, 1, nt)) .* MUUB_MIN <= muUB_aux <= (1 - repmat(v_s, 1, nt)) .* MUUB_MAX ):          'Linearized complementarity lower and upper bound (5i)'];
opts   = sdpsettings('solver', 'Gurobi', 'verbose', 0);
diag   = optimize(constB, -objB, opts);

% Check validity of Upper bound of Dual variables
betaUB_max_cort  = all(value(betaUB) < BETAUB_MAX, 'all');
gammaLB_max_cort = all(value(gammaLB) < GAMMALB_MAX, 'all');
gammaUB_max_cort = all(value(gammaUB) < GAMMAUB_MAX, 'all');
muUB_max_cort    = all(value(muUB) < MUUB_MAX, 'all');
if ~all([betaUB_max_cort, gammaLB_max_cort, gammaUB_max_cort, muUB_max_cort])
    % ref: https://doi.org/10.1287/opre.2019.1944
    disp(['Assumptions on the upper bound of dual variables are wrong, ' ...
    'but we can do little about this since chosing a valid bound is NP-hard']);
end

u_g_vau    = value(u_g);
v_s_vau    = value(v_s);
p_gt_vau   = value(p_gt);
p_st_vau   = value(p_st);
d_t_vau    = value(d_t);
lambda_vau = value(lambda);

resB.TolCost  = ( sum(sum(repmat(CG,1 ,nt) .* p_gt_vau)) + CS * sum(D - d_t_vau) + sum(IG .* u_g_vau) + sum(IS .* v_s_vau) ) * 365; % total costs per year (M$)
resB.InveCost = ( sum(IG .* u_g_vau) + sum(IS .* v_s_vau) ) * 365;                                                                  % investment costs per unit of generating units per year (M$)
resB.OperCost = ( sum(sum(repmat(CG, 1, nt) .* p_gt_vau)) + CS * sum(D - d_t_vau) ) * 365;                                          % operating costs per year (M$)
resB.LoadShed = sum(D - d_t_vau) / sum(D) * 1e2;                                                                                    % load shedding (%)
resB.ThInv    = sum(u_g_vau(1 : th) .* PG(1 : th)) * 1e3;                                                                           % thermal capacity (MW)
resB.RenInv   = sum(u_g_vau(th + (1 : so)) .* PG(th + (1 : so))) * 1e3;                                                             % solar investment (MW)
resB.SoInv    = sum(v_s_vau .* PS) * 1e3;                                                                                           % storage investment (MW)
resB.Profit   = value(objB) * 365;                                                                                                  % power producer profit per year (M$)
resB.AvePrc   = sum(lambda_vau .* d_t_vau) / sum(d_t_vau) * 1e3;                                                                    % average price ($/MWh)
