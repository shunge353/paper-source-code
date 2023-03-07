%% Solving linear bilevel problems
%    minimize   a'*x + b'*y
%    subject to ci'*x + di'*y <= ei
%               minimize   p'*x + q'*y
%               subject to rj'*x + sj'*y <= tj: lambdaj
% using big-M.
%
% single-level reformulation
%    minimize   a'*x + b'*y
%    subject to ci'*x + di'*y <= ei
%               rj'*x + sj'*y <= tj
%               q + sum_j{lambdaj*sj} = 0
%               lambdaj >= 0
%               lambdaj*(rj'*x + sj'*y - tj) = 0
%
% big-M reformulation
%    minimize   a'*x + b'*y
%    subject to ci'*x + di'*y <= ei
%               rj'*x + sj'*y <= tj
%               q + sum_j{lambdaj*sj} = 0
%               lambdaj >= 0
%               lambdaj <= uj*MDj
%               -rj'*x -sj'*y + tj <= (1 - uj)*MPj
%               uj is binary
%

%% Example
clear; close all; clc;

n = 2; % number of outer-level variables
m = 3; % number of lower-level variables
consto = 2; % number of outer-level constraints
consti = 6; % number of inner-level constraints

% Problem data
a = [-8; -4]; b = [4; -40; -4];
c1 = [-1; 0]; d1 = [0; 0; 0]; e1 = 0;
c2 = [0; -1]; d2 = [0; 0; 0]; e2 = 0;

p = [1; 2]; q = [1; 1; 2];
r1 = [0; 0]; s1 = [-1; 0; 0]; t1 = 0;
r2 = [0; 0]; s2 = [0; -1; 0]; t2 = 0;
r3 = [0; 0]; s3 = [0; 0; -1]; t3 = 0;
r4 = [0; 0]; s4 = [-1; 1; 1]; t4 = 10;
r5 = [2; 0]; s5 = [-1; 2; -0.5]; t5 = 10;
r6 = [0; 2]; s6 = [2; -1; -0.5]; t6 = 9.7;

MP = 10*ones(consti,1);
MD = 10*ones(consti,1);

% Solving lin bilevel problems using 'solvebilevel' command
x = sdpvar(n,1);
y = sdpvar(m,1);

OO = a'*x + b'*y;
CO = [c1'*x + d1'*y <= e1, ...
    c2'*x + d2'*y <= e2];

OI = p'*x + q'*y;
CI = [r1'*x + s1'*y <= t1, ...
    r2'*x + s2'*y <= t2, ...
    r3'*x + s3'*y <= t3, ...
    r4'*x + s4'*y <= t4, ...
    r5'*x + s5'*y <= t5, ...
    r6'*x + s6'*y <= t6];

opts = sdpsettings('verbose',0);
solvebilevel(CO,OO,CI,OI,y,opts);

% Solving lin bilevel problems using big-M reformulation
x = sdpvar(n,1);
y = sdpvar(m,1);
lambda = sdpvar(6,1);
u = binvar(6,1);

O = a'*x + b'*y;
C = [[c1'; c2']*x + [d1'; d2']*y <= [e1; e2], ...
    [r1'; r2'; r3'; r4'; r5'; r6']*x + [s1'; s2'; s3'; s4'; s5'; s6']*y <= [t1; t2; t3; t4; t5; t6], ...
    q + [s1, s2, s3, s4, s5, s6]*lambda == 0, ...
    lambda >= 0, ...
    lambda <= u.*MD, ...
    -[r1'; r2'; r3'; r4'; r5'; r6']*x - [s1'; s2'; s3'; s4'; s5'; s6']*y + [t1; t2; t3; t4; t5; t6] <= (1 - u).*MP];

opts = sdpsettings('solver','gurobi','verbose',0);
diag = optimize(C, O, opts);
comple = value(lambda).*([r1';r2';r3';r4';r5';r6']*value(x)+[s1';s2';s3';s4';s5';s6']*value(y)-[t1;t2;t3;t4;t5;t6]);

if (diag.problem == 0) && (isequal(comple,zeros(consti,1)))
    disp('Problems have sucessfully solved');
end

corrD = MD > value(lambda);
corrP = MP > -[r1';r2';r3';r4';r5';r6']*value(x)-[s1';s2';s3';s4';s5';s6']*value(y)+[t1;t2;t3;t4;t5;t6];

if all(corrP) && all(corrD)
    disp('We have chosen a correct big-M');
elseif ~all(corrP)
    disp('New constraints binding, turn to choose a larger MP');
elseif ~all(corrD)
    disp('New constraints binding, turn to choose a larger MD');
end

%% Counterexample
clear; close all; clc;

sdpvar x y

% Solving lin bilevel problems using 'solvebilevel' command
OO = -x - y;
CO = [0 <= x <= 2];

OI = y;
CI = [y >= 0, x - 0.01*y <= 1];

opts = sdpsettings('verbose',0);
solvebilevel(CO,OO,CI,OI,y,opts);

% Solving lin bilevel problems using big-M reformulation
MP1 = 200; MP2 = 200;
MD1 = 200; MD2 = 200;

sdpvar x y
sdpvar lambda1 lambda2
binvar u1 u2

O = -x - y;
C = [0 <= x <= 2, ...
    y >= 0, ...
    x - 0.01*y <= 1, ...
    1 - lambda1 - 0.01*lambda2 == 0, ...
    [lambda1, lambda2] >= 0, ...
    lambda1 <= u1*MD1, ...
    y <= (1 - u1)*MP1, ...
    lambda2 <= u2*MD2, ...
    -x + 0.01*y + 1 <= (1 - u2)*MP2];

opts = sdpsettings('solver', 'gurobi', 'verbose', 2);
diag = optimize(C, O, opts);
