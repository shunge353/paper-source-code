% -------------------- nested formulation of 2S-RLO --------------------
%   min_{y} c'*y + max_{u in U} min_{x in F(y,u)} b'*x
%   s.t.    A*y >= d, y >= 0
%   where   F(y,u) = {x: G*x >= h - E*y - M*u, x >= 0}, U is a polyhedron
%           uncertainty set for the uncertain parameter u.
%
clear; close all; clc;

%% Data
m = 3;              % number of facilities
n = 3;              % number of customers

c = [400; 414; 326; 18; 25; 20];
b = [22; 33; 24; 33; 23; 30; 20; 25; 27];
A = [800 0   0  -1  0  0;
     0   800 0   0 -1  0;
     0   0   800 0  0 -1];
d = [0; 0; 0];
G = [-1 -1 -1  0  0  0  0  0  0;
      0  0  0 -1 -1 -1  0  0  0;
      0  0  0  0  0  0 -1 -1 -1;
      1  0  0  1  0  0  1  0  0;
      0  1  0  0  1  0  0  1  0;
      0  0  1  0  0  1  0  0  1];
h = [0; 0; 0; 206; 274; 220];
E = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0];
M = [0    0    0;
     0    0    0;
     0    0    0;
     -40  0    0;
     0  -40    0;
     0    0  -40];

tol = 1e-2;         % tolerance of optimality

%% Variables
y = binvar(m,1);
z = sdpvar(m,1);
y = [y;z];                      % first-stage decision variables
x = sdpvar(m*n,100);            % second-stage decision variables
eta = sdpvar(1,1);              % auxiliary variable for the value function
g = sdpvar(n,1);                % uncertain parameter
pi = sdpvar(size(G,1),1);       % dual variable for the second-stage constraints
v = binvar(size(G,1),1);        % auxiliary variable for the comple. slack. cond.
w = binvar(size(G,2),1);        % auxiliary variable for the comple. slack. cond.
u = binvar(size(G,1),1);        % auxiliary variable for the comple. slack. cond.
s = sdpvar(size(G,1),1);        % slack variable for feasibility test

%% CCG
% ------------------------------ master problem ------------------------------
%   min_    {y,eta,x^{l} for l in {1,...,k}} c'*y + eta             (1)
%   s.t.    A*y >= d, y >= 0                                        (2)
%           eta >= L                                                (3)
%           eta >= b'*x^{l}, l = 1,...k                             (4)
%           E*y + G*x^{l} >= h - M*u^{l}, l = 1,...k                (5)
%           x^{l} >= 0, l = 1,...k                                  (6)
%   where   x^{l}, l = 1,...,k are the generated second-stage variables until iter. k,
%           L is an initial lower bound on value function, u^{l}, l = 1,...,k are the
%           detected extreme point or extreme ray of the polyhedron uncertainty set U
%           until iter. k.
%           (2)-(3) are the first-stage constraints, (4)-(6) are the optimality cuts
%           assoc. with new created variable x^{l}, (5)-(6) are the feasibility cuts
%           assoc. with new created variable x^{l}.
% 
% ------------------------------ subproblem ------------------------------
%   max_    {x,s,pi,u} e'*s                                         (7)
%   s.t.    0 <= G*x + s - h + E*y^{k+1} + M*u comple, pi >= 0      (8)
%           0 <= x comple. -G'*pi >= 0                              (9)
%           0 <= s comple. e - pi >= 0                              (10)
%           u in U                                                  (11)
%   where   (7)-(11) are the feasibility subproblem to detect a worst case scenario
%           u in U such that the second stage feasible set is empty, i.e. F(y,u) null,
%           (8)-(10) are the KKT conditions of the modified second-stage problem
%           for a fixed first-stage variable y^{k+1} and a realized scenario u in U,
%           (11) constrains the uncertain parameter u in the predefined U.
% 
%   max_    {x,pi,u} b'*x                                           (12)
%   s.t.    0 <= G*x - h + E*y^{k+1} + M*u comple, pi >= 0          (13)
%           0 <= x comple. b - G'*pi >= 0                           (14)
%           u in U                                                  (15)
%   where   (12)-(15) are the optimality subproblem to detect a optimal scenario
%           u in U such that the second-stage optimal value is finite, i.e.
%           Q(y^{k+1}) < +inf, (13)-(14) are the KKT conditions of the original
%           second-stage problem.
%
LB = -inf; UB = +inf; iter = 1; bigM = 1e3;
disp('CCG');

constMP = [A*y >= d, y >= 0, eta >= -1e3];
objMP = c'*y + eta;
opts = sdpsettings('solver', 'gurobi', 'verbose', 0);

while abs(UB - LB) > tol
    disp(['iters:', num2str(iter)]);
    diag = optimize(constMP, objMP, opts);
    LB = value(objMP);
    
    objSP = sum(s);
    constSP = [];
    constSP = [constSP, G*x(:,iter) + s - h + E*value(y) + M*g >= 0, x(:,iter) >= 0, s >= 0];
    constSP = [constSP, -G'*pi >= 0, 1 >= pi >= 0];
    constSP = [constSP, G*x(:,iter) + s - h + E*value(y) + M*g <= bigM*v, pi <= bigM*(1 - v)];
    constSP = [constSP, x(:,iter) <= bigM*w, -G'*pi <= bigM*(1 - w)];
    constSP = [constSP, s <= bigM*u, 1 - pi <= bigM*(1 - u)];
    constSP = [constSP, 0 <= g <= 1, g(1) + g(2) + g(3) <= 1.8, g(1) + g(2) <= 1.2];
    diag = optimize(constSP, -objSP, opts);
    
    % check the bigM value
    corrtbigM = all(value([G*x(:,iter) + s - h + E*value(y) + M*g; pi; x(:,iter); -G'*pi; s; 1 - pi]) < bigM);
    if corrtbigM ~= 1
        disp('invalid bigM');
        break;
    end
    
    % If the second stage problem is infeasible under the worst case scenario
    if sum(value(s)) > tol
        constMP = [constMP, x(:,iter+1) >= 0, G*x(:,iter+1) - h + E*y + M*value(g) >= 0];
        iter = iter+1;
        continue;
    % else the second stage problem is feasible
    else
        objSP = b'*x(:,iter);
        constSP = [];
        constSP = [constSP, G*x(:,iter) - h + E*value(y) + M*g >= 0, x(:,iter) >= 0];
        constSP = [constSP, b - G'*pi >= 0, pi >= 0];
        constSP = [constSP, G*x(:,iter) - h + E*value(y) + M*g <= bigM*v, pi <= bigM*(1 - v)];
        constSP = [constSP, x(:,iter) <= bigM*w, b - G'*pi <= bigM*(1 - w)];
        constSP = [constSP, 0 <= g <= 1, g(1) + g(2) + g(3) <= 1.8, g(1) + g(2) <= 1.2];
        diag = optimize(constSP, -objSP, opts);
        
        % check the bigM value
        corrtbigM = all(value([G*x(:,iter) - h + E*value(y) + M*g; pi; x(:,iter); b - G'*pi]) < bigM);
        if corrtbigM ~= 1
            disp('invalid bigM');
            break;
        end
        
        UB = min(UB, c'*value(y) + value(objSP));
        constMP = [constMP, x(:,iter+1) >= 0, eta >= b'*x(:,iter+1), G*x(:,iter+1) - h + E*y + M*value(g) >= 0];
        iter = iter+1;
    end
end

%% Benders
% ------------------------------ master problem ------------------------------
%   min_    {y,eta} c'*y + eta                                      (16)
%   s.t.    A*y >= d, y >= 0                                        (17)
%           eta >= L                                                (18)
%           ets >= pi^{l}'*(h - E*y - M*u^{l}), l = 1,...,k         (19)
%           0 >= pi^{l}'*(h - E*y - M*u^{l}), l = 1,...,k           (20)
%   where   L is an initial lower bound on value function, pi^{l}, l = 1,...,k
%           are the identified extrem point or extrem ray of the second-stage
%           dual polyhedron, u^{l}, l = 1,...,k are the identified extreme point
%           or extreme ray of the polyhedron uncertainty set U.
%           (17)-(18) are the fist-stage constraints, (19) is the optimality cuts
%           (20) is the feasibility cuts.
%
% ------------------------------ subproblem ------------------------------
%   same as (7)-(11) and (12)-(15).
%
x = sdpvar(m*n,1);
LB = -inf; UB = +inf; iter = 1; bigM = 1e3;
disp('Benders');

constMP = [A*y >= d, y >= 0, eta >= -1e3];
objMP = c'*y + eta;
opts = sdpsettings('solver', 'gurobi', 'verbose', 0);

while abs(UB - LB) > tol
    disp(['iters:', num2str(iter)]);
    diag = optimize(constMP, objMP, opts);
    LB = value(objMP);
    
    objSP = sum(s);
    constSP = [];
    constSP = [constSP, G*x + s - h + E*value(y) + M*g >= 0, x >= 0, s >= 0];
    constSP = [constSP, -G'*pi >= 0, 1 >= pi >= 0];
    constSP = [constSP, G*x + s - h + E*value(y) + M*g <= bigM*v, pi <= bigM*(1 - v)];
    constSP = [constSP, x <= bigM*w, -G'*pi <= bigM*(1 - w)];
    constSP = [constSP, s <= bigM*u, 1 - pi <= bigM*(1 - u)];
    constSP = [constSP, 0 <= g <= 1, g(1) + g(2) + g(3) <= 1.8, g(1) + g(2) <= 1.2];
    diag = optimize(constSP, -objSP, opts);
    
    % check the bigM value
    corrtbigM = all(value([G*x + s - h + E*value(y) + M*g; pi; x; -G'*pi; s; 1 - pi]) < bigM);
    if corrtbigM ~= 1
        disp('invalid bigM');
        break;
    end
    
    % If the second stage problem is infeasible under the worst case scenario
    if sum(value(s)) > tol
        constMP = [constMP, 0 >= value(pi)'*(h - E*y - M*value(g))];
        iter = iter+1;
        continue;
    % else the second stage problem is feasible
    else
        objSP = b'*x;
        constSP = [];
        constSP = [constSP, G*x - h + E*value(y) + M*g >= 0, x >= 0];
        constSP = [constSP, b - G'*pi >= 0, pi >= 0];
        constSP = [constSP, G*x - h + E*value(y) + M*g <= bigM*v, pi <= bigM*(1 - v)];
        constSP = [constSP, x <= bigM*w, b - G'*pi <= bigM*(1 - w)];
        constSP = [constSP, 0 <= g <= 1, g(1) + g(2) + g(3) <= 1.8, g(1) + g(2) <= 1.2];
        diag = optimize(constSP, -objSP, opts);
        
        % check the bigM value
        corrtbigM = all(value([G*x - h + E*value(y) + M*g; pi; x; b - G'*pi]) < bigM);
        if corrtbigM ~= 1
            disp('invalid bigM');
            break;
        end
        
        UB = min(UB, c'*value(y) + value(objSP));
        constMP = [constMP, eta >= value(pi)'*(h - E*y - M*value(g))];
        iter = iter+1;
    end
end

disp('done');