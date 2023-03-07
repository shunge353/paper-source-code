% Test System data for the Equilibria in an Oligopolistic Electricity Pool with Stepwise Offer Curves

%% system MVA base
baseMVA = 100;

%% branch data
% The index and meaning of each column of the branch matrix is given below:
%   1   "from" bus number
%   2   "to" bus number
%   3   susceptance (p.u.)
%   4   transmission capacity (MW)

%   fbus    tbus    b   rateA
branch = [
    1   2   5   30;
    2   1   5   30;
];

% change unit of branch susceptance to siemens
branch(:, 3) = baseMVA * branch(:, 3);

%% generator data
% The index and meaning of each column of the gen matrix is given below:
%   1   bus number
%   2   capacity of block 1 (MW)
%   3   capacity of block 2 (MW)
%   4   marginal cost of block 1 ($/MWh)
%   5   marginal cost of block 2 ($/MWh)

%   bus   PG1max   PG2max   lambdaG1   lambdaG2
gen = [
    1   25   25   15   19;
    2   25   25   17   21;
];

%% demand data
% The index and meaning of each column of the load matrix is given below:
%   1   bus number
%   2   capacity of block 1 (MW)
%   3   capacity of block 2 (MW)
%   4   marginal utility of block 1 ($/MWh)
%   5   marginal utility of block 2 ($/MWh)

%   bus   PD1max   PD2max   lambdaD1   lambdaD2
load = [
    1   30   20   20   17;
    2   30   20   23   18;
];
