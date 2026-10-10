clear all
global X SAIDAS

load result_opt_subsea2_v01

X = [];
SAIDAS = [];
clc

% Initial Estimate
xopt0 = [xopt];

% Lower and Upper Limits
LB = [0.5;0.5;0.5;0.5;0.1;0.1;0.1;0.1];
UB = [20;20;20;20;4;4;4;4];

xopt=fmincon(@objective_function_v01,xopt0,[],[],[],[],LB,UB,[],opts, parametersOpt, parameters);

%xopt = xopt0;
[fob,saidasOpt,saidas,x] = objective_function_v01(xopt, parametersOpt, parameters);
save result_opt_subsea_v01_2