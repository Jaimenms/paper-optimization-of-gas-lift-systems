clear all
global X SAIDAS

load result_opt_subsea2_v01_4

X = [];
SAIDAS = [];
clc

% Initial Estimate
xopt0 = [xopt];

xopt=fmincon(@objective_function_v02,xopt0,[],[],[],[],LB,UB,[],opts, parametersOpt, parameters);

%xopt = xopt0;
[fob,saidasOpt,saidas,x] = objective_function_v02(xopt, parametersOpt, parameters);
save result_opt_subsea2_v01_5