close all
clear all
clc
load result_opt_subsea_v01_2

global X SAIDAS

X =[];
SAIDAS =[];


% Initial Estimate
xopt0 = xopt;

% Prices
parametersOpt.gasPrice = 5.34.*0.035315; % USD/m³ wellhead @ US in 2007
parametersOpt.oilPrice = 51.15/4; % USD/Brent barrel

% FOB Selection
parametersOpt.injGasCost  = 0;
parametersOpt.fixedCost   = 0;
parametersOpt.fuelGasCost = 1;

opts = optimset('Display','iter','OutputFcn',@stopfunction);
xopt=fmincon(@objective_function_v02,xopt0,[],[],[],[],LB,UB,[],opts, parametersOpt, parameters);
[fob,saidasOpt,saidas,x] = objective_function_v02(xopt, parametersOpt, parameters);
save result_opt_subsea_v41_2