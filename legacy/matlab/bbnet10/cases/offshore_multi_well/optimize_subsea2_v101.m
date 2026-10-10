close all
clear all
global X SAIDAS

load result_opt_subsea2_v02

X = [];
SAIDAS = [];
clc

% Initial Estimate
xopt0 = [xopt(1:8);2;2;2;2];

% Lower and Upper Limits
LB = [0.5;0.5;0.5;0.5;0.1;0.1;0.1;0.1;1.0;1.0;1.0;1.0];
UB = [20;20;20;20;4;4;4;4;4;4;4;4];


% Objective Function Parameters
parametersOpt.pwell = [220;220;220;220]; % bar
parametersOpt.IP    = [1;2.5;1;0.5]; % (kg/s)/bar
parametersOpt.showDisplay = 1;
parametersOpt.ix_WellFlowRate = [1;2;3;4];
parametersOpt.ix_GasFlowRate  = [5;6;7;8];
parametersOpt.ix_PipeDiameter = [9;10;11;12];
parametersOpt.iy_PipeDiameter{1} = [parameters.ix_GasPipes(1) parameters.ix_GasPipes(2)];
parametersOpt.iy_PipeDiameter{2} = [parameters.ix_GasPipes(3) parameters.ix_GasPipes(4)];
parametersOpt.iy_PipeDiameter{3} = [parameters.ix_GasPipes(5) parameters.ix_GasPipes(6)];
parametersOpt.iy_PipeDiameter{4} = [parameters.ix_GasPipes(7) parameters.ix_GasPipes(8)];

% Compressor Parameters
parametersOpt.comp.n = 3;
parametersOpt.comp.gama = 1.4;
parametersOpt.comp.T1 = 300; % K
parametersOpt.comp.ef = 0.8;
parametersOpt.comp.P1 = 5; % bar

% Pipes Parameters
parametersOpt.pipe.Fw=1.8;

parametersOpt.injGasCost  = 0;
parametersOpt.fixedCost   = 1;
parametersOpt.fuelGasCost = 1;

opts = optimset('Display','iter','OutputFcn',@stopfunction);
%[fob,saidasOpt,saidas,xout] = objective_function_v03(xopt0, parametersOpt, parameters);
xopt=fmincon(@objective_function_v03,xopt0,[],[],[],[],LB,UB,...
    @(xopti,parametersOpti, parametersi)objective_nonlcon_v01(xopti,parametersOpti, parametersi),...
    opts, parametersOpt, parameters);

%xopt = xopt0;
[fob,saidasOpt,saidas,x] = objective_function_v03(xopt, parametersOpt, parameters);
save result_opt_subsea2_v101