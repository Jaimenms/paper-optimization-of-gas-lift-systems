clear all
close all

load result_opt_subsea_v01_2

global X SAIDAS


X = [];
SAIDAS = [];

% Initial Estimate
xopt0 = [xopt(1);xopt(2);2];

% Case Specifications
case_03
parametersOpt.options = options;

% Lower and Upper Limits
LB = [5;1;1.5];
UB = [20;4;5];
parametersOpt.n_GasPipePerWell = 2;

% Objective Function Parameters
parametersOpt.pwell = 220; % bar
parametersOpt.showDisplay = 1;
parametersOpt.ix_WellFlowRate = 1;
parametersOpt.ix_GasFlowRate  = 2;
parametersOpt.ix_PipeDiameter = 3;

% Indices da tubulações a serem associadas às variaveis de decisão
parametersOpt.iy_PipeDiameter{1} = [parameters.ix_GasPipes(1) parameters.ix_GasPipes(2)];

% Compressor Parameters
parametersOpt.comp.n = 3;
parametersOpt.comp.gama = parameters.gasPolytrophicCoef;
parametersOpt.comp.T1 = 300; % K
parametersOpt.comp.ef = 0.8;
parametersOpt.comp.P1 = 5; % bar

% Pipes Parameters
parametersOpt.pipe.c_L=100; %USD/m
parametersOpt.pipe.Fc=5.0;
parametersOpt.pipe.Fw=1.8;

% Prices
parametersOpt.gasPrice = 5.34.*0.035315; % USD/m³ wellhead @ US in 2007
% Gas Price
% http://www.eia.doe.gov/emeu/aer/txt/ptb0607.html
parametersOpt.oilPrice = 51.15; % USD/Brent barrel
% Petroleum Price
% http://www.eia.doe.gov/emeu/international/Crude1.xls
% 02/04/2009 - The New York Times

% FOB Selection
parametersOpt.injGasCost  = 0;
parametersOpt.fixedCost   = 1;
parametersOpt.fuelGasCost = 1;

opts = optimset('Display','iter','OutputFcn',@stopfunction);

xopt=fmincon(@objective_function_v02,xopt0,[],[],[],[],LB,UB,[],opts, parametersOpt, parameters);

[fob,saidasOpt,saidas,x] = objective_function_v02(xopt, parametersOpt, parameters);
save result_opt_subsea_v101