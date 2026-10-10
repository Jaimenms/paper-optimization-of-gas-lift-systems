close all
clear all
global X SAIDAS
clc

% Initial Estimate
qoleo = 10;
qgas = 2;
xopt0 = [qoleo;qgas];

% Case Specifications
case_03
parametersOpt.options = options;

% Lower and Upper Limits
LB = [5;1];
UB = [20;4];

% Objective Function Parameters
parametersOpt.pwell = 220; % bar
parametersOpt.IP    = 1; % (kg/s)/bar
parametersOpt.showDisplay = 1;
parametersOpt.ix_WellFlowRate = 1;
parametersOpt.ix_GasFlowRate  = 2;
parametersOpt.ix_PipeDiameter = [];

% Compressor Parameters
parametersOpt.comp.n = 3;
parametersOpt.comp.gama = 1.4;
parametersOpt.comp.T1 = 300; % K
parametersOpt.comp.ef = 0.8;
parametersOpt.comp.P1 = 5; % bar

% Pipes Parameters
parametersOpt.pipe.c_L=20; %USD/m
parametersOpt.pipe.Fc=0.5;
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
parametersOpt.fixedCost   = 0;
parametersOpt.fuelGasCost = 0;

opts = optimset('Display','iter','OutputFcn',@stopfunction);
%[fob,saidasOpt,saidas,xout] = objective_function_v03(xopt0, parametersOpt, parameters);
xopt=fmincon(@objective_function_v03,xopt0,[],[],[],[],LB,UB,...
    @(xopti,parametersOpti, parametersi)objective_nonlcon_v01(xopti,parametersOpti, parametersi),...
    opts, parametersOpt, parameters);

%xopt = xopt0;
[fob,saidasOpt,saidas,x] = objective_function_v03(xopt, parametersOpt, parameters);
save result_opt_subsea_v01