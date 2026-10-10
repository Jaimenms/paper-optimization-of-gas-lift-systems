clear all
global X SAIDAS
clc

% Initial Estimate
qoleo1 = 3;
qoleo2 = 3;
qoleo3 = 3;
qoleo4 = 3;
qgas1 = 1;
qgas2 = 1;
qgas3 = 1;
qgas4 = 1;

% qoleo1 = 14;
% qoleo2 = 11;
% qoleo3 = 8;
% qoleo4 = 6;
% qgas1 = 1.5;
% qgas2 = 1.75;
% qgas3 = 1.9;
% qgas4 = 2;
xopt0 = [qoleo1;qoleo2;qoleo3;qoleo4;qgas1;qgas2;qgas3;qgas4];

% Case Specifications
case_04
parametersOpt.options = options;

% Lower and Upper Limits
LB = [0.2;0.2;0.2;0.2;0.2;0.2;0.2;0.2];
UB = [20;20;20;20;10;10;10;10];

% Objective Function Parameters
parametersOpt.pwell = [240;240;240;240]; % bar
parametersOpt.IP    = [1;2.5;1;0.5]; % (kg/s)/bar
parametersOpt.showDisplay = 1;
parametersOpt.ix_WellFlowRate = [1;2;3;4];
parametersOpt.ix_GasFlowRate  = [5;6;7;8];
parametersOpt.ix_PipeDiameter = [];
parametersOpt.iy_PipeDiameter = [];

% Compressor Parameters
parametersOpt.comp.n = 3;
parametersOpt.comp.gama = 1.4;
parametersOpt.comp.T1 = 300; % K
parametersOpt.comp.ef = 0.8;
parametersOpt.comp.P1 = 5; % bar

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
save result_opt_subsea2_v01