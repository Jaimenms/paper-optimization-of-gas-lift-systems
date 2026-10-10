function [F,J,saidas]=model_network(x,parameters,rigor)
% Function that determines the residuos of the system of algebraic
% non-linear equations for a two-phase flow pipe network model based on
% Beggs and Brill correlation (1973) via rigorous momentum balance.
%
% Input Variables:
% x -> state variables vector
% parameters -> struct
%  
% Output Variables:
% F -> Residuous vector
persistent zrsSN zrsNS zrsNN unsNN
global dF3dqG dF3dqL

%% (1) Initial Calculations

% Getting the specification vector 
E    = parameters.E;
e    = parameters.e;
M    = parameters.M;
M1 = parameters.M1;
M2 = parameters.M2;
N = parameters.N;
S = parameters.S;

if isempty(zrsSN)
    % Starting zero and identity matrix
    zrsSN = zeros(S,N);
    zrsNS = zeros(N,S);
    zrsNN = zeros(N,N);
    unsNN = eye(N);
end
    
%% (2) Isolating variables
qL  =x(parameters.ixqL,:);
qG  =x(parameters.ixqG,:);
PI  =x(parameters.ixPI,:);
PO  =x(parameters.ixPO,:);
PW  =x(parameters.ixPW,:);
wL  =x(parameters.ixwL,:);
wG  =x(parameters.ixwG,:);

%% (3) Model - part 1 - Mass Balances and Specification Equations
% Balanço de massa
F1 = M*qL+wL;
F2 = M*qG+wG;

% Specifications
F5 = E*x-e;

%% (4) Integrating pressure drop along the pipes
if nargout==3
    [F3,F41,F42,saidas] = momentum_balances(PI,PO,PW,qL,qG,parameters,rigor);
else
    [F3,F41,F42] = momentum_balances(PI,PO,PW,qL,qG,parameters,rigor);
end

%% (5) Generating the residuos vector (model)
F =[F1;F2;F3;F41;F42;F5];

%% (6) Calculating Jacobian (analytical part)
% Liquid Mass balance part
dF1dqL = M;
dF1dqG = zrsNS;
dF1dPI = zrsNS;
dF1dPO = zrsNS;
dF1dPW = zrsNN;
dF1dwL = unsNN;
dF1dwG = zrsNN;

% Gas Mass balance part
dF2dqL = zrsNS;
dF2dqG = M;
dF2dPI = zrsNS;
dF2dPO = zrsNS;
dF2dPW = zrsNN;
dF2dwL = zrsNN;
dF2dwG = unsNN;

% Momentum balance at pipes
dF3dPW = zrsSN;
dF3dwL = zrsSN;
dF3dwG = zrsSN;

% Momentum balance at nodes
dF41dwL = zrsSN;
dF42dwL = zrsSN;
dF41dwG = zrsSN;
dF42dwG = zrsSN;

dF41dPW = M2';
dF42dPW = M1';
    
    
%% (7) Calculating Jacobian (numerical part)
    
% Deviations for numerical Jacobian part
qL_ast = qL.*(1+1e-3)+1e-6;
qG_ast = qG.*(1+1e-3)+1e-6;
PI_ast = PI.*(1+1e-3)+1e-6;
PO_ast = PO.*(1+1e-3)+1e-6;

% Calculating inverse of deviations
aux_qL = 1./(qL_ast-qL);
aux_qG = 1./(qG_ast-qG);
aux_PI = 1./(PI_ast-PI);
aux_PO = 1./(PO_ast-PO);

% Recalculating momentum balances for numerical Jacobian
[F3_1,F41_1,F42_1] = momentum_balances(PI    ,PO    ,PW,qL_ast,qG    ,parameters,rigor);
[F3_2,F41_2,F42_2] = momentum_balances(PI    ,PO    ,PW,qL    ,qG_ast,parameters,rigor);
[F3_3,F41_3,F42_3] = momentum_balances(PI_ast,PO    ,PW,qL    ,qG    ,parameters,rigor);
[F3_4,F41_4,F42_4] = momentum_balances(PI    ,PO_ast,PW,qL    ,qG    ,parameters,rigor);

% Determining numerical derivates
dF3dqL = diag((F3_1-F3).*aux_qL);
dF3dqG = diag((F3_2-F3).*aux_qG);
dF3dPI = diag((F3_3-F3).*aux_PI);
dF3dPO = diag((F3_4-F3).*aux_PO);

% Determining numerical derivates
dF41dqL = diag((F41_1-F41).*aux_qL);
dF41dqG = diag((F41_2-F41).*aux_qG);
dF41dPI = diag((F41_3-F41).*aux_PI);
dF41dPO = diag((F41_4-F41).*aux_PO);
    
% Determining numerical derivates
dF42dqL = diag((F42_1-F42).*aux_qL);
dF42dqG = diag((F42_2-F42).*aux_qG);
dF42dPI = diag((F42_3-F42).*aux_PI);
dF42dPO = diag((F42_4-F42).*aux_PO);

    
%% (8) Constructing Jacobian

J = [ dF1dqL dF1dqG dF1dPI dF1dPO dF1dPW dF1dwL dF1dwG; ...
      dF2dqL dF2dqG dF2dPI dF2dPO dF2dPW dF2dwL dF2dwG; ...
      dF3dqL dF3dqG dF3dPI dF3dPO dF3dPW dF3dwL dF3dwG; ...
      dF41dqL dF41dqG dF41dPI dF41dPO dF41dPW dF41dwL dF41dwG; ...
      dF42dqL dF42dqG dF42dPI dF42dPO dF42dPW dF42dwL dF42dwG; ...
      E];  

function [F3,F41,F42,saidas] = momentum_balances(PI,PO,PW,qL,qG,parameters,rigor)
global vL vG densG
persistent opts

densL= parameters.liquidDensity;
PR   = parameters.PR;
MinP = parameters.MinP;
AI = parameters.AI;
AO = parameters.AO;
M = parameters.M;
M1 = parameters.M1;
M2 = parameters.M2;
N = parameters.S;
S = parameters.S;

% Alocating variables
vLI = zeros(S,1);
vGI = vLI;
vLO = vLI;
vGO = vLI;
densGI = vLI;
densGO = vLI;

if isempty(opts)
    % Determining odeset parameters
    opts = odeset('initialstep',1e-2,...
    'abstol',parameters.OdeAbsTol,...
    'reltol',parameters.OdeRelTol,...
    'maxstep',parameters.MaxStep);
end

% Integrating pressure drop in respect to x
M_L = (0:0.01:1)';
M_dP= [];

switch rigor
    case 0
        for j=1:S
            [dp] = pipe_model(0.5,0,(PO(j)+PI(j))./2.*PR,qL(j),qG(j),parameters,j,1);
            dP(j,1) = dp./PR;
            % Obtaining velocities and density
            vLI(j,1)=vL(:,  1);
            vGI(j,1)=vG(:,  1);
            densGI(j,1)=densG(:,  1);
            vLO(j,1)=vL(:,1);
            vGO(j,1)=vG(:,1);
            densGO(j,1)=densG(:,1);
            saidas = [];
            M_dP(:,j) = linspace(0,dp,length(M_L));
        end
    case 1

        for j=1:S

            % Determining Pressure Drop - from the begining to the end
            [l,d_pressure] = ode45(@pipe_model,[0 1],0,opts,PI(j).*PR,qL(j),qG(j),parameters,j,1);

            % Obtaining velocities and density at the begining
            vLI(j,1)=vL(:,  1);
            vGI(j,1)=vG(:,  1);
            densGI(j,1)=densG(:,  1);

            % Filtering low pressure
            pressure = max(PI(j).*PR+interp1(l,d_pressure,M_L),MinP);
            d_pressure = pressure-PI(j).*PR;


            % Determining Pressure Drop - from the end to the begining
            [l2,d_pressure2] = ode45(@pipe_model,[0 1],0,opts,PO(j).*PR,qL(j),qG(j),parameters,j,-1);

            % Obtaining velocities and density at the begining
            vLO(j,1)=vL(:,1);
            vGO(j,1)=vG(:,1);
            densGO(j,1)=densG(:,1);

            % Filtering low pressure
            pressure2 = max(PO(j).*PR-interp1(l2,d_pressure2,M_L),MinP);
            d_pressure2 = PO(j).*PR-pressure2;

            % Generating the matrix of state variables along the pipes
            M_dP=[M_dP  d_pressure];
            
            % Calculating mean pressure drop
            dP(j,1) = (d_pressure(end,1)+d_pressure2(end,1))./2./PR;

        end
end

% Calculating Saidas
[d_pressure,saidas] = pipe_model(M_L,M_dP,PI.*PR,qL,qG,parameters,1:S,1);        


% Momentum balance at the pipes
F3 = PO - PI - dP;

% Momentum balance at nodes
F41 = M2'*PW-PO-(abs(vLO)./(2.*densL.^2./AO).*qL + abs(vGO)./(2.*densGO.^2./AO).*qG )./PR; 
F42 = M1'*PW-PI-(abs(vLI)./(2.*densL.^2./AI).*qL + abs(vGI)./(2.*densGI.^2./AI).*qG )./PR; 