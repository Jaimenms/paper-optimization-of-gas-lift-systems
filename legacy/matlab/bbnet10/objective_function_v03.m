function [fob,saidasOpt,saidas,xout] = objective_function_v02(xopt, parametersOpt, parameters)
global x saidasOpt X SAIDAS 
persistent fob_last status_last

%% 1 - Initial Calculations

% Starting initial value of objective function
if isempty(fob_last),fob_last = inf; end
if isempty(status_last),status_last.errorFlag=0; end

showDisplay = parametersOpt.showDisplay;

%% 2 - Flow Model Simulation

% Executing Simulation
[x,saidas,status] = objective_simulation(xopt,x,parametersOpt,parameters);
xout = x;

%% 3-Fixed Cost Estimation
[C_c,C_d,W_c,P2] = calculate_compressor(x, parametersOpt, parameters);

if parametersOpt.fixedCost

    % Calculating Costs
    [C_p, L_p, D_p] = calculate_pipes(xopt,parametersOpt, parameters);

    %C_c -> Fixed-capital cost with compressors
    %C_p -> Fixed-capital cost with gas pipes
    ISBL = sum(C_c)+sum(C_d)+sum(C_p); %kUSD
    OSBL = 3.*ISBL;

    %A_TCC -> Annual Fixed-capital Cost
    CCF = 1/3;
    A_TCC = CCF*OSBL; %1e6kUSD/year

else
    C_c =0;
    C_p = 0;
    C_d = 0;
    ISBL = 0;
    OSBL = 0;
    A_TCC = 0;
    D_p = 0;
    L_p = 0;
end


%% 4-Production Estimation

%R  -> Rate of Production
barrel2m3 = 0.1589873;
m_o = x(parameters.ix_qLwell,1)./parameters.liquidDensity./barrel2m3; % barrel/s
m2_o = x(parameters.ix_qLwell,1); % kg/s
R = m_o.*(3600.*24.*365)./1e6; %1e6 barrel/year


%% 4-Operational Cost Estimation
[f_g_lift_calc,m_g_lift_calc] = calculate_injectedGas(x,parameters);
if parametersOpt.injGasCost
    f_g_lift=f_g_lift_calc;
    m_g_lift=m_g_lift_calc;
else
    f_g_lift = 0;
    m_g_lift = 0;
end

if parametersOpt.fuelGasCost
    [f_g_turb] = calculate_fuelGas(x,W_c,parameters);
else
    f_g_turb = 0;
end

%% 5-Profit Estimation

%c_S    -> Petroleum Sales Price per Unit
c_S = parametersOpt.oilPrice; % USD/Brent barrel
%c_VE_g -> Variable Expense with Natural Gas Cost per Energy Unit
c_VE_g = parametersOpt.gasPrice; % USD/m³

in  = sum(R.*c_S); %1e6USD/year
out = ((sum(f_g_lift)+f_g_turb).*c_VE_g); %1e6USD/year

A_NP = in - out - A_TCC; %1e6USD/year

fob = -A_NP;

%% 6-Penalties
[fob,penalty_D] = divergency_penalty(fob,fob_last,status,status_last);
fob_last = fob;
status_last = status;

pdownhole = x(parameters.ix_Pwell).*parameters.PR./1e5;
pwell     = parametersOpt.pwell;
q_viaPI   = parametersOpt.IP.*(pwell-pdownhole);


%% 7-Penalties

if showDisplay
    disp(' ')
    disp('################################################')
	disp(['Fob:                       ' num2str(fob) ' 1e6 USD/year'])
	disp(['Profit:                    ' num2str(A_NP) ' 1e6 USD/year'])
	disp(['InValue:                   ' num2str(in) ' 1e6 USD/year'])
	disp(['OutValue:                  ' num2str(out) ' 1e6 USD/year'])
	disp(['Annual Total-capital Cost: ' num2str(A_TCC) ' 1e6 USD/year'])
	disp(['Annual Injected Gas Cost:  ' num2str(f_g_lift(:)'.*c_VE_g) ' 1e6 USD/year'])
	disp(['Annual Fuel Gas Cost:      ' num2str(f_g_turb.*c_VE_g) ' 1e6 USD/year'])
	disp(['ISBL:                      ' num2str(ISBL) ' 1e6 USD'])
	disp(['OSBL:                      ' num2str(OSBL) ' 1e6 USD'])
    
    disp(['Oil Flow Rate:        ' num2str(R(:)') ' 1e6 barrel/year'])
    disp(['Oil Flow Rate:        ' num2str(m2_o(:)') ' kg/s'])
    disp(['Oil Flow Rate via IP: ' num2str(q_viaPI(:)') ' kg/s'])
    disp(['Gas Lift Flow Rate:   ' num2str(f_g_lift_calc(:)') ' 1e6 m³/year'])
    disp(['Gas Lift Flow Rate:   ' num2str(m_g_lift_calc(:)') ' kg/s'])
    disp(['Fuel Gas Flow Rate:   ' num2str(f_g_turb(:)') ' 1e6 m³/year'])
    
	disp(['Compressor Power:          ' num2str(W_c) ' kW'])
	disp(['Discharge Pressure:        ' num2str(P2) ' bar'])
	disp(['Downhole Pressure:         ' num2str(pdownhole(:)') ' bar'])
	disp(['Well Pressure:             ' num2str(pwell(:)') ' bar'])
    
    disp(['Compressor-capital Cost: ' num2str(sum(C_c)) ' 1e6 USD/year'])
    disp(['Turbine-capital Cost:    ' num2str(sum(C_d)) ' 1e6 USD/year'])
    disp(['Pipes-capital Cost:      ' num2str(sum(C_p)) ' 1e6 USD/year'])
    disp(['Pipes length:   ' num2str(L_p(:)') ' m'])
    disp(['Pipes diameter: ' num2str(D_p(:)') ' m'])

    disp(['Divergence Pennalty: ' num2str(penalty_D) ' 1e6 USD/year'])
    disp(['Model Residue:       ' num2str(status.residue ) ' (-)'])
    disp('################################################')
    disp(' ')
    
end

saidasOpt.Fob=fob;
saidasOpt.Profit=A_NP;
saidasOpt.InValue=in;
saidasOpt.OutValue=out;
saidasOpt.ISBL=ISBL;
saidasOpt.OSBL=OSBL;
saidasOpt.AnnualTotalCapitalCost=A_TCC;
saidasOpt.AnnualInjectedGasCost=f_g_lift(:).*c_VE_g;
saidasOpt.AnnualFuelGasCost=f_g_turb.*c_VE_g;

saidasOpt.OilFlowRate=R;
saidasOpt.m2_o=m2_o;

saidasOpt.GasLiftFlowRate=f_g_lift_calc;
saidasOpt.GasLiftFlowRate2=m_g_lift_calc;
saidasOpt.FuelGasFlowRate=f_g_turb;

saidasOpt.CompressorPower=W_c;
saidasOpt.DischargePressure=P2;
saidasOpt.DownholePressure=pdownhole;
saidasOpt.WellPressure=pwell;

saidasOpt.CompressorCapitalCost=sum(C_c);
saidasOpt.TurbineCapitalCost=sum(C_d);
saidasOpt.PipesCapitalCost=sum(C_p);
saidasOpt.Pipeslength=L_p(:);
saidasOpt.Pipesdiameter=D_p(:);

saidasOpt.DivergencePennalty = penalty_D;
saidasOpt.Residue = status.residue;

outputfunction;


X = [X x];
n = length(SAIDAS);
SAIDAS{n+1}=saidasOpt;



%% Subfunctions
function [C_c,C_d,W_c,P2] = calculate_compressor(x, parametersOpt, parameters)
% Function that calculates the compression cost


% (1) Determining Power Consumption

% n-> number of stages
% gama -> polytrophic coefficient
% MW -> gas molecular weight (kg/kmol)
% R -> universal gas constant (8.314 kJ/(kmol K))
% T1-> initial absolute temperature (K)
% ef-> compressor efficiency
% P1-> initial absolute pressure (bar)
% P2-> final absolute pressure (bar)
% Q -> gas mass flow rate (kg/s)
% w -> compression work (kJ/kmol)
% W_c -> compression power (kW)

R = 8.314;

if isfield(parameters,'MW')
    MW   = parameters.MW;
else
    MW = 25;
end    

if isfield(parametersOpt,'comp')
    n    = parametersOpt.comp.n;
    gama = parametersOpt.comp.gama;
    T1   = parametersOpt.comp.T1;
    ef   = parametersOpt.comp.ef;
    P1   = parametersOpt.comp.P1;
else
    n = 3;
    gama = 1.4;
    T1 = 300;
    ef = 0.8;
    P1 = 15;
    MW = 25;
end

P2  = max(x(parameters.ix_Pgaslift).*parameters.PR./1e5); % bar
Q   = sum(x(parameters.ix_qGgaslift)); % kg/s

w = n.*gama.*R.*T1./ef./(gama-1).*((P2./P1).^((gama-1)./(n.*gama))-1);
W_c = w.*(Q./MW); %kJ/s = kW


% (2) Determining Fixed Cost

% Fc -> Correction factor according to compressor type
% bhp -> power (brake horsepower)
bhp = 1.34102209.*W_c;

C_c = 14.3e-3.*bhp.^0.62; %1e6USD
C_d = 7.19e-3.*bhp.^0.81; %1e6USD

function [C_p,L_p,D_p] = calculate_pipes(xopt,parametersOpt,parameters)
% Function that determines the cost of gas pipes

% L -> pipe length (m)
% D -> pipe diameter (m)
% c_L -> cost per length (USD/m)
% c_D -> cost per diameter (USD/m(diameter)/m(length)
% Fc -> correction factor due to diameter
% Fw -> construction factor

if ~isfield(parametersOpt,'iy_PipeDiameter')
    parametersOpt.iy_PipeDiameter=[];
end

nx = length(parametersOpt.ix_PipeDiameter);
if nx>0
    D_p = zeros(nx,1);
    L_p = D_p;
    for i=1:nx
        D_p(i,1) = xopt(parametersOpt.ix_PipeDiameter(i));
        L_p(i,1)=0;
        aux = parametersOpt.iy_PipeDiameter{i};
        ny = length(aux);
        for j=1:ny
            L_p(i,1) = L_p(i,1)+parameters.pipe_length(aux(j),1);
        end
    end

%     Estimated model
%     c_L = 30.*D_p+40;
%     C_FC_p = sum(L_p.*c_L)./1e6; %1e6USD
%     C_WC_p = Fw.*C_FC_p; %1e6USD
    
    % McCoy model
    L_p_km = L_p./1000;
    D_nps = D_p;
    % Considering SW of USA
    a0 =  3.112;  a4 =  0.000;  a6 =  0.901;  a7 =  1.590;
    b = 10.^(a0+a4);
    C_p_materials = sum(b.*(L_p_km.^a6).*(D_nps.*a7));

    a0 =  4.487;  a4 =  -0.216;  a6 =  0.820;  a7 =  0.940;
    b = 10.^(a0+a4);
    C_p_labor = sum(b.*(L_p_km.^a6).*(D_nps.*a7));

    a0 =  3.950;  a4 =   0.000;  a6 =  1.049;  a7 =  0.403;
    b = 10.^(a0+a4);
    C_p_row = sum(b.*(L_p_km.^a6).*(D_nps.*a7));

    a0 =  4.390;  a4 =   0.000;  a6 =  0.783;  a7 =  0.791;
    b = 10.^(a0+a4);
    C_p_misc = sum(b.*(L_p_km.^a6).*(D_nps.*a7));

    % Disconsidering Right-of-way costs
    C_p = sum(C_p_materials+C_p_labor+C_p_misc)./1e6; %1e6USD
    
else
    D_p = 0;
    L_p = 0;
    C_p = 0;
end






function [fob2,penalty_D] = divergency_penalty(fob,fob_last,status,status_last)
% Penalty due to model divergency
% if status.errorFlag&&status_last.errorFlag==0
%     penalty_D = fob_last + 1e3 - fob;
%     fob2 = fob_last + 1e3;
% else
%     fob2 = fob;
%     penalty_D=0;
% end

if status.errorFlag
    fob2 = 1e6;
    penalty_D = fob2 - fob;
    
else
    fob2 = fob;
    penalty_D=0;
end


function [f_g_turb,m_g_turb] = calculate_fuelGas(x,W_c,parameters)
% Function that calculates gas consumptions as fuel gas
%f_g_turb -> gas consumption at compressor driver

f_g_turb = W_c.*8.45e-3; %1e6m³/year

function [f_g_lift, m_g_lift] = calculate_injectedGas(x,parameters)
% Function that calculates gas consumptions as injected gas

%f_g_lift -> gas consumption for gas lift
if isfield(parameters,'gasNormalDensity')
    gasNormalDensity = parameters.gasNormalDensity;
else
    gasNormalDensity = 1.2; %kg/m³
end


m_g_lift = x(parameters.ix_qGgaslift); % kg/s
f_g_lift = m_g_lift./gasNormalDensity.*(3600.*24.*365)./1e6; % 1e6m³/year