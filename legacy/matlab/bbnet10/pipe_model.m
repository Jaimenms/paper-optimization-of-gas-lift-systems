function [dPdx,outputs]  = pipe_model(position, d_pressure, pressure, liquidMassRate, gasMassRate,  parameters, ixPipe, way)
% Function that determines the pressure gradient as a function of the pipe
% pressure (pressure), the liquid mass flow rate (liquidMassRate), the gas mass
% flow rate (gasMassRate), a set of parameters presented a parameter struct
% (parameters), the pipe index (ixPipe) and the way described with a flag
% (way). The model adopted is the Beggs and Brill correlation (1973).
%
% It is important to highlight that this function is able to execute a
% forward integration (from the begining to the end of the pipe) or a
% backward integration (from the end to the begining). For this
% functionality a flag variable way is used with value 1 indicating forward
% integration and equal to -1 indicating backward integration.
%
% The pressure pressure corresponds to the value at an extremity and depends on
% the direction of integration. For a foward integration the value of pressure
% corresponds to the pressure at the inlet of the pipe. For backward pressure
% corresponds to the pressure at the outlet of the pipe.
%
% Input Variables:
% position -> axial position along the pipe (is equal to zero at
% the pipe inlet and equal to 1 at the pipe outlet)
% pressure -> pressure in Pa
% liquidMassRate -> mass flow rate of the liquid phase in kg/s
% gasMassRate -> mass flow rate of the gas phase in kg/s
% parameters -> struct of model parameters
% ixPipe -> pipe index vector
% way -> direction flag (+1 mean forward and -1 means backward)
%  
% Output Variables:
% dPdx -> gradient of pressure in respect to normalized axial coordinate
% outputs -> struct with model outputs
%
% By Jaime Souza at April 3rd 2009.


%% (1) Initial calculations

% Defining global variables ( necessary for the pipe network model )
global vL vG densG

% Determining sizes
nt = length(position); % number of abcissas
nx = length(ixPipe); % number of pipes

% Transposing states for final output of the ode solver
if nt>1
    position=position(:)';
    d_pressure=d_pressure';
end

% Resizing data
uns=ones(1,nt);
qL = liquidMassRate(:)*uns;
qG = gasMassRate(:)*uns;
pressure  = pressure(:)*uns;

% Recalculating the position according to the integration direction
switch way
    case 1
        position2 = position;
        pressure = pressure + d_pressure;
    case -1
        position2 = 1- position;
        pressure = pressure - d_pressure;
end

% Filtering pressures smaller than MinP
MinP = parameters.MinP./parameters.PR;
ix = pressure<MinP;
pressure(ix) = MinP;

% Unpacking parameters
gasPolytrophicCoef     = parameters.gasPolytrophicCoef;
referencePressure      = parameters.referencePressure;
liquidDensity          = parameters.liquidDensity;
gasReferenceDensity    = parameters.gasReferenceDensity;
L(:,1)                 = parameters.pipe_length(ixPipe);


%% (2) Interpolating Pipe Geometry Information
for i=1:nx
    auxL = parameters.pipe_L{ixPipe(i)};
    rough(i,:)  =interp1(auxL,parameters.pipe_rough{ixPipe(i)},position2.*L(i));
    teta(i,:)   =interp1(auxL,parameters.pipe_teta{ixPipe(i)},position2.*L(i));
    D(i,:)      =interp1(auxL,parameters.pipe_D{ixPipe(i)},position2.*L(i));
    factor(i,:) =(1+parameters.pipe_Leq(ixPipe(i))./L(i)).*ones(1,nt);
end
% Calculating Area
Area = pi.*D.^2./4;

%% (6) Determining closure informations using Beggs and Brill

% Packing parameters
parameters2                  = parameters;
parameters2.diameter         = D;
parameters2.roughness        = rough;
parameters2.inclinationAngle = teta;

% Calculating densities
densL=liquidDensity;
densG=gasReferenceDensity.*(pressure./referencePressure).^...
    (1/gasPolytrophicCoef);
Dp=densG./gasPolytrophicCoef./pressure;

% Calculating superficial velocities
vSL = qL./Area./densL;
vSG = qG./Area./densG;

% Calculating the closure relations
[tension,densS,vL,vG,HL,HG]  = ...
  closure_relations_bb(pressure, vSL, vSG,  parameters2);

%% (6.1) Calculating numerical derivates

% Determining velocities deviations for numerical derivate calculation
vSG_ast = vSG.*(1+1e-6)+1e-9;

% Determining DvsG numerically
[tension_ast,densS_ast,vL_ast,vG_ast,HL_ast]  = ...
  closure_relations_bb(pressure, vSL, vSG_ast,  parameters2);

DvsG = (HL_ast-HL)./(vSG_ast-vSG);

%% (7) Momentum Balance Calculation

% Constants definition
g    = 9.80665; % m/s²

% Determining compressibility coeficients
aL = (Dp.*vSG./densG).*((abs(vL).*DvsG  )./(HL+1e-6));
aG =-(Dp.*vSG./densG).*((abs(vG).*DvsG./(HG+1e-6))+1);

% Calculating dPdx (where x is dimentionless)
dPdL = -(Area.*densS.*g.*sin(teta)+pi.*D.*tension.*factor)./...
    (Area+aL.*qL+aG.*qG);

% Calculating dPdL (where L is the pipe length)
dPdx = diag(L)*dPdL;

% Filtering pressure drop to forbid lower pressures
ix = pressure==MinP&dPdx<0;
dPdx(ix) = 0;

% Generating output
outputs.dPdx = dPdx;
outputs.Area = Area;
outputs.densS = densS;
outputs.HL = HL;
outputs.HG = HG;
outputs.teta = teta;
outputs.pressure = pressure;
outputs.vSL = vSL;
outputs.vSG = vSG;
outputs.qL = qL;
outputs.qG = qG;
outputs.densG = densG;
outputs.tension = tension;

function [tension,densS,vL,vG,HL,HG,outputs]  = closure_relations_bb(pressure, vSL, vSG,  parameters)
% Function that determines shear stress (tension), the mixture density 
% (densS), the phases velocities and the hold-ups according to
% Beggs and Brill correlation (1973).
%
% Input Variables:
% pressure -> pressure in Pa
% vSL -> Liquid superficial velocity (m/s)
% vSG -> Gas superficial velocity (m/s)
% parameters -> struc com os seguintes campos
%   .diameter -> pipe diameters in m
%   .roughness -> pipe roughness in m
%   .inclinationAngle -> inclination angle in relation to horizon plane
%   (rad)
%   .gasPolytrophicCoef -> gas polytrophic coefficient (-)
%   .interfatialTension -> interfatial tension (Pa.m)
%   .liquidViscosity -> liquid viscosity (Pa.s)
%   .gasViscosity -> gas viscosity (Pa.s)
%   .referencePressure -> reference pressure (Pa)
%   .liquidDensity -> liquid density
%   .gasReferenceDensity -> gas reference density
%  
%
% Output Variables:
% tension -> shear stress in Pa
% densS -> mixture density
% vL -> liquid velocity
% vG -> gas velocity
% HL -> liquid hold-up
% HG -> gas hold-up

%% (1) Initial Calculations

% Constants definition
g    = 9.80665; % m/s²

% Filtering pressures lower then 1 bar
pressure = max(1e5,pressure);

% Unpacking parameters
diameter               = parameters.diameter;
roughness              = parameters.roughness;
inclinationAngle       = parameters.inclinationAngle;
gasPolytrophicCoef     = parameters.gasPolytrophicCoef;
interfatialTension     = parameters.interfatialTension;
liquidViscosity        = parameters.liquidViscosity;
gasViscosity           = parameters.gasViscosity;
referencePressure      = parameters.referencePressure;
liquidDensity          = parameters.liquidDensity;
gasReferenceDensity    = parameters.gasReferenceDensity;

% Direction
direction = zeros(size(vSL));
direction(vSL>0)  =  1;
direction(vSL<=0) = -1;

vSL = abs(vSL);
vSG = abs(vSG);

% Cross section area calculation
A = pi.*diameter.^2./4;

% Densities Calculation
densL=liquidDensity;
densG=gasReferenceDensity.*(pressure./referencePressure).^...
    (1/gasPolytrophicCoef);

% Mixture velocity calculation
vm = vSL + vSG;


%% (2) Calculating non-slip variables

% Non-slip hold-up 
lambdaL = vSL ./ (vm+1e-9);
lambdaG = vSG ./ (vm+1e-9);

ix = lambdaL>1;
lambdaL(ix) = 1;
lambdaG(ix) = 0;

ix = lambdaG>1;
lambdaG(ix) = 1;
lambdaL(ix) = 0;

% Mixture properties
densN = densL.*lambdaL + densG.*lambdaG;
viscN = liquidViscosity.*lambdaL + gasViscosity.*lambdaG;

% Dimensionless numbers
NFR = vm.^2./(g.*diameter);% Froude Number
NLV = vSL.*(densL./interfatialTension./g).^(1/4); % Liquid Velocity Number
NRE = densN.*vm.*diameter./viscN; % Mixture Reynolds number

%% (3) Flow Pattern Identification
[fl,L1,L2,L3] =  pattern_identification(lambdaL,NFR);

%% (4) Liquid Hold-up for horizontal flow
[HL0] =  calculate_HL0(lambdaL,NFR);

%% (5) Inclination Correction factor
[qsi] =  calculate_qsi(lambdaL,NFR,NLV,inclinationAngle,direction);

%% (6) Liquid Hold-up for any inclination and phases velocities
[HL,HG] = calculate_HL(NFR,fl,HL0,qsi,L2,L3);
ix1 = find(HG<1e-6);
HG(ix1)=1e-6;
HL(ix1)=1-1e-6;
ix2 = find(HL<1e-6);
HL(ix2)=1e-6;
HG(ix2)=1-1e-6;

vL = direction.*vSL./HL;
vG = direction.*vSG./HG;

%% (7) Mixture density calculation
densS = HL.*densL+HG.*densG;

%% (7) Single fase Friction Factor Calculation
A1=(2.457*log(1./(0.27*roughness./diameter+(7./NRE).^0.9))).^16;
B1=(37530./NRE).^16;
f=4.*(2.*((8./NRE).^12+1./(A1+B1).^1.5).^(1./12));

%% (8) Beggs and Brill Friction Factor Calculation
k = log(lambdaL./(HL.^2));
S = log(2.2.*exp(k)-1.2);
ix = find(k<=0|k>=0.1823);
S(ix) = (-0.0523.*k(ix).^(-1)+3.1820-0.8725.*k(ix)+0.01853.*k(ix).^3).^(-1);
S(ix1) = 0; S(ix2) = 0;
ftp = f .*exp(S);

%% (9) Beggs and Brill Shear Stress Calculation

tension = direction.*(ftp.*densS.*vm.^2)./8;

outputs.pressure  = pressure;
outputs.HL = HL;
outputs.HG = HG;
outputs.NFR = NFR;
outputs.NLV = NLV;
outputs.NRE = NRE;
outputs.lambdaL = lambdaL;
outputs.lambdaG = lambdaG;
outputs.densL = densL;
outputs.densG = densG;
outputs.densS = densS;
outputs.fl = fl;
outputs.vm = vm;


%% Subfunction
function [fl,L1,L2,L3,reg] =  pattern_identification(lambdaL,NFR)
% Function that determines the flow pattern according to Beggs and Brill
% correlation (1973)
%
% Input Variables:
% lambdaL -> non-slip liquid hold-up
% NFR -> Liquid Froude number
%
% Output Variables:
% fl -> struct variable with the following fields:
%   .Tra -> is near one when transition pattern, otherwise is near zero.
%   .Est -> is near one when segregatted pattern, otherwise is near zero.
%   .Int -> is near one when Intermittent pattern, otherwise is near zero.
%   .Dis -> is near one when Distributed pattern, otherwise is near zero.

% (Sub1 - 1) Flow Pattern identification limits
L1 = 316      .*lambdaL.^(+0.3020);
L2 = 0.0009252.*lambdaL.^(-2.4684);
L3 = 0.10     .*lambdaL.^(-1.4516);
L4 = 0.5      .*lambdaL.^(-6.7380);

% (Sub1 - 2) Sigmoidal-if tests
factorExp = 1e6; 
testA = sigmoidal(factorExp,lambdaL,0.01);
testB = sigmoidal(factorExp,lambdaL,0.4);
testC = sigmoidal(factorExp,L1,NFR);
testD = sigmoidal(factorExp,L2,NFR);
testE = sigmoidal(factorExp,L3,NFR);
testF = sigmoidal(factorExp,L4,NFR);

% (Sub1 - 3) Pattern Flags generation

% Segregated
fl.Est = (1-testA).*testC + testA.*testD;
% Transition
fl.Tra = testA.*testE.*(1-testD);
% Intermittent
fl.Int = (1-testB).*testA.*testC.*(1-testE) + testB.*testF.*(1-testE);
% Distributed
fl.Dis = (1-testB).*(1-testC) + testB.*(1-testF);

% (Sub1 - 4) Flag filtering
fl.Tra=max(fl.Tra,0);fl.Tra=min(fl.Tra,1);
fl.Est=max(fl.Est,0);fl.Est=min(fl.Est,1);
fl.Int=max(fl.Int,0);fl.Int=min(fl.Int,1);
fl.Dis=max(fl.Dis,0);fl.Dis=min(fl.Dis,1);

% (Sub1 - 5) Pattern Identifier
reg = 1.*fl.Est + 2.*fl.Tra + 3.*fl.Int + 4.*fl.Dis;

function [F] = sigmoidal(factorExp,x1,x2)
% Sigmoidal function that eliminates discontinuities during if-tests:
% if x1 = x2, then F = 0.5
% elseif x1 > X2, then F tends to 1
% elseif x1 < x2, then F tends to 0
%
% Input variables:
% x1, x2 -> coordinates
%
% Output variables:
% F -> function result
%
% Parameter:
% fatorExp -> factor that smooths the discontinuity

% (Sub2 - 1) Calculation
F= 1./(1+exp(-factorExp.*(x1-x2)));


function [HL0] =  calculate_HL0(lambdaL,NFR)
% Function that determines the liquid hold-up for horizontal flow
% according to Beggs and Brill correlation (1973)
%
% Input Variables:
% lambdaL -> non-slip liquid hold-up
% NFR -> Liquid Froude number
%
% Output Variables:
% HL0 -> Struct with liquid hold-up for horizontal flow for different
% patterns
%   .Est -> horizontal hold-up for Segregated pattern
%   .Int -> horizontal hold-up for Intermittent pattern
%   .Dis -> horizontal hold-up for Distributed pattern

% (Sub3 - 1) Calculating liquid hold-up
a = 0.98; b = 0.4846; c = 0.0868;
HL0.Est = max(a.*(lambdaL.^b)./(NFR.^c),lambdaL);
a = 0.845; b = 0.5351; c = 0.0173;
HL0.Int = max(a.*(lambdaL.^b)./(NFR.^c),lambdaL);
a = 1.065; b = 0.5824; c = 0.0609;
HL0.Dis = max(a.*(lambdaL.^b)./(NFR.^c),lambdaL);

% (Sub3 - 2) Liquid hold-up filtering
HL0.Est=max(HL0.Est,lambdaL);HL0.Est=min(HL0.Est,1);
HL0.Int=max(HL0.Int,lambdaL);HL0.Int=min(HL0.Int,1);
HL0.Dis=max(HL0.Dis,lambdaL);HL0.Dis=min(HL0.Dis,1);


function [qsi] =  calculate_qsi(lambdaL,NFR,NLV,inclinationAngle,direction)
% Function that determines the inclination correction factor
% according to Beggs and Brill correlation (1973)
%
% Input Variables:
% lambdaL -> non-slip liquid hold-up
% NFR -> Liquid Froude number
% NLV -> Liquid velocity number
% inclinationAngle -> pipe inclination in relation to horizon plane
% direction -> flow direction (+1 meand the same direction of pipe and -1
% means opposite direction)
%
% Output Variables:
% qsi -> Struct with inclination correction factors for different
% patterns
%   .Est -> horizontal hold-up for Segregated pattern
%   .Int -> horizontal hold-up for Intermittent pattern
%   .Dis -> horizontal hold-up for Distributed pattern

% (Sub4 - 1) Determining if uphill or downhill
flagAsc = (1 + sign(direction).*sign(inclinationAngle))./2;
flagDesc = 1-flagAsc;

% (Sub4 - 2) Determining correction for downhill
Cd = 4.70; Ce = -0.3692; Cf = 0.1244; Cg = -0.5056;
C = (1-lambdaL).*log(Cd.*(lambdaL.^Ce).*(NLV.^Cf).*(NFR.^Cg));
C = max(C,0);
qsiDesc = 1 + C.*(sin(1.8.*inclinationAngle)-...
    0.333.*(sin(1.8.*inclinationAngle)).^3);

% (Sub4 - 3) Determining correction for segregated
Cd = 0.011; Ce = -3.768; Cf = 3.539; Cg = -1.614;
C = (1-lambdaL).*log(Cd.*(lambdaL.^Ce).*(NLV.^Cf).*(NFR.^Cg));
C = max(C,0);
qsiAsc = 1 + C.*(sin(1.8.*inclinationAngle)-...
    0.333.*(sin(1.8.*inclinationAngle)).^3);
qsi.Est = qsiAsc.*flagAsc + qsiDesc.*flagDesc;

% (Sub4 - 4) Determining correction for intermittent
Cd = 2.96; Ce = 0.305; Cf = -0.4473; Cg = 0.0978;
C = (1-lambdaL).*log(Cd.*(lambdaL.^Ce).*(NLV.^Cf).*(NFR.^Cg));
C = max(C,0);
qsiAsc = 1 + C.*(sin(1.8.*inclinationAngle)-...
    0.333.*(sin(1.8.*inclinationAngle)).^3);
qsi.Int = qsiAsc.*flagAsc + qsiDesc.*flagDesc;

% (Sub4 - 5) Determining correction for distributed
qsiAsc = 1;
qsi.Dis = qsiAsc.*flagAsc + qsiDesc.*flagDesc;


function [HL,HG] = calculate_HL(NFR,fl,HL0,qsi,L2,L3)
% Function that determines the liquid and the gas hold up for any
% inclination according to Beggs and Brill correlation (1973)
%
% Input Variables:
% lambdaL -> non-slip liquid hold-up
% NFR -> Liquid Froude number
% NLV -> Liquid velocity number
% fl  -> struct with flow pattern flags
% HL0 -> struct with horizontal hold-up
% qsi -> struct with inclination correction factors
% L2, L3 -> pattern test factors
%
% Output Variables:
% HL,HG -> Structs with liquid and gas hold-up for any inclination
%   .Est -> horizontal hold-up for Segregated pattern
%   .Int -> horizontal hold-up for Intermittent pattern
%   .Dis -> horizontal hold-up for Distributed pattern

% (Sub5 -1) Calculating liquid hold-up for different flow patterns
HL.Est=HL0.Est.*qsi.Est;
HL.Int=HL0.Int.*qsi.Int;
HL.Dis=HL0.Dis.*qsi.Dis;

% (Sub5 -2) Interpolating liquid hold-up for transition
fi = (L3-NFR)./(L3-L2);
HL.Tra = fi.*HL.Est + (1-fi).*HL.Int;

% (Sub5 -3) Calculating liquid hold-up
HL = HL.Est.*fl.Est + HL.Int.*fl.Int + HL.Dis.*fl.Dis + HL.Tra.*fl.Tra;

% (Sub5 -4) Filtering
HL = max(0,HL);
HL = min(1,HL);

% (Sub5 -5) Calculating gas hold-up
HG = 1-HL;