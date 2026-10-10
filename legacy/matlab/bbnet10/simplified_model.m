function [x2]=simplified_model(parameters)
% Function that determines the residuos of the system of algebraic
% linear equations for a two-phase laminar flow pipe network model based
% on a very simplified laminar model. This model is applied for initial
% estimative of a rigorous model.
%
% Input Variables:
% parameters -> struct
%  
% Output Variables:
% x -> state variables vector


%% (1) Initial Calculations

% Unpacking parameters
M    = parameters.M;
E    = parameters.E;
e    = parameters.e;
PR   = parameters.PR;
L    = parameters.pipe_length;
viscL= parameters.liquidViscosity;
viscG= parameters.gasViscosity;
densL= parameters.liquidDensity;
densG= parameters.gasReferenceDensity;
z    = parameters.z;

if isfield(parameters,'HL0'),
    HL0 = parameters.HL0;
else
    HL0 = 0.5;
end

% Calculating dimensions
[N,S] = size(M);

% Calculating number of scenarios
nt = size(e,2);

%% (2) Evaluating pipe limits
% Obtaining cross section area at pipes extremities
for i=1:S
    D(i,1) = mean(parameters.pipe_D{i});
end


%% (3) Special incidence matrix
% Generating byproduct incidence matrix
M1 = M<0;
M2 = M>0;

%% (4) Model

% Considering Laminar Flow
FL=diag(128.*HL0.*viscL.*L./(D.^4.*pi))./densL./PR;
FG=diag(128.*(1-HL0).*viscG.*L./(D.^4.*pi))./densG./PR;

% Calculating Coefficient Matrix
 A = [[M zeros(N,S) zeros(N,S) zeros(N,S) zeros(N) eye(N) zeros(N)]; ... % Liquid mass balance
      [zeros(N,S) M zeros(N,S) zeros(N,S) zeros(N) zeros(N) eye(N)]; ... % Gas mass balance
      [-FL -FG eye(S) -eye(S) zeros(S,N) zeros(S,N) zeros(S,N)]; ... % Pipe Momentum Balance
      [zeros(S) zeros(S) zeros(S) -eye(S) M2' zeros(S,N) zeros(S,N)]; ... % Pipe exit balance
      [zeros(S) zeros(S) -eye(S) zeros(S) M1' zeros(S,N) zeros(S,N)]; ... % Pipe entrance balance
      E];
 
% Determining static head
g = 9.81;
d = (diag(densG.*g./PR.*(1-HL0))+diag(densL.*g./PR.*HL0))*(M'*z);

% Determining non-homogeneus terms
c = [zeros(N,nt);zeros(N,nt);d;zeros(2*S,nt);e];
 
% Solving the linear system
x2=A\c;
