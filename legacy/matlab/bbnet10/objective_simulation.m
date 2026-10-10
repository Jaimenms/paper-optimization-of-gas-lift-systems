function [x2,saidas,status] = objective_simulation(xopt,x1,parametersOpt,parameters)

% Colecting Flow Rate
WellFlowRate = xopt(parametersOpt.ix_WellFlowRate,:);
GasFlowRate = xopt(parametersOpt.ix_GasFlowRate,:);

% Setting Flow Rate
parameters.e(parameters.ie_qLwell) = WellFlowRate;
parameters.e(parameters.ie_qGgaslift) = GasFlowRate;

% Selecting Diameters
PipeDiameter = xopt(parametersOpt.ix_PipeDiameter,:);

if ~isfield(parametersOpt,'iy_PipeDiameter')
    parametersOpt.iy_PipeDiameter=[];
end
for i=1:length(parametersOpt.ix_PipeDiameter)
    D = xopt(parametersOpt.ix_PipeDiameter(i)).*0.0254;
    aux = parametersOpt.iy_PipeDiameter{i};
    for j=1:length(aux)
        parameters.pipe_meanD(aux(j),1) = D;
        parameters.pipe_D{aux(j)} = ones(size(parameters.pipe_D{aux(j)})).*D;
    end
end

%% Part 1 -Initial Estimate
options = parametersOpt.options;
xin0=simplified_model(parameters);
if isempty(x1)||sum(isnan(x1)),     x1=xin0; end

%% Part 2 - Finalization

[xout,saidas,status]=simulate(@rigorous_model,x1,options,parameters);

if status.errorFlag
    [xout,saidas,status]=simulate(@nonrigorous_model,x1,options,parameters);
    [xout,saidas,status]=simulate(@rigorous_model,xout,options,parameters);
end

if status.errorFlag,
    x2 = x1;
else
    x2 = xout;
end