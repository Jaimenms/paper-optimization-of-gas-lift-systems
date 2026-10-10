function [c,ceq,GC,GCeq] = objective_nonlcon_v01(xopt, parametersOpt, parameters)
global x dF3dqG dF3dqL
persistent path iPipe_well


% Determining dimentions
Np  = length(parametersOpt.ix_PipeDiameter);
Neq = length(iPipe_well);
S0  = size(parameters.p,1);
[N,S] = size(parameters.M);
Nw  = length(parameters.ixwells);

% Searchingg paths
if isempty(path)
    path = search_path(parameters.M);
    iPipe_well = S0+(6:6:(6*Nw))-4;
end

% Colecting states
qLwell    = x(parameters.ix_qLwell);
pdownhole = x(parameters.ix_Pwell).*parameters.PR./1e5;
pwell     = parametersOpt.pwell;

% Colecting Production index
IP = parametersOpt.IP;

% Determining total derivatives
dF3dqL_tot = zeros(size(dF3dqL,1),1);
dF3dqG_tot = zeros(size(dF3dqG,1),1);
for i=1:S
    path_i = path{i};
    for j = 1:length(path_i)
        dF3dqL_tot(i,1)=dF3dqL_tot(i,1)+dF3dqL(path_i(j),path_i(j));
        dF3dqG_tot(i,1)=dF3dqG_tot(i,1)+dF3dqG(path_i(j),path_i(j));
    end
end
dpdownhole_dqLwell = dF3dqL_tot(iPipe_well,:);
dpdownhole_dqGwell = dF3dqG_tot(iPipe_well,:);


% Determining Jacobian terms
GCeq_dqLwell = -diag(dpdownhole_dqLwell) - (diag(1./IP))*eye(Nw);
GCeq_dqGwell = -diag(dpdownhole_dqGwell);
GCeq_dD      = zeros(Neq,Np); 

% Calculating function outputs
c =[];
ceq = pwell-pdownhole-diag(1./IP)*qLwell;
GC =[];
GCeq = [GCeq_dqLwell GCeq_dqGwell GCeq_dD];


function path = search_path(M)
[N,S] = size(M);
path = [];
for j=1:S
    iy = j;
    path{j}=j;
    while ~isempty(iy)
        ix = find(M(:,iy)==1);
        iy = find(M(ix,:)==-1);
        path{j}=[path{j} iy];
    end
end   