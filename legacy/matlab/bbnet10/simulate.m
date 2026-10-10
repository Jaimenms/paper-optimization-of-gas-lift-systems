function [X,outputs,status]=simulate(FUN,X1,options,parameters)
% Simulate - Non-linear algebraic equation solver via Newton-Raphson.
% This algorithm is adusted for pipe flow simulations where unfeasible
% values of pressure are protected using a lower pass filter with set equal
% to the minimum possible pressure.
% 
% Input variables
%   FUN -> Function handle that generates the residue vector and the
%   Jacobian matrix. The sintaxe of this function is:
%   [F,J,outputs]=FUN(X1,parameters), where F is the residous vector,
%   J is the Jacobian Matrix and outputs is the function outuput struct.
%   X1  -> Vector of initial estimations of the state variables.
%   options -> struct containing convergence parameters:
%     .Tol -> norm convergence criteria limit;
%     .maxIter -> maximum number of iteractions.
%
% Output variables
%   X -> model solution vector
%   outputs -> strct containing model outputs
%   status -> struct containing the simulation results and status
%
% By Jaime Souza at April 3rd 2009.


%% (1) Initial calculations

% Defining index vectors of pressures at the state vector
global ixPI ixPO ixPW

% Defining persistent information
persistent start options2

% Starting Variables
cont=0;
dX = 0;
t=cputime;
errorFlag = 0;
convCriteria = 0;

% Stating Parameters
if isempty(start)
    options2 = options;
    start = 1;
    if isempty(options2)||isempty(options2.TolFun) , options2.TolFun = 1e-9; end
    if isempty(options2)||isempty(options2.MaxIter), options2.MaxIter = 15; end
    if isempty(options2)||isempty(options2.Display), options2.Display = 'iter'; end
    if isempty(options2)||isempty(options2.TolX)   , options2.TolX = 1e-6; end
end

% Unpacking parameters
Tol=options2.TolFun;
maxIter=options2.MaxIter;


%% (2) Initial residuous vector calculation
% Executing initial simulation
showDisplay= strcmp(options2.Display,'iter');
showEnd = strcmp(options2.Display,'iter')|strcmp(options2.Display,'final');
if showDisplay,     disp('Initiating stationary simulation.'), end
X = X1;
% ix = find(X(ixPI)<0);
% if ~isempty(ix), if showDisplay,     disp(['Error: Negative pressure at pipe inlet number: ' num2str(ix(:)')]),end
%     errorFlag = 4;
% end
% ix = find(X(ixPO)<0);
% if ~isempty(ix), if showDisplay,     disp(['Error: Negative pressure at pipe outlet number: ' num2str(ix(:)')]),end
%     errorFlag = 5;
% end
% ix = find(X(ixPW)<0);
% if ~isempty(ix), if showDisplay,     disp(['Error: Negative pressure at node: ' num2str(ix(:)')]),end
%     errorFlag = 6;
% end
% [F,dF,outputs]=feval(FUN,X,parameters);

% Calculating residue
%Res = norm(F);
Res = inf;
% if showDisplay,     
%     conv0 = max(0,min(100,round(log10(Res)./log10(Tol).*100)));
%     disp(['it.:  ' num2str(0) '  log(norm) :  ' num2str(round(log10(Res))) '  convergence : ' num2str(conv0) ])
% end

%% (3) Newton-Raphson

% Starting variables
while Res>Tol
    
    [F,dF,outputs]=feval(FUN,X,parameters);
    Res = norm(F);
    dX = dF\F;
    X=X-dX;

    
    % Error treatment
    if ~isreal(X),if showDisplay,     disp('Error: Crashed due to complex number!'),end
        errorFlag = -1;
        break
    end

    if cont==maxIter, if showDisplay,     disp('Error: Crashed due to excess of iterection!'),end
        errorFlag = 1;
        break
    end
    
    if sum(sum(isnan(X))), if showDisplay,     disp('Error: Nan reaults!'),end
        errorFlag = 2;
        break
    end
    
    if rank(dF)~=size(dF,1), if showDisplay,     disp('Error: Singular Jacobian!'),end
        errorFlag = 3;
        break
    end

    % Newton-Raphson step calculation
%    dX = dF\F;
        
    if norm(dX)<options2.TolX
        convCriteria = 2;
        break
    end
    
    % Generating new state vector
%    X=X-dX;
    
    % Forcing lower and upper bounds
    ix=X<parameters.x_LB;
    X(ix)=parameters.x_LB(ix);
    ix=X>parameters.x_UB;
    X(ix)=parameters.x_UB(ix);
    
    ix = find(X(ixPI)<0);
    if ~isempty(ix), if showDisplay,     disp(['Error: Negative pressure at pipe inlet number: ' num2str(ix(:)')]),end
        errorFlag = 4;
        break
    end
    ix = find(X(ixPO)<0);
    if ~isempty(ix), if showDisplay,     disp(['Error: Negative pressure at pipe outlet number: ' num2str(ix(:)')]),end
        errorFlag = 5;
        break
    end
    ix = find(X(ixPW)<0);
    if ~isempty(ix), if showDisplay,     disp(['Error: Negative pressure at node: ' num2str(ix(:)')]),end
        errorFlag = 6;
        break
    end
    
    % Recalculating residuous and Jacobian
%    [F,dF,outputs]=feval(FUN,X,parameters);
    
    % Calculating norm of residuous
%    Res =norm(F);
    cont=cont+1;    
    
    % Displaying iteractions
    
    if showDisplay,         
        conv0 = max(0,min(100,round(log10(Res)./log10(Tol).*100)));
        disp(['it.:  ' num2str(cont) '  log(norm) :  ' num2str(round(log10(Res))) '  convergence : ' num2str(conv0) ])
    end

end

if Res<Tol, convCriteria = 1; end

% Generating status struct
status.convCriteria = convCriteria;
status.errorFlag = errorFlag;
status.residue = Res;
status.time = cputime-t;
status.variables = length(X);
status.iteractions = cont;
status.direction = dX;
status.jacobian = dF;

% End
if showEnd,     
    switch convCriteria
        case 0
            disp('Finished due convergence error.')
        case 1
            disp('Converged due to minimum residue.')
        case 2
            disp('Converged due to minimum step.')
    end
    disp(['Number of Iteractions: ' num2str(status.iteractions)])
    disp(['Residue: ' num2str(status.residue)])
end