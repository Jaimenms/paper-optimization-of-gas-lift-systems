function stop = stopfunction(x,optimValues,state, parametersOpt, parameters)
% OPTIMPLOTRESNORM Plot value of the norm of residuals at each iteration.
%
%   STOP = OPTIMPLOTRESNORM(X,OPTIMVALUES,STATE) plots OPTIMVALUES.resnorm.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTRESNORM as the plot
%   function
%     options = optimset('PlotFcns',@optimplot);
%
%   Pass the options into an optimization problem to view the plot
%     lsqnonlin(@(x) sin(3*x),[1 4],[],[],options);

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:46:32 $

global STATUS_OUTPUT_BUTTON


stop = STATUS_OUTPUT_BUTTON;