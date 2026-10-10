function [F,J,outputs]=rigorous_model(x,parameters)
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
% J -> Jacobian Matrix
% outputs -> model output struct
%
% See also: MODEL_NETWORK, RIGOROUS_MODEL, SIMPLIFIED_MODEL
%
% By Jaime Souza at April 3rd 2009.

rigor = 1;
[F,J,outputs]=model_network(x,parameters,rigor);