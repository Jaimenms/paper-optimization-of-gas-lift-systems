function [F,J,outputs]=nonrigorous_model(x,parameters)
% Function that determines the residuos of the system of algebraic
% non-linear equations for a two-phase flow pipe network model based on
% Beggs and Brill correlation (1973) via non-rigorous momentum balance.
% In this approach the pressure gradient is calculated using the mean
% pressure of each pipe, with geometry equal to the geometry at the half
% position (at length equal to 50% of the total length of the pipe).
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

rigor = 0;
[F,J,outputs]=model_network(x,parameters,rigor);