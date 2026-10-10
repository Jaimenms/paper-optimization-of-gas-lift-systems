% Case 4
% Objective: One Well with mass sources only
% Network Topology
%    1     (1)      2      (2)      3       (3)            4
%    O------------->O-------------->O--------------------->o
%    ^              ^     
%    |(4)           |(10)
%    |              |
%    X well 1/2     X well 3/4

%% (1) Fluid Specification

% Fluids Parameters
parameters.interfatialTension    = 0.071;  % SCALAR -> Pa.m
parameters.liquidViscosity       = 0.012;  % SCALAR -> Pa.s
parameters.gasViscosity          = 8.0e-5; % SCALAR -> Pa.s
parameters.gasPolytrophicCoef    = 1.25;   % SCALAR -> -
parameters.liquidDensity         = 850;    % SCALAR -> kg/m³
parameters.gasReferenceDensity   = 80;     % SCALAR -> kg/m³
parameters.referencePressure     = 1e7;    % SCALAR -> Pa

%% (2) Simulation specification
% Simulation Parameters
parameters.PR                    = 1e7;    % SCALAR -> Pa
parameters.MinP                  = 1e5;    % SCALAR -> Pa
parameters.MinqL                 = 0;    % SCALAR -> Pa
parameters.MinqG                 = 0;    % SCALAR -> Pa
parameters.MaxP                  = inf;    % SCALAR -> Pa
parameters.MaxqL                 = inf;    % SCALAR -> Pa
parameters.MaxqG                 = inf;    % SCALAR -> Pa

% Convergence parameters
parameters.MaxStep               = 0.5;   % SCALAR -> -
parameters.OdeAbsTol             = 1e-6;   % SCALAR -> -
parameters.OdeRelTol             = 1e-2;      % SCALAR -> -


parameters.HL0 = [0.5;0.5;0.5;0.5;0.5;1.0;0.;0.;0.;0.5;0.5;1.0;0.;0.;0.;0.5;0.5;1.0;0.;0.;0.;0.5;0.5;1.0;0.;0.;0.]; % Initial stimatives for simplified simulation

%% (3) Geometry (excluding wells)
% Determining non-well elements

% Net
parameters.p(1,:) = [1 2];          % VECTOR 1 x 2 ->[(from node) (to node)]
parameters.p(2,:) = [2 3];          % VECTOR 1 x 2 ->[(from node) (to node)]
parameters.p(3,:) = [3 4];          % VECTOR 1 x 2 ->[(from node) (to node)]

% Nodes
parameters.x = [    0;    0;    0;    0;];       % VECTOR Nnodes x 1 -> m -> coordinates
parameters.y = [    0; 1500; 4000; 4000;];       % VECTOR Nnodes x 1 -> m -> coordinates
parameters.z = [-1400;-1400;-1500;    0;];       % VECTOR Nnodes x 1 -> m -> coordinates

% Pipes
parameters.pipe_x{1} = [0 0];              % CELL ARRAY of VECTOR -> m -> coordinates
parameters.pipe_y{1} = [0 1500];           % CELL ARRAY of VECTOR -> m -> coordinates
parameters.pipe_z{1} = [-1400 -1400];      % CELL ARRAY of VECTOR -> m -> coordinates
parameters.pipe_D{1} = [6 6].*0.0254;      % CELL ARRAY of VECTOR -> m -> diameters
parameters.pipe_rough{1} = [100 100].*1e-6;% CELL ARRAY of VECTOR -> m -> roughness
parameters.pipe_Leq(1,1) = 0;              % CELL ARRAY of VECTOR -> m -> lenght of accidents

parameters.pipe_x{2} = [0 0 0];
parameters.pipe_y{2} = [1500 3000 4000];
parameters.pipe_z{2} = [-1400 -1400 -1500];
parameters.pipe_D{2} = [6 6 6].*0.0254;
parameters.pipe_rough{2} = 100.*[1 1 1].*1e-6;
parameters.pipe_Leq(2,1) = 0;

parameters.pipe_x{3} = [0 0];
parameters.pipe_y{3} = [4000 4000];
parameters.pipe_z{3} = [-1500 0];
parameters.pipe_D{3} = [6 6].*0.0254;
parameters.pipe_rough{3} = 100.*[1 1].*1e-6;
parameters.pipe_Leq(3,1) = 0;

%% (4) Geometry (well elements)

% index of node interconnected to wells
parameters.ixwells = [1 1 2 2];                       % VECTOR nwells x 1 -> (-)

% wellhead coordinates
parameters.wellhead_x = [  400 400 400 400];                % VECTOR Nwells x 1 -> m -> coordinates
parameters.wellhead_y = [ -300 300 1200 1800];                % VECTOR Nwells x 1 -> m -> coordinates
parameters.wellhead_z = [-1400 -1400 -1400 -1400];                % VECTOR Nwells x 1 -> m -> coordinates

% well geometry
parameters.annulusDiameter  = [10 10 10 10 ].*0.0254;     % VECTOR Nwells x 1 -> m
parameters.annulusRough     = [100 100 100 100 ].*1e-6;      % VECTOR Nwells x 1 -> m
parameters.annulusLength    = [1710 1710 1710 1710 ];           % VECTOR Nwells x 1 -> m %1710
parameters.wellboreDiameter = [4 4 4 4 ].*0.0254;      % VECTOR Nwells x 1 -> m
parameters.wellboreRough    = [100 100 100 100 ].*1e-6;      % VECTOR Nwells x 1 -> m
parameters.wellboreLength   = [2048 2048 2048 2048 ];           % VECTOR Nwells x 1 -> m %2048
parameters.gasPipeDiameter  = [2 2 2 2 ].*0.0254;      % VECTOR Nwells x 1 -> m
parameters.gasPipeRough     = 0.*[1 1 1 1].*1e-6;      % VECTOR Nwells x 1 -> m

parameters.gasPipe1_x{1} = [0 0 0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_y{1} = [4000 3000 1500 0];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_z{1} = [-1500 -1400 -1400 -1400]+5; % CELL ARRAY of VECTOR -> m -> coordinates

parameters.gasPipe2_x{1} = [0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_y{1} = [4000 4000];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_z{1} = [15 -1500]+5; % CELL ARRAY of VECTOR -> m -> coordinates

parameters.gasPipe1_x{2} = [0 0 0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_y{2} = [4000 3000 1500 0];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_z{2} = [-1500 -1400 -1400 -1400]+5; % CELL ARRAY of VECTOR -> m -> coordinates

parameters.gasPipe2_x{2} = [0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_y{2} = [4000 4000];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_z{2} = [15 -1500]+5; % CELL ARRAY of VECTOR -> m -> coordinates

parameters.gasPipe1_x{3} = [0 0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_y{3} = [4000 3000 1500];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_z{3} = [-1500 -1400 -1400]+5; % CELL ARRAY of VECTOR -> m -> coordinates

parameters.gasPipe2_x{3} = [0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_y{3} = [4000 4000];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_z{3} = [15 -1500]+5; % CELL ARRAY of VECTOR -> m -> coordinates

parameters.gasPipe1_x{4} = [0 0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_y{4} = [4000 3000 1500];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe1_z{4} = [-1500 -1400 -1400]+5; % CELL ARRAY of VECTOR -> m -> coordinates

parameters.gasPipe2_x{4} = [0 0];                    % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_y{4} = [4000 4000];        % CELL ARRAY of VECTOR -> m -> coordinates
parameters.gasPipe2_z{4} = [15 -1500]+5; % CELL ARRAY of VECTOR -> m -> coordinates

% Specifications
% value equal to nan means that this information is not specified

% NODAL
espec.nodePressure = [nan;nan;nan;0.15]; % VECTOR Nwells x 1 -> m
espec.nodewL       = [0;0;0;nan];      % VECTOR Nwells x 1 -> m
espec.nodewG       = [0;0;0;nan];      % VECTOR Nwells x 1 -> m

% WELL
espec.gasInletPressure = [nan nan nan nan ];      % VECTOR Nwells x 1 -> (-) Pa/PR
espec.gasInletwL       = 1e-6.*[1 1 1 1];     % VECTOR Nwells x 1 -> kg/s
espec.gasInletwG       = [qgas1 qgas2 qgas3 qgas4];      % VECTOR Nwells x 1 -> kg/s
espec.downholePressure = [nan nan nan nan ];      % VECTOR Nwells x 1 -> (-) Pa/PR
espec.downholewL  = [qoleo1 qoleo2 qoleo3 qoleo4];        % VECTOR Nwells x 1 -> kg/s
espec.downholewG  = 1e-6.*[1 1 1 1];      % VECTOR Nwells x 1 -> kg/s

%% (6) Initial calculations
[parameters] = initial_calculations(parameters,espec);

%% (7) Simulation

options = optimset('Display','iter','TolFun',1e-9,'MaxIter',15,'TolX',1e-20);

% Simulating
%[x1]=simplified_model(parameters);
%[x2,saidas,status]=simulate(@nonrigorous_model,x1,options,parameters);
%[   x2,saidas,status]=simulate(@rigorous_model,x2,options,parameters);