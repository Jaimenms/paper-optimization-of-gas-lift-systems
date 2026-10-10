function [parameters] = initial_calculations(parameters,espec,tiniL)
% Function that prepares the data and organizing this at the model 
% parameters struct. This function is also able to include the wells at the
% pipe network geometry and to generate the specification matrix.
%
% The following inputs are necessary:
%
% parameters -> struc of parameters with the following fields (all units
%    according to SI):
%    .p -> matrix of connectivity, where each line corresponds to a network
%    pipe. Has two column: the first indicate the vertice of origin and the
%    second the vertice of destination.
%    .x -> column vector containing the vertice x-coordinate
%    .y -> column vector containing the vertice y-coordinate
%    .z -> column vector containing the vertice z-coordinate
%    .pipe_x{i} -> cell array of the i-th pipe with the column vector
%    of x-coordinates from the begining to the end of the pipe.
%    .pipe_y{i} -> cell array of the i-th pipe with the column vector
%    of y-coordinates from the begining to the end of the pipe.
%    .pipe_z{i} -> cell array of the i-th pipe with the column vector
%    of z-coordinates from the begining to the end of the pipe.
%    .pipe_D{i} -> cell array of the i-th pipe with the column vector
%    of diameters from the begining to the end of the pipe.
%    .pipe_rough{i} -> cell array of the i-th pipe with the column vector
%    of rougthness from the begining to the end of the pipe.
%    .pipe_Leq -> vector containing the equivalent length of accidents at
%    the pipes.
%
% espec -> 
%   .nodePressure -> column vector with pressure specifications at the
%   nodes. If the i-th node does not have a pressure specification the
%   i-th element vector shall receive nan.
%   .nodewL -> column vector with liquid mass flor rate specifications at
%   the nodes. If the i-th node does not have this specification the
%   i-th element vector shall receive nan. Positive value means that the
%   fluid is entering the net and negative the fluid is leaving the net.
%   .nodewG -> sabe of .nodewL but related to gas.
%
% The following inputs are also necessary if network has petroleum wells:
% parameters -> struc of parameters with the following fields (all units
% according to SI):
%   .ixwells -> vector of index of the node that is connected to the
%       well head
%   .wellhead_x -> vector of x-coordinate of the wellhead
%   .wellhead_y -> vector of y-coordinate of the wellhead
%   .wellhead_z -> vector of z-coordinate of the wellhead
%   .annulusDiameter -> vector of annulus diameter
%   .annulusRough -> vector of annulus roughness
%   .annulusLength -> vector of annulus length
%   .wellboreDiameter -> vetor of wellbore diameter
%   .wellboreRough -> vector of wellbore roughness
%   .wellboreLength -> vetor of wellbore length
%   .gasPipeDiameter ->  vector of gas pipe diamter
%   .gasPipeRough -> vetor of gas pipe roughness
%   .gasPipe1_x{1} -> cell array of vector of x-coodinates of the first gas
%   pipe.
%   .gasPipe1_y{1} -> cell array of vector of y-coodinates of the first gas
%   pipe.
%   .gasPipe1_z{1} -> cell array of vector of z-coodinates of the first gas
%   pipe.
%   .gasPipe2_x{1} -> cell array of vector of x-coodinates of a second gas
%   pipe.
%   .gasPipe2_y{1} -> cell array of vector of y-coodinates of a second gas
%   pipe.
%   .gasPipe2_z{1} -> cell array of vector of z-coodinates of a second gas
%   pipe.
%
% By Jaime Souza at April 3rd 2009.


% Crating incidence matrix
parameters = creating_M(parameters);

% Generating Specification Data
parameters = creating_E(parameters,espec);

% Adjusting Network for Including Petroleum Wells
parameters = creating_wells(parameters,espec);

% Treating pipes gemometry vectors
if ~exist('tiniL'),tiniL = 1;;end
parameters = interpolate_angle(parameters,tiniL);

% Calculating areas at the limits
parameters = creating_AreaAtLimits(parameters);

% Generating variable index
parameters = creating_index(parameters);

% Generating variables limits
parameters = creating_VariableLimits(parameters);


%% Subfunctions
function parameters = creating_VariableLimits(parameters)
parameters.x_LB(parameters.ixPI,1) = parameters.MinP./parameters.PR;
parameters.x_LB(parameters.ixPO,1) = parameters.MinP./parameters.PR;
parameters.x_LB(parameters.ixPW,1) = parameters.MinP./parameters.PR;
parameters.x_UB(parameters.ixPI,1) = parameters.MaxP./parameters.PR;
parameters.x_UB(parameters.ixPO,1) = parameters.MaxP./parameters.PR;
parameters.x_UB(parameters.ixPW,1) = parameters.MaxP./parameters.PR;
parameters.x_LB(parameters.ixqL,1) = parameters.MinqL;
parameters.x_UB(parameters.ixqL,1) = parameters.MaxqL;
parameters.x_LB(parameters.ixqG,1) = parameters.MinqG;
parameters.x_UB(parameters.ixqG,1) = parameters.MaxqG;
parameters.x_LB(parameters.ixwL,1) = -abs(parameters.MaxqL);
parameters.x_UB(parameters.ixwL,1) = abs(parameters.MaxqL);
parameters.x_LB(parameters.ixwG,1) = -abs(parameters.MaxqG);
parameters.x_UB(parameters.ixwG,1) = abs(parameters.MaxqG);


function parameters = creating_index(parameters)
% Determining index vectors
N = parameters.N;
S = parameters.S;
parameters.ixqL = 0*N+0*S+[1:S];
parameters.ixqG = 0*N+1*S+[1:S];
parameters.ixPI = 0*N+2*S+[1:S];
parameters.ixPO = 0*N+3*S+[1:S];
parameters.ixPW = 0*N+4*S+[1:N];
parameters.ixwL = 1*N+4*S+[1:N];
parameters.ixwG = 2*N+4*S+[1:N];


function parameters = creating_AreaAtLimits(parameters)
% Obtaining diameters and cross section area at pipes extremities
[N,S] = size(parameters.M);
zrs = zeros(S,1);
DI = zrs;
DO = zrs;
for i=1:S
    auxL = parameters.pipe_L{i};
    DI(i,1)=interp1(auxL,parameters.pipe_D{i},auxL(1));
    DO(i,1)=interp1(auxL,parameters.pipe_D{i},auxL(end));
end
parameters.AO = pi.*DO.^2./4;
parameters.AI = pi.*DI.^2./4;



function parameters = creating_M(parameters);
% Function that creates the incidence matrix based on the connectivity
% matrix

% Creating incidence matrix
N0 = max(parameters.p(:));
S0 = size(parameters.p,1);
M = zeros(N0,S0);
for i=1:S0,
    M(parameters.p(i,:),i)=[-1;1];
end
[N,S] = size(M);
M1 = M<0;
M2 = M>0;
parameters.M = M;
parameters.M1 = M1;
parameters.M2 = M2;
parameters.N = N;
parameters.S = S;

function parameters = creating_E(parameters,espec)
% Function that determines the specification vector and matrix

% Sizing the network
[N,S] = size(parameters.M);

% Alocating
E = zeros(N+S,4*S+3*N);
e = zeros(N+S,1);

% Creating E and e
aux = 1;
for i=1:N
    if ~isnan(espec.nodePressure(i))
        E(aux,4*S+i) = 1;
        e(aux,1) = espec.nodePressure(i);
        aux = aux+1;
    end
    if ~isnan(espec.nodewL(i))
        E(aux,4*S+N+i) = 1;
        e(aux,1) = espec.nodewL(i);
        aux = aux+1;
    end
    if ~isnan(espec.nodewG(i))
        E(aux,4*S+2*N+i) = 1;
        e(aux,1) = espec.nodewG(i);
        aux = aux+1;
    end
end

% Packing E and e
parameters.E = E;
parameters.e = e;




function parameters = creating_wells(parameters,espec)
% Function that modifies all the network parameters to consider the
% petroleum wells.

if isfield(parameters,'ixwells')
    if ~isempty(parameters.ixwells)
        % Unpacking parameters
        p = parameters.p;
        ixwells = parameters.ixwells;
        wellhead_x = parameters.wellhead_x;
        gasPipe1_x = parameters.gasPipe1_x;
        gasPipe2_x=parameters.gasPipe2_x;
        wellhead_y = parameters.wellhead_y;
        gasPipe1_y = parameters.gasPipe1_y;
        gasPipe2_y = parameters.gasPipe2_y;
        wellhead_z = parameters.wellhead_z;
        gasPipe1_z = parameters.gasPipe1_z;
        gasPipe2_z = parameters.gasPipe2_z;
        wellboreDiameter = parameters.wellboreDiameter;
        wellboreRough = parameters.wellboreRough;
        wellboreLength = parameters.wellboreLength;
        annulusLength = parameters.annulusLength;
        annulusDiameter = parameters.annulusDiameter;
        annulusRough = parameters.annulusRough;
        gasPipeDiameter = parameters.gasPipeDiameter;
        gasPipeRough = parameters.gasPipeRough;


        % Changing incidence matrix
        N0 = max(p(:));
        S0 = size(p,1);
        nwells = length(ixwells);
        N = N0 + 6*nwells;
        S = S0 + 6*nwells;
        for i=1:nwells
            p(S0+6*(i-1)+1,:) = [N0+6*(i-1)+1 ixwells(i)];
            p(S0+6*(i-1)+2,:) = [N0+6*(i-1)+2 N0+6*(i-1)+1];
            p(S0+6*(i-1)+3,:) = [N0+6*(i-1)+3 N0+6*(i-1)+2];
            p(S0+6*(i-1)+4,:) = [N0+6*(i-1)+4 N0+6*(i-1)+2];
            p(S0+6*(i-1)+5,:) = [N0+6*(i-1)+5 N0+6*(i-1)+4];
            p(S0+6*(i-1)+6,:) = [N0+6*(i-1)+6 N0+6*(i-1)+5];
        end
        M = zeros(N,S);for i=1:S, M(p(i,:),i)=[-1;1];end
        parameters.M                     = M;


        % Node Elements
        for i=1:nwells
            gasSupply1_xi = gasPipe1_x{1};
            gasSupply1_yi = gasPipe1_y{1};
            gasSupply1_zi = gasPipe1_z{1};

            gasSupply2_xi = gasPipe2_x{1};
            gasSupply2_yi = gasPipe2_y{1};
            gasSupply2_zi = gasPipe2_z{1};


            aux = wellhead_x([i; i; i; i]);
            parameters.x = [parameters.x; aux(:); gasSupply1_xi(1); gasSupply2_xi(1)];
            aux = wellhead_y([i; i; i; i]);
            parameters.y = [parameters.y; aux(:); gasSupply1_yi(1); gasSupply2_yi(1)];
            parameters.z = [parameters.z; wellhead_z(i);...
                                          wellhead_z(i)-annulusLength(i);...
                                          wellhead_z(i)-wellboreLength(i);...
                                          wellhead_z(i); gasSupply1_zi(1); gasSupply2_zi(1)];
        end

        % Pipe Elements
        m = S0;
        for i=1:nwells
            % Generating wellhead to production node
            parameters.pipe_x{m+6*(i-1)+1} = [wellhead_x(i) parameters.x(ixwells(i))];
            parameters.pipe_y{m+6*(i-1)+1} = [wellhead_y(i) parameters.y(ixwells(i))];
            parameters.pipe_z{m+6*(i-1)+1} = [wellhead_z(i) parameters.z(ixwells(i))];
            parameters.pipe_D{m+6*(i-1)+1} = [1 1].*wellboreDiameter(i);
            parameters.pipe_rough{m+6*(i-1)+1} = [1 1].*wellboreRough(i);
            parameters.pipe_Leq(m+6*(i-1)+1,1)  = 0;
            parameters.HL0(m+6*(i-1)+1,1)  = 0.5;

            % Generating wellbore - part 1
            parameters.pipe_x{m+6*(i-1)+2} = [wellhead_x(i) wellhead_x(i)];
            parameters.pipe_y{m+6*(i-1)+2} = [wellhead_y(i) wellhead_y(i)];
            parameters.pipe_z{m+6*(i-1)+2} = [wellhead_z(i)-annulusLength(i) wellhead_z(i)];
            parameters.pipe_D{m+6*(i-1)+2} = [1 1].*wellboreDiameter(i);
            parameters.pipe_rough{m+6*(i-1)+2} = [1 1].*wellboreRough(i);
            parameters.pipe_Leq(m+6*(i-1)+2,1)  = 0;
            parameters.HL0(m+6*(i-1)+2,1)  = 0.5;

            % Generating wellbore - part 2
            parameters.pipe_x{m+6*(i-1)+3} = [wellhead_x(i) wellhead_x(i)];
            parameters.pipe_y{m+6*(i-1)+3} = [wellhead_y(i) wellhead_y(i)];
            parameters.pipe_z{m+6*(i-1)+3} = [wellhead_z(i)-wellboreLength(i) wellhead_z(i)-annulusLength(i) ];
            parameters.pipe_D{m+6*(i-1)+3} = [1 1].*wellboreDiameter(i);
            parameters.pipe_rough{m+6*(i-1)+3} = [1 1].*wellboreRough(i);
            parameters.pipe_Leq(m+6*(i-1)+3,1)  = 0;
            parameters.HL0(m+6*(i-1)+3,1)  = 0.5;

            % Generating annulus
            boreArea = pi.*annulusDiameter(i).^2./4;
            boreS    = pi.*annulusDiameter(i);
            extArea  = pi.*wellboreDiameter(i).^2./4;
            extS     = pi.*wellboreDiameter(i);
            Deq = 4*(boreArea+extArea)./(boreS+extS);
            parameters.pipe_x{m+6*(i-1)+4} = [wellhead_x(i) wellhead_x(i)];
            parameters.pipe_y{m+6*(i-1)+4} = [wellhead_y(i) wellhead_y(i)];
            parameters.pipe_z{m+6*(i-1)+4} = [wellhead_z(i)-annulusLength(i) wellhead_z(i)];
            parameters.pipe_D{m+6*(i-1)+4} = [1 1].*Deq;
            parameters.pipe_rough{m+6*(i-1)+4} = [1 1].*annulusRough(i);
            parameters.pipe_Leq(m+6*(i-1)+4,1)  = 0;
            parameters.HL0(m+6*(i-1)+4,1)  = 0;

            % Generating gas pipe
            parameters.pipe_x{m+6*(i-1)+5} = [gasPipe1_x{i} wellhead_x(i)];
            parameters.pipe_y{m+6*(i-1)+5} = [gasPipe1_y{i} wellhead_y(i)];
            parameters.pipe_z{m+6*(i-1)+5} = [gasPipe1_z{i} wellhead_z(i)];
            parameters.pipe_D{m+6*(i-1)+5} = [ones(1,1+length(gasPipe1_x{i}))].*gasPipeDiameter(i);
            parameters.pipe_rough{m+6*(i-1)+5} = [ones(1,1+length(gasPipe1_x{i}))].*gasPipeRough(i);
            parameters.pipe_Leq(m+6*(i-1)+5,1)  = 0;
            parameters.HL0(m+6*(i-1)+5,1)  = 0;
            parameters.ix_GasPipes(2*(i-1)+1,1) = m+6*(i-1)+5;

            % Generating gas pipe
            parameters.pipe_x{m+6*(i-1)+6} = gasPipe2_x{i};
            parameters.pipe_y{m+6*(i-1)+6} = gasPipe2_y{i};
            parameters.pipe_z{m+6*(i-1)+6} = gasPipe2_z{i};
            parameters.pipe_D{m+6*(i-1)+6} = ones(1,length(gasPipe2_x{i})).*gasPipeDiameter(i);
            parameters.pipe_rough{m+6*(i-1)+6} = ones(1,length(gasPipe2_x{i})).*gasPipeRough(i);
            parameters.pipe_Leq(m+6*(i-1)+6,1)  = 0;
            parameters.HL0(m+6*(i-1)+6,1)  = 0;
            parameters.ix_GasPipes(2*(i-1)+2,1) = m+6*(i-1)+6;
        end

        % Modifying Specification
        E = zeros(N+S,4*S+3*N);
        e = zeros(N+S,1);

        E(1:(N0+S0),[(1:S0) S+(1:S0) 2*S+(1:S0) 3*S+(1:S0) 4*S+(1:N0) 4*S+N+(1:N0) 4*S+2*N+(1:N0)]) = parameters.E;
        e(1:(N0+S0),1) = parameters.e;

        ie_Pwell    = [];
        ie_qLwell    = [];
        ie_qGwell    = [];
        ie_Pgaslift    = [];
        ie_qLgaslift = [];
        ie_qGgaslift = [];
        ix_Pwell    = [];
        ix_qLwell    = [];
        ix_qGwell    = [];
        ix_Pgaslift    = [];
        ix_qLgaslift = [];
        ix_qGgaslift = [];

        aux = N0+S0+1;
        for i=1:nwells
            % Creating connection nodes
            E(aux,4*S+N+N0+6*(i-1)+1) = 1; % n node connection (Liquid)
            e(aux,1) = 0;
            aux = aux+1;
            E(aux,4*S+2*N+N0+6*(i-1)+1) = 1; % n node connection (Gas)
            e(aux,1) = 0;
            aux = aux+1;

            E(aux,4*S+N+N0+6*(i-1)+2) = 1; % n+1 node connection (Liquid)
            e(aux,1) = 0;
            aux = aux+1;
            E(aux,4*S+2*N+N0+6*(i-1)+2) = 1; % n+1 node connection (Gas)
            e(aux,1) = 0;
            aux = aux+1;

            E(aux,4*S+N+N0+6*(i-1)+4) = 1; % n+3 node connection (Liquid)
            e(aux,1) = 0;
            aux = aux+1;
            E(aux,4*S+2*N+N0+6*(i-1)+4) = 1; % n+3 node connection (Gas)
            e(aux,1) = 0;
            aux = aux+1;

            E(aux,4*S+N+N0+6*(i-1)+5) = 1; % n+5 node connection (Liquid)
            e(aux,1) = 0.;
            aux = aux+1;
            E(aux,4*S+2*N+N0+6*(i-1)+5) = 1; % n+5 node connection (Gas)
            e(aux,1) = 0.;
            aux = aux+1;

            % Creating source nodes
            if ~isnan(espec.downholewL(i))
                E(aux,4*S+N+N0+6*(i-1)+3) = 1; % n+2 node connection (Liquid)
                e(aux,1) = espec.downholewL(i);
                ie_qLwell = [ie_qLwell;aux];
                aux = aux+1;
            end
            
            if ~isnan(espec.downholewG(i))
                E(aux,4*S+2*N+N0+6*(i-1)+3) = 1; % n+2 node connection (Gas)
                e(aux,1) = espec.downholewG(i);
                ie_qGwell = [ie_qGwell;aux];
                aux = aux+1;
            end

            if ~isnan(espec.gasInletwL(i))
                E(aux,4*S+N+N0+6*(i-1)+6) = 1; % n+4 node connection (Liquid)
                e(aux,1) = espec.gasInletwL(i);
                ie_qLgaslift = [ie_qLgaslift;aux];
                aux = aux+1;
            end
            if ~isnan(espec.gasInletwL(i))
                E(aux,4*S+2*N+N0+6*(i-1)+6) = 1; % n+4 node connection (Gas)
                e(aux,1) = espec.gasInletwG(i);
                ie_qGgaslift = [ie_qGgaslift;aux];
                aux = aux+1;
            end

            % Creating fixed pressure nodes
            if ~isnan(espec.downholePressure(i))
                E(aux,4*S+N0+6*(i-1)+3) = 1; % n+2 node connection (Liquid)
                e(aux,1) = espec.downholePressure(i);
                ie_Pwell = [ie_Pwell;aux];
                aux = aux+1;
            end
            if ~isnan(espec.gasInletPressure(i))
                E(aux,4*S+N0+6*(i-1)+5) = 1; % n+2 node connection (Liquid)
                e(aux,1) = espec.gasInletPressure(i);
                ie_Pgaslift = [ie_Pgaslift;aux];
                aux = aux+1;
            end
            
            ix_Pwell = [ix_Pwell;4*S+N0+6*(i-1)+3];
            ix_Pgaslift = [ix_Pgaslift;4*S+N0+6*(i-1)+5];
            ix_qLwell = [ix_qLwell;4*S+N+N0+6*(i-1)+3];
            ix_qGwell = [ix_qGwell;4*S+2*N+N0+6*(i-1)+3];
            ix_qLgaslift = [ix_qLgaslift;4*S+N+N0+6*(i-1)+6];
            ix_qGgaslift = [ix_qGgaslift;4*S+2*N+N0+6*(i-1)+6];

        end

        % Packing E and e
        parameters.E = E;
        parameters.e = e;

        % Packing Decision Variable Indexes
        parameters.ie_Pwell    = ie_Pwell;
        parameters.ie_qLwell    = ie_qLwell;
        parameters.ie_qGwell    = ie_qGwell;
        parameters.ie_Pgaslift    = ie_Pgaslift;
        parameters.ie_qLgaslift = ie_qLgaslift;
        parameters.ie_qGgaslift = ie_qGgaslift;
        parameters.ix_Pwell    = ix_Pwell;
        parameters.ix_qLwell    = ix_qLwell;
        parameters.ix_qGwell    = ix_qGwell;
        parameters.ix_Pgaslift    = ix_Pgaslift;
        parameters.ix_qLgaslift = ix_qLgaslift ;
        parameters.ix_qGgaslift = ix_qGgaslift ;
        

        % Final calculations
        [N,S] = size(parameters.M);
        M1 = parameters.M<0;
        M2 = parameters.M>0;

        parameters.M = M;
        parameters.M1 = M1;
        parameters.M2 = M2;
        parameters.N = N;
        parameters.S = S;


    end
end

function [parameters] = interpolate_angle(parameters,tinyl)
% Function that prepared the geometry vectors of the pipes

% Determining the number of pipes
ni = length(parameters.pipe_x);


parameters.pipe_x0 = parameters.pipe_x;
parameters.pipe_y0 = parameters.pipe_y;
parameters.pipe_z0 = parameters.pipe_z;
parameters.pipe_D0 = parameters.pipe_D;
parameters.pipe_rough0 = parameters.pipe_rough;


for i=1:ni % For each pipe
    
    % Capturing coordinates
    x = parameters.pipe_x{i};
    y = parameters.pipe_y{i};
    z = parameters.pipe_z{i};
    D = parameters.pipe_D{i};
    rough = parameters.pipe_rough{i};

    % Calculating differences
    dx = (x(2:end)-x(1:end-1));
    dy = (y(2:end)-y(1:end-1));
    dz = (z(2:end)-z(1:end-1));
    mD = (D(2:end)+D(1:end-1))./2;
    mrough = (rough(2:end)+rough(1:end-1))./2;
    
    % Calculating lengths
    l = sqrt(dx.^2+dy.^2+dz.^2);
    
    % Calculating angles
    teta = asin(dz./l);

    % Preparing axial position vector
    l = [0 cumsum(l)];
    parameters.pipe_L0{i} = l;
    
    % Determining number of elements
    nl = length(x); % number of elements of data vector
    nl2 = (nl-1)*2; % number of elements of inclination vector
    
    % Creating special vectors
    ix = (1:nl)'*[1 1];ix = ix';ix = ix(:)';ix(2)= []; ix(end-1)=[];
    ix2= (1:(nl-1))'*[1 1];ix2 = ix2';ix2 = ix2(:)';
    ls = [ones(nl2/2,1) -ones(nl2/2,1)]; ls = ls'; ls = ls(:); ls(1)=0; ls(end)=0; ls = ls';

    parameters.pipe_L{i}       = l(ix)+ ls.*tinyl;
    parameters.pipe_teta{i}    = teta(ix2);
    parameters.pipe_D{i}       = mD(ix2);
    parameters.pipe_rough{i}   = mrough(ix2);
    parameters.pipe_length(i,1)     = max(parameters.pipe_L{i});
end

for i=1:ni
    aux = linspace(0,parameters.pipe_length(i,1),100);
    parameters.pipe_meanD(i,1) = mean(interp1(parameters.pipe_L{i},parameters.pipe_D{i},aux));
end