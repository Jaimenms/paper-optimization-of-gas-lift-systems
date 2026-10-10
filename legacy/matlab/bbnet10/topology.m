function aux=topology(parameters,ix_verts,ix_well,ix_lines,ix_well_lines,seafloor,axis_value)

%figure

if ~exist('seafloor'),seafloor = -1500; end
if ~exist('ix_lines'),n = length(parameters.pipe_x0);ix_lines = 1:n; end
n = length(ix_lines);
for i=1:n
    a=plot3(parameters.pipe_x0{ix_lines(i)},...
            parameters.pipe_y0{ix_lines(i)},...
            parameters.pipe_z0{ix_lines(i)});
    set(a,'Linewidth',4)
    hold on
end

if ~exist('ix_well_lines'),ix_well_lines = []; end
n = length(ix_well_lines);
b=plot3(nan,nan,nan,'b-.');
set(b,'Linewidth',4)
hold on
for i=1:n
    b=plot3(parameters.pipe_x0{ix_well_lines(i)},...
            parameters.pipe_y0{ix_well_lines(i)},...
            parameters.pipe_z0{ix_well_lines(i)},'b-.');
    set(b,'Linewidth',4)
    hold on
end

if ~exist('ix_verts'),n = length(parameters.x); ix_verts = 1:n; end
if ~exist('ix_well'),ix_well = []; end
c = plot3(parameters.x(ix_verts),parameters.y(ix_verts),parameters.z(ix_verts),'ok');
set(c,'MarkerSize',10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
d = plot3(parameters.x(ix_well),parameters.y(ix_well),parameters.z(ix_well),'ok');
set(d,'MarkerSize',10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])


axis equal
if exist('axis_value'), axis(axis_value), end
grid on
aux = axis;
e(1) = patch([aux(1) aux(1) aux(2) aux(2)],[aux(3) aux(4) aux(4) aux(3)],[seafloor seafloor seafloor seafloor],[0.75    0.7    0.55]);
e(2) = patch([aux(1) aux(1) aux(1) aux(1)],[aux(3) aux(4) aux(4) aux(3)],[aux(5) aux(5) seafloor seafloor ],[0.75    0.7    0.55]);
e(3) = patch([aux(2) aux(2) aux(2) aux(2)],[aux(3) aux(4) aux(4) aux(3)],[aux(5) aux(5) seafloor seafloor],[0.75    0.7    0.55]);
e(4) = patch([aux(1) aux(1) aux(2) aux(2)],[aux(3) aux(3) aux(3) aux(3)],[seafloor aux(5) aux(5) seafloor],[0.75    0.7    0.55]);
e(5) = patch([aux(1) aux(1) aux(2) aux(2)],[aux(4) aux(4) aux(4) aux(4)],[seafloor aux(5) aux(5) seafloor],[0.75    0.7    0.55]);

e(6) = patch([aux(1) aux(1) aux(2) aux(2)],[aux(3) aux(4) aux(4) aux(3)],[0 0 0 0],[0 1 1]);
e(7) = patch([aux(1) aux(1) aux(1) aux(1)],[aux(3) aux(4) aux(4) aux(3)],[0 0 seafloor seafloor],[0 1 1]);
e(8) = patch([aux(1) aux(1) aux(2) aux(2)],[aux(4) aux(4) aux(4) aux(4)],[seafloor 0 0 seafloor],[0 1 1]);
e(9) = patch([aux(2) aux(2) aux(2) aux(2)],[aux(3) aux(4) aux(4) aux(3)],[0 0 seafloor seafloor],[0 1 1]);
e(10) = patch([aux(1) aux(1) aux(2) aux(2)],[aux(3) aux(3) aux(3) aux(3)],[seafloor 0 0 seafloor],[0 1 1]);

set(e,'FaceAlpha',0.25)

xlabel('x (m)','Fontsize',12)
ylabel('y (m)','Fontsize',12)
zlabel('z (m)','Fontsize',12)
%title('Oil/Gas network','Fontsize',14)
set(gca,'View',[-67.5 50])


d=plot3(nan,nan,nan,'r--');    set(d,'Linewidth',4)

h = legend([a;d;b;c],'Oil/Gas Lines','Gas Lines','Wellbores','Nodes','Location','North');
set(h,'Fontsize',10)
