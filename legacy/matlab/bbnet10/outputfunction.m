function outputfunction(x, optimValues, parametersOpt, parameters)
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

global saidasOpt HAND STATUS_OUTPUT_BUTTON

if ~isempty(HAND)
    if ishandle(HAND(100))
        STATUS_OUTPUT_BUTTON = get(HAND(100),'value');
    end
else
    STATUS_OUTPUT_BUTTON = 0;
end

aux = findobj('Tag','majorFigure');

if isempty(HAND)||~ishandle(HAND(1))

    HAND(101) = figure;
    set(gcf,'units','normalized','Tag','majorFigure',...
        'position',[0.2 0.1 0.55 0.75],'Color',[1 1 1],...
        'Menubar','none','NumberTitle','off','Name','Optimization');

    HAND(100) =  uicontrol(gcf,'style','togglebutton','string','STOP',...
    'FontSize',6,'units','normalized','position',[0.01 0.01 0.1 0.05],'fontsize',12);
                
    HAND(1) = subplot(3,3,1);
    generate_optimplot('optimplot1','FOB Value: ',saidasOpt.Fob,'FOB',saidasOpt.Fob)

    HAND(2) = subplot(3,3,2);
    generate_optimplot('optimplot2','Profit: ',saidasOpt.Profit,'10_6 USD',saidasOpt.Profit)

    HAND(3) = subplot(3,3,3);
    generate_optimplot('optimplot3','Pennalty: ',saidasOpt.DivergencePennalty,'10_6 USD',aux)

    HAND(4) = subplot(3,3,4);
    generate_optimplot('optimplot5','Compressor Power: ',saidasOpt.CompressorPower,'kW',saidasOpt.CompressorPower)
    
    HAND(5) = subplot(3,3,5);
    generate_optimplot('optimplot6','Capital Cost: ',saidasOpt.AnnualTotalCapitalCost,'10_6 USD',saidasOpt.AnnualTotalCapitalCost)

    HAND(6) = subplot(3,3,6);
    aux = sum(saidasOpt.GasLiftFlowRate)+saidasOpt.FuelGasFlowRate;
    generate_optimplot('optimplot7','Gas Consumption: ',sum(saidasOpt.GasLiftFlowRate),'10_6 m³/year',aux)
    hold on
    plot8 = plot(0,saidasOpt.FuelGasFlowRate,'kd-', ...
        'MarkerFaceColor',[0 0 1],'linewidth',2);
    xlabel('Iteration','interp','none','fontsize',8);
    set(plot8,'Tag','optimplot8');
    a = legend('Gas Lift','Fuel Gas','Location','SouthOutside'); set(a,'Fontsize',8)
    
    HAND(7) = subplot(3,3,7);
    generate_optimplot('optimplot9','Oil Production: ',sum(saidasOpt.OilFlowRate),'10_6 barrel/year',sum(saidasOpt.OilFlowRate))

    HAND(8) = subplot(3,3,8);
    generate_optimplot('optimplot10','Model Residue: ',saidasOpt.Residue,'norm(x)',saidasOpt.Residue)

    set(HAND(1:8),'Fontsize',6)

    drawnow
else
                
    %figure(HAND(101))

    aux2 = sum(saidasOpt.GasLiftFlowRate)+saidasOpt.FuelGasFlowRate;                

    actualize_optimplot(HAND(1),'optimplot1',  saidasOpt.Fob,'FOB Value: ',saidasOpt.Fob);
    actualize_optimplot(HAND(2),'optimplot2',  saidasOpt.Profit,'Profit: ',saidasOpt.Profit);
    actualize_optimplot(HAND(3),'optimplot3',  saidasOpt.DivergencePennalty,'Pennalty: ',saidasOpt.DivergencePennalty);
    actualize_optimplot(HAND(4),'optimplot5',  saidasOpt.CompressorPower,'Compressor Power: ',saidasOpt.CompressorPower);
    actualize_optimplot(HAND(5),'optimplot6',  saidasOpt.AnnualTotalCapitalCost,'Capital Cost: ',saidasOpt.AnnualTotalCapitalCost);
    actualize_optimplot(HAND(6),'optimplot7',  sum(saidasOpt.GasLiftFlowRate),'Gas Consumption: ',aux2);
    actualize_optimplot(HAND(6),'optimplot8',  saidasOpt.FuelGasFlowRate,'Gas Consumption: ',aux2);
    actualize_optimplot(HAND(7),'optimplot9',  sum(saidasOpt.OilFlowRate),'Oil Production: ',sum(saidasOpt.OilFlowRate));
    actualize_optimplot(HAND(8),'optimplot10',  saidasOpt.Residue,'Model Residue: ',saidasOpt.Residue);

end

%figure(HAND(101))
drawnow

function actualize_optimplot(h,tag,value,title_aux,title_value)
plot1 = findobj('Tag',tag);

X = get(plot1,'Xdata');
Y = get(plot1,'Ydata');
if isempty(X), iter = 1;
else, iter = X(end)+1; end

newX = [X iter];
newY = [Y value];

nx = length(newX);
if nx>50
    newX = newX(nx-49:end);
    newY = newY(nx-49:end);
end    
set(plot1,'Xdata',newX, 'Ydata',newY);
set(get(h,'Title'),'String',[title_aux,num2str(title_value)]);



function generate_optimplot(tag,title_aux,value,label_aux,title_value)
    plot1 = plot(0,value,'kd-', ...
        'MarkerFaceColor',[1 0 1],'linewidth',2);
    xlabel('Iteration','interp','none','fontsize',8);
    set(plot1,'Tag',tag);
    ylabel(label_aux,'interp','none','fontsize',10);
    title([title_aux,num2str(title_value)],'interp','none','fontsize',10);