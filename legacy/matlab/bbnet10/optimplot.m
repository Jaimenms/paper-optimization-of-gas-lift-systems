function stop = optimplot(x,optimValues,state)
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

persistent plotavailable
global SAIDAS
global HAND

stop = false;
if ~isempty(SAIDAS)
    saidasOpt = SAIDAS;
end

if iter

switch state
     case 'init'
         plotavailable = true;
         HAND(101) = figure;
         set(gcf,'units','normalized',...
             'position',[0.2 0.1 0.55 0.75],'Color',[1 1 1],...
             'Menubar','none','NumberTitle','off','Name','Optimization');

         HAND(100) =  uicontrol(gcf,'style','togglebutton','string','STOP',...
        'FontSize',6,'units','normalized','position',[0.01 0.01 0.1 0.05],'fontsize',12);

    case 'iter'
        if plotavailable
            if optimValues.iteration == 0
                % The 'iter' case is  called during the zeroth iteration,
                % but it now has values that were empty during the 'init' case
                
            
                if ishandle(HAND(100)),stop = get(HAND(100),'value');end
                
                figure(HAND(101))
                
                HAND(1) = subplot(3,2,1);
                plot1 = plot(0,saidasOpt.Fob,'kd-', ...
                    'MarkerFaceColor',[1 0 1],'linewidth',2);
                xlabel('Iteration','interp','none','fontsize',8);
                set(plot1,'Tag','optimplot1');
                ylabel('FOB','interp','none','fontsize',10);
                title(['FOB Value: ',num2str(saidasOpt.Fob)],'interp','none','fontsize',10);
                
                HAND(2) = subplot(3,2,2);
                plot2 = plot(0,saidasOpt.Profit,'kd-', ...
                    'MarkerFaceColor',[1 0 1],'linewidth',2);
                xlabel('Iteration','interp','none','fontsize',8);
                set(plot2,'Tag','optimplot2');
                ylabel('10_6 USD','interp','none','fontsize',10);
                title(['Profit: ',num2str(saidasOpt.Profit)],'interp','none','fontsize',10);
                
                HAND(3) = subplot(3,2,3);
                plot3 = plot(0,saidasOpt.PressurePennalty,'kd-', ...
                    'MarkerFaceColor',[1 0 1],'linewidth',2);
                xlabel('Iteration','interp','none','fontsize',8);
                set(plot3,'Tag','optimplot3');
                ylabel('10_6 USD','interp','none','fontsize',10);
                title(['Pennalty: ',num2str(saidasOpt.PressurePennalty+saidasOpt.DivergencePennalty)],'interp','none','fontsize',10);
                hold on
                plot4 = plot(0,saidasOpt.DivergencePennalty,'kd-', ...
                    'MarkerFaceColor',[0 0 1],'linewidth',2);
                set(plot4,'Tag','optimplot4');
                legend('Pressure Pennalty','Divergence Pennalty','Fontsize',8)
                
                HAND(4) = subplot(3,2,4);
                plot5 = plot(0,saidasOpt.CompressorPower,'kd-', ...
                    'MarkerFaceColor',[1 0 1],'linewidth',2);
                xlabel('Iteration','interp','none','fontsize',8);
                set(plot5,'Tag','optimplot5');
                ylabel('kW','interp','none','fontsize',10);
                title(['Compressor Power: ',num2str(saidasOpt.CompressorPower)],'interp','none','fontsize',10);
                
                HAND(5) = subplot(3,2,5);
                plot6 = plot(0,saidasOpt.AnnualTotalCapitalCost,'kd-', ...
                    'MarkerFaceColor',[1 0 1],'linewidth',2);
                xlabel('Iteration','interp','none','fontsize',8);
                set(plot6,'Tag','optimplot6');
                ylabel('10_6 USD','interp','none','fontsize',10);
                title(['Capital Cost: ',num2str(saidasOpt.AnnualTotalCapitalCost)],'interp','none','fontsize',10);

                HAND(6) = subplot(3,2,6);
                plot7 = plot(0,saidasOpt.GasLiftFlowRate,'kd-', ...
                    'MarkerFaceColor',[1 0 1],'linewidth',2);
                xlabel('Iteration','interp','none','fontsize',8);
                set(plot7,'Tag','optimplot7');
                ylabel('10_6 m³/year','interp','none','fontsize',10);
                title(['Gas Consumption: ',num2str(saidasOpt.GasLiftFlowRate+saidasOpt.FuelGasFlowRate)],'interp','none','fontsize',10);
                hold on
                plot8 = plot(0,saidasOpt.FuelGasFlowRate,'kd-', ...
                    'MarkerFaceColor',[0 0 1],'linewidth',2);
                xlabel('Iteration','interp','none','fontsize',8);
                set(plot8,'Tag','optimplot8');

                set(HAND(1:6),'Fontsize',6)

                drawnow
            else
                
                aux1 = saidasOpt.PressurePennalty+saidasOpt.DivergencePennalty;
                aux2 = saidasOpt.GasLiftFlowRate+saidasOpt.FuelGasFlowRate;                
                
                actualize_optimplot(HAND(1),'optimplot1',  saidasOpt.Fob,'FOB Value: ',saidasOpt.Fob);
                actualize_optimplot(HAND(2),'optimplot2',  saidasOpt.Profit,'Profit: ',saidasOpt.Profit);
                actualize_optimplot(HAND(3),'optimplot3',  saidasOpt.PressurePennalty,'Pennalty: ',aux1);
                actualize_optimplot(HAND(3),'optimplot4',  saidasOpt.DivergencePennalty,'Pennalty: ',aux1);
                actualize_optimplot(HAND(4),'optimplot5',  saidasOpt.CompressorPower,'Compressor Power: ',saidasOpt.CompressorPower);
                actualize_optimplot(HAND(5),'optimplot6',  saidasOpt.AnnualTotalCapitalCost,'Capital Cost: ',saidasOpt.AnnualTotalCapitalCost);
                actualize_optimplot(HAND(6),'optimplot7',  saidasOpt.GasLiftFlowRate,'Gas Consumption: ',aux2);
                actualize_optimplot(HAND(6),'optimplot8',  saidasOpt.FuelGasFlowRate,'Gas Consumption: ',aux2);
                
                drawnow
            end
        end
end



function actualize_optimplot(h,tag,  value,title_aux,title_value)
plot1 = findobj('Tag',tag);
iter = length(get(plot1,'Xdata'))+1;
newX = [get(plot1,'Xdata') iter];
newY = [get(plot1,'Ydata') value];
set(plot1,'Xdata',newX, 'Ydata',newY);
set(get(h,'Title'),'String',[title_aux,num2str(title_value)]);

