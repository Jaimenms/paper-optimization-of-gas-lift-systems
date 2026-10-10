%% Gerando Primeira Figura

arquivos = {'result_opt_subsea_v41','result_opt_subsea_v31','result_opt_subsea_v21','result_opt_subsea_v01'};

for k=1:length(arquivos)
    
    clear SAIDAS SAIDAS2 X
    load(arquivos{k})

    names = fieldnames(SAIDAS{1});

    % Vetorizandos saídas da otimização
    SAIDAS2 = SAIDAS{1};
    for i=2:length(SAIDAS)
        for j=1:length(names)
            SAIDAS2.(names{j}) = [SAIDAS2.(names{j}) SAIDAS{i}.(names{j})];
        end
    end

    % Determinando melhor iteração
    [minfob,ix_minfob] = min(SAIDAS2.Fob);

    Pgaslift  = X(parameters.ix_Pgaslift,:).*parameters.PR./1e5;
    Pwell     = X(parameters.ix_Pwell,:).*parameters.PR./1e5;
    qLwell    = X(parameters.ix_qLwell,:);
    qGgaslift = X(parameters.ix_qGgaslift,:);

    QG(k) = qGgaslift(end);
    QL(k) = qLwell(end);
    PW(k) = Pwell(end);
    FG(k) = SAIDAS2.OutValue(end);
    PR(k) = SAIDAS2.Profit(end);
    CP(k) = SAIDAS2.CompressorPower(end);
    PG(k) = Pgaslift(end);
    
    if k==4
        parametersOpt.gasPrice = 0;
    end
    RATIO(k) = parametersOpt.gasPrice./parametersOpt.oilPrice;
    GasPrice(k)=parametersOpt.gasPrice;
    OilPrice(k)=parametersOpt.oilPrice;
    
end

disp(['Vazão de Óleo (kg/s) :  ' num2str(QL)])
disp(['Vazão de Gás  (kg/s) :  ' num2str(QG)])
disp(['RATIO ((USD/barrel)/(USD/Nm³)) :  ' num2str(RATIO)])
disp(['Well Pressure (bar):  '  num2str(PW)])
disp(['Custo Gás Comb (10^6 USD/year) :  ' num2str(FG)])

figure(1)
subplot(3,1,1)
plot(RATIO,QL,'-ok','MarkerSize',6,'linewidth',2);
xlabel('Cost Ratio ((USD/Nm³)/(USD/barrel))','Fontsize',12)
ylabel('Produced Oil (kg/s)','Fontsize',12)
title('Maximum profit operational conditions','Fontsize',14)
subplot(3,1,2)
plot(RATIO,QG,'-ok','MarkerSize',6,'linewidth',2);
xlabel('Cost Ratio ((USD/Nm³)/(USD/barrel))','Fontsize',12)
ylabel('Injected Gas (kg/s)','Fontsize',12)

subplot(3,1,3)
plot(RATIO,PR,'-ok','MarkerSize',6,'linewidth',2);
xlabel('Cost Ratio ((USD/Nm³)/(USD/barrel))','Fontsize',12)
ylabel('Profit (million USD/year)','Fontsize',12)


