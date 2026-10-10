%% Gerando Primeira Figura

arquivos = {'result_opt_subsea_v01_5','result_opt_subsea_v01','result_opt_subsea_v01_2','result_opt_subsea_v01_3','result_opt_subsea_v01_4'};

for k=1:length(arquivos)
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
end

PI = [0.5 1 2 5 10];
disp(['Vazão de Óleo (kg/s) :  ' num2str(QL)])
disp(['Vazão de Gás  (kg/s) :  ' num2str(QG)])
disp(['Pressão poço  (bar)  :  ' num2str(PW)])
disp(['Índice ((kg/s)/bar)  :  ' num2str(PI)])


figure(1)
subplot(3,1,1)
plot(PI,QL,'-ok','MarkerSize',6,'linewidth',2);
xlabel('Productivity Index ((kg/s)/bar)','Fontsize',12)
ylabel('Produced Oil (kg/s)','Fontsize',12)
title('Maximum production operational conditions','Fontsize',14)

subplot(3,1,2)
plot(PI,QG,'-ok','MarkerSize',6,'linewidth',2);
xlabel('Productivity Index ((kg/s)/bar)','Fontsize',12)
ylabel('Injected Gas (kg/s)','Fontsize',12)

subplot(3,1,3)
plot(PI,PW,'-ok','MarkerSize',6,'linewidth',2);
xlabel('Productivity Index ((kg/s)/bar)','Fontsize',12)
ylabel('Well Pressure (bar)','Fontsize',12)


load results_Subsea_1_2


%% Montando Superfícies (Pressão no Poço e Pressão de Gás na Plataforma)
P=[];
for ii=1:nx
    for jj=1:ny
        x_aux = X{ii,jj};
        P(ii,jj) = x_aux(parameters.ix_Pwell).*parameters.PR./1e5;
    end
end
NY = ny;

%% Gerando Primeira Figura
figure(2)
colormap ([0 0 0])
[C,h] = contour(QGAS(:,1:NY),QOLEO(:,1:NY),P(:,1:NY));
set(h,'linewidth',2)
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])
xlabel('Injected Gas (kg/s)','Fontsize',12)
ylabel('Produced Oil (kg/s)','Fontsize',12)
title('Maximum production locus','Fontsize',14)

hold on
a=plot(QG,QL,'--or','MarkerSize',6,'linewidth',2);
legend([h;a],'Constant Well Pressure curves','Optimum locus')