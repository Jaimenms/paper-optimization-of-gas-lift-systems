clear all

load results_Subsea_1_2


%% Montando Superfícies (Pressão no Poço e Pressão de Gás na Plataforma)
P=[];
for ii=1:nx
    for jj=1:ny
        x_aux = X{ii,jj};
        P(ii,jj) = x_aux(parameters.ix_Pwell).*parameters.PR./1e5;
        P2(ii,jj) = x_aux(parameters.ix_Pgaslift).*parameters.PR./1e5;
    end
end
NY = ny;

%% Gerando Primeira Figura
figure(1)
colormap ([0 0 0])
[C,h] = contour(QGAS(:,1:NY),QOLEO(:,1:NY),P(:,1:NY),[220 220]);
set(h,'linewidth',2)
text_handle = clabel(C,h,'labelspacing',300);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])
xlabel('Node 2 - Gas Mass Flow Rate (kg/s)','Fontsize',12)
ylabel('Node 1 - Liquid Mass Flow Rate (kg/s)','Fontsize',12)
title('Node 1 pressure (bar)','Fontsize',14)


load result_opt_subsea_v01

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

hold on
a = plot(qGgaslift,qLwell,'ok','MarkerSize',6);
b = plot(qGgaslift(ix_minfob),qLwell(ix_minfob),'r*','MarkerSize',10);
legend([h;a;b],'Simulation Results','Iteration values','Optimum point')

disp(['Vazão de Óleo (kg/s) :  ' num2str(qLwell(ix_minfob))])
disp(['Vazão de Gás  (kg/s) :  ' num2str(qGgaslift(ix_minfob))])
disp(['Pressão poço  (bar)  :  ' num2str(Pwell(ix_minfob))])
disp(['Numero de Iterações  :  ' num2str(length(SAIDAS))])