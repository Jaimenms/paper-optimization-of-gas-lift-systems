clear all
close all
clc
qgas = 0.5;
qoleo = 0.5;
case_03

[x2]=simplified_model(parameters);
%[x2,saidas,status]=simulate(@nonrigorous_model,x1,options,parameters);
[   x2,saidas,status]=simulate(@rigorous_model,x2,options,parameters);

vetor_qgas = [0.5:0.5:4.0];
vetor_qoleo = [0.5 1:1:20];

[QGAS,QOLEO] = meshgrid(vetor_qgas,vetor_qoleo);
[nx,ny] = size(QGAS);

% Executando Simulações
for jj=1:ny
    for ii=1:nx
        % Caminhando sem saltos
        if round(jj/2)==(jj/2)
            kk = nx-ii+1;
        else
            kk = ii;
        end
        
        % Atualizando vazões
        qgas = QGAS(kk,jj);
        qoleo= QOLEO(kk,jj);
        disp('__________________________________')
        disp(['Vazao de Gas (' num2str(jj) ') -> ' num2str(qgas)])
        disp(['Vazao de Oleo (' num2str(kk) ') -> ' num2str(qoleo)])
        parameters.e(parameters.ie_qLwell) = qoleo;
        parameters.e(parameters.ie_qGgaslift) = qgas;
        
        % Simulando rigorosamente
        [xout,saidas,status]=simulate(@rigorous_model,x2,options,parameters);
        if status.errorFlag==0
            x = xout;
            x2 = x;
        else
           [x1]=simplified_model(parameters);
           [xout,saidas,status]=simulate(@rigorous_model,x1,options,parameters);
        end
        if status.errorFlag==0
            x = xout;
            x2 = x;
        else
            x = nan.*xout;
            
        end
        
            
        % Guardando resultados
        SAIDAS{kk,jj} = saidas;
        X{kk,jj} = x;
        disp('__________________________________')

    end
end

% Salcando Resultados
save results_Subsea_1_2
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
[C,h] = contour(QGAS(:,1:NY),QOLEO(:,1:NY),P(:,1:NY));
set(h,'linewidth',2)
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])
xlabel('Node 2 - Gas Mass Flow Rate (kg/s)','Fontsize',12)
ylabel('Node 1 - Liquid Mass Flow Rate (kg/s)','Fontsize',12)
title('Node 1 pressure (bar)','Fontsize',14)

figure(2)
colormap ([0 0 0])
PI = QOLEO(:,1:NY)./(220-P(:,1:NY));
[C,h] = contour(QGAS(:,1:NY),QOLEO(:,1:NY),PI,[0.5 1 2 5 10]);
set(h,'linewidth',2)
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])
xlabel('Node 2 - Gas Mass Flow Rate (kg/s)','Fontsize',12)
ylabel('Node 1 - Liquid Mass Flow Rate (kg/s)','Fontsize',12)
title('Production Index (bar)','Fontsize',14)

figure(3)
colormap ([0 0 0])
[C,h] = contour(QGAS(:,1:NY),QOLEO(:,1:NY),P2(:,1:NY));
set(h,'linewidth',2)
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
    'Edgecolor',[.7 .7 .7])
xlabel('Node 2 - Gas Mass Flow Rate (kg/s)','Fontsize',12)
ylabel('Node 1 - Liquid Mass Flow Rate (kg/s)','Fontsize',12)
title('Gas Lift Inlet Pressure (bar)','Fontsize',14)