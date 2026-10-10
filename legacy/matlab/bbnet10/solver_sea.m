function [X,flag]=solver_sea(FUN,X1,options,parametros)
% Solver de sistema de equações algebricas
%
% Onde options contem os parametros do solver, FUN corresponde a função a
% ser resolvida e t,X1 e parametros são variáveis de entrada de FUN.
global X2 F2

% Tol -> Tolerância do método (norma do vetor de residuos)
Tol=options.Tol;
% Numero máximo de iterações
maxIter=options.maxIter;

% Calculando valores iniciais
[F0]=feval(FUN,X1,parametros);

% Partindo variaveis
cont=0;
flag=1;
T0=cputime;
Res=norm(F0);

disp(['it.:     ' num2str(cont) '    norma :  ' num2str(norm(F0))])
while Res>Tol

    if ~isreal(X1)
        disp('Crashed due to complex number')
        break
    end

    if cont==maxIter
        disp('Crashed due to excess of iterection')
        break
    end
    
    dF=jacobiana(options,F0,FUN,X1,parametros);
    passo=dF\F0;
    %[fator,Res] = fminbnd(@(fator) veros(fator,passo,FUN,X1,parametros),0,1);
    fator=1;
    X2=X1-fator.*passo;
    [F2]=feval(FUN,X2,parametros);
    Res =norm(F2);

    X1 = X2;
    F0 = F2;
    cont=cont+1;    
    
    disp(['it.:     ' num2str(cont) '    norma :  ' num2str(norm(F0))])
    

end
X=X1;

if sum(sum(isnan(X))),disp('Não Convergiu'), end


function [Res] = veros(fator,passo,FUN,X1,parametros)
global X2 F2

X2=X1-fator.*passo;   
[F2]=feval(FUN,X2,parametros);
Res =norm(F2);

function [dF]=jacobiana(options,F0,FUN,X1,parametros)
% Função para cálculo da matriz Jacobiana.

% Delta e Delta2 -> Parametros para calculo da Jacobiana numérica
Delta=options.Delta;
Delta2=options.Delta2;
% F0 -> vetor de resíduos no X0 médio da Jacobiana

dF=[];
for j=1:length(X1)
    X2=X1;
    X2(j)=(1+Delta)*X1(j)+Delta2;
    dF=[dF (feval(FUN,X2,parametros)-F0)./(X2(j)-X1(j))];
end