%Trabalho modelacao
%Usando o metodo de Euler explicito para o modelo SIR
function SIR_eulerExp(beta,gama,iniS,iniI,iniR,b,h)

a=0;            %Inicio no dia 0
%h=step size   
n=(b-a)/h;      %numero de partições
x=a:h:b;        %Ponto de descritizacao
beta=beta/(iniS+iniI+iniR);  

%s: suscetiveis
ys=zeros(1,n+1);
ys(1)=iniS;
%i: infetados
yi=zeros(1,n+1);
yi(1)=iniI;
%r: recuperados
yr=zeros(1,n+1);
yr(1)=iniR;

for(i=1:n)
    ys(i+1)=ys(i)+h*sirS(ys(i),yi(i),beta);
    yi(i+1)=yi(i)+h*sirI(ys(i),yi(i),beta,gama);
    yr(i+1)=yr(i)+h*sirR(yi(i),gama);
end

%Maximo de infetados
fprintf('Maximo de infetados = %d, dia = %d\n',round(max(yi)),round(find(yi==max(yi))*h));

%Plot
plot(x,ys,'Color',[1 0.7 0],'DisplayName','S'); hold on;
plot(x,yi,'Color',[0 0.7 1],'DisplayName','I');
plot(x,yr,'Color',[1 0 0],'DisplayName','R'); hold off;
xlabel('Tempo(dias)');
legend('show');
end