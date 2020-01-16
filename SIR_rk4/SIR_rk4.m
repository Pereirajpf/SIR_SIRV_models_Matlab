%Trabalho modelacao
%Usando o metodo de Runge-Kuta de ordem 4 para o modelo SIR
function SIR_rk4(beta,gama,iniS,iniI,iniR,b,h)

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
%K1
s_k1=sirS(ys(i),yi(i),beta);
i_k1=sirI(ys(i),yi(i),beta,gama);
r_k1=sirR(yi(i),gama);
%K2
s_k2=sirS(ys(i)+(0.5*s_k1*h),yi(i)+(0.5*i_k1*h),beta);
i_k2=sirI(ys(i)+(0.5*s_k1*h),yi(i)+(0.5*i_k1*h),beta,gama);
r_k2=sirR(yi(i)+(0.5*i_k1*h),gama);
%K3
s_k3=sirS(ys(i)+(0.5*s_k2*h),yi(i)+(0.5*i_k2*h),beta);
i_k3=sirI(ys(i)+(0.5*s_k2*h),yi(i)+(0.5*i_k2*h),beta,gama);
r_k3=sirR(yi(i)+(0.5*i_k2*h),gama);
%K4
s_k4=sirS(ys(i)+(s_k3*h),yi(i)+(i_k3*h),beta);
i_k4=sirI(ys(i)+(s_k3*h),yi(i)+(i_k3*h),beta,gama);
r_k4=sirR(yi(i)+(i_k3*h),gama);
%y(i+1)
ys(i+1)=ys(i)+((h/6)*(s_k1+(2*s_k2)+(2*s_k3)+s_k4));
yi(i+1)=yi(i)+((h/6)*(i_k1+(2*i_k2)+(2*i_k3)+i_k4));
yr(i+1)=yr(i)+((h/6)*(r_k1+(2*r_k2)+(2*r_k3)+r_k4));
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