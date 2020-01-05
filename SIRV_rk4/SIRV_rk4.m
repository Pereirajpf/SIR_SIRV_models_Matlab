%Trabalho modelacao
%Usando o metodo de Runge-Kuta de ordem 4 para o modelo SIRV
function SIRV_rk4(beta,gama,mi,ro,iniS,iniI,iniR,iniV,b)

a=0;            %Inicio no dia 0
h=0.01;         %step size   
n=(b-a)/h;      %numero de partições
x=a:h:b;        %Ponto de descritizacao
beta=beta/(iniS+iniI+iniR+iniV);  

%s: suscetiveis
ys=zeros(1,n+1);
ys(1)=iniS;
%i: infetados
yi=zeros(1,n+1);
yi(1)=iniI;
%r: recuperados
yr=zeros(1,n+1);
yr(1)=iniR;
%v: vacinados
yv=zeros(1,n+1);
yv(1)=iniV;

for i=1:n
%K1
s_k1=sirvS(ys(i),yi(i),beta,ro,yr(i),mi);
i_k1=sirvI(ys(i),yi(i),beta,gama);
r_k1=sirvR(yi(i),yr(i),gama,mi);
v_k1=sirvV(ys(i),ro);
%K2
s_k2=sirvS(ys(i)+(0.5*s_k1*h),yi(i)+(0.5*i_k1*h),beta,ro,yr(i)+(0.5*r_k1*h),mi);
i_k2=sirvI(ys(i)+(0.5*s_k1*h),yi(i)+(0.5*i_k1*h),beta,gama);
r_k2=sirvR(yi(i)+(0.5*i_k1*h),yr(i)+(0.5*r_k1*h),gama,mi);
v_k2=sirvV(ys(i)+(0.5*s_k1*h),ro);
%K3
s_k3=sirvS(ys(i)+(0.5*s_k2*h),yi(i)+(0.5*i_k2*h),beta,ro,yr(i)+(0.5*r_k2*h),mi);
i_k3=sirvI(ys(i)+(0.5*s_k2*h),yi(i)+(0.5*i_k2*h),beta,gama);
r_k3=sirvR(yi(i)+(0.5*i_k2*h),yr(i)+(0.5*r_k2*h),gama,mi);
v_k3=sirvV(ys(i)+(0.5*s_k2*h),ro);
%K4
s_k4=sirvS(ys(i)+(s_k3*h),yi(i)+(i_k3*h),beta,ro,yr(i)+(r_k3*h),mi);
i_k4=sirvI(ys(i)+(s_k3*h),yi(i)+(i_k3*h),beta,gama);
r_k4=sirvR(yi(i)+(i_k3*h),yr(i)+(r_k3*h),gama,mi);
v_k4=sirvV(ys(i)+(s_k3*h),ro);
%y(i+1)
ys(i+1)=ys(i)+((h/6)*(s_k1+(2*s_k2)+(2*s_k3)+s_k4));
yi(i+1)=yi(i)+((h/6)*(i_k1+(2*i_k2)+(2*i_k3)+i_k4));
yr(i+1)=yr(i)+((h/6)*(r_k1+(2*r_k2)+(2*r_k3)+r_k4));
yv(i+1)=yv(i)+((h/6)*(v_k1+(2*v_k2)+(2*v_k3)+v_k4));
end

%Maximo de infetados
fprintf('Maximo de infetados = %d, dia = %d\n',round(max(yi)),round(find(yi==max(yi))*0.01));

%Plot
plot(x,ys,'Color',[1 0.7 0],'DisplayName','S'); hold on;
plot(x,yi,'Color',[0 0.7 1],'DisplayName','I');
plot(x,yr,'Color',[1 0 0],'DisplayName','R'); 
plot(x,yv,'Color',[0 1 0],'DisplayName','V'); hold off;
xlabel('Tempo(dias)');
legend('show');
end