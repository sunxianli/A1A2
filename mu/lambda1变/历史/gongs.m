%时变的A1A2
clc
clear
%参数设置
N=1000;
lamda1=[0.3,0.6,0.9];
mun=length(lamda1);
lamda2=0.8;
delta1=0.5;
delta2=0.5;
gama1=2;%衰减系数，当其为0时，知道信息的个体是完全免疫的。在这里若等于1，
gama2=0.5;
ma=5;%AD网络每次连接的边数
mb=5;
ter=20;% mu的取值个数
stp=200;%时间长度
%定义矩阵
 beta_MMCA=zeros(mun,ter);
 AActivity=1:N;
 BActivity=1:N;
 PG=zeros(MMCA_rep,N); 
 PA1=zeros(MMCA_rep,N);
 PA2=zeros(MMCA_rep,N); 
PG_ave=zeros(1,N); 
PA1_ave=zeros(1,N); 
PA2_ave=zeros(1,N);       
ter2=300;%β的取值个数
 for l1=1:mun
for xun = 1:ter
   mu=xun/ter;
    for xun2= 1:ter2
      beta_R=xun2/700;
        for rep = 1:MMCA_rep
            PA2S=0.01*ones(1,N);
            PUS=0.98*ones(1,N);
            PA1G=0.01*ones(1,N);
            PA1S=0*ones(1,N);

            PUS_UPDATE=zeros(1,N);
            PA2S_UPDATE=zeros(1,N);
            PA1G_UPDATE=zeros(1,N);
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);
            rA2=zeros(1,N);
            qA1=zeros(1,N);
            qA2=zeros(1,N);
            qU=zeros(1,N);

            RA1=zeros(N,N);
            RA2=zeros(N,N);
            QA1=zeros(N,N);
            QA2=zeros(N,N);
            QU=zeros(N,N);
            %生成邻接矩阵
            temAct=rand(1,N);
            for i = 1:N
                AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
            end
            temAct=rand(1,N);
            for i = 1:N
                BActivity(i)= yita*minActivity*(1-(temAct(i))^(1/(Exponent_B-1)));
            end
%计算马氏链
            for t = 1:stp
                for i =1:N
                    for j =1:N
                       RA1(j,i)=1-(AActivity(i)+AActivity(j))*(PA1G(1,j)+PA1S(1,j))*lamda1(l1)*ma/N;
                        RA2(j,i)=1-(AActivity(i)+AActivity(j))*PA2S(1,j)*lamda2*ma/N;
                        QA1(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*gama1*mb/N;
                        QA2(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*gama2*mb/N;
                        QU(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*mb/N;
                    end
                    tempprodRA1=cumprod(RA1(:,i));
                    rA1(1,i)=tempprodRA1(N);
                    tempprodRA2=cumprod(RA2(:,i));
                    rA2(1,i)=tempprodRA2(N);
                    tempprodQA1=cumprod(QA1(:,i));
                    qA1(1,i)=tempprodQA1(N);
                    tempprodQA2=cumprod(QA2(:,i));
                    qA2(1,i)=tempprodQA2(N); 
                    tempprodQU=cumprod(QU(:,i));
                    qU(1,i)=tempprodQU(N);
 PUS_UPDATE(1,i)=PA1S(1,i)*delta1*qU(1,i)+PA2S(1,i)*delta2*qU(1,i)+PA1G(1,i)*delta1*mu+PUS(1,i)* rA1(1,i)*rA2(1,i)*qU(1,i);
 PA1G_UPDATE(1,i)=PA1S(1,i)*(delta1*(1-qU(1,i))+(1-delta1)*(1-qA1(1,i)))+PA2S(1,i)*(delta2*(1-qU(1,i))+(1-delta2)*(1-qA2(1,i)))...
                 +PUS(1,i)*(rA1(1,i)*rA2(1,i)*(1-qU(1,i))+(1-rA1(1,i)-(1-rA1(1,i))*(1-rA2(1,i)))*(1-qA1(1,i))+(1-rA2(1,i))*(1-qA2(1,i)))...
                 +PA1G(1,i)*(1-mu);
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1)*qA1(1,i)+PA1G(1,i)*(1-delta1)*mu+PUS(1,i)*(1-rA1(1,i)-(1-rA1(1,i))*(1-rA2(1,i)))*qA1(1,i);  
 PA2S_UPDATE(1,i)=PA2S(1,i)*(1-delta2)*qA2(1,i)+PUS(1,i)*(1-rA2(1,i))*qA2(1,i); 

                end
                PUS=PUS_UPDATE;
                PA1G=PA1G_UPDATE;
                PA1S= PA1S_UPDATE;
                PA2S=PA2S_UPDATE;
            end
   PA1(rep,:)=PA1S+PA1G;
   PA2(rep,:)=PA2S;
   PG(rep,:)=PA1G;
        end
  PG_ave=mean(PG,1);
  PA1_ave=mean(PA1,1);
  PA2_ave=mean(PA2,1);
 rhoG=sum(PG_ave)/N;
 rhoA1=sum(PA1_ave)/N;
 rhoA2=sum(PA2_ave)/N;
if rhoG<0.001
    Theta_b2A1=0;
    Theta_b2A2=0;
    Theta_bA1=0;
    Theta_bA2=0;
    Ba2=0;
    for i =1:N
        Theta_bA1=Theta_bA1+PA1_ave(1,i)*BActivity(i);
        Theta_bA2=Theta_bA2+PA2_ave(1,i)*BActivity(i);
        Theta_b2A1=Theta_b2A1+PA1_ave(1,i)*(BActivity(i)^2);
        Theta_b2A2=Theta_b2A2+PA2_ave(1,i)*(BActivity(i)^2);
        Ba2=Ba2+BActivity(i)^2;
    end
    Theta_bA1=Theta_bA1/N;
    Theta_bA2=Theta_bA2/N;
    Theta_b2A2=Theta_b2A2/N;
    Theta_b2A1=Theta_b2A1/N;
    Ba=sum(BActivity)/N;
    Ba2=Ba2/N;
    x=(1-(1-gama1)* rhoA1-(1-gama2)* rhoA2)*(Ba2-(1-gama1)*Theta_b2A1-(1-gama2)*Theta_b2A2);
    EH=sqrt(x)+Ba-(1-gama1)*Theta_bA1-(1-gama2)*Theta_bA2;
    beta_c=mu/(mb*EH);%阈值的表达式，但是只有在β大于0.001时，他的上个值才是真的阈值
 else 
 beta_MMCA(l1,xun)= beta_c; 
    break;
end 
end
end
end
hold on;
 box on;
 grid off;
set(gca,'Fontsize',15);
xzhou=(1:ter)/ter;
plot(xzhou,beta_MMCA(1,:)','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,beta_MMCA(2,:)','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
plot(xzhou,beta_MMCA(3,:),'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\mu$','FontSize',15);ylabel('$\beta_c$','FontSize',15);   % 坐标轴解释
h=legend('$lamda_1$=0.3','$lamda_1$=0.6','$lamda_1$=0.9');
set(h,'Interpreter','latex','FontSize',15)%,





