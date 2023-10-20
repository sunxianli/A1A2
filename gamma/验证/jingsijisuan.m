%时变的A1A2
clc
clear
%参数设置
N=1000;
lamda1=0.8;
delta1=0.05;
mu=0.05;
gama1=0;%衰减系数，当其为0时，知道信息的个体是完全免疫的。在这里若等于1，
yita=10;%缩放因子
minActivity=0.001;%下界
Exponent_A=2.1;%AD网络活跃指数
%  Exponent_B=2.1;
 ma=4;%AD网络每次连接的边数
mb=4;
ter=40;% mu的取值个数
stp=150;%时间长度
%定义矩阵
 MMCA_rep=1;
beta=zeros(1,ter);%不同gamma下的阈值情况
 beta_MMCA=zeros(MMCA_rep,1);
   temAct=rand(1,N);
   temAct2=rand(1,N);
for xun = 1:ter
 Exponent_B=(xun+13)/7;

  %每个Exponent_B都会产生随机数，因此这里需要多次仿真
        for rep = 1:MMCA_rep
%产生活跃度矩阵，每一个gamma都会产生新的活跃度矩阵，并且每一次仿真都会生成随机数，所以应该在同一个gamma不同仿真下重置
AActivity=1:N;
BActivity=1:N;

   for i = 1:N
      AActivity(i)=yita*power(((power(1,-Exponent_A+1)-power(minActivity,-Exponent_A+1))*temAct(i)+power(minActivity,-Exponent_A+1)),1/(-Exponent_A+1));
end
   
   for i = 1:N 
 BActivity(i)=yita*power(((power(1,-Exponent_B+1)-power(minActivity,-Exponent_B+1))*temAct2(i)+power(minActivity,-Exponent_B+1)),1/(-Exponent_B+1));
end
   %设置初始值，开始迭代
            PA1S=0.5*ones(1,N);
            
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);

            RA1=zeros(N,N);      
%计算马氏链
            for t = 1:stp
                for i =1:N
                    for j =1:N
                        RA1(j,i)=1-(AActivity(i)+AActivity(j))*PA1S(1,j)*lamda1*ma/N;
                    end
                    tempprodRA1=cumprod(RA1(:,i));
                    rA1(1,i)=tempprodRA1(N);
          
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1)+(1-PA1S(1,i))*(1-rA1(1,i));  
        end
                PA1S= PA1S_UPDATE;
            end
%记录每一次仿真的取值
% PA1(rep,:)=PA1S;  
PA1_new=PA1S; 
rhoA1=sum(PA1_new)/N;

%再每一次仿真下计算阈值
    Theta_b2A1=0;
    Theta_bA1=0;
    Ba2=0;
    for i =1:N
        Theta_bA1=Theta_bA1+PA1_new(i)*BActivity(i);
        Theta_b2A1=Theta_b2A1+PA1_new(i)*(BActivity(i)^2);
        Ba2=Ba2+BActivity(i)^2;
    end
    Theta_bA1=Theta_bA1/N;
    Theta_b2A1=Theta_b2A1/N;
    Ba=sum(BActivity)/N;
    Ba2=Ba2/N;
    x=(1-(1-gama1)* rhoA1)*(Ba2-(1-gama1)*Theta_b2A1);
    EH=sqrt(x)+Ba-(1-gama1)*Theta_bA1;
    beta_MMCA(rep,xun)=mu/mb/EH;        
        end
%同一个gamma下多次仿真取均值
end
beta=mean(beta_MMCA);
% rhoA1=sum(PA1_ave)/N;
xzhou=1:ter;
plot(xzhou,beta);
hold on;




