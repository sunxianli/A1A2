%非时变的A1A2
clc
clear
%参数设置
load('AActivity.mat')
load('BActivity.mat')
N=1000;
delta1=0.8;
gama1=0.3;%衰减系数，当其为0时，知道信息的个体是完全免疫的。在这里若等于1，
ter=50;% mu的取值个数
MMCA_rep=1;%仿真次数
stp=150;%时间长度
%定义矩阵
mu=0.6;
bata_c=zeros(1,ter);
PA1=zeros(MMCA_rep,N);
PA1_ave=zeros(1,N); 
for xun = 1:ter
  lamda1=(xun-1)/ter;
         for rep = 1:MMCA_rep
            PA1S=0.02*ones(1,N);
            
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);
            
            RA1=zeros(N,N);
%计算马氏链
            for t = 1:stp
                for i =1:N
                    for j =1:N
                        RA1(j,i)=1-AActivity(i,j)*PA1S(1,j)*lamda1;
                    end
                    tempprodRA1=cumprod(RA1(:,i));
                    rA1(1,i)=tempprodRA1(N);
                   
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1)+(1-PA1S(1,i))*(1-rA1(1,i));  

                end
                PA1S= PA1S_UPDATE;
            end
%在不同的mu下的取值
             PA1(rep,:)=PA1S;          
         end

%每个mu下每个β下的均值以及密度
          PA1_ave=mean(PA1,1);
       
         %画函数图
     h=zeros(N,N);
for i=1:N
    for j=1:N
h(i,j)=(1-(1-gama1)*PA1_ave(1,j))*BActivity(i,j);
    end
end
[m1,n1]=eig(h);
tez=max(diag(n1));
zan=mu/tez;
    bata_c(xun)=zan;  
end
hold on;
yzhou=(1-1:ter-1)/ter;
xlabel('\mu');
ylabel('proportion');
plot(bata_c,yzhou);


