%定义参数
clc
clear
% warning off
N=1000;
stp=300;
kb=1;
lamda1=0.8;
lamda2=0.8;
delta1=0.5;
delta2=0.5;
csmd_G=0.01*N;
csmd_A1=0.01*N;
csmd_A2=0.01*N;
% yita=10;
% minActivity=0.001;%下界
% Exponent_A=2.1;%AD网络活跃指数
% Exponent_B=2.1;
gama1=1.2;
gama2=0.8;
mu=0.5;
ma=5;
mb=5;
beta_R=0.2875;
%开始迭代
MC_rep=1;%仿真次数
BETA_A1_AVER=zeros(1,stp);
BETA_A2_AVER=zeros(1,stp);
BETA_U_AVER=zeros(1,stp);
BETA_G_AVER=zeros(1,stp);
neighbora=1:N;%产生随机重联矩阵的存放地
neighborb=1:N;
 load('BActivity.mat')
 load('AActivity.mat')
  BETA_A1=zeros(MC_rep,stp);%每次仿真的结果
  BETA_A2=zeros(MC_rep,stp);
  BETA_U=zeros(MC_rep,stp);
  BETA_G=zeros(MC_rep,stp);
 for rep = 1:MC_rep%仿真次数
   %定义初始值已经系统中的状态规则，不存在UG,A2G,当上层是U,下层必然是S;下层是G,上面必然是A1
    x=zeros(1,N);
m=zeros(1,N);
neighborA1=randi([1,N],1,2*csmd_G);
 for nA1 = 1:csmd_G;
          x(neighborA1(nA1))=1;%下层G对应上层A1
          m(neighborA1(nA1))=1;
 end
 for nA2 =csmd_G+1:2*csmd_G;%%挑选出上层的A2，下层自动对应0，即S
          m(neighborA1(nA2))=2;
 end
   %开始按照每个步长进行仿真，这里时间是stp=40
    for t = 1:stp
      n=1:N;
      n(:)=0;
      y=1:N;
      y(:)=0;
          mediA1=zeros(1,N);
          mediA2=zeros(1,N);
%m和x分别为上下层的状态矩阵，m由0，1，2组成；x由0，1组成；m和x对应的是t时刻，n和y对应t+1时刻的上下层
      A=zeros(N,N);%邻接矩阵A上层
      B=zeros(N,N);%下层邻接矩阵B
      for i = 1:N;
        temPro=rand(1,1);
        if temPro < AActivity(i);
            neighbora=randi([1,N],1,ma);
            for na = 1:ma;
                A(i,neighbora(na))=1;
                A(neighbora(na),i)=1;
            end
        end
        temPro=rand(1,1);
        if temPro < BActivity(i);
            neighborb=randi([1,N],1,mb);
            for nb = 1:mb
                B(i,neighborb(nb))=1;
                B(neighborb(nb),i)=1;
            end
        end
      end
%%开始进行上下层的耦合传播
     for i = 1:N
          %上层传播
        if m(i)==0 ;%(U,系统只中存在US，即此时x必定为0，US可变为A1S,A2S）
            for j = 1:N;
                 p1=rand(1,1);
                 p2=rand(1,1);
            if (A(j,i)==1)&&(m(j)==1)&&(p1<lamda1);
                   mediA1(i)=1;
              elseif (A(j,i)==1)&&(m(j)==2)&&(p2<lamda2);
                   mediA2(i)=2;
              end
            end
              if mediA2(i)==2%(mediA2(i)==2)&&(mediA1(i)==1);
                   n(i)=2;
              elseif  (mediA2(i)==0)&&(mediA1(i)==1);
                   n(i)=1;
              end  
        elseif m(i)==1;%A1，系统中存在A1S,A1G，A1S可变为US，A1G可变UG 
            if x(i)==0;%A1S,变成US
                p1=rand(1,1);
                if p1<delta1;
                    n(i)=0;%US
                else
                    n(i)=1;
                end
            else %x(i)==1，A1G变成A1G或者UG
                 p1=rand(1,1);
                 if p1<delta1;
                     n(i)=0;
                 else
                      n(i)=1;
                 end
            end
        else m(i)==2;%A2，只存在A2S,即此时x必然为0，A2S可变为US
            p2=rand(1,1);
            if p2<delta2;
              n(i)=0;  
            else
                n(i)=2;
            end
        end
         %下层传播
        if x(i)==0;%S状态,系统中存在US,A1S,A2S，A1S会变成A1G；US先变成UG,再变成A1G;A2S先变成A2G,再变成A1G
          if n(i)==1;%A1S，先处理复杂的A1S，再处理不会改变的
               for j = 1:N
                   p1=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p1<gama1*beta_R);
                    y(i)=1;
                end
               end
          elseif n(i)==0;%US变成A1G
              for j = 1:N
                   p0=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p0<beta_R);
                    y(i)=1;
                    n(i)=1;
                end
              end
          else n(i)==2;%A2S，先变成A2G,再变成A1G
                for j = 1:N
                   p2=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p2<gama2*beta_R);
                    y(i)=1;
                    n(i)=1;
                end
              end
          end
        else  x(i)==1;%G状态，系统中只存在A1G,还有中间状态UG,A1G可变为A1S,UG可变为US和A1G
           if n(i)==1;%A1G
           p1=rand(1,1);
              if p1<mu;
                y(i)=0;   %A1G免疫，变成A1S
            else
                y(i)=1;    
              end
           elseif n(i)==0; %UG
                  p0=rand(1,1);
                  if p0<mu;%US
                      y(i)=0; 
                  else
                       y(i)=1; %A1G
                       n(i)=1;
                  end
           end
        end
      end
      %更新m，x矩阵
      m=n;
      x=y;  
    A1=length(find(m(:)==1)); %计算A1的数量
    A2=length(find(m(:)==2)); %计算A2的数量        
    A3=length(find(m(:)==0)); %计算A2的数量    
   
   BETA_A1(rep,t)=sum(A1)/N;%第l次仿真中，A1的密度
    BETA_A2(rep,t)=sum(A2)/N;%第l次仿真中，A2的密度
    BETA_U(rep,t)=sum(A3)/N;%第l次仿真中，A1的密度
    BETA_G(rep,t)=sum(x)/N;%第l次仿真中，G的密度 
    end
 end
 %统计每层中的不同节点的数量，x层中只有0，1，所以很好计算；m中有0，1，2要复杂一点，先计算m中的节点数
% %计算多次仿真的平均值
for av = 1:stp;
  BETA_A1_AVER(av)=sum(BETA_A1(:,av))/MC_rep;
  BETA_A2_AVER(av)=sum(BETA_A2(:,av))/MC_rep;
  BETA_G_AVER(av)=sum(BETA_G(:,av))/MC_rep;
  BETA_U_AVER(av)=sum(BETA_U(:,av))/MC_rep;
end
% xzhou=(1:termi)/termi;
xzhou=0:stp-1;
hold on;
 box on;
 grid off;
 set(gca,'Fontsize',15);
 plot(xzhou,BETA_G_AVER','-o','color',[77/256 133/256 189/256]);
plot(xzhou,BETA_A1_AVER','-^','color',[247/256 144/256 61/256]);
plot(xzhou,BETA_A2_AVER,'-v','color',[89/256 169/256 90/256]);
% plot(xzhou,BETA_G_AVER,'b-.',xzhou,BETA_A1_AVER,'r-.',xzhou,BETA_A2_AVER,'g-.');
% set(gca,'DefaultAxesFontsize',5);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('t','FontSize',15);ylabel('proportion','FontSize',15);   % 坐标轴解释
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,

