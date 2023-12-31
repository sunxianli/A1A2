%定义参数
clc
clear
N=1000;
stp=30;
kb=1;
lamda1=0.6;
lamda2=0.6;
delta1=0.3;
delta2=0.3;
csmd_G=0.01;
csmd_A1=0.01;
csmd_A2=0.02;
gama1=1;
gama2=0.3;
mu=0.5;
yita=1000;
minActivity=0.001;
Exponent_A=2.5;
Exponent_B=2.5;
ma=5;
mb=5;
%开始迭代
termi=20;%β的取值个数2
MC_rep=20;%仿真次数
BETA_A1_NEW=zeros(MC_rep,termi);
BETA_A2_NEW=zeros(MC_rep,termi);
BETA_G_NEW=zeros(MC_rep,termi);
neighbora=1:N;%产生随机重联矩阵的存放地
neighborb=1:N;

BETA_A1_AVER =1:termi;%均值矩阵
BETA_A2_AVER =1:termi;
BETA_G_AVER =1:termi;
%开始仿真循环
parfor rep = 1:MC_rep%仿真次数
    AActivity=1:N;%活跃度矩阵
BActivity=1:N;
  BETA_A1=1:termi;%每次仿真的结果
  BETA_A2=1:termi;
  BETA_G=1:termi;
  %开始不同的β值循环
  for l = 1:termi
   disp(['第',num2str(rep),'次仿真','β为',num2str(l)])
   beta_R=l/40;
   %定义初始值已经系统中的状态规则，不存在UG,A2G,当上层是U,下层必然是S;下层是G,上面必然是A1
      x=rand(1,N);
    for i = 1:N;
        if x(i)<csmd_G;
            x(i)=1;
        else
            x(i)=0;
        end
    end
    m=rand(1,N);
    for i = 1:N
      if m(i)<csmd_A1
        m(i)=1;
      elseif (m(i)>csmd_A1)&&(m(i)<csmd_A2)
        m(i)=2;
      else
          m(i)=0; 
      end
      if x(i)==1
        m(i)=1;
      end
       if m(i)==0
        x(i)=0;
       end
    end
    %开始定义时变网络中每个节点的活跃程度
   temAct=rand(1,N);
   for i = 1:N
      AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
   end
   %开始定义时变网络中每个节点的活跃程度
   temAct2=rand(1,N);
   for i = 1:N 
%       BActivity(i)=1-AActivity(i);
       BActivity(i)= yita*minActivity*(1-(temAct2(i))^(1/(Exponent_B-1)));
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
%产生邻接矩阵
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
        if(temPro < BActivity(i));
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
    end
 %统计每层中的不同节点的数量，x层中只有0，1，所以很好计算；m中有0，1，2要复杂一点，先计算m中的节点数
 %for i=1:N
A1=length(find(m(:)==1)); %计算A1的数量
A2=length(find(m(:)==2)); %计算A2的数量             
    BETA_A1(l)=sum(A1)/N;%第l次仿真中，A1的密度

    BETA_A2(l)=sum(A2)/N;%第l次仿真中，A2的密度
    BETA_G(l)=sum(x)/N;%第l次仿真中，G的密度
  end
%储存每次仿真的结果
  BETA_A1_NEW(rep,:)=BETA_A1;
  BETA_A2_NEW(rep,:)=BETA_A2;
  BETA_G_NEW(rep,:)=BETA_G;
end
%计算多次仿真的平均值
for av = 1:termi
  BETA_A1_AVER(av)=sum(BETA_A1_NEW(:,av))/MC_rep;
  BETA_A2_AVER(av)=sum(BETA_A2_NEW(:,av))/MC_rep;
  BETA_G_AVER(av)=sum(BETA_G_NEW(:,av))/MC_rep;
end
xzhou=(1:termi)/40;
hold on;
 box on;
 grid off;
 set(gca,'Fontsize',15);
plot(xzhou,BETA_G_AVER','-o','color',[77/256 133/256 189/256]);
plot(xzhou,BETA_A1_AVER','-^','color',[247/256 144/256 61/256]);
plot(xzhou,BETA_A2_AVER,'-v','color',[89/256 169/256 90/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',15);ylabel('proportion','FontSize',15);   % 坐标轴解释
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,
