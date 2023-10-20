%参数设置
N=1000;
lamda1=0.8;
lamda2=0.8;
delta1=0.5;
delta2=0.5;
mu=0.5;
gama1=2;%衰减系数，当其为0时，知道信息的个体是完全免疫的。在这里若等于1，
gama2=0.5;
ma=5;%AD网络每次连接的边数
mb=5;
MMCA_rep=1;%仿真次数
stp=100;%时间长度
%定义矩阵
 load('BActivity.mat')
 load('AActivity.mat')
rhoG=0.01*ones(MMCA_rep,stp);%节点G的密度
rhoA1=0.01*ones(MMCA_rep,stp);
rhoA2=0.01*ones(MMCA_rep,stp);
% rhoG=0.01*ones(MMCA_rep,stp);%节点G的密度
% rhoA1=0.51*ones(MMCA_rep,stp);
% rhoA2=0.2*ones(MMCA_rep,stp);

rhoG_AVE=zeros(1,stp);%节点G的密度
rhoA1_AVE=zeros(1,stp);
rhoA2_AVE=zeros(1,stp);
%%师兄推导
% x0=(1-minActivity^(Exponent_A-1))*rand(1,N);
% y0=(1-minActivity^(Exponent_B-1))*rand(1,N);
% AActivity=yita*minActivity*(1-x0).^(1/(1-Exponent_A));
% BActivity=yita*minActivity*(1-y0).^(1/(1-Exponent_B));
%%网上结论
%     temAct=rand(1,N);
%    for i = 1:N
% %       AActivity(i)=yita*power(((power(1,-Exponent_A+1)-power(minActivity,-Exponent_A+1))*temAct(i)+power(minActivity,-Exponent_A+1)),1/(-Exponent_A+1));
% end
% %    temAct2=rand(1,N);
%    for i = 1:N 
% %  BActivity(i)=yita*power(((power(1,-Exponent_B+1)-power(minActivity,-Exponent_B+1))*temAct2(i)+power(minActivity,-Exponent_B+1)),1/(-Exponent_B+1));
% end
% %%蔡老师结论
%     temAct=rand(1,N);
%    for i = 1:N
%       AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
%    end
%    temAct2=rand(1,N);
%    for i = 1:N 
%   BActivity(i)= yita*minActivity*(1-(temAct2(i))^(1/(Exponent_B-1)));
%    end
 beta_R=0.2875;
        for rep = 1:MMCA_rep
            PA2S=0.01*ones(1,N);
            PUS=0.98*ones(1,N);
            PA1G=0.01*ones(1,N);
            PA1S=0*ones(1,N);

%   PA2S=0.2*ones(1,N);
%             PUS=0.29*ones(1,N);
%             PA1G=0.01*ones(1,N);
%             PA1S=0.5*ones(1,N);
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
% %计算马氏链
            for t = 1:stp-1
                for i =1:N
                    for j =1:N
                        RA1(j,i)=1-(AActivity(i)+AActivity(j))*(PA1G(1,j)+PA1S(1,j))*lamda1*ma/N;
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
                PA2S=PA2S_UPDATE;%求A2的数量
            PG=PA1G;%求G的数量
            PA1=PA1S+PA1G;%求A1的数量
            PA2=PA2S;%求A2的数量
          rhoG(rep,t+1)=sum(PG)/N;
          rhoA1(rep,t+1)=sum(PA1)/N;
          rhoA2(rep,t+1)=sum(PA2)/N;
            end
        end
rhoG_AVE(1,:)=sum(rhoG,1)/MMCA_rep;%节点G的密度
rhoA1_AVE(1,:)=sum(rhoA1,1)/MMCA_rep;
rhoA2_AVE(1,:)=sum(rhoA2,1)/MMCA_rep;
xzhou=1:stp;
hold on;
 box on;
 grid off;
plot(xzhou,rhoG_AVE','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,rhoA1_AVE','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
plot(xzhou,rhoA2_AVE,'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
% save('AActivity.mat','AActivity')