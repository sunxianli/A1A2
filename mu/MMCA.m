%ʱ���A1A2
clc
clear
%��������
N=1000;
lamda1=0.2;
lamda2=0.2;
delta1=0.3;
delta2=0.3;
% mu=0.5;
gama1=0;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
gama2=0.5;
yita=1000;%��������
minActivity=0.001;%�½�
Exponent_A=2.5;%AD�����Ծָ��
Exponent_B=2.5;
ma=5;%AD����ÿ�����ӵı���
mb=5;
ter=20;% mu��ȡֵ����
MMCA_rep=100;%�������
stp=100;%ʱ�䳤��
%�������
AActivity=1:N;
BActivity=1:N;
% rhoG=zeros(1,ter);%�ڵ�G���ܶ�
% rhoA1=zeros(1,ter);
% rhoA2=zeros(1,ter);
 beta_R=0.5;
 beta_MMCA=zeros(1,ter);
 PA2=zeros(MMCA_rep,N); 
PA1=zeros(MMCA_rep,N);
PA2_ave=zeros(1,N); 
PA1_ave=zeros(1,N);
for xun = 1:ter
   mu=xun/ter;
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
            %�����ڽӾ���
            temAct=rand(1,N);
            for i = 1:N
                AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
            end
            temAct=rand(1,N);
            for i = 1:N
                BActivity(i)= yita*minActivity*(1-(temAct(i))^(1/(Exponent_B-1)));
            end
%����������
            for t = 1:stp
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
                PA2S=PA2S_UPDATE;
            end
%�ڲ�ͬ��mu�µ�ȡֵ(���ǻ�������ͼ��
%             PG=PA1G;
%             PA1=PA1S+PA1G;
%����ֵ��ͼ��ʱ����ϵͳ������̬ʱ��G��ֵΪ0.����A1=A1S
 PA2(rep,:)=PA2S;
 PA1(rep,:)=PA1S;
        end
          PA2_ave=mean(PA2,1); 
          PA1_ave=mean(PA1,1);
%           rhoG(xun)=sum(PG)/N;%��ͬmu�µ�ֵ
          rhoA1=sum(PA1_ave)/N;
          rhoA2=sum(PA2_ave)/N;
%������ͼ
    Theta_b2A1=0;
    Theta_b2A2=0;
    Theta_bA1=0;
    Theta_bA2=0;
    Ba2=0;
    for i =1:N
        Theta_bA1=Theta_bA1+PA1(1,i)*BActivity(i);
        Theta_bA2=Theta_bA2+PA2(1,i)*BActivity(i);
        Theta_b2A1=Theta_b2A1+PA1(1,i)*(BActivity(i)^2);
        Theta_b2A2=Theta_b2A2+PA2(1,i)*(BActivity(i)^2);
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
    beta_MMCA(xun)=mu/mb/EH;        
end
hold on;
xzhou=(1:ter)/ter;
xlabel('mu');
ylabel('\beta_c');
plot(xzhou, beta_MMCA);




