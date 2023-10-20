%ʱ���A1A2
clc
clear
%��������
N=500;
lamda1=0.2;
MMCA_rep=10;
delta1=0.05;
mu=0.05;
gama1=0;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
yita=1000;%��������
minActivity=0.001;%�½�
% Exponent_A=2.1;%AD�����Ծָ��
Exponent_B=3.1;
ma=4;%AD����ÿ�����ӵı���
mb=4;
ter=20;% mu��ȡֵ����
stp=250;%ʱ�䳤��
%�������

rhoG=zeros(1,ter);%�ڵ�G���ܶ�
rhoA1=zeros(1,ter);
 beta_R=0.2;
 beta_MMCA=zeros(1,ter);
parfor xun = 1:ter
    AActivity=1:N;
BActivity=1:N;
%  beta_R=xun/ter;
 Exponent_A=(xun+13)/7;
        for rep = 1:MMCA_rep
            PUS=0.4*ones(1,N);
            PA1G=0.1*ones(1,N);
            PA1S=0.5*ones(1,N);

            PUS_UPDATE=zeros(1,N);
            PA1G_UPDATE=zeros(1,N);
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);
            qA1=zeros(1,N);
            qU=zeros(1,N);

            RA1=zeros(N,N);
            QA1=zeros(N,N);
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
                        QA1(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*gama1*mb/N;
                        QU(j,i)=1-(BActivity(i)+BActivity(j))*PA1G(1,j)*beta_R*mb/N;
                    end
                    tempprodRA1=cumprod(RA1(:,i));
                    rA1(1,i)=tempprodRA1(N);
                    tempprodQA1=cumprod(QA1(:,i));
                    qA1(1,i)=tempprodQA1(N);
                    tempprodQU=cumprod(QU(:,i));
                    qU(1,i)=tempprodQU(N);
 PUS_UPDATE(1,i)=PA1S(1,i)*delta1*qU(1,i)+PA1G(1,i)*delta1*mu+PUS(1,i)*rA1(1,i)*qU(1,i);
 PA1G_UPDATE(1,i)=PA1S(1,i)*(delta1*(1-qU(1,i))+(1-delta1)*(1-qA1(1,i)))+PUS(1,i)*((1-rA1(1,i))*(1-qA1(1,i))+rA1(1,i)*(1-qU(1,i)))+PA1G(1,i)*(1-mu);
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1)*qA1(1,i)+PA1G(1,i)*(1-delta1)*mu+PUS(1,i)*(1-rA1(1,i))*qA1(1,i);  

                end
                PUS=PUS_UPDATE;
                PA1G=PA1G_UPDATE;
                PA1S= PA1S_UPDATE;
            end
%�ڲ�ͬ��mu�µ�ȡֵ
            PG=PA1G;
            PA1=PA1S+PA1G;


%           rhoG(xun)=sum(PG)/N;%��ͬmu�µ�ֵ
%           rhoA1(xun)=sum(PA1)/N;
          %��η���ȡֵ
         rhoG(xun)=rhoG(xun)+sum(PG)/N/MMCA_rep;
          rhoA1(xun)=rhoA1(xun)+sum(PA1)/N/MMCA_rep;
      
        end
%������ͼ
%     Theta_b2A1=0;
%     Theta_bA1=0;
%     Ba2=0;
%     for i =1:N
%         Theta_bA1=Theta_bA1+PA1(1,i)*BActivity(i);
%         Theta_b2A1=Theta_b2A1+PA1(1,i)*(BActivity(i)^2);
%         Ba2=Ba2+BActivity(i)^2;
%     end
%     Theta_bA1=Theta_bA1/N;
%     Theta_b2A1=Theta_b2A1/N;
%     Ba=sum(BActivity)/N;
%     Ba2=Ba2/N;
%     x=(1-(1-gama1)* rhoA1(xun))*(Ba2-(1-gama1)*Theta_b2A1);
%     EH=sqrt(x)+Ba-(1-gama1)*Theta_bA1;
%     beta_MMCA(xun)=mu/mb/EH;        
end

% xlabel('lamda');
% ylabel('\beta_c');
xzhou=(1+13:ter+13)/7;
plot(xzhou,rhoG);
hold on;



