%ʱ���A1A2
clc
clear
%��������
N=1000;
lamda1=0.3;
lamda2=0.3;
delta1=0.3;
delta2=0.3;
gama1=2;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
gama2=0.5;
yita=1000;%��������
minActivity=0.001;%�½�
Exponent_A=2.5;%AD�����Ծָ��
Exponent_B=2.5;
ma=5;%AD����ÿ�����ӵı���
ter=20;% mu��ȡֵ����
MMCA_rep=1;%�������
stp=200;%ʱ�䳤��
%�������
%  beta_MMCA=zeros(mun,ter);
 bata_c=zeros(1,ter);
% ter2=50;
AActivity=1:N;
BActivity=1:N;
PA1=zeros(MMCA_rep,N); 
PA1_ave=zeros(1,N); 
PA2=zeros(MMCA_rep,N); 
PA2_ave=zeros(1,N);
   temAct=rand(1,N);
   mu=0.5;
    
for xun = 1:ter
   mb=xun;
   disp(['mb=' num2str(mb) ])
        for rep = 1:MMCA_rep
            PA2S=0.1*ones(1,N);
            PA1S=0.01*ones(1,N);

            PA2S_UPDATE=zeros(1,N);
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);
            rA2=zeros(1,N);

            RA1=zeros(N,N);
            RA2=zeros(N,N);

            %�����ڽӾ���
                 for i = 1:N
                AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
            end
            temAct=rand(1,N);
            for i = 1:N
                BActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_B-1)));
            end
%����������
            for t = 1:stp
                for i =1:N
                    for j =1:N 
                        RA1(j,i)=1-(AActivity(i)+AActivity(j))*PA1S(1,j)*lamda1*ma/N;
                        RA2(j,i)=1-(AActivity(i)+AActivity(j))*PA2S(1,j)*lamda2*ma/N;
                    end
                    tempprodRA1=cumprod(RA1(:,i));
                    rA1(1,i)=tempprodRA1(N);
                    tempprodRA2=cumprod(RA2(:,i));
                    rA2(1,i)=tempprodRA2(N);
                    
 PA1S_UPDATE(1,i)=PA1S(1,i)*(1-delta1)+(1-PA1S(1,i)-PA2S(1,i))*(rA2(1,i)-rA1(1,i)*rA2(1,i));  
 PA2S_UPDATE(1,i)=PA2S(1,i)*(1-delta2)+(1-PA1S(1,i)-PA2S(1,i))*(1-rA2(1,i)); 

                end
                PA1S= PA1S_UPDATE;
                PA2S=PA2S_UPDATE;
            end

%       
%����ֵ��ͼ��ʱ����ϵͳ������̬ʱ��G��ֵΪ0.����A1=A1S
  PA1(rep,:)=PA1S;
  PA2(rep,:)=PA2S;
        end
 PA1_ave=mean(PA1,1);
 PA2_ave=mean(PA2,1);

 rhoA1=sum(PA1_ave)/N;
 rhoA2=sum(PA2_ave)/N;%�����ĳ�����µľ�ֵ
 %����ֵ
 %if rhoG<0.001
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
%      zan=mu/(mb*EH);%��ֵ�ı��ʽ������ֻ���ڦ´���0.001ʱ�������ϸ�ֵ���������ֵ
     bata_c(xun)= mu/(mb*EH);  
end
hold on;
 box on;
 grid off;
set(gca,'Fontsize',15);
xzhou=1:ter;
plot(xzhou, bata_c);





