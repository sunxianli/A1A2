%ʱ���A1A2
clc
clear
%��������
N=1000;
lamda1=0.8;
delta1=0.05;
mu=0.05;
gama1=0;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
yita=10;%��������
minActivity=0.001;%�½�
Exponent_A=2.1;%AD�����Ծָ��
%  Exponent_B=2.1;
 ma=4;%AD����ÿ�����ӵı���
mb=4;
ter=40;% mu��ȡֵ����
stp=150;%ʱ�䳤��
%�������
 MMCA_rep=1;
beta=zeros(1,ter);%��ͬgamma�µ���ֵ���
 beta_MMCA=zeros(MMCA_rep,1);
   temAct=rand(1,N);
   temAct2=rand(1,N);
for xun = 1:ter
 Exponent_B=(xun+13)/7;

  %ÿ��Exponent_B�����������������������Ҫ��η���
        for rep = 1:MMCA_rep
%������Ծ�Ⱦ���ÿһ��gamma��������µĻ�Ծ�Ⱦ��󣬲���ÿһ�η��涼�����������������Ӧ����ͬһ��gamma��ͬ����������
AActivity=1:N;
BActivity=1:N;

   for i = 1:N
      AActivity(i)=yita*power(((power(1,-Exponent_A+1)-power(minActivity,-Exponent_A+1))*temAct(i)+power(minActivity,-Exponent_A+1)),1/(-Exponent_A+1));
end
   
   for i = 1:N 
 BActivity(i)=yita*power(((power(1,-Exponent_B+1)-power(minActivity,-Exponent_B+1))*temAct2(i)+power(minActivity,-Exponent_B+1)),1/(-Exponent_B+1));
end
   %���ó�ʼֵ����ʼ����
            PA1S=0.5*ones(1,N);
            
            PA1S_UPDATE=zeros(1,N);
           
            rA1=zeros(1,N);

            RA1=zeros(N,N);      
%����������
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
%��¼ÿһ�η����ȡֵ
% PA1(rep,:)=PA1S;  
PA1_new=PA1S; 
rhoA1=sum(PA1_new)/N;

%��ÿһ�η����¼�����ֵ
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
%ͬһ��gamma�¶�η���ȡ��ֵ
end
beta=mean(beta_MMCA);
% rhoA1=sum(PA1_ave)/N;
xzhou=1:ter;
plot(xzhou,beta);
hold on;




