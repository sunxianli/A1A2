%��ʱ���A1A2
clc
clear
%��������
load('AActivity.mat')
load('BActivity.mat')
N=1000;
delta1=0.8;
gama1=0.3;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
ter=50;% mu��ȡֵ����
MMCA_rep=1;%�������
stp=150;%ʱ�䳤��
%�������
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
%����������
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
%�ڲ�ͬ��mu�µ�ȡֵ
             PA1(rep,:)=PA1S;          
         end

%ÿ��mu��ÿ�����µľ�ֵ�Լ��ܶ�
          PA1_ave=mean(PA1,1);
       
         %������ͼ
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


