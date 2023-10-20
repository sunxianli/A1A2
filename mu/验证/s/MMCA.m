%��ʱ���A1A2
clc
clear
%��������
load('AActivity.mat')
load('BActivity.mat')
N=1000;
lamda1=0.2;
lamda2=0.2;
delta1=0.3;
delta2=0.3;
gama1=0;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
gama2=0.5;
ter=10;% mu��ȡֵ����
MMCA_rep=1;%�������
stp=50;%ʱ�䳤��
%�������
% PA2=zeros(MMCA_rep,N); 
% PA1=zeros(MMCA_rep,N);
% PA2_ave=zeros(1,N); 
% PA1_ave=zeros(1,N);
% beta_R=0.7;
ter2=50;
bata_c=zeros(1,ter);
  PG=zeros(MMCA_rep,N);
PG_ave=zeros(1,N); 
for xun = 1:ter
  mu=xun/ter
  for xun2= 1:ter2
      beta_R=xun2/200
         for rep = 1:MMCA_rep
             rep
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
%����������
            for t = 1:stp
                for i =1:N
                    for j =1:N
                        RA1(j,i)=1-AActivity(i,j)*(PA1G(1,j)+PA1S(1,j))*lamda1;
                        RA2(j,i)=1-AActivity(i,j)*PA2S(1,j)*lamda2;
                        QA1(j,i)=1-BActivity(i,j)*PA1G(1,j)*beta_R*gama1;
                        QA2(j,i)=1-BActivity(i,j)*PA1G(1,j)*beta_R*gama2;
                        QU(j,i)=1-BActivity(i,j)*PA1G(1,j)*beta_R;
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
%�ڲ�ͬ��mu�µ�ȡֵ
             PG(rep,:)=PA1G;

%             PA2(rep,:)=PA2S;
%             PA1(rep,:)=PA1S;
         end
%           PA2_ave=mean(PA2,1); 
%           PA1_ave=mean(PA1,1);
%ÿ��mu��ÿ�����µľ�ֵ�Լ��ܶ�
          PG_ave=mean(PG,1);
 rhoG=sum(PG_ave)/N;
%������ͼ
% h=zeros(N,N);
% for i=1:N
%     for j=1:N
% h(i,j)=(1-(1-gama1)*PA1_ave(1,j)-(1-gama2)*PA2_ave(1,j))* BActivity(i,j);
%     end
% end
% [m1,n1]=eig(h);
% tez=max(diag(n1));
% bata_c(xun)=mu/tez;
%�ж�G���ܶ��Ƿ����0
if  rhoG>0.001
    bata_c(xun)=beta_R;
    break;
end
  end%xun2��Ӧ��
end
hold on;
xzhou=(1:ter)/ter;
xlabel('\mu');
ylabel('proportion');
plot(xzhou,bata_c,'-o');
% plot(xzhou,rhoG','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
% plot(xzhou,rhoA1','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
% plot(xzhou,rhoA2,'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);


