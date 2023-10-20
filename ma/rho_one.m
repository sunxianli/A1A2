clc
clear
%��������
N=1000;
lamda1=0.3;
lamda2=0.7;
delta1=0.5;
delta2=0.5;
mu=0.5;
gama1=2;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
gama2=0.5;
% ma=5;%AD����ÿ�����ӵı���
mb=5;
ter=20;% mu��ȡֵ����
stp=100;%ʱ�䳤��
mu=0.5;
yita=1;%��������
minActivity=0.1;%�½�
%�������
MMCA_rep=1;
% rhoG_num=zeros(num,ter);
  beta_R=0.15;

  %��������������󣬱������仯
    load AActivity ;
  load BActivity ;
rhoG=zeros(MMCA_rep,ter);%�ڵ�G���ܶ�
rhoA1=zeros(MMCA_rep,ter);
rhoA2=zeros(MMCA_rep,ter);
 rhoG_ave=zeros(1,ter);
    for rep = 1:MMCA_rep
for xun = 1:ter
ma=xun;
%ÿ��Exponent_B�����������������������Ҫ��η���
%������Ծ�Ⱦ���ÿһ��gamma��������µĻ�Ծ�Ⱦ��󣬲���ÿһ�η��涼�����������������Ӧ����ͬһ��gamma��ͬ����������

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

            PG=PA1G;
            PA1=PA1S+PA1G;
            PA2=PA2S;

          rhoG(rep,xun)=sum(PG)/N;
          rhoA1(rep,xun)=sum(PA1)/N;
          rhoA2(rep,xun)=sum(PA2)/N;
end
    end
rhoG_ave(1,:)=mean(rhoG,1); 

   
xzhou=((1:ter)+13)/7;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
plot(xzhou,rhoG_ave(1,:)','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
% plot(xzhou,rhoG_num(2,:)','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
% plot(xzhou,rhoG_num(3,:),'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
% set(gcf,'DefaultTextInterpreter','latex');
% xlabel('$\gamma_1$','FontSize',15);ylabel('$\beta_c$','FontSize',15);   % 鍧愭爣杞磋В閲�?
% h=legend('$\lambda_1=0.8,\lambda_2=0.5$','$$\lambda_1=0.5,\lambda_2=0.5$','$$\lambda_1=0.5,\lambda_2=0.8$');
% set(h,'Interpreter','latex','FontSize',15)%,
% save('rhoG_ave.mat','rhoG_ave')
% saveas(gcf, 'save.fig')
% 
% 
% 
% 
% 
% 