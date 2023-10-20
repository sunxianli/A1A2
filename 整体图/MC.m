%�������
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
% minActivity=0.001;%�½�
% Exponent_A=2.1;%AD�����Ծָ��
% Exponent_B=2.1;
gama1=1.2;
gama2=0.8;
mu=0.5;
ma=5;
mb=5;
beta_R=0.2875;
%��ʼ����
MC_rep=1;%�������
BETA_A1_AVER=zeros(1,stp);
BETA_A2_AVER=zeros(1,stp);
BETA_U_AVER=zeros(1,stp);
BETA_G_AVER=zeros(1,stp);
neighbora=1:N;%���������������Ĵ�ŵ�
neighborb=1:N;
 load('BActivity.mat')
 load('AActivity.mat')
  BETA_A1=zeros(MC_rep,stp);%ÿ�η���Ľ��
  BETA_A2=zeros(MC_rep,stp);
  BETA_U=zeros(MC_rep,stp);
  BETA_G=zeros(MC_rep,stp);
 for rep = 1:MC_rep%�������
   %�����ʼֵ�Ѿ�ϵͳ�е�״̬���򣬲�����UG,A2G,���ϲ���U,�²��Ȼ��S;�²���G,�����Ȼ��A1
    x=zeros(1,N);
m=zeros(1,N);
neighborA1=randi([1,N],1,2*csmd_G);
 for nA1 = 1:csmd_G;
          x(neighborA1(nA1))=1;%�²�G��Ӧ�ϲ�A1
          m(neighborA1(nA1))=1;
 end
 for nA2 =csmd_G+1:2*csmd_G;%%��ѡ���ϲ��A2���²��Զ���Ӧ0����S
          m(neighborA1(nA2))=2;
 end
   %��ʼ����ÿ���������з��棬����ʱ����stp=40
    for t = 1:stp
      n=1:N;
      n(:)=0;
      y=1:N;
      y(:)=0;
          mediA1=zeros(1,N);
          mediA2=zeros(1,N);
%m��x�ֱ�Ϊ���²��״̬����m��0��1��2��ɣ�x��0��1��ɣ�m��x��Ӧ����tʱ�̣�n��y��Ӧt+1ʱ�̵����²�
      A=zeros(N,N);%�ڽӾ���A�ϲ�
      B=zeros(N,N);%�²��ڽӾ���B
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
%%��ʼ�������²����ϴ���
     for i = 1:N
          %�ϲ㴫��
        if m(i)==0 ;%(U,ϵͳֻ�д���US������ʱx�ض�Ϊ0��US�ɱ�ΪA1S,A2S��
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
        elseif m(i)==1;%A1��ϵͳ�д���A1S,A1G��A1S�ɱ�ΪUS��A1G�ɱ�UG 
            if x(i)==0;%A1S,���US
                p1=rand(1,1);
                if p1<delta1;
                    n(i)=0;%US
                else
                    n(i)=1;
                end
            else %x(i)==1��A1G���A1G����UG
                 p1=rand(1,1);
                 if p1<delta1;
                     n(i)=0;
                 else
                      n(i)=1;
                 end
            end
        else m(i)==2;%A2��ֻ����A2S,����ʱx��ȻΪ0��A2S�ɱ�ΪUS
            p2=rand(1,1);
            if p2<delta2;
              n(i)=0;  
            else
                n(i)=2;
            end
        end
         %�²㴫��
        if x(i)==0;%S״̬,ϵͳ�д���US,A1S,A2S��A1S����A1G��US�ȱ��UG,�ٱ��A1G;A2S�ȱ��A2G,�ٱ��A1G
          if n(i)==1;%A1S���ȴ����ӵ�A1S���ٴ�����ı��
               for j = 1:N
                   p1=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p1<gama1*beta_R);
                    y(i)=1;
                end
               end
          elseif n(i)==0;%US���A1G
              for j = 1:N
                   p0=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p0<beta_R);
                    y(i)=1;
                    n(i)=1;
                end
              end
          else n(i)==2;%A2S���ȱ��A2G,�ٱ��A1G
                for j = 1:N
                   p2=rand(1,1);
                if (B(j,i)==1)&&(x(j)==1)&&(p2<gama2*beta_R);
                    y(i)=1;
                    n(i)=1;
                end
              end
          end
        else  x(i)==1;%G״̬��ϵͳ��ֻ����A1G,�����м�״̬UG,A1G�ɱ�ΪA1S,UG�ɱ�ΪUS��A1G
           if n(i)==1;%A1G
           p1=rand(1,1);
              if p1<mu;
                y(i)=0;   %A1G���ߣ����A1S
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
      %����m��x����
      m=n;
      x=y;  
    A1=length(find(m(:)==1)); %����A1������
    A2=length(find(m(:)==2)); %����A2������        
    A3=length(find(m(:)==0)); %����A2������    
   
   BETA_A1(rep,t)=sum(A1)/N;%��l�η����У�A1���ܶ�
    BETA_A2(rep,t)=sum(A2)/N;%��l�η����У�A2���ܶ�
    BETA_U(rep,t)=sum(A3)/N;%��l�η����У�A1���ܶ�
    BETA_G(rep,t)=sum(x)/N;%��l�η����У�G���ܶ� 
    end
 end
 %ͳ��ÿ���еĲ�ͬ�ڵ��������x����ֻ��0��1�����Ժܺü��㣻m����0��1��2Ҫ����һ�㣬�ȼ���m�еĽڵ���
% %�����η����ƽ��ֵ
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
xlabel('t','FontSize',15);ylabel('proportion','FontSize',15);   % ���������
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,

