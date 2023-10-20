%�������
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
gama1=2;
gama2=0.5;
mu=0.5;
yita=10;
minActivity=0.001;
Exponent_A=2.5;
Exponent_B=2.5;
ma=5;
mb=5;
%��ʼ����
termi=20;%�µ�ȡֵ����2
MC_rep=1;%�������
BETA_A1_NEW=zeros(MC_rep,termi);
BETA_A2_NEW=zeros(MC_rep,termi);
BETA_G_NEW=zeros(MC_rep,termi);
neighbora=1:N;%���������������Ĵ�ŵ�
neighborb=1:N;

BETA_A1_AVER =1:termi;%��ֵ����
BETA_A2_AVER =1:termi;
BETA_G_AVER =1:termi;
AActivity=1:N;%��Ծ�Ⱦ���
BActivity=1:N;
  BETA_A1=1:termi;%ÿ�η���Ľ��
  BETA_A2=1:termi;
  BETA_G=1:termi;
temAct=rand(1,N);
   for i = 1:N
%       AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
       AActivity(i)=yita*power(((power(1,-Exponent_A+1)-power(minActivity,-Exponent_A+1))*temAct(i)+power(minActivity,-Exponent_A+1)),1/(-Exponent_A+1));
   end
   %��ʼ����ʱ��������ÿ���ڵ�Ļ�Ծ�̶�
   temAct2=rand(1,N);
   for i = 1:N 
       BActivity(i)=yita*power(((power(1,-Exponent_B+1)-power(minActivity,-Exponent_B+1))*temAct(i)+power(minActivity,-Exponent_B+1)),1/(-Exponent_B+1));
%        BActivity(i)= yita*minActivity*(1-(temAct2(i))^(1/(Exponent_B-1)));
   end
%��ʼ����ѭ��
for rep = 1:MC_rep%�������
    
  %��ʼ��ͬ�Ħ�ֵѭ��
  for l = 1:termi
   disp(['��',num2str(rep),'�η���','��Ϊ',num2str(l)])
   beta_R=l/40;
   %�����ʼֵ�Ѿ�ϵͳ�е�״̬���򣬲�����UG,A2G,���ϲ���U,�²��Ȼ��S;�²���G,�����Ȼ��A1
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
    %��ʼ����ʱ��������ÿ���ڵ�Ļ�Ծ�̶�
   
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
%�����ڽӾ���
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
    end
 %ͳ��ÿ���еĲ�ͬ�ڵ��������x����ֻ��0��1�����Ժܺü��㣻m����0��1��2Ҫ����һ�㣬�ȼ���m�еĽڵ���
 %for i=1:N
A1=length(find(m(:)==1)); %����A1������
A2=length(find(m(:)==2)); %����A2������             
    BETA_A1(l)=sum(A1)/N;%��l�η����У�A1���ܶ�

    BETA_A2(l)=sum(A2)/N;%��l�η����У�A2���ܶ�
    BETA_G(l)=sum(x)/N;%��l�η����У�G���ܶ�
  end
%����ÿ�η���Ľ��
  BETA_A1_NEW(rep,:)=BETA_A1;
  BETA_A2_NEW(rep,:)=BETA_A2;
  BETA_G_NEW(rep,:)=BETA_G;
end
%�����η����ƽ��ֵ
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
xlabel('$\beta$','FontSize',15);ylabel('proportion','FontSize',15);   % ���������
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,
