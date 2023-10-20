clc
clear
N=1000;
csmd_I=0.08;
csmd_N=0.02;
lamda=0.10;
delta=0.64;
gama=0;
miu=0.48;
yita=1000;
minActivity=0.001;
Exponent_A=2.5;
Exponent_B=2.5;
AActivity=1:N;
BActivity=1:N;
ma=5;
mb=5;
mx=1 ;
ks=3;
kb=1;
ter=50;
MMCA=5;
stp=50;
% BETA1=zeros(mx,ter);
% BETA2=zeros(mx,ter);
        for l = 1:ter
            beta_R = l*kb/100;
            PRI=0.06*ones(1,N);
            PRS=0.92*ones(1,N);
            PNI=0.02*ones(1,N);
            PNS=0.00*ones(1,N);

            PRI_UPDATE=zeros(1,N);
            PRS_UPDATE=zeros(1,N);
            PNI_UPDATE=zeros(1,N);
            PNS_UPDATE=zeros(1,N);
            PR=1:N;
            PR(:)=0;
            r=zeros(1,N);
            q=zeros(1,N);

            R=zeros(N,N);
            Q=zeros(N,N);

            temAct=rand(1,N);
            for i = 1:N
                AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
            end
            temAct=rand(1,N);
            for i = 1:N
                BActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_B-1)));
            end

            for t = 1:stp
                for i =1:N
                    for j =1:N
                        R(j,i)=1-(AActivity(i)+AActivity(j))*(PRI(1,j)+PRS(1,j))*delta*ma/N;
                        Q(j,i)=1-(BActivity(i)+BActivity(j))*(PRI(1,j)+PNI(1,j))*beta_R*mb/N;
                    end
                    tempprod=cumprod(R(:,i));
                    r(1,i)=tempprod(N);
                    tempprod=cumprod(Q(:,i));
                    q(1,i)=tempprod(N);


                PNI_UPDATE(1,i)=PNI(1,i)*r(1,i)*(1-miu*gama)+PRI(1,i)*lamda*(1-miu*gama)+PRS(1,i)*lamda*(1-q(1,i));

                PRI_UPDATE(1,i)=PNI(1,i)*(1-r(1,i))*(1-miu)+PRI(1,i)*(1-lamda)*(1-miu)+PRS(1,i)*(1-lamda)*(1-q(1,i));

                PRS_UPDATE(1,i)=PNI(1,i)*((1-r(1,i))*miu+r(1,i)*miu*gama)+PRI(1,i)*(miu*(1-lamda)+lamda*miu*gama)+PRS(1,i)*(q(1,i));
                end
                PNI=PNI_UPDATE;
                PRI=PRI_UPDATE;
                PNS=0;
                PRS=PRS_UPDATE;
            end

            PN=PNI+PNS;
            PR=PRI;
            PI=PRI+PNI;

            BETA1=sum(PR)/N;
            BETA2=sum(PI)/N;
 if BETA2<0.001
    Theta_b2=0;
    Theta_b=0;
    Ba2=0;
    for i =1:N
        Theta_b=Theta_b+PR(1,i)*BActivity(i);
        Theta_b2=Theta_b2+PR(1,i)*(BActivity(i)^2);
        Ba2=Ba2+BActivity(i)^2;
    end
    Theta_b=Theta_b/N;
    Theta_b2=Theta_b2/N;
    Ba=sum(BActivity)/N;
    Ba2=Ba2/N;
    EH=sqrt((Ba2-Theta_b2)*(1-BETA1))+Ba-Theta_b;
    beta_MMCA=miu/mb/EH;
%    beta_c=mu/(mb*EH);%阈值的表达式，但是只有在β大于0.001时，他的上个值才是真的阈值
 else 
bata_c(l)= beta_MMCA;
    break;
end 

   end
        
% xzhou=(1:mx)*ks+2;
% plot(xzhou,beta_rec(1,:),'b-',xzhou,pic_rec2(1,:),'bo');





