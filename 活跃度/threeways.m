%��������
N=1000;
yita=10;%��������
minActivity=0.001;%�½�
Exponent_A=2.1;%AD�����Ծָ��
Exponent_B=2.1;
%�������
AActivity=1:N;
BActivity=1:N;
%%ʦ���Ƶ�
% x0=(1-minActivity^(Exponent_A-1))*rand(1,N);
% y0=(1-minActivity^(Exponent_B-1))*rand(1,N);
% AActivity=yita*minActivity*(1-x0).^(1/(1-Exponent_A));
% BActivity=yita*minActivity*(1-y0).^(1/(1-Exponent_B));
% %���Ͻ���
%     temAct=rand(1,N);
%    for i = 1:N
%       AActivity(i)=yita*power(((power(1,-Exponent_A+1)-power(minActivity,-Exponent_A+1))*temAct(i)+power(minActivity,-Exponent_A+1)),1/(-Exponent_A+1));
% end
%    temAct2=rand(1,N);
%    for i = 1:N 
%  BActivity(i)=yita*power(((power(1,-Exponent_B+1)-power(minActivity,-Exponent_B+1))*temAct2(i)+power(minActivity,-Exponent_B+1)),1/(-Exponent_B+1));
% end
% %%����ʦ����
    temAct=rand(1,N);
   for i = 1:N
      AActivity(i)=yita*minActivity*(1-(temAct(i))^(1/(Exponent_A-1)));
   end
   temAct2=rand(1,N);
   for i = 1:N 
  BActivity(i)= yita*minActivity*(1-(temAct2(i))^(1/(Exponent_B-1)));
   end
save('AActivity200.mat','AActivity');
save('BActivity200.mat','BActivity');