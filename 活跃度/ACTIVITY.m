%ADN�������ɻ�Ծ��
%��������
N=1000;
yita=1;%��������
minActivity=0.1;%�½�
Exponent_A=2.1;%AD�����Ծָ��
Exponent_B=2.1;
%�������
AActivity=1:N;
BActivity=1:N;
    temAct=rand(1,N);
   for i = 1:N
      AActivity(i)=yita*power(((power(1,-Exponent_A+1)-power(minActivity,-Exponent_A+1))*temAct(i)+power(minActivity,-Exponent_A+1)),1/(-Exponent_A+1));
end
   temAct2=rand(1,N);
   for i = 1:N 
 BActivity(i)=yita*power(((power(1,-Exponent_B+1)-power(minActivity,-Exponent_B+1))*temAct2(i)+power(minActivity,-Exponent_B+1)),1/(-Exponent_B+1));
   end
save('AActivity.mat','AActivity')
save('BActivity.mat','BActivity')