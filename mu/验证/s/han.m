load('PA2.mat')
load('PA1.mat')
load('BActivity.mat')
N=1000;
h=zeros(N,N);
gama1=0;%˥��ϵ��������Ϊ0ʱ��֪����Ϣ�ĸ�������ȫ���ߵġ�������������1��
gama2=0.5;
for j=1:N
    for i=1:N
h(j,i)=(1-(1-gama1)*PA1(1,i)-(1-gama2)*PA2(1,i))* BActivity(j,i);
    end
end
[m1,n1]=eig(h);
tez=max(diag(n1));
% bata_c=mu/tez;
% end
% hold on;
% xzhou=(1:20)/20;
% xlabel('\mu');
% ylabel('proportion');
% plot(xzhou,bata_c);