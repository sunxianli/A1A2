load('PA2.mat')
load('PA1.mat')
load('BActivity.mat')
N=1000;
h=zeros(N,N);
gama1=0;%衰减系数，当其为0时，知道信息的个体是完全免疫的。在这里若等于1，
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