function main
global m n x y %����ȫ�ֱ���
clc, n=9; t=0:2*pi/n:2*pi; 
m=nchoosek(n,2); %����n���ڵ����ȫͼ�ı���m
x=cos(t); y=sin(t); 
axis([-1.1,1.1,-1.1,1.1])
subplot(1,2,1),plot(x,y,'o','Color','k')
subplot(1,2,2),hold on
[i1,j1]=myfun(0.1)  %���ݸ���0.1,���ú����������ߵĵ�ַ,����ͼ
figure, subplot(1,2,1), hold on
[i2,j2]=myfun(0.15)  %���ݸ���0.15,���ú����������ߵĵ�ַ,����ͼ
subplot(1,2,2), hold on
[i3,j3]=myfun(0.25)  %���ݸ���0.25,���ú����������ߵĵ�ַ,����ͼ
function [i,j]=myfun(p);%�ú������ݸ����ĸ���p,�������ߵĽڵ�i��j
global m n x y %����ȫ�ֱ���
z=rand(1,m);  %����m�������
ind1=(z<=p); %��z��С�ڵ���p�����������Ӧ�ĵ�ַ��������
ind2=squareform(ind1); %��0-1����ת�����ڽӾ���
[i,j]=find(ind2); %��ߵĽڵ���
plot(x,y,'o','Color','k')
for k=1:length(i)
  line([x(i(k)),x(j(k))],[y(i(k)),y(j(k))],'Color','k')
end
