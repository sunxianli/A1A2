function main
global m n x y %定义全局变量
clc, n=9; t=0:2*pi/n:2*pi; 
m=nchoosek(n,2); %计算n个节点的完全图的边数m
x=cos(t); y=sin(t); 
axis([-1.1,1.1,-1.1,1.1])
subplot(1,2,1),plot(x,y,'o','Color','k')
subplot(1,2,2),hold on
[i1,j1]=myfun(0.1)  %根据概率0.1,调用函数计算连边的地址,并画图
figure, subplot(1,2,1), hold on
[i2,j2]=myfun(0.15)  %根据概率0.15,调用函数计算连边的地址,并画图
subplot(1,2,2), hold on
[i3,j3]=myfun(0.25)  %根据概率0.25,调用函数计算连边的地址,并画图
function [i,j]=myfun(p);%该函数根据给定的概率p,计算连边的节点i与j
global m n x y %定义全局变量
z=rand(1,m);  %生成m个随机数
ind1=(z<=p); %找z中小于等于p的随机数，对应的地址将来连边
ind2=squareform(ind1); %把0-1向量转换成邻接矩阵
[i,j]=find(ind2); %求边的节点编号
plot(x,y,'o','Color','k')
for k=1:length(i)
  line([x(i(k)),x(j(k))],[y(i(k)),y(j(k))],'Color','k')
end
