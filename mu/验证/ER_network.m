N = 1000  %100个点
p = 0.005  %点与点之间以0.1的概率形成连边
%可以修改上方的参数，得到不同的模型
 
position=zeros(N,2);    %点位置信息position,一共有N组数据，每组数据有2个信息
adj = zeros(N,N);  %创建邻接矩阵，初始化邻接矩阵全零
for m=1:N           %给每个点安排位置，围成一个圆
    position(m,1)=cos(m/N*2*pi);
    position(m,2)=sin(m/N*2*pi);
end
 
figure('name','ER随机图');
hold on;
plot(position(:,1),position(:,2),'d')
 
for m=1:N
    for n=m+1:N
        if(rand(1,1)<p)  %以0.1的概率生成边
            adj(m,n)=1;  %这里两句给邻接表赋值
            adj(n,m)=1;  
        end
    end
end
for m = 1:N
    for n = m+1:N
        if(adj(m,n)==1)  %如果有边就画出来
            plot(position([m,n],1),position([m,n],2));
        end
    end
end
% hold off;
 