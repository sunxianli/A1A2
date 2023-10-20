N = 1000  %100����
p = 0.005  %�����֮����0.1�ĸ����γ�����
%�����޸��Ϸ��Ĳ������õ���ͬ��ģ��
 
position=zeros(N,2);    %��λ����Ϣposition,һ����N�����ݣ�ÿ��������2����Ϣ
adj = zeros(N,N);  %�����ڽӾ��󣬳�ʼ���ڽӾ���ȫ��
for m=1:N           %��ÿ���㰲��λ�ã�Χ��һ��Բ
    position(m,1)=cos(m/N*2*pi);
    position(m,2)=sin(m/N*2*pi);
end
 
figure('name','ER���ͼ');
hold on;
plot(position(:,1),position(:,2),'d')
 
for m=1:N
    for n=m+1:N
        if(rand(1,1)<p)  %��0.1�ĸ������ɱ�
            adj(m,n)=1;  %����������ڽӱ�ֵ
            adj(n,m)=1;  
        end
    end
end
for m = 1:N
    for n = m+1:N
        if(adj(m,n)==1)  %����б߾ͻ�����
            plot(position([m,n],1),position([m,n],2));
        end
    end
end
% hold off;
 