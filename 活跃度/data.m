%%�鿴���ݵķֲ�
% ����һ���������
data = randn(1000, 1);
% ����ֱ��ͼ
histogram(data, 'Normalization', 'probability'); % �����ʹ�һ��
title('���ݷֲ�ֱ��ͼ');
xlabel('��ֵ');
ylabel('����');
%����һ��
% data = randn(1000, 1);
% histogram(data, 'Normalization', 'count'); % ����һ���������ʾ���ݵ������
% title('���ݷֲ�ֱ��ͼ');
% xlabel('��ֵ');
% ylabel('����');
%%�鿴�Ƿ��Ǿ��ȷֲ�
data = rand(1, 100);  % ���������ǴӾ��ȷֲ���������ɵ�
histogram(data, 'Normalization', 'probability');
