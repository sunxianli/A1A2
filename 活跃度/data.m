%%查看数据的分布
% 创建一组随机数据
data = randn(1000, 1);
% 绘制直方图
histogram(data, 'Normalization', 'probability'); % 按概率归一化
title('数据分布直方图');
xlabel('数值');
ylabel('概率');
%不归一化
% data = randn(1000, 1);
% histogram(data, 'Normalization', 'count'); % 不归一化，纵轴表示数据点的数量
% title('数据分布直方图');
% xlabel('数值');
% ylabel('概率');
%%查看是否是均匀分布
data = rand(1, 100);  % 假设数据是从均匀分布中随机生成的
histogram(data, 'Normalization', 'probability');
