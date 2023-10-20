% 导入数据（假设数据存储在名为data的变量中）
% data = ...
 load('BActivity.mat')
data =  BActivity;
% 设置显著性水平（通常为0.05）
alpha = 0.05;
% 数据排序
data = sort(data, 'descend'); % 从大到小排序
% 计算累积分布函数（CDF）
n = length(data);
cdf = (1:n) / n;

% 估计幂律分布参数
[power_law_params, ~, ~] = plfit(data);

% 生成理论幂律分布的CDF
x_min = min(data);
alpha_hat = power_law_params(1);
theory_cdf = 1 - (data / x_min).^(1 - alpha_hat);

% 使用 Kolmogorov-Smirnov 检验比较观察CDF和理论CDF
[h, p] = kstest2(cdf, theory_cdf);

% 输出假设检验结果
if p < alpha
    fprintf('数据不符合幂律分布 (p = %.4f)\n', p);
else
    fprintf('数据符合幂律分布 (p = %.4f)\n', p);
end

% 可视化结果
figure;
subplot(2, 1, 1);
loglog(data, cdf, 'b.', 'MarkerSize', 10);
title('数据的双对数坐标图');
xlabel('数据值');
ylabel('CDF');

subplot(2, 1, 2);
loglog(data, theory_cdf, 'r-', 'LineWidth', 2);
title('观察CDF与理论CDF的比较');
xlabel('数据值');
ylabel('CDF');

legend('观察CDF', '理论CDF', 'Location', 'SouthEast');

