% �������ݣ��������ݴ洢����Ϊdata�ı����У�
% data = ...
 load('BActivity.mat')
data =  BActivity;
% ����������ˮƽ��ͨ��Ϊ0.05��
alpha = 0.05;
% ��������
data = sort(data, 'descend'); % �Ӵ�С����
% �����ۻ��ֲ�������CDF��
n = length(data);
cdf = (1:n) / n;

% �������ɷֲ�����
[power_law_params, ~, ~] = plfit(data);

% �����������ɷֲ���CDF
x_min = min(data);
alpha_hat = power_law_params(1);
theory_cdf = 1 - (data / x_min).^(1 - alpha_hat);

% ʹ�� Kolmogorov-Smirnov ����ȽϹ۲�CDF������CDF
[h, p] = kstest2(cdf, theory_cdf);

% ������������
if p < alpha
    fprintf('���ݲ��������ɷֲ� (p = %.4f)\n', p);
else
    fprintf('���ݷ������ɷֲ� (p = %.4f)\n', p);
end

% ���ӻ����
figure;
subplot(2, 1, 1);
loglog(data, cdf, 'b.', 'MarkerSize', 10);
title('���ݵ�˫��������ͼ');
xlabel('����ֵ');
ylabel('CDF');

subplot(2, 1, 2);
loglog(data, theory_cdf, 'r-', 'LineWidth', 2);
title('�۲�CDF������CDF�ıȽ�');
xlabel('����ֵ');
ylabel('CDF');

legend('�۲�CDF', '����CDF', 'Location', 'SouthEast');

