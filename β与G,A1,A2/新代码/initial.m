%%��ʼֵ������
excluded_numbers = [5, 10, 15]; % ��������Ҫ�ų���Щ��
lower_limit = 1;
upper_limit = 20;
num_samples = 5; % ������Ҫ����5�������

% �����ų�ĳЩ����������
random_numbers = zeros(1, num_samples);
count = 0;
while count < num_samples
    temp_rand = randi([lower_limit, upper_limit], 1, 1); % ����һ���������
    if ~ismember(temp_rand, excluded_numbers) % ����������Ƿ����ų��б���
        count = count + 1;
        random_numbers(count) = temp_rand;
    end
end

disp(random_numbers);
