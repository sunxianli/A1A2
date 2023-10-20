%%初始值的生成
excluded_numbers = [5, 10, 15]; % 假设你想要排除这些数
lower_limit = 1;
upper_limit = 20;
num_samples = 5; % 假设你要生成5个随机数

% 生成排除某些数后的随机数
random_numbers = zeros(1, num_samples);
count = 0;
while count < num_samples
    temp_rand = randi([lower_limit, upper_limit], 1, 1); % 生成一个随机整数
    if ~ismember(temp_rand, excluded_numbers) % 检查该随机数是否在排除列表中
        count = count + 1;
        random_numbers(count) = temp_rand;
    end
end

disp(random_numbers);
