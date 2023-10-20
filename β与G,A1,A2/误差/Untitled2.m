%从图像读取数据
% fig = openfig('MMCA.fig');
% lines = findobj(fig, 'type', 'line');
% xdata = get(lines, 'XData');
% ydata = get(lines, 'YData');
% ha = get(gcf,'Figure 1');  % 获取当前的图形的子对象：Axes坐标轴对象
% hl = get(ga,'Figure 1')    % 获取坐标轴的子对象：Line对象
% xdata = get(hl,'XData');
% ydata = get(hl,'YData');

handle = findobj(gca,'Type','line');%获取曲线的handle，如果图中有多条曲线，handle为一个数组
xdata = get(handle,'XData');%将handle的x数据赋给xdata一维数组
ydata = get(handle,'YData');%将handle的y数据赋给ydata一维数组