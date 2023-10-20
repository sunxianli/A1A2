
handle = findobj(gca,'Type','line');%获取曲线的handle，如果图中有多条曲线，handle为一个数组
xdata = get(handle,'XData');%将handle的x数据赋给xdata一维数组
ydata = get(handle,'YData');%将handle的y数据赋给ydata一维数组