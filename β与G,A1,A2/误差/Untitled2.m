%��ͼ���ȡ����
% fig = openfig('MMCA.fig');
% lines = findobj(fig, 'type', 'line');
% xdata = get(lines, 'XData');
% ydata = get(lines, 'YData');
% ha = get(gcf,'Figure 1');  % ��ȡ��ǰ��ͼ�ε��Ӷ���Axes���������
% hl = get(ga,'Figure 1')    % ��ȡ��������Ӷ���Line����
% xdata = get(hl,'XData');
% ydata = get(hl,'YData');

handle = findobj(gca,'Type','line');%��ȡ���ߵ�handle�����ͼ���ж������ߣ�handleΪһ������
xdata = get(handle,'XData');%��handle��x���ݸ���xdataһά����
ydata = get(handle,'YData');%��handle��y���ݸ���ydataһά����