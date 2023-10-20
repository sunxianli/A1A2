load rhoG
x=(1:40)/40;
y=(1:40)/40;
image(x,y,rhoG,'CDataMapping','scaled');
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\lambda_1$');
ylabel('$\lambda_2$');