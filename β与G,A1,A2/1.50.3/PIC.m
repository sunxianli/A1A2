load rhoG
load rhoA1
load rhoA2
ter=40;
xzhou=(1:ter)/80;
hold on;
box on;
grid off;
set(gca,'Fontsize',15);
plot(xzhou,rhoG','-o','color',[77/256 133/256 189/256],'MarkerFaceColor',[77/256 133/256 189/256]);
plot(xzhou,rhoA1','-^','color',[247/256 144/256 61/256],'MarkerFaceColor',[247/256 144/256 61/256]);
plot(xzhou,rhoA2,'-v','color',[89/256 169/256 90/256],'MarkerFaceColor',[89/256 169/256 90/256]);
set(gcf,'DefaultTextInterpreter','latex');
xlabel('$\beta$','FontSize',15);ylabel('proportion','FontSize',15);   % 鍧愭爣杞磋В閲�?
h=legend('$\rho^{G}$','$\rho^{A_1}$','$\rho^{A_2}$');
set(h,'Interpreter','latex','FontSize',15)%,
save('rhoG.mat','rhoG')
save('rhoA1.mat','rhoA1')
save('rhoA2.mat','rhoA2')
saveas(gcf, 'save.fig')