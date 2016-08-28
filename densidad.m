load out_press_225p_v4_g0.5.txt

minimo=0;
maximo=20;
delta=21;

xx=linspace(minimo,maximo,n);
yy=linspace(minimo,maximo,m);
contourf(xx,yy,N);
colorbar;
set(colorbar,'fontsize',17);
set(gca,'FontSize',17)

xlabel('Position x (m)') % x-axis label
ylabel('Position y (m)') % y-axis label
axis([15 20 0 20])

