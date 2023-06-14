
load xe.dat
load ye.dat
load xc.dat
load yc.dat
load xs.dat
load ys.dat
load Analytic_u.dat
load Analytic_v.dat
load Analytic_p.dat
load Exact_Uxsys.dat
load Exact_Vxsys.dat
load Numer_Uxsys.dat
load Numer_Vxsys.dat
load Uvelocity.dat
load Vvelocity.dat

ns = length(Exact_Uxsys);

figure(1)
hold on
mesh(xe(2:end-1),yc,Analytic_u)
plot3(xs(1:ns),ys(1:ns),Numer_Uxsys,'.r')
title('velocity u')
hold off

figure(2)
hold on
mesh(xc,ye(2:end-1),Analytic_v)
plot3(xs(1:ns),ys(1:ns),Numer_Vxsys,'.r')
title('velocity v')
hold off

figure(3)
mesh(xc,yc,Analytic_p)
title('Pressure p')

figure(4)
Error_Uxsys = abs(Exact_Uxsys - Numer_Uxsys);
Error_Vxsys = abs(Exact_Vxsys - Numer_Vxsys);

subplot(1,2,1)
plot(1:ns, Error_Uxsys,'.r')
subplot(1,2,2)
plot(1:ns, Error_Vxsys,'.r')

%hold on
%object    = plot(xs,ys,'*r');
%objectOut = plot(xsout,ysout,'ob');
%objectIn  = plot(xsin,ysin,'.k');
%objects   = [object,objectOut,objectIn];
%legend(objects,'Interface','Out Points',' In points')
%hold off