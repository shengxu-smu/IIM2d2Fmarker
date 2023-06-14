load xs.dat
load ys.dat
load xsout.dat
load ysout.dat
load xsin.dat
load ysin.dat

hold on
object    = plot(xs,ys,'*r');
objectOut = plot(xsout,ysout,'ob');
objectIn  = plot(xsin,ysin,'.k');
objects   = [object,objectOut,objectIn];
legend(objects,'Interface','Out Points',' In points')
hold off