load gx.dat;
load gy.dat;
ip=2:4;


figure(1)
plot(gx(:,1),-2*gx(:,1),'-k',gy(:,1),gy(:,ip),'-o')
legend('Exact u','u','v','p')
title('u')

figure(2)
plot(gx(:,1),2*gx(:,1),'-k',gx(:,1),gx(:,ip),'-o')
legend('Exact v','u','v','p')
title('v')



