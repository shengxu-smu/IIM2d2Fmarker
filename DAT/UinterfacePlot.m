
load E_Uxsys.dat
load E_Vxsys.dat
load N_Uxsys.dat
load N_Vxsys.dat

ns = length(E_Uxsys);

Error_Uxsys = abs(E_Uxsys - N_Uxsys);
Error_Vxsys = abs(E_Vxsys - N_Vxsys);

subplot(2,3,1)
plot(1:ns, E_Uxsys,'.r')
title('Exact U interface')

subplot(2,3,2)
plot(1:ns, N_Uxsys,'.r')
title('Numerical')

subplot(2,3,3)
plot(1:ns, Error_Uxsys,'.r')
title('Error')

subplot(2,3,4)
plot(1:ns, E_Vxsys,'.r')
title('Exact V interface')

subplot(2,3,5)
plot(1:ns, N_Vxsys,'.r')
title('Numerical')

subplot(2,3,6)
plot(1:ns, Error_Vxsys,'.r')
title('Error')
