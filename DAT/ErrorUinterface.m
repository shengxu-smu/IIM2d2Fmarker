%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Calculating the errors of the interface velocity              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load E_Uxsys.dat
load E_Vxsys.dat
load N_Uxsys.dat
load N_Vxsys.dat

ns = length(E_Uxsys);
Error_Uxsys = abs(E_Uxsys - N_Uxsys);
Error_Vxsys = abs(E_Vxsys - N_Vxsys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1)
plot(1:ns,E_Uxsys,1:ns,N_Uxsys,'.r')
title('U interface')
legend('Analytical','Numerical',4)

subplot(2,2,2)
plot(Error_Uxsys,'.r')
title('Error')

subplot(2,2,3)
plot(1:ns,E_Vxsys,1:ns,N_Vxsys,'.r')
title('V interface')
legend('Analytical','Numerical',4)

subplot(2,2,4)
plot(Error_Vxsys,'.r')
title('Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%