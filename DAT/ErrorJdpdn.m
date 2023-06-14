%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Calculating the errors of the jump condition [dp/dn]              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ExactJdpdn.dat
load AproxJdpdn.dat

err = abs(ExactJdpdn - AproxJdpdn);
disp('____________________________________________')
disp('              Error [dp/dn]                 ')
disp('____________________________________________')
disp('                                            ')
disp(norm(err))
disp('____________________________________________')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,1,1)
plot(ExactJdpdn)
title('Exact [dp/dn]')

subplot(3,1,2)
plot(AproxJdpdn)
title('Approx')

subplot(3,1,3)
plot(err)
title('Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
