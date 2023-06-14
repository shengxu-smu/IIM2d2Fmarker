%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Checking the erros of the Principal Jumps                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exact jumps data

load Ejc_dudn.dat
load Ejc_dvdn.dat
load Ejc_p.dat
load Ejc_1rhodpdn.dat
load Ejc_Lp.dat

% Numerical jumps data


load jc_dudn.dat
load jc_dvdn.dat
load jc_p.dat
load jc_1rhodpdn.dat
load jc_Lp.dat

ns = length(Ejc_dudn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors

err1 = abs(jc_dudn - Ejc_dudn);
err2 = abs(jc_dvdn - Ejc_dvdn);
err3 = abs(jc_p    - Ejc_p);
err4 = abs(jc_1rhodpdn - Ejc_1rhodpdn);
err5 = abs(jc_Lp - Ejc_Lp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure(1)
subplot(5,3,1)
plot(1:ns, jc_dudn,'.r')
title('Approx [du/dn]')

subplot(5,3,2)
plot(1:ns, Ejc_dudn,'.b')
title('Exact [du/dn]')

subplot(5,3,3)
plot(1:ns, err1,'.r')
title('Error [du/dn]')

%----------

subplot(5,3,4)
plot(1:ns, jc_dvdn,'.r')
title('Approx [dv/dn]')

subplot(5,3,5)
plot(1:ns, Ejc_dvdn,'.b')
title('Exact')

subplot(5,3,6)
plot(1:ns, err2,'.r')
title('Error')

%----------

subplot(5,3,7)
plot(1:ns, jc_p,'.r')
title('Approx [p]')

subplot(5,3,8)
plot(1:ns, Ejc_p,'.b')
title('Exact')

subplot(5,3,9)
plot(1:ns, err3,'.r')
title('Error')

%----------

subplot(5,3,10)
plot(1:ns, jc_1rhodpdn,'.r')
title('Approx [1/rho*dp/dn]')

subplot(5,3,11)
plot(1:ns, Ejc_1rhodpdn,'.b')
title('Exact')

subplot(5,3,12)
plot(1:ns, err4,'.r')
title('Error')

%----------

subplot(5,3,13)
plot(1:ns, jc_Lp,'.r')
title('Approx [\Deltap]')

subplot(5,3,14)
plot(1:ns, Ejc_Lp,'.b')
title('Exact')

subplot(5,3,15)
plot(1:ns, err5,'.r')
title('Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%