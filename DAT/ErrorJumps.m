%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Checking the erros of the Principal Jumps                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exact jumps data

load Ejc_u.dat
load Ejc_dudn.dat
load Ejc_Lu.dat
load Ejc_v.dat
load Ejc_dvdn.dat
load Ejc_Lv.dat
load Ejc_p.dat
load Ejc_1rhodpdn.dat
load Ejc_Lp.dat

% Numerical jumps data

load Njc_u.dat
load Njc_dudn.dat
load Njc_Lu.dat
load Njc_v.dat
load Njc_dvdn.dat
load Njc_Lv.dat
load Njc_p.dat
load Njc_1rhodpdn.dat
load Njc_Lp.dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors

err1u = abs(Njc_u - Ejc_u);
err2u = abs(Njc_dudn - Ejc_dudn);
err3u = abs(Njc_Lu - Ejc_Lu);
err1v = abs(Njc_v - Ejc_v);
err2v = abs(Njc_dvdn - Ejc_dvdn);
err3v = abs(Njc_Lv - Ejc_Lv);
err1p = abs(Njc_p - Ejc_p);
err2p = abs(Njc_1rhodpdn - Ejc_1rhodpdn);
err3p = abs(Njc_Lp - Ejc_Lp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

ns = length(Ejc_dudn);

figure(1)

subplot(3,3,1)
plot(1:ns, Ejc_u)
title('Exact [u]')
subplot(3,3,2)
plot(1:ns, Njc_u,'.r')
title('Approx')
subplot(3,3,3)
plot(1:ns, err1u,'.r')
title('Error')

subplot(3,3,4)
plot(1:ns, Ejc_dudn)
title('Exact [du/dn]')
subplot(3,3,5)
plot(1:ns, Njc_dudn,'.r')
title('Approx')
subplot(3,3,6)
plot(1:ns, err2u,'.r')
title('Error')

subplot(3,3,7)
plot(1:ns, Ejc_Lu)
title('Exact [Lu]')
subplot(3,3,8)
plot(1:ns, Njc_Lu,'.r')
title('Approx')
subplot(3,3,9)
plot(1:ns, err3u,'.r')
title('Error')

%-------------------------------------------------------------------------

figure(2)

subplot(3,3,1)
plot(1:ns, Ejc_v)
title('Exact [v]')
subplot(3,3,2)
plot(1:ns, Njc_v,'.r')
title('Approx')
subplot(3,3,3)
plot(1:ns, err1v,'.r')
title('Error')

subplot(3,3,4)
plot(1:ns, Ejc_dvdn)
title('Exact [dv/dn]')
subplot(3,3,5)
plot(1:ns, Njc_dvdn,'.r')
title('Approx')
subplot(3,3,6)
plot(1:ns, err2v,'.r')
title('Error')

subplot(3,3,7)
plot(1:ns, Ejc_Lv)
title('Exact [Lv]')
subplot(3,3,8)
plot(1:ns, Njc_Lv,'.r')
title('Approx')
subplot(3,3,9)
plot(1:ns, err3v,'.r')
title('Error')

%-------------------------------------------------------------------------

figure(3)

subplot(3,3,1)
plot(1:ns, Ejc_p)
title('Exact [p]')
subplot(3,3,2)
plot(1:ns, Njc_p,'.r')
title('Approx')
subplot(3,3,3)
plot(1:ns, err1p,'.r')
title('Error')

subplot(3,3,4)
plot(1:ns, Ejc_1rhodpdn)
title('Exact [1/rho*dp/dn]')
subplot(3,3,5)
plot(1:ns, Njc_1rhodpdn,'.r')
title('Approx')
subplot(3,3,6)
plot(1:ns, err2p,'.r')
title('Error')

subplot(3,3,7)
plot(1:ns, Ejc_Lp)
title('Exact [Lp]')
subplot(3,3,8)
plot(1:ns, Njc_Lp,'.r')
title('Approx')
subplot(3,3,9)
plot(1:ns, err3p,'.r')
title('Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%