%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Checking the erros of the Principal Jumps                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exact jumps data

load E_dudnp.dat
load E_dudnm.dat
load E_dvdnp.dat
load E_dvdnm.dat

% Numerical jumps data

load N_dudnp.dat
load N_dudnm.dat
load N_dvdnp.dat
load N_dvdnm.dat

ns = length(E_dudnp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors

err1 = abs(E_dudnp - N_dudnp);
err2 = abs(E_dudnm - N_dudnm);
err3 = abs(E_dvdnp - N_dvdnp);
err4 = abs(E_dvdnm - N_dvdnm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure(1)
subplot(4,3,1)
plot(1:ns, N_dudnp,'.r')
title('Approx (du/dn)+ ')

subplot(4,3,2)
plot(1:ns, E_dudnp,'.b')
title('Exact')

subplot(4,3,3)
plot(1:ns, err1,'.r')
title('Error')

%----------

subplot(4,3,4)
plot(1:ns, N_dudnm,'.r')
title('Approx (du/dn)- ')

subplot(4,3,5)
plot(1:ns, E_dudnm,'.b')
title('Exact')

subplot(4,3,6)
plot(1:ns, err2,'.r')
title('Error')

%----------

subplot(4,3,7)
plot(1:ns, N_dvdnp,'.r')
title('Approx (dv/dn)+ ')

subplot(4,3,8)
plot(1:ns, E_dvdnp,'.b')
title('Exact')

subplot(4,3,9)
plot(1:ns, err3,'.r')
title('Error')

%----------

subplot(4,3,10)
plot(1:ns, N_dvdnm,'.r')
title('Approx (dv/dn)- ')

subplot(4,3,11)
plot(1:ns, E_dvdnm,'.b')
title('Exact')

subplot(4,3,12)
plot(1:ns, err4,'.r')
title('Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%