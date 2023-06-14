%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Checking the erros of the Principal Jumps                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exact jumps data

load Edudnp.dat
load Edudnm.dat
load Edvdnp.dat
load Edvdnm.dat

% Numerical jumps data

load Ndudnp.dat
load Ndudnm.dat
load Ndvdnp.dat
load Ndvdnm.dat

ns = length(Edudnp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors

err1 = abs(Edudnp - Ndudnp);
err2 = abs(Edudnm - Ndudnm);
err3 = abs(Edvdnp - Ndvdnp);
err4 = abs(Edvdnm - Ndvdnm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure(1)
subplot(4,3,1)
plot(1:ns, Ndudnp,'.r')
title('Approx (du/dn)+ ')

subplot(4,3,2)
plot(1:ns, Edudnp,'.b')
title('Exact')

subplot(4,3,3)
plot(1:ns, err1,'.r')
title('Error')

%----------

subplot(4,3,4)
plot(1:ns, Ndudnm,'.r')
title('Approx (du/dn)- ')

subplot(4,3,5)
plot(1:ns, Edudnm,'.b')
title('Exact')

subplot(4,3,6)
plot(1:ns, err2,'.r')
title('Error')

%----------

subplot(4,3,7)
plot(1:ns, Ndvdnp,'.r')
title('Approx (dv/dn)+ ')

subplot(4,3,8)
plot(1:ns, Edvdnp,'.b')
title('Exact')

subplot(4,3,9)
plot(1:ns, err3,'.r')
title('Error')

%----------

subplot(4,3,10)
plot(1:ns, Ndvdnm,'.r')
title('Approx (dv/dn)- ')

subplot(4,3,11)
plot(1:ns, Edvdnm,'.b')
title('Exact')

subplot(4,3,12)
plot(1:ns, err4,'.r')
title('Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%