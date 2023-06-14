%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Calculating the errors of the pressure                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load xc.dat
load yc.dat
load ap.dat
load p.dat

m = length(xc);
c1 = 6.8548e+04;
c2 = 159870;
err = zeros(m,m);
for i=1:m
    for j=1:m
        r = sqrt(xc(i)*xc(i)+yc(j)*yc(j));
        if r > 0.51       
           err(i,j) = (ap(i,j)-p(i,j)-2.8); 
        elseif r < 0.49
            err(i,j) = (ap(i,j)-p(i,j)-1.9);
        else
            err(i,j) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,3,1)
mesh(ap)
title('Exact')

subplot(1,3,2)
mesh(p)
title('Approx')

subplot(1,3,3)
mesh(err)
title('Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
