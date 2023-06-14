%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Calculating the errors of the pressure                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load xc.dat
load yc.dat
load ap.dat
load p.dat


% Errors __________________________________________________________________

m = length(xc);
errI = zeros(m,m);
errE = zeros(m,m);
for i=1:m
    for j=1:m
        r = sqrt(xc(i)*xc(i)+yc(j)*yc(j));
        if r > 0.5
          errE(i,j) = abs(ap(i,j)-p(i,j));
           
        else
          errI(i,j) = abs(ap(i,j)-p(i,j));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1)
mesh(ap)
title('Exact')

subplot(2,2,2)
mesh(p)
title('Approx')

subplot(2,2,3)
mesh(errI)
title('Interior Error')

subplot(2,2,4)
mesh(errE)
title('Exterior Error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
