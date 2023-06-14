%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Calculating the errors of the pressure                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load xc.dat
load yc.dat
load ap.dat
load p.dat
load SolExample2.dat

p = SolExample2;
m = length(xc);
cE = 2.65234;
cI = 2.15245;
errI = zeros(m,m);
errE = zeros(m,m);
for i=1:m
    for j=1:m
        r = sqrt(xc(i)*xc(i)+yc(j)*yc(j));
        if r > 0.52       
           errE(i,j) = (ap(i,j)-p(i,j)-cE); 
        elseif r < 0.48
            errI(i,j) = (ap(i,j)-p(i,j)-cI);
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
