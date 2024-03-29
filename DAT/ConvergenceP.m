%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            ERROR ANALISIS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Error Data                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     _______________________________________________________
%    |         P <-- U_exact               U <- P(one_step)  | 
%    |     N   p(one_step)-  p(one_step)   u(one_step)       | 
%    |_______________________________________________________|

 Sol1 = [ 40   2.9954e-004   3.1074e-004   1.0328e-004    
          80   1.3665e-005   2.7153e-005   2.8246e-005 
         160   6.3160e-006   1.5546e-005   7.3870e-006 
         320   3.1805e-006   5.2207e-006   1.8891e-006               
         640   1.7351e-006   1.9470e-006   4.7762e-007       
        1280   7.3760e-007   3.1248e-007   1.2007e-007];
%     _______________________________________________________
%    |    P <-- U(one_step)             U <- P(two_step)     |
%    |    p(two_step)-  p(two_step)+    u(two_step)          |
%    |_______________________________________________________|

 Sol2 = [ 5.4220e-003    5.4239e-003     1.2352e-004   
          6.0167e-004    4.1211e-003     3.3628e-005          
          3.0832e-004    2.0589e-003     8.7242e-006  %<-------
          2.7818e-004    7.8969e-004     2.3183e-006
          2.4626e-004    6.1702e-004     5.8408e-007
          1.2156e-004    3.4663e-004     1.4933e-007];
          
%  Tis are the silulation of p and u after eleminanting d

  Sol3 = [   1.2395e-003  1.9202e-004
             5.4395e-004  5.2888e-005
             2.7306e-004  1.4201e-005
             1.4190e-004  3.7037e-006  
             7.0139e-005  9.4572e-007
             3.5686e-005  2.3894e-007];
           
%     ______________________________________
%    |       P(infty) <-- U(infty)          |
%    |       p(t=infty)   u(t=infty)        |
%    |______________________________________|
                   
 Sol_inf = [ 3.6435e-003  1.8800e-004 ];
              
M   = Sol1(:,1);       
Err = Sol2(:,2);   
Dr =2./M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Plotting                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logDr   = log(Dr);
logErr  = log(Err);
logErr1 = 1*(logDr(end)-logDr(1))+logErr(1);  % First order
logErr2 = 2*(logDr(end)-logDr(1))+logErr(1);  % Second order
logErr3 = 3*(logDr(end)-logDr(1))+logErr(1);  % Third order
hold on
num = plot(logDr,logErr,'*k');
fir = plot([logDr(1),logDr(end)],[logErr(1),logErr1],':r');
sec = plot([logDr(1),logDr(end)],[logErr(1),logErr2],'--b');
thi = plot([logDr(1),logDr(end)],[logErr(1),logErr3],'-g');
plot_group = [num,fir,sec,thi];
legend(plot_group,'Numerical error','1st-order','2nd-order','3rd-order',2)
title('Error Analysis')
xlabel('log(\Deltax)')
ylabel('log(Error)')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Ratio                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format short g

for i=1:length(M)-1
    order(i,1) = abs(log(Err(1)/Err(i+1))/log(M(1)/M(i+1)));
end
disp('______________________________________________')
disp('                    Order                     ')
disp('______________________________________________')
disp('            M        Ratio     ')
disp([M,[0;order]])
disp('______________________________________________')

format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


