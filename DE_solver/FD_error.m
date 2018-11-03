%% Preliminaries
clear all
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);

%% Specify grid spacing(???)
h = [0.001,0.005,0.01,0.05,0.1];

%% Specify u and dudx as inline functions  (http://www.mathworks.com/help/matlab/ref/inline.html)
u          = inline('exp(-x.^2)');
dudx    = inline('-2.*x.*exp(-x.^2)');

%% Compute the FD approximations (??.)
Dplus   = (u(1+h) - u(1))./h;
Dminus = (u(1) - u(1-h))./h;
Dzero   = (u(1+h) - u(1-h))./(2*h);
Dthree  = (2*u(1+h) + 3*u(1) - 6*u(1-h) + u(1-2*h))./(6*h);

%% Compute the exact value of dudx
duexact = dudx(1);

%% Plot the error on a log-log scale (http://www.mathworks.com/help/matlab/ref/loglog.html?searchHighlight=loglog)
figure(1);
clf;
loglog(h, abs(Dplus-duexact),'-ob','LineWidth',1.5,'MarkerSize',10);
hold on; %%?????????
loglog(h, abs(Dminus-duexact),'-xr','LineWidth',1.5,'MarkerSize',10);
loglog(h, abs(Dzero-duexact),'-sg','LineWidth',1.5,'MarkerSize',10);
loglog(h, abs(Dthree-duexact),'-vm','LineWidth',1.5,'MarkerSize',10);
hold off; %%?????
axis([0.0007 0.11 10^(-10) 10^(-1)]); %%axis([xmin xmax ymin ymax])
legend('D_{+}','D_{-}','D_0','D_3',4);
title('Error in FD methods to approximate du/dx');
xlabel('h');
grid on;

%% Find the slopes of the lines ( polyfit-finds the coefficients of a polynomial p(x) of degree n that fits the data)
 [p1,s]=polyfit(log(h), log(abs(Dplus-duexact)),1);
 [p2,s]=polyfit(log(h), log(abs(Dminus-duexact)),1);
 [p3,s]=polyfit(log(h), log(abs(Dzero-duexact)),1);
 [p4,s]=polyfit(log(h), log(abs(Dthree-duexact)),1);

 disp(['The slopes Dplus is ', num2str(p1(1))])
 disp(['The slopes Dminus is ', num2str(p2(1))])
 disp(['The slopes Dzero is ', num2str(p3(1))])
 disp(['The slopes Dthree is ', num2str(p4(1))])
 