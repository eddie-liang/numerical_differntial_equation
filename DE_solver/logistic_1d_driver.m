%% logistic_1d_driver.m
%
% This script will solve the Logistic-Diffusion Eqn:
%
% du/dt = nu*d^2u/dx^2 - c*du/dx + r*u*(1-u/K)
%
% Uses second-order centre differencing in space
%
% units: [nu] = m^2/s
%        [c]  = m/s
% 
% Numerical stability requires:
%  1) dt < dx^2/max(nu)     (with diffusion)
%  2) dx < dx/max(c)        (with advection)

% Define biological parameters
global r k nu dx c bcx
r = 0.05;
k = 100;
nu= 0.3;
bcx='periodic';

% Define temporal variables
t0 = 0;
tf = 30;
dt = 0.5;
M  = length(t0:dt:tf);

% Define grid
L  = 30;
N  = 60;
dx = L/N;
x  = [0:dx:L-dx]';
c  = 0+0*x;

% Initial Condition
%ic = 0.5*(1+sin(2*pi*x/L));
ic = 0*x;
ic(30:31) = 1;

% Compute numerical solutions
[U1,t1]=RK4(@logistic_1d,t0,tf,dt,ic); 

% Plot the solution
figure(1); clf;
for ii=1:1:M
   plot(x,U1(:,ii),'-ob');
   axis([0 L 0 1.05]);
   title(['ii = ', num2str(ii)]);
   drawnow;

   
   
   pause(0.1);
end

figure(2); clf;
surf(t1,x,U1);
xlabel('time');
ylabel('space');
title('Logistic-Diffusion Eqn');

figure(3); clf;
contourf(x,t1,U1');
colorbar;
caxis([0 0.5]);
colorbar;
ylabel('time');
xlabel('space');
title('Logistic-Diffusion Eqn');

figure(4); clf;
plot(x,U1(:,1),'-b','LineWidth',2);
hold on;
plot(x,U1(:,101),'-r','LineWidth',2);
plot(x,U1(:,201),'-g','LineWidth',2);
plot(x,U1(:,301),'-k','LineWidth',2);
plot(x,U1(:,401),'-m','LineWidth',2);
hold off;
legend('t=0 (days)','t=100 (days)','t=200 (days)','t=300 (days)','t=400 (days)');
xlabel('space');
title('Aphid Model');

figure(5); clf;
plot(t1,U1(30,:),'-b','LineWidth',2);
hold on;
plot(t1,U1(1,:),'-r','LineWidth',2);
hold off;
xlabel('time');
legend('Middle','Edge');
title('Aphid Model');





