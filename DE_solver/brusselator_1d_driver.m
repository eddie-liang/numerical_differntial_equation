%% brusselator_1d_driver.m
%
% This script solves the Lotka-Volterra-Advection-Diffusion Eqn:
%
% du/dt =  a*u   - b*u*v + nu(1)*d^2u/dx^2 - c(:,1)*du/dx
% dv/dt = -c*u*v - d*v   + nu(2)*d^2v/dx^2 - c(:,2)*dv/dx
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
global nu dx c bcx N

nu= 1e-4*[0; 0];
bcx='periodic';

% Define temporal variables
t0 = 0;
tf = 10;
dt = 0.1;
M  = length(t0:dt:tf);

% Define grid
L  = 1;
N  = 100;
dx = L/N;
x  = [0:dx:L-dx]';
c  = [0e0+0*x, 0*x];

% Initial Condition
ic = [1+sin(2*pi*x); 3+0*x];
 
% Compute numerical solutions
[U1,t1]=RK4(@brusselator_1d,t0,tf,dt,ic); 
u1 = U1(1:N,:);
v1 = U1(N+1:2*N,:);

% Plot the solution
figure(3); clf;
surf(t1,x,u1);
xlabel('time');
ylabel('space');
title('prey');

figure(4); clf;
surf(t1,x,v1);
xlabel('time');
ylabel('space');
title('predator');

for ii=1:1:M
   figure(1); clf;
   plot(x,u1(:,ii),'-ob');
   hold on;
   plot(x,v1(:,ii),'--r');
   hold off;
   axis([0 L 0 6]);
   legend('X1','X2');
   drawnow;
   pause(0.1);
end






