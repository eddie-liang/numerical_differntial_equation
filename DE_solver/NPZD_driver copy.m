%% NPZD_driver.m
%
% Nutrent-Phytoplankton-Zooplankton-Detritus (NPZD) model
% From section 2.9 of "A practical guide to ecological modelling"

% Define biological parameters
global parms 

parms.mumax=2;
parms.kn=1;
parms.mp=0.1;
parms.rfa=0.7;
parms.mz=0.2;
parms.po=0.2;
parms.z0=0.2;
parms.no=1.6;
parms.a=1;
parms. imax=1.5;
parms. b=1;
parms. c=1.5;




% Define temporal parameters
t0 = 0;
tf = 70;
dt = 0.1;
ic = [0.2; 0.2; 1.6];

% Compute numerical solutions
[U1,t1]=RK4(@NPZD,t0,tf,dt,ic);

% Plot the solution
figure(1); clf;
subplot(2,2,1);
plot(t1, U1(1,:), 'LineWidth',4)
xlabel('time')
title('PHYTO');

subplot(2,2,2);
plot(t1, U1(2,:), 'LineWidth',4)
xlabel('time')
title('ZOO');

subplot(2,2,4);
plot(t1, U1(3,:), 'LineWidth',4)
xlabel('time')
title('DIN');
