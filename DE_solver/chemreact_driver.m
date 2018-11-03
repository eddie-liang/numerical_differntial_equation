%% chemreact_driver.m
%
% Chemical-Reaction
% 
% Sovles a system of ODEs that describe a chemical reaction
% 
%         k1->
% E + D <----> I
%       <-k2
%
%         k3->				
% I + F  ----> E + G
%
% taken from section 2.4 from 
% "A practical guide to ecological modelling"

% Define parameters
global parms 

parms.k1=0.01/24;
parms.k2=0.1/24;
parms.k3=0.1/24;

t0 = 0;
tf = 300;
dt = 0.5;
ic = [100; 10; 1; 1; 0]; 

% Compute numerical solutions
[U1,t1]=RK4(@chemreact,t0,tf,dt,ic);

% Plot the solution
figure(1); clf;
subplot(2,3,1);
plot(t1,U1(1,:),'-b','LineWidth',4);
xlabel('time');
title('D');
subplot(2,3,2);
plot(t1,U1(2,:),'-b','LineWidth',4)
xlabel('time');
title('I (RK4)');
subplot(2,3,3);
plot(t1,U1(3,:),'-b','LineWidth',4)
xlabel('time');
title('E');
subplot(2,3,4);
plot(t1,U1(4,:),'-b','LineWidth',4)
xlabel('time');
title('F');
subplot(2,3,5);
plot(t1,U1(5,:),'-b','LineWidth',4)
xlabel('time');
title('G');





