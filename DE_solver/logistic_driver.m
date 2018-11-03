% Define parameters

global r k

r  = 1;
k  = 1;

t0 = 0;
tf = 10;
dt = 0.1;
ic = [0.1];

% Compute numerical solution using Euler Method
[U1,t1]=RK4(@logistic,t0,tf,dt,ic);


% Plot the solution
figure(1); clf;
plot(t1,U1(1,:),'-b','LineWidth',4);


