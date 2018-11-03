%% linearode_driver.m
%
% This code will numerically integrate the linear ODE
%    du/dt = a*u

% Publish using the following:
% publish('linearode_driver.m','format','html','imageFormat','jpg')
% Define parameters
global a;
a  = 4;
t0 = 0;
tf = 1;
dv = 0.1./2.^[0:6];
ic = 1;

N  = length(dv);
tm = zeros(1,N);
er = zeros(1,N);

% Loop over different time steps to find error
for ii = 1:N

  dt = dv(ii);

  % Compute numerical solution using Euler Method
  tic
  [U1,t1]=ForwardEuler(@linearode,t0,tf,dt,ic);
  tm(ii)=toc;

  % Define exact solution
  u1 = ic*exp(a*[t0:dt:tf]);
  er(ii) = max(abs(u1-U1));

  % Display computing time and error and plot soln
  output(tm(ii), er(ii), dt, t1, u1, U1);

end

% Print error
figure(3); clf;
loglog(dv,er,'-or','LineWidth',4);
xlabel('log(time)');
ylabel('log(error)');

% Compute slope of global error
[p,s]=polyfit(log(dv), log(er), 1);
disp(['Best fit for the global error is ', num2str(p(1))]);
disp(' ');




