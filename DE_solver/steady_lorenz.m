%% Stability of the Lorenz Model
%
clear all

%% Fix spatial parameters
p    = 10;             % Prandtl number
b    = 8/3;            % geometry of domain
rvec = 0:0.1:5;        % Rayleigh number

%% Set up vectors for maximum growth rates
D1 = 0*rvec;
D2 = 0*rvec;
D3 = 0*rvec;

%% Loop over different values of r
cnt = 1;
for r = rvec

    % Pick out solution for a given r: steady solution 1
    X  = 0;
    Y  = 0;
    Z  = 0;

    % Compute Jacobian
    J = [ - p,  p,  0;...
          r-Z, -1, -X;...
            Y,  X, -b];
         
    % Find the largest real part of the eigenvalues
    [~,D]=eig(J);
    D1(cnt) = max(real(diag(D)));

    % Pick out solution for a given r: steady solution 2
    X  = sqrt(b*(r-1));
    Y  = X;
    Z  = r-1;

    % Compute Jacobian
    J = [ - p,  p,  0;...
          r-Z, -1, -X;...
            Y,  X, -b];

     % Find the largest real part of the eigenvalues
     [~,D]=eig(J);
     D2(cnt) = max(real(diag(D)));


    % Pick out solution for a given r: steady solution 3
    X  = -sqrt(b*(r-1));
    Y  =  X;
    Z  =  r-1;

    % Compute Jacobian
    J = [ - p,  p,  0;...
          r-Z, -1, -X;...
            Y,  X, -b];

     % Find the largest real part of the eigenvalues
     [~,D]=eig(J);
     D3(cnt) = max(real(diag(D)));

     cnt = cnt+1;
end

figure(1);clf;
plot(rvec, D1,'-g');
hold on;
plot(rvec, D2,'-b');
plot(rvec, D3,'--r');
plot(rvec,0*rvec,'-k');
plot([1 1],[-1.5 3.5],'-k');
hold off;
legend('zero','positive','negative',2);
xlabel('r')
title('Maximum Growth Rates');
grid on;
