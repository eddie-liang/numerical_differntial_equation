%% convection.m
%
% Solve the 2D Nonlinear Vorticity Equation with convection
% This is to compliment our study of the Lorenz attractor.
%
% Geometry: periodic in x and solid walls in z
%
% Evolution Eqns:	
%	q_t = - u q_x - v q_y + g*alpha*theta_x + Dissipation
%	T_t = - u T_x + v T_y + con*v			+ Dissipaiton
%
% Diagnostic Eqn:	q = psi_xx + psi_yy
%
% Numerical Method:
% 1) FFT to compute the derivatives in spectral space
% 2) Adams-Bashforth for Advection and Temperature terms
% 3) Crank-Nicholson for the Hyperviscosity 
%
%     Requires scripts:
%        dudx.m   - compute first derivative using spectral method
%        du2dx2.m - compute second derivative using spectral method

clear all

global nu nu4 ikxx ikzz kxx2 kzz2 
global Nz kLap kiLap khyp
global g alpha dToH dt

%% Grid Parameters
Lx  = 1.0;
Lz  = 1.0;
Nx  = 32;
Nz  = Nx;
dx  = 2*Lx/Nx;
dz  = 2*Lz/Nz;

%% Temporal Parameters
t0  = 0.0;
tf  = 50.0;
dt  = 1e-2;
npt = 50;

%% Physical Parameters
g 		= 9.81;                         % m/s^2
alpha 	= 1e-2;                         % /Kelvin
rho0	= 1024;                         % kg/m^3
dT		= 1;                            % K
nu		= 1e-4;                         % m^2/s
nu4     = 1e1*(min(dx,dz)/pi)^4/dt;     % m^4/s
kappa	= 1e-5;                         % m^2/s
dToH	= dT/Lz;                        % K/m

%% Define Grid
x   = [-Lx+dx/2:dx:Lx-dx/2]';
z   = [-Lz+dz/2:dz:Lz-dz/2]';
[xx,zz]=meshgrid(x,z);

t   = t0:dt:tf;
Nt  = length(t);
Ntp = round(Nt/npt);
tp  = t(1:npt:end);

%% Define wavenumber (frequency)
kx  =     pi/Lx*[0:Nx/2 -Nx/2+1:-1]';
kz  = 0.5*pi/Lz*[0:Nz   -Nz+1:-1]';
[kxx, kzz] = meshgrid(kx,kz);
ikxx = 1i*kxx;
kxx2 = kxx.^2;
ikzz = 1i*kzz;
kzz2 = kzz.^2;
kLap = kxx2 + kzz2;
kiLap= 1./kLap;      kiLap(1,1) = 0;
khyp = kLap.^2;

%% Initial Conditions
psi0  = 1e-4*rand(Nz,Nx);
T0    = 1e-4*rand(Nz,Nx);
q0    = Lap(psi0);

%% Create power spectrum arrays to store data
nkx = 4; 
nky = 4;
pow_spec1 = zeros(nkx,nky,Ntp+1);
pow_spec2 = zeros(nkx,nky,Ntp+1);

% Compute power spectrum
pext = [psi0; -psi0(end:-1:1,:)];
flux = abs(fft2(pext));
pow_spec1(:,:,1) = flux(1:nky,1:nkx);
text = [T0; -T0(end:-1:1,:)];
flux = abs(fft2(text));
pow_spec2(:,:,1) = flux(1:nky,1:nkx);

%% Time-Step the solution

% Use an Euler step and Crank-Nicholson for dissipation
[psi,q,T,NL1nm,NL2nm] = euler_step(psi0,q0,T0);

% Use an AB2 step and Crank-Nicholson for dissipation
[psi,q,T,NL1n, NL2n ] = AB2_step(psi,q,T,NL1nm,NL2nm);

cnt = 2;
for ii=3:Nt

    % Use an AB3 step and Crank-Nicholson for dissipation
    [psi,q,T,NL1nm,NL2nm,NL1n,NL2n] = AB3_step(psi,q,T,NL1nm,NL2nm,NL1n,NL2n);

    if mod(ii-1,npt)==0

       % Plot temperature
       figure(1); clf;
       contourf(x,z,psi);
       colorbar;
       title(['Streamfunction at t = ',num2str(t(ii))]);
       drawnow;

       % Plot streamfunction
       figure(2); clf;
       contourf(x,z,T);
       colorbar;
       title(['Temperature at t = ',num2str(t(ii))]);
       drawnow;
       
       % Compute power spectrum
       pext = [psi; -psi(end:-1:1,:)];
       flux = abs(fft2(pext));
       pow_spec1(:,:,cnt) = flux(1:nky,1:nkx);
       text = [T; -T(end:-1:1,:)];
       flux = abs(fft2(text));
       pow_spec2(:,:,cnt) = flux(1:nky,1:nkx);

       cnt = cnt+1;
    end
end

