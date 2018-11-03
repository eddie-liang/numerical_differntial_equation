
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
npt = 333;

%% Physical Parameters
g 		= 9.81;                         % m/s^2
alpha 	= 1e-2;                         % /Kelvin
rho0	= 1024;                         % kg/m^3
dT		= 1000;                            % K
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
Nps = (Nt-1)/2;
Ntp = round(Nt/npt);
tp  = t(1:npt:end);
k   = [0:Nps/2 -Nps/2+1:-1]*2*pi/(dt*Nps);

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
flux1 = abs(fft2(pext));
pow_spec1(:,:,1) = flux1(1:nky,1:nkx);
text = [T0; -T0(end:-1:1,:)];
flux2 = abs(fft2(text));
pow_spec2(:,:,1) = flux2(1:nky,1:nkx);

%plot
figure(2); clf;
    subplot(3,1,1);
    semilogy(k(1:Nps/2), flux1(1:Nps/2),'-ob');
    xlim([0 k(Nps/2)]);
    title('X');
    subplot(3,1,2);
    semilogy(k(1:Nps/2), flux2(1:Nps/2),'-ob');
    xlim([0 k(Nps/2)]);
    %ylim([1 1e6]);
    title('Y');
    
    pause;