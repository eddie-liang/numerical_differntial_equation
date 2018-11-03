function flux = chemreact(u,t);  

global parms

% Define parameters
k1 = parms.k1;
k2 = parms.k2;
k3 = parms.k3;

% Define Varaibles
D = u(1);
I = u(2);
E = u(3);
F = u(4);
G = u(5);

% Define flux
flux = [-k1*E*D + k2*I;...
         k1*E*D - k2*I - k3*I*F;...
	-k1*E*D + k2*I + k3*I*F;...
	               - k3*I*F;...
	                 k3*I*F];

end
