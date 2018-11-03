function flux = brusselator_1d(uv,t) 

global nu dx c bcx N

X1 = uv(1:N);
X2 = uv(N+1:2*N);

reaction  = [1 + X1.^2.*X2 - 4*X1;...
	       - X1.^2.*X2 + 3*X1];
diffusion = [nu(1)*d2udx2(X1,dx,bcx);...
             nu(2)*d2udx2(X2,dx,bcx)];
%advection = [- c(:,1).*upwind(  u,dx,bcx);...
%             - c(:,2).*upwind(  v,dx,bcx)];
advection = [- c(:,1).*dudx(  X1,dx,bcx);...
             - c(:,2).*dudx(  X2,dx,bcx)];

flux = reaction + diffusion + advection;

end
