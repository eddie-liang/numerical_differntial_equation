function flux = logistic_1d(u,t); 

global r k nu dx c bcx

reaction  = r*u.*(1-u/k);                
diffusion = nu*d2udx2(u,dx,bcx);         
advection = - c.*dudx(  u,dx,bcx);       

flux = reaction + diffusion + advection;
end
