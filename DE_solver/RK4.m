function [u,t] = RK4(f, a, b, dt, u0);

t=a:dt:b;
N=length(t);
u=zeros(length(u0), length(t));
u(:,1) = u0;

for I=1:(N-1)
    Y1 = u(:,I);
    k1 = f(Y1,t(I));
    Y2 = u(:,I) + 0.5*dt*k1;
    k2 = f(Y2,t(I)+dt/2);
    Y3 = u(:,I) + 0.5*dt*k2;
    k3 = f(Y3,t(I)+dt/2);
    Y4 = u(:,I) + dt*k3;
    k4 = f(Y4,t(I+1));
    u(:,I+1) = u(:,I) + dt/6*(k1+2*k2+2*k3+k4);
end

end