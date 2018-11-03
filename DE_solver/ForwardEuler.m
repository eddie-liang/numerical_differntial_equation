function [u,t] = ForwardEuler(f,a,b,dt,u0)

t       = a:dt:b;
N       = length(t);
u       = zeros(length(u0), N);
u(:,1)  = reshape(u0,length(u0),1);

for I=1:(N-1)
    u(:,I+1) = u(:,I) + dt*f(u(:,I), t(I));
end

end