function [flux]=linearode(U,t)

global a

%% Compute's the flux of the linear ode
%  du/dt = a*u

flux = a*U;

end


