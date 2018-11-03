function [psi, q, T, NL1nm, NL2nm] = euler_step(psi0, q0, T0)

global g alpha dToH nu nu4 dt kLap khyp Nz 

    u   = -dudz(psi0);
    v   =  dudx(psi0);
    q_x =  dudx(q0);
    q_z =  dudz(q0);
    T_x =  dudx(T0);
    T_z =  dudz(T0);

    NL1nm = - u.*q_x - v.*q_z + g*alpha*dudx(T0);
    Dis   = nu*Lap(q0) - nu4*Hyper(q0);
    rhs   = q0 + dt*(NL1nm + 0.5*Dis);
    rhsext= [rhs; -rhs(end:-1:1,:)];
    rhshat= fft2(rhsext);
    qhat  = rhshat./(1-0.5*dt*(nu*kLap - nu4*khyp));
    q     = real(ifft2(qhat));
    q     = q(1:Nz,:);
    psi   = invLap(q);

    NL2nm = - u.*T_x - v.*T_z + dToH*dudx(psi0);
    Dis   = nu*Lap(T0) - nu4*Hyper(T0);
    rhs   = T0 + dt*(NL2nm + 0.5*Dis);
    rhsext= [rhs; -rhs(end:-1:1,:)];
    rhshat= fft2(rhsext);
    That  = rhshat./(1-0.5*dt*(nu*kLap - nu4*khyp));
    T     = real(ifft2(That));
    T     = T(1:Nz,:);

end

