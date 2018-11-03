function [psi, q, T, NL1n, NL2n] = AB2_step(psi0, q0, T0, NL1nm, NL2nm)

global g alpha dToH nu nu4 dt kLap khyp Nz 

    u   = -dudz(psi0);
    v   =  dudx(psi0);
    q_x =  dudx(q0);
    q_z =  dudz(q0);
    T_x =  dudx(T0);
    T_z =  dudz(T0);

    NL1n  = - u.*q_x - v.*q_z + g*alpha*dudx(T0);
    Dis   = nu*Lap(q0) - nu4*Hyper(q0);
    rhs   = q0 + 0.5*dt*(3*NL1n - NL1nm + Dis);
    rhsext= [rhs; -rhs(end:-1:1,:)];
    rhshat= fft2(rhsext);
    qhat  = rhshat./(1-0.5*dt*(nu*kLap - nu4*khyp));
    q     = real(ifft2(qhat));
    q     = q(1:Nz,:);
    psi   = invLap(q);

    NL2n  = - u.*T_x - v.*T_z + dToH*dudx(psi0);
    Dis   = nu*Lap(T0) - nu4*Hyper(T0);
    rhs   = T0 + 0.5*dt*(3*NL2n - NL2nm + Dis);
    rhsext= [rhs; -rhs(end:-1:1,:)];
    rhshat= fft2(rhsext);
    That  = rhshat./(1-0.5*dt*(nu*kLap - nu4*khyp));
    T     = real(ifft2(That));
    T     = T(1:Nz,:);

end

