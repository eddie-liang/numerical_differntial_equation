function flux = advect(psi, q)

  u   = -dudz(psi);
  v   =  dudx(psi);
  q_x =  dudx(q);
  q_z =  dudz(q);

  flux = u.*q_x + v.*q_z;
 
end