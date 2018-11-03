function flux = invLap(u)

global Nz kiLap

  uext = [u; -u(end:-1:1,:)];
  flux = real(ifft2(-kiLap.*fft2(uext)));
  flux = flux(1:Nz,:);
 
end