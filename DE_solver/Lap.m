function flux = Lap(u)

global Nz kLap

  uext = [u; -u(end:-1:1,:)];
  flux = real(ifft2(-kLap.*fft2(uext)));
  flux = flux(1:Nz,:);

end