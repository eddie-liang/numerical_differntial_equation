function flux = Hyper(u)

global Nz khyp

  uext = [u; -u(end:-1:1,:)];
  flux = real(ifft2(khyp.*fft2(uext)));
  flux = flux(1:Nz,:);

end