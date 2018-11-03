function flux = dudz(u)

global Nz ikzz 

  uext = [u; -u(end:-1:1,:)];
  flux = real(ifft(ikzz.*fft(uext,[],1),[],1));
  flux = flux(1:Nz,:);

end