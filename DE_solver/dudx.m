function flux = dudx(u)

global Nz ikxx 

  uext = [u; -u(end:-1:1,:)];
  flux = real(ifft(ikxx.*fft(uext,[],2),[],2));
  flux = flux(1:Nz,:);
  
end