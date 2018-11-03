function flux = d2udx2(u)

global kxx2 

  flux = real(ifft(-kxx2.*fft(u,[],2),[],2));
 
end