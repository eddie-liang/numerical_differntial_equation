function flux = d2udz2(u)

global kzz2 

  flux = real(ifft(-kzz2.*fft(u,[],1),[],1));
 
end