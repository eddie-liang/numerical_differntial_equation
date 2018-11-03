function flux = logistic(N,t); 

global r k

flux = r*N*(1-N/k);

end
