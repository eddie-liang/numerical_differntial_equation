function flux = NPZD(u,t);  

global parms

% Define parameters
mumax            = parms.mumax;
kn                = parms.kn;
mp                = parms.mp;
rfa           = parms.rfa;
mz            = parms.mz;
po              = parms.po;
z0        = parms.z0;
no        = parms.no;
a   = parms.a;
imax            = parms. imax;
b            = parms. b;
c           = parms. c;

% Define Varaibles
PHYTO    = u(1);
ZOO      = u(2);
DIN      = u(3);

% Define the flux
      
Nuptake        = (DIN*mumax*PHYTO)/(kn+DIN);
Grazing        = ZOO*(1-exp(-1*rfa*PHYTO))*imax;



dPHYTO    = Nuptake - Grazing-mp*PHYTO;
dZOO      = rfa*(1-exp(-1*rfa*PHYTO))*imax*ZOO-mz*ZOO;
dDIN      = Nuptake*(-1)+(1-rfa)*ZOO*(1-exp(-1*rfa*PHYTO))*imax;


flux = [dPHYTO; dZOO; dDIN];

end
