module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral

 !Define quantities to be preserved between recontouring and evolution:

 !Physical fields:
double precision:: ee(ng,nt),see(ng,nt),eep(ng,nt)
double precision:: dd(ng,nt),sdd(ng,nt),ddp(ng,nt)
double precision:: tt(ng,nt),stt(ng,nt),ttp(ng,nt)
double precision:: qq(ng,nt),qc(ng,nt)
double precision:: zz(ng,nt),pp(ng,nt)
double precision:: sqs(ng,nt),sqd(ng,nt)
double precision:: hhe(ng,nt),hhb(ng,nt),hhp(ng,nt) 
double precision:: topo(ng,nt) 

 !Logical to indicate presence or not of bottom topography:
logical:: topogr

end module common
