module variables

 !Include all variables which may change during the course of a simulation.

implicit none

 !Time-stepping variables:
double precision:: t,dt,tdsave
double precision:: tmax,tsavori,tmaxori
integer:: iref,nper,idump

 !Logicals:
logical:: ramping


end module
