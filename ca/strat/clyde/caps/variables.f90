module variables

 !Include all variables which may change during the course of a simulation.

 !Time & time step:
double precision:: t,dt,hfdt,qudt

 !Sea and river uniform flow speeds:
double precision:: usea,vsea,uriv,vriv

 !Record for direct access writes:
integer:: igrids

end module
