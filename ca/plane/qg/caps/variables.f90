module variables

 !Include all variables which may change during the course of a simulation.

 !Time, time step & fractional time steps:
double precision:: t,dt,dt2,dt3,dt6

 !Twist parameter for timing of surgery:
double precision:: twist

 !Counters for writing direct-access data (gridded and contours):
integer:: igrec,icrec

end module variables
