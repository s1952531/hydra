module variables

use parameters

 !Include all variables which may change during the course of a simulation.

 !Point mass positions:
double precision:: x(n),y(n)

 !Point mass velocities:
double precision:: u(n),v(n)

 !Point mass accelerations:
double precision:: a(n),b(n)

 !A vector used in the calculation of velocities:
double precision:: unit(0:n)

 !Time, time step & their fractions used in RK4 scheme:
double precision:: t,dt,dt2,dt3,dt6

 !Record index for writing direct access data:
integer:: loop

end module
