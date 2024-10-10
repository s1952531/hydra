module variables

use parameters

 !Include all variables which may change during the course of a simulation.

 !Point vortex positions:
double precision:: x(n),y(n),z(n)

 !Point vortex velocities:
double precision:: u(n),v(n),w(n)

 !A vector used in the calculation of velocities:
double precision:: unit(0:n)

 !Time, time step & their fractions used in RK4 scheme:
double precision:: t,dt,dt2,dt3,dt6

end module
