LibODESolve \
\
A C++ library to solve ordinary differential equations \
\
The following components are currently available. \
&ensp;IntTrapezoid (module + OpenMP) \
&ensp;&ensp;Trapezoidal rule for 1D numerical integration \
&ensp;ODESolve::LibRK (class) \
&ensp;&ensp;Runge-Kutta method for initial condition problems \
&ensp;ODESolve::LibRX (class) \
&ensp;&ensp;Relaxation method for boundary condition problems \
&ensp;&ensp;(based on Numerical Recipes in C) \
\
It also contains the following examples to demonstrate library usage. \
&ensp;Calculation of pi \
&ensp;&ensp;by means of trapezoidal rule \
&ensp;&ensp;for 1-dimensional integration \
&ensp;OrbitKepler \
&ensp;&ensp;calculates Kepler orbits \
&ensp;&ensp;by solving Newtonian dynamics with Runge-Kutta method \
&ensp;GeodesicSphere \
&ensp;&ensp;finds geodesics (shortest path between two points) \
&ensp;&ensp;on a sphere \
\
This library can be built with cmake. \
One can build at a subdirectory with the following commands. \
&ensp;$ mkdir [subdirectory name] \
&ensp;$ cd [subdirectory name] \
&ensp;$ cmake [directory for the LibODESolve local repository] \
&ensp;$ cmake --build .
