LibODESolve \
\
A C++ library to solve ordinary differential equations \
\
Currently, the following components are available. \
&ensp;IntTrapezoid (module + OpenMP) \
&ensp;&ensp;Trapezoidal rule for 1D numerical integration \
&ensp;ODESolve::LibRK (class) \
&ensp;&ensp;Runge-Kutta method for initial condition problems \
&ensp;ODESolve::LibRX (class) \
&ensp;&ensp;Relaxation method for boundary condition problems \
&ensp;&ensp;(based on Numerical Recipes in C) \
\
This library can be built with cmake. \
One can build at a subdirectory with the following commands. \
&ensp;$ mkdir [subdirectory name] \
&ensp;$ cd [subdirectory name] \
&ensp;$ cmake [directory for the LibODESolve local repository] \
&ensp;$ cmake --build .
