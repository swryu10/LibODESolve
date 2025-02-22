# About
A C++ library to solve ordinary differential equations

# Ingredients
The following components are currently available.
* **IntTrapezoid** (module + OpenMP) \
Trapezoidal rule for 1-dimensional numerical integration
* **ODESolve::LibRK** (class) \
Runge-Kutta method for initial condition problems
* **ODESolve::LibRX** (class) \
Relaxation method for boundary condition problems (based on Numerical Recipes in C)

# Examples
It also contains the following examples to demonstrate library usage.
Please find documentations in doc subdirectory for more detail.
* Calculation of &pi; \
by means of trapezoidal rule for 1-dimensional integration
* **OrbitKepler** \
calculates Kepler orbits by solving Newtonian dynamics with Runge-Kutta method
* **GeodesicSphere** \
finds geodesics (shortest path between two points) on a sphere by solving geodesic equation with relaxation method

# Build
This library can be built with cmake. \
One can build at a subdirectory with the following commands. \
&ensp;$ mkdir [subdirectory name] \
&ensp;$ cd [subdirectory name] \
&ensp;$ cmake [directory for the LibODESolve local repository] \
&ensp;$ cmake --build .
