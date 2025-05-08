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
* **OrbitKepler** (class) \
calculates Kepler orbits by solving Newtonian dynamics with Runge-Kutta method
* **GeodesicSphere** (class) \
finds geodesics (shortest path between two points) on a sphere by solving geodesic equation with relaxation method

Note that, in the case of running the geodesic example, Python modules `plotly` and `airportsdata` are required to visualize and obtain latitude and longitude from IATA codes (for the origin and destination airports).
In a **Linux/UNIX** system, one can install `plotly` and `airportsdata` in a Python virtual environment (venv) by running the following commands.
```
$ python3 -m venv [directory for venv]
$ source [directory for venv]/bin/activate
$ python3 -m pip install plotly
$ python3 -m pip install airportsdata
```

# Build
This library can be built with **cmake**. \
In a **Linux/UNIX** system, one can build at a subdirectory with the following commands.
```
$ mkdir [subdirectory name]
$ cd [subdirectory name]
$ cmake [directory for the LibODESolve local repository]
$ cmake --build .
```
