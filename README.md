# fbbc-lib
C++ library for numerical FBBC-detector simulation (for NICA project).

## This library contains of 3 files:
- **fbbc-lib.h** which contains function prototypes and in which our
FBBC-class members and ParticleFBBC structure are defined;

- **fbbc-lib.cpp** in which the methods of your analysis are implemented;

- **triangle_distribution.h** unneccessary suport file, which provide us with triangle random distribution,
based on uniform random distribution and Monte-Carlo method calculations; anyway,
this file is unneccessary for **fbbc-lib** is working.

## Support Functions description

- **double RZtoEta(double r, double z)** - 
calculates pseudorapidity Eta using formulas theta=atan(r/z) and Eta=-log(tan(theta/2)),
where (r,phi,z) is given particle 
coordinate (in cilindrical coordinates) relatively interaction point (0,0,0).

- **double GetRandomNumber(const double first=0, const double last=1)** - 
returns a random number, generated with uniform distribution in the interval [first, last].
If no arguments, given interval equals [0, 1].

## ParticleFBBC structure
In short, this structure defines the way that FBBC "see" particles.
It is supposed, that every particle in a single collision can be fully described (from a point of
view of a FBBC detector facility) with:
- index number Id,
- energy E,
- momenta P,
- its longitudinal component Pz,
- polar angle Phi,
- Z coordinate  of interaction point.

With a knowledge of E, P, Pz, Phi and Z we can calculate coordinates of particle contact with MCP plates
and time of this contact. Index number Id provide us an ability go back from ParticleFBBC to initial
particle format (e.g. ROOT branch).

## FBBCDetector class
This is the main library class...
