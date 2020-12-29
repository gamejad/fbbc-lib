# fbbc-lib
C++ library for numerical FBBC-detector simulation (for NICA project).

## This library contains of 3 files:
- **fbbc-lib.h** which contains function prototypes and in which our
FBBC-class members and ParticleFBBC structure are defined;

- **fbbc-lib.cpp** in which the methods of your analysis are implemented;

- **triangle_distribution.h** which provide us with triangle random distribution,
based on uniform random distribution and Monte-Carlo method calculations;

## Functions and methods description

### Support functions

- **double RZtoEta(double r, double z)** - 
calculates pseudorapidity using formulas $\theta=\atan{\frac{r}{z}}$ and $\eta=-\log{\tan{\frac{\theta}{2}}}$,
where $(r,\varphi,z)$ is given particle 
coordinate (in cilindrical coordinates) relatively interaction point $(0,0,0)$

- **double GetRandomNumber(const double first=0, const double last=1)** - 
returns a random number, generated with uniform distribution in the interval [first, last].
If no arguments, given interval equals [0, 1].

## ParticleFBBC structure
...

## FBBCDetector class
...
