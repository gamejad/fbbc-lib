# fbbc-lib
C++ library for numerical FBBC-detector simulation (for NICA project).

## This library contains of 3 files:
- fbbc-lib.h which contains function prototypes and in which our
FBBC-class members and ParticleFBBC structure are defined;

- fbbc-lib.cpp in which the methods of your analysis are implemented;

- triangle_distribution.h which provide us with triangle random distribution,
based on uniform random distribution and Monte-Carlo method calculations;

## Functions and methods description
### double RZtoEta(double r, double z)
...

### double GetRandomNumber(const double first, const double last);
...
