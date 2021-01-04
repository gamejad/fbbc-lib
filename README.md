# fbbc-lib
C++ library for numerical FBBC-detector simulations (for NICA project). 

To meet the challenges of the fast monitoring of the beam-beam collisions of the
high intensity NICA beams, the compact MCP-based Fast Beam-Beam Collisions (FBBC) 
detector with high timing properties was proposed in
(A.A. Baldin, G.A. Feofilov, P. Har'yuzov, F.F. Valiev, NIMA, v.958, 1 (2020), 162154).
The FBBC is based on the application of Microchannel plate detectors (Fig.1) for the
precise determination of arrival times of charged particles produced in the nucleus-nucleus collisions at NICA.
![alt text](https://github.com/vsandul/fbbc-lib/blob/master/pictures/fbbc.png)
We supposed FBBC as a two sets of MCPs placed on the left and right sides  symmetrically
from interaction point (Fig.2). The number of MCPs on the single side is optional.
We suppose every single MCP as a ring with some inner and outer diameters.
MCP separated on some amount of sectors, which one of them 
is connected to read-out electronic channel (Fig.1).

## This library contains of 3 files:
- **fbbc-lib.h** which contains function prototypes and in which our
FBBC-class members and ParticleFBBC structure are defined;
- **fbbc-lib.cpp** in which the methods of your analysis are implemented;
- **triangle_distribution.h** - unneccessary support file, which provide us with triangle random distribution,
based on uniform random distribution and Monte-Carlo method calculations; anyway,
this file is unneccessary for **fbbc-lib** is working.

## Support Functions description

- **double RZtoEta(double r, double z)** - 
calculates pseudorapidity η using formulas θ=atan(r/z) and η=-log(tan(θ/2)),
where (r,φ,z) is given particle 
coordinate (in cilindrical coordinates) relatively interaction point (0,0,0).
- **double GetRandomNumber(const double first=0, const double last=1)** - 
returns a random number, generated with uniform distribution in the interval [first, last].
If no arguments, given interval equals [0, 1].

## ParticleFBBC structure
In short, this structure defines the way that FBBC "see" particles.
It is supposed, that every particle in a single collision can be fully described (from a point of
view of a FBBC detector facility) with:
- index number **Id**,
- energy **E**,
- momentum **P**,
- its longitudinal component **Pz**,
- polar angle **Phi**,
- **Z** coordinate  of interaction point.

With a knowledge of E, P, Pz, Phi and Z we can calculate coordinates of particle contact with MCP plates
and time of this contact. Index number Id provide us an ability go back from ParticleFBBC to initial
particle format (e.g. ROOT branch).

## FBBCDetector class
This is the main library class, which simulates the work of FBBC detector facility. At the input it accepts 
a vector of particles in ParticleFBBC format (vector<ParticleFBBC>). At the output is returns a vector of vectors
of pairs <time, Id>, where Id is particle index number Id (see ParticleFBBC structure) and time is the time of particle
contact with the i-th MCP plate of FBBC.

### Private variables and methods of FBBCDetector class
- **vector<double> plates_distances** is the vector of distances of MCP plates from the detector center;
- **size_t plates_number** is the number of plates (size of plates_distances vector);
- **double r_in** is inner (small) radius of every MCP plate;
- **double r_out** is outer (big) radius of every MCP plate;
- **size_t rad_sec_num** is the number of radial sector separations;
- **size_t ang_sec_num** is the number of angular sector separations;
- **double efficiency** is the efficiency of every MCP plate (a probability
to registrate particle which passed through MCP);
- **double time_prec** is the time precision of every sector of every MCP plate;
- **vector<ParticleFBBC> particles** is the vector of particles in considered event 
(converted to ParticleFBBC format);

- **vector<PartTime> PassThrowMCP(const size_t mcp_num)** is the method which "simulates" passage of particles through 
MCP with a nubmer #mcp_num and returns vector of times of particles passed and registrated through the plate;
- **vector<vector<PartTime>> PassThrowDetector()** return the vector of the results of PassThrowMCP method 
calcuations for every MCP plate.

### Public methods of FBBCDetector class
- **const size_t GetPlatesNumber() const** returns plates_number;
- **const double GetEfficiency() const** returns efficiency;
- **const double GetTimePrecision() const** returns time_prec;
- **const double GetRInner() const** returns r_in;
- **const double GetROuter() const** returns r_out;
- **const size_t GetRadialSecNumber() const** returns rad_sec_num;
- **const size_t GetAngularSecNumber() const** returns ang_sec_num;
- **const vector<double> GetPlatesDistances() const** returns plates_distances vector;
- **const vector<vector<double>> GetPlatesPseudorapidity() const** converts plates_distances, r_in, 
r_out and Z coordinate of event into pseudorapidity interavals (using function RZtoEta) covered by 
every MCP plate;

- **void SetParticlesFBBC(const vector<ParticleFBBC> parts)** initializes "vector<ParticleFBBC> particles" from "parts" vector;
- **vector<vector<PartTime>> GetOutputVector()** returns the result of PassThrowDetector() evaluation.


## You can find an example of **fbbc-lib** usage in "example script" folder.
