//FBBC-lib...
#pragma once

#include <vector>
#include <set>
#include <random>
#include <cmath>
#include <algorithm>

const double PI = 3.14159265359;
const double SPEED_OF_LIGHT = 0.299792458; // mm/ps

double RZtoEta();
double GetRandomNumber(const double first, const double last);
//........................................................................
struct ParticleFBBC{
	int Id;
	double E; // GeV
	double P; // GeV/c
	double Pz; // GeV/c
	double Phi; // radian
	double Z; // mm
};

struct PartTime{
   double Time;
   int Id;
   size_t Rad_sec;
   size_t Ang_sec;
};
//........................................................................

class FBBCDetector{
public:
	FBBCDetector(std::vector<double> mcp_dists, double rin, double rout, size_t rad_num, size_t ang_num, double eff, double t_pr):
		plates_distances(mcp_dists),
		plates_number(mcp_dists.size()),
		r_in(rin),
		r_out(rout),
		rad_sec_num(rad_num),
		ang_sec_num(ang_num),
		efficiency(eff),
		time_prec(t_pr)
		{
			sort(plates_distances.begin(), plates_distances.end());
		}

		const size_t GetPlatesNumber() const;
		const double GetEfficiency() const;
		const double GetTimePrecision() const;
		const double GetRInner() const;
		const double GetROuter() const;
		const size_t GetRadialSecNumber() const;
		const size_t GetAngularSecNumber() const;
		const std::vector<double> GetPlatesDistances() const;
		const std::vector<std::vector<double>> GetPlatesPseudorapidity() const;

		void SetParticlesFBBC(const std::vector<ParticleFBBC> parts);
		std::vector<std::vector<PartTime>> GetOutputVector();

private:
	std::vector<double> plates_distances; // MCPs' distances from center (sorted)
	size_t plates_number; //number of MCPs
	double r_in; // mm, inner radius
	double r_out; //mm, outer radius
	size_t rad_sec_num;
	size_t ang_sec_num;
	double efficiency; // efficiency of MCP, 1 = 100%
	double time_prec; // ps,

	std::vector<ParticleFBBC> particles; // particles in event

	std::vector<PartTime> PassThrowMCP(const size_t mcp_num); //return vector<{time, PartID}>
	std::vector<std::vector<PartTime>> PassThrowDetector(); // return vector<PassThrowMCP>

};
