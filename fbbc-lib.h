//FBBC-lib...
#pragma once

#include <vector>
#include <map>
#include <random>
#include <cmath>
#include <algorithm>

using namespace std;

const double PI = 3.14159265359;
const double c = 0.299792458; // mm/ps

double RZtoEta();
//........................................................................
struct Particle{
	int id;
	double E;
	double P;
	double Pz;
	double Phi;

	double Z;
};

//........................................................................
class MCPlateSector {
public:
	MCPlateSector(double z, pair<double,double> phi_pair, pair<double,double> r_pair, double eff, double t_prec):
		Z_coord(z),
		phi_limit(phi_pair),
		r_limit(r_pair),
		efficiency(eff),
		time_prec(t_prec) {}

	void AddParticle(const Particle part);
	vector<double> GetTimesOutput() const;
	const pair<double,double> GetRLimits() const;
	const pair<double,double> GetPhiLimits() const;
	const pair<double, double> GetEtaLimits() const;

private:
	double Z_coord;
	pair<double,double> phi_limit; //{phi1, phi2}
	pair<double,double> r_limit;   // {r1, r2}
	vector<double> part_times;
	vector<Particle> particles;
	double efficiency = 1; // 1=100%
	double time_prec = 0; //ps
};

//........................................................................
class MCPlate {
public:
	MCPlate(double z, double rin, double rout, double rad_num, double phi_num, double eff, double t_prec):
		Z_coord(z),
		Rin(rin),
		Rout(rout),
		Rad_sec_num(rad_num),
		Phi_sec_num(phi_num),
		efficiency(eff),
		time_prec(t_prec) {
			double delta_r = (rout-rin)/rad_num;
			double delta_phi = 2*PI/phi_num;
			mc_sectors.resize(rad_num);
			for(int i = 0; i < rad_num; i++){
				pair<double,double> rads = {rin + i*delta_r, rin+(i+1)*delta_r};
				for(int j = 0; j < phi_num; j++){
					pair<double,double> phis = {j*delta_phi, (j+1)*delta_phi};
					MCPlateSector plate_sec(z, rads, phis, efficiency, time_prec);
					mc_sectors[i].push_back(plate_sec);
				}
			}
	}

	void AddParticle(const Particle part);
	vector<vector<vector<double>>> GetSectorsTimesOutput();
	const vector<double> GetTimesOutput();
	const pair<double, double> GetEtaLimits() const;
	const pair<double,double> GetRLimits() const;
	const double GetTimePrecision() const;

private:
	vector<vector<MCPlateSector>> mc_sectors;
	double Z_coord;
	double Rin;
	double Rout; // inner and outter radiuses, in mm
	size_t Rad_sec_num = 1;
	size_t Phi_sec_num = 1;
  vector<Particle> particles;
	double efficiency = 1; // 1=100%
	double time_prec = 0; //ps
};

//........................................................................
class DetectorFacility {
public:
    DetectorFacility(vector<double> plt_coords, double rin, double rout, double rad_num, double phi_num, double eff, double t_prec) {
            plate_num = plt_coords.size();
						for(int i = 0; i < plate_num; i++){
							MCPlate plate(plt_coords.at(i), rin, rout, rad_num, phi_num, eff, t_prec);
							mcplates.insert({plt_coords.at(i), plate});
						}
        }

		void SetParticles(const vector<Particle> parts);
		vector<vector<double>> GetPlatesTimes() ;
		vector<double> GetPlatesDistances() const;
		vector<pair<double,double>> GetPlatesPseudorapidities() const;
		vector<pair<double,double>> GetPlatesRadiuses() const;
		vector<double> GetTimePrecisions() const;


private:
	map<double, MCPlate > mcplates;
	int plate_num = 0;
	vector<Particle> particles;
};
