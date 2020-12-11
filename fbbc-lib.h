//FBBC-lib...
#pragma once

#include <vector>
#include <pair>
#include <random>
#include <cmath>

#define PI = 3.14159265359;

double RZtoEta();

class MCSector {
public:
	MCSector(pair<double> phi_pair, pair<double> r_pair, double eff, double t_prec):
		phi_limit(phi_pair),
		r_limit(r_pair),
		efficiency(eff),
		time_prec(t_prec) {}

	void SetParticleTimes(const vector<double> times);

	vector<double> GetParticleTimes() const;
private:
	pair<double> phi_limit; //{phi1, phi2}
	pair<double> r_limit;   // {r1, r2}
	double efficiency = 1; // 100%
	double time_prec = 0; //ps
	vector<double> part_times;
}


class MCPlate {
public:
	MCPlate(double rin, double rout, double rad_num, double phi_num, double eff, double t_prec):
		Rin(rin),
		Rout(rout),
		Rad_sec_num(rad_num),
		Phi_sec_num(phi_num) {
		double delta_r = (rout-rin)/rad_num;
		double delta_phi = 2*PI/phi_num;
		mc_sectors.resize(rad_num);
		for(int i = 0; i < rad_num; i++){
			pair<double> rads = {rin + i*delta_r, rin+(i+1)*delta_r};
			for(int j = 0; j < phi_num; j++){
				pair<double> phis = {j*delta_phi, (j+1)*delta_phi};
				MCSector plate(rads, phis, eff, t_prec);
				mc_sectors[i].push_back(plate);
			}
		}
	}

private:
	vector<vector<MCSector>> mc_sectors;
	double Rin, Rout; // inner and outter radiuses, in mm
	size_t Rad_sec_num = 1;
	size_t Phi_sec_num = 1;
}

class DetectorFacility {
public:


private:
	vector<pair<MCPlate, double> > mcplates;
	int plate_num = 0;
	double Z_coordinate = 0;
}




