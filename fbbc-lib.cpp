//fbbc-lib...
#include "fbbc-lib.h"


double RZtoEta(double r, double z){
	double theta = atan(r/z);
	return -log(tan(theta/2));
}
//--------------------------------------

void MCPlateSector::SetParticleTimes(const vector<double> times){
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	for(const auto t : times){
		double number = distribution(generator);
		if(number <= efficiency){
			part_times.push_back(t);
		} else continue;
	}
	part_times.sort();
}

vector<double> MCPlateSector::GetParticleTimes() const{
	vector<double> result;
	result.push_back(part_times.front());
	for(const auto t : part_times){
		if(t > part_times.back()+time_prec)
			result.push_back(t);
	}
	return result;
}
//.............................................................
void MCPlate::AddParticle(const Particle part){
		double pseudorapidity = atanh(part.Pz/part.P);
		if (pseudorapidity > RZtoEta(Rout, Z_coord-part.Z) && pseudorapidity < RZtoEta(Rin, Z_coord-part.Z)){
			particles.push_back(part);
		}
}
//..............................................................

void DetectorFacility::SetParticles(const vector<Particle> parts){
	particles = parts;
	for(const auto part : parts){
		for(auto [z, plate] : mcplates){
			plate.AddParticle();
		}
	}
}
