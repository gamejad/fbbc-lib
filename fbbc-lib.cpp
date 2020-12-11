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
//#
vector<double> MCPlateSector::GetParticleTimesOutput() const{
	vector<double> result;
	result.push_back(part_times.front());
	for(const auto t : part_times){
		if(t > part_times.back()+time_prec)
			result.push_back(t);
	}
	return result;
}
//#
void AddParticleTime(const double time){
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	double number = distribution(generator);
	if(number <= efficiency){
		part_times.push_back(time);
	}
}
//.............................................................
void MCPlate::AddParticle(const Particle part){
		double pseudorapidity = atanh(part.Pz/part.P);
		if (pseudorapidity > RZtoEta(Rout, Z_coord-part.Z) && pseudorapidity < RZtoEta(Rin, Z_coord-part.Z)){
			particles.push_back(part);
		}
}
//#
vector<vector<double>> GetSectorsTimesOutput() {
	for(const auto part : particles){
		double pseudorapidity = atanh(part.Pz/part.P);
		double detect_time = (part.P/part.Pz)*(Z_coord-part.Z))/(c*c);
		bool is_matched = false;
		for(const auto r : mc_sectors){
			if(is_matched) break;
			for(const auto sec : r){
				if (pseudorapidity > RZtoEta(..., Z_coord-part.Z) && pseudorapidity < RZtoEta(..., Z_coord-part.Z)
						&& part.Phi > ... & part.Phi < ....){
					sec.AddParticleTime(detect_time);
					is_matched = true;
					break;
				}
			}
		}
	}

	vector<vector<double>> result(Rad_sec_num);
	for(size_t i = 0; i < result.size(); i++){
		result[i].resize(Phi_sec_num);
		for(size_t j = 0; j < result[i].size; j++){
			result[i][j] = mc_sectors[i][j].GetParticleTimesOutput();
		}
	}

	return result;
}
//#
const vector<double> GetTimesOutput(){
	auto v = GetSectorsTimesOutput();
	vector<double> output;
	for(auto vi : v){
		for(auto vij: vi){
			output.push_back(move(vij));
		}
	}
	return output;
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
