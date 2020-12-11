//fbbc-lib...
#include "fbbc-lib.h"


double RZtoEta(double r, double z){
	double theta = atan(r/z);
	return -log(tan(theta/2));
}
//--------------------------------------

void MCSector::SetParticleTimes(const vector<double> times){
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

vector<double> MCSector::GetParticleTimes() const{
	vector<double> result;
	result.push_back(part_times.front());
	for(const auto t : part_times){
		if(t > part_times.back()+time_prec)
			result.push_back(t);
	}
	return result;
}



