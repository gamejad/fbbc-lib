//fbbc-lib...
#include "fbbc-lib.h"


double RZtoEta(double r, double z){
	double theta = atan(r/z);
	return -log(tan(theta/2));
}
//--------------------------------------

vector<double> MCPlateSector::GetTimesOutput() const{
	vector<double> result;
	result.push_back(part_times.front());
	for(const auto t : part_times){
		if(t > part_times.back()+time_prec)
			result.push_back(t);
	}
	return result;
}
//#
void MCPlateSector::AddParticle(const Particle part){
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	double number = distribution(generator);
	if(number <= efficiency){
		particles.push_back(part);
		double t = ((part.P/part.Pz)*(Z_coord-part.Z))/(c*c);
		part_times.push_back(t);
	}
}
//#
const pair<double, double> MCPlateSector::GetRLimits() const {
	return r_limit;
}
//#
const pair<double, double> MCPlateSector::GetPhiLimits() const{
	return phi_limit;
}
//#
const pair<double, double> MCPlateSector::GetEtaLimits() const {
	return {RZtoEta(r_limit.second, Z_coord), RZtoEta(r_limit.first, Z_coord)};
}
//.............................................................
void MCPlate::AddParticle(const Particle part){
		double pseudorapidity = atanh(part.Pz/part.P);
		if (pseudorapidity > RZtoEta(Rout, Z_coord-part.Z) && pseudorapidity < RZtoEta(Rin, Z_coord-part.Z)){
			particles.push_back(part);
		}
		//double detect_time = (part.P/part.Pz)*(Z_coord-part.Z))/(c*c);
		bool is_matched = false;
		for(auto r : mc_sectors){
			if(is_matched) break;
			for(auto sec : r){
				if (pseudorapidity > RZtoEta(min(sec.GetRLimits().first, sec.GetRLimits().second), Z_coord-part.Z) &&
						pseudorapidity < RZtoEta(max(sec.GetRLimits().first, sec.GetRLimits().second), Z_coord-part.Z) &&
						part.Phi > min(sec.GetPhiLimits().first, sec.GetPhiLimits().second) &&
						part.Phi < max(sec.GetPhiLimits().first, sec.GetPhiLimits().second)){
					sec.AddParticle(part);
					is_matched = true;
					break;
				}
			}
		}
}
//#
vector<vector<vector<double>>> MCPlate::GetSectorsTimesOutput() {
	vector<vector<vector<double>>> result(Rad_sec_num);
	for(size_t i = 0; i < result.size(); i++){
		result[i].resize(Phi_sec_num);
		for(size_t j = 0; j < result[i].size(); j++){
			result[i][j] = mc_sectors[i][j].GetTimesOutput();
		}
	}
	return result;
}
//#
const vector<double> MCPlate::GetTimesOutput(){
	auto v = GetSectorsTimesOutput();
	vector<double> output;
	for(auto vi : v){
		for(auto vij: vi){
			for(auto i : vij){
				output.push_back(move(i));
			}
		}
	}
	return output;
}
//#
const pair<double, double> MCPlate::GetEtaLimits() const {
	return {RZtoEta(Rout, Z_coord), RZtoEta(Rin, Z_coord)};
}
//#
const pair<double, double> MCPlate::GetRLimits() const {
	return {Rin, Rout};
}
//#
const double MCPlate::GetTimePrecision() const{
	return time_prec;
}
//..............................................................

void DetectorFacility::SetParticles(const vector<Particle> parts){
	particles = parts;
	for(const auto part : parts){
		for(auto [z, plate] : mcplates){
			plate.AddParticle(part);
		}
	}
}
//#
vector<vector<double>> DetectorFacility::GetPlatesTimes() {
	vector<vector<double>> result;
	for( auto [z, plate] : mcplates){
		result.push_back(plate.GetTimesOutput());
	}
	return result;
}
//#
vector<double> DetectorFacility::GetPlatesDistances() const{
	vector<double> result;
	for(const auto [z, plate] : mcplates){
		result.push_back(z);
	}
	return result;
}
//#
vector<pair<double,double>> DetectorFacility::GetPlatesPseudorapidities() const{
	vector<pair<double,double>> result;
	for(const auto [z, plate] : mcplates){
		result.push_back(plate.GetEtaLimits());
	}
	return result;
}
//#
vector<pair<double,double>> DetectorFacility::GetPlatesRadiuses() const{
	vector<pair<double,double>> result;
	for(const auto [z, plate] : mcplates){
		result.push_back(plate.GetRLimits());
	}
	return result;
}
//#
vector<double> DetectorFacility::GetTimePrecisions() const{
	vector<double> result;
	for(const auto [z, plate] : mcplates){
		result.push_back(plate.GetTimePrecision());
	}
	return result;
}
