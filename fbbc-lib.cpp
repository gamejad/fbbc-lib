//fbbc-lib...
#include "fbbc-lib.h"
#include <exception>
#include <list>
#include <iostream>


double RZtoEta(double r, double z){
	double theta = atan(r/z);
	if (theta < 0) theta += PI;
	return -log(tan(theta/2));
}
////
double GetRandomNumber(const double first=0, const double last=1){
	std::random_device rd; // rd is a random seed
	std::default_random_engine generator(rd());
	std::uniform_real_distribution<double> distribution(first, last);
	return distribution(generator);
}
//--------------------------------------

const size_t FBBCDetector::GetPlatesNumber() const {
	return plates_number;
}
////
const double FBBCDetector::GetEfficiency() const{
	return efficiency;
}
////
const double FBBCDetector::GetTimePrecision() const{
	return time_prec;
}
////
const double FBBCDetector::GetRInner() const{
	return r_in;
}
////
const double FBBCDetector::GetROuter() const{
	return r_out;
}
////
const size_t FBBCDetector::GetRadialSecNumber() const{
	return rad_sec_num;
}
////
const size_t FBBCDetector::GetAngularSecNumber() const{
	return ang_sec_num;
}
////
const vector<double> FBBCDetector::GetPlatesDistances() const{
	return plates_distances;
}
////
const vector<vector<double>> FBBCDetector::GetPlatesPseudorapidity() const{
	vector<vector<double>> result;
	for(size_t i = 0; i < plates_number; i++){
		result.push_back(
			{min(RZtoEta(r_in, plates_distances.at(i)), RZtoEta(r_out, plates_distances.at(i))),
			 max(RZtoEta(r_in, plates_distances.at(i)), RZtoEta(r_out, plates_distances.at(i)))}
		);
	}
	return result;
}
////
void FBBCDetector::SetParticlesFBBC(const vector<ParticleFBBC> parts){
	particles = parts;
}
////
vector<PartTime> FBBCDetector::PassThrowMCP(const size_t mcp_num){
	const double c = 0.299792458; // mm/ps
	vector<vector<vector<PartTime>>> sector_counts(rad_sec_num, vector<vector<PartTime>>(ang_sec_num, vector<PartTime>())) ;

	try{
		for(const auto part : particles)
		{
			double detect_dist  = plates_distances.at(mcp_num) - part.Z;
			double eta_part = atanh(part.Pz/part.P);
			double DR = (r_out-r_in)/rad_sec_num;
			double DPHI = 2*PI/ang_sec_num;
			if (eta_part > max(RZtoEta(r_in, detect_dist), RZtoEta(r_out, detect_dist))
			 || eta_part < min(RZtoEta(r_in, detect_dist), RZtoEta(r_out, detect_dist)) ) continue;
			for(size_t r = 0; r < rad_sec_num; r++)
			{
				if (eta_part > max(RZtoEta(r_in+r*DR, detect_dist), RZtoEta(r_in+(r+1)*DR, detect_dist))
				|| eta_part < min(RZtoEta(r_in+r*DR, detect_dist), RZtoEta(r_in+(r+1)*DR, detect_dist))  ) continue;
				for (size_t phi = 0; phi < ang_sec_num; phi++)
				{
					if (part.Phi > -PI+(phi+1)*DPHI || part.Phi < -PI+phi*DPHI ) continue;
					double number = GetRandomNumber(0, 1);
					if (number > efficiency) continue;
					double detect_time = (part.E/part.Pz)*detect_dist/(c*c);
					sector_counts[r][phi].push_back({detect_time, part.Id});
				}
			}
		}

		vector<vector<vector<PartTime>>> sector_counts_reduced(rad_sec_num, vector<vector<PartTime>>(ang_sec_num, vector<PartTime>())) ;

		for(size_t r = 0; r < rad_sec_num; r++)
		{
			for (size_t phi = 0; phi < ang_sec_num; phi++)
			{
				sort(sector_counts[r][phi].begin(), sector_counts[r][phi].end());
				for(auto p_t : sector_counts[r][phi]){
					if(sector_counts_reduced[r][phi].size() == 0){
						sector_counts_reduced[r][phi].push_back(move(p_t));
					} else {
						if(sector_counts_reduced[r][phi].back().first < p_t.first + time_prec){
							sector_counts_reduced[r][phi].push_back(move(p_t));
						} else continue;
					}
				}
			}
		}

		vector<PartTime> result;
		for(size_t r = 0; r < rad_sec_num; r++)
		{
			for (size_t phi = 0; phi < ang_sec_num; phi++)
			{
				for(auto p_t : sector_counts_reduced[r][phi]){
					result.push_back(p_t);
				}
			}
		}

		return result;

	} catch (const out_of_range& e) {
			cerr << "Error! There is no MCP with number " << mcp_num << endl;
	}
}
////
vector<vector<PartTime>> FBBCDetector::PassThrowDetector(){
	vector<vector<PartTime>> result;
	for(size_t i = 0; i < plates_number; i++){
		result.push_back(PassThrowMCP(i)); // проблема с пустыми детекторами - решить
	}
	return result;
}
////
vector<vector<PartTime>> FBBCDetector::GetOutputVector(){
	return PassThrowDetector();
}
