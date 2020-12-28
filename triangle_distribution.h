#pragma once

#include <random>

double TriangleRandom(const double xmin, const double xmax, const double x0){
    if(xmin>x0 ||xmin > xmax || x0 > xmax){
        cerr <<  "Error! (xmin < x0 < xmax) is not satisfied!" << endl;
    }
    std::random_device rd; // rd is a random seed
	std::default_random_engine generator(rd());
    double l = xmax-xmin;
    double h = 2/l;
    
	std::uniform_real_distribution<double> distributionX(xmin, xmax);
    std::uniform_real_distribution<double> distributionY(0, h);
	double ksiX = distributionX(generator);
    double ksiY = distributionY(generator);
    
    if(ksiX <= x0 && ksiY <= h*(ksiX-xmin)/(x0-xmin)){
        return ksiX;
    } else if(ksiX > x0 && ksiY <= h*(ksiX-xmax)/(x0-xmax)){
        return ksiX;
    } else {
        return TriangleRandom(xmin, xmax, x0);
    }
    
}
