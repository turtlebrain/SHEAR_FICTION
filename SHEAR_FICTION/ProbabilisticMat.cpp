#include "ProbabilisticMat.h"
#include "PlasticMat.h"
#include <iostream>
#include "Tensors.h"
#include <vector>
#include <math.h>
using namespace std;

#define PI 3.14159265358979

vector<double> ProbabilisticMat::GetInitialDistribution(eDistribution dist, const double mean, const double std_dev, const int vec_size) const
{
	vector<double> pdf(vec_size, 0);
	vector<double> r_prop(vec_size, 0);
	double r_min = mean - 5 * std_dev;
	double r_max = mean + 5 * std_dev;
	double r_delta = (r_max - r_min) / ((double)vec_size - 1);

	r_prop[0] = r_min;
	for (int i = 1; i < vec_size; i++) {
		r_prop[i] = r_prop[i - 1] + r_delta;
	}

	switch (dist) {
		case eDistribution::NORMAL:
			vector<double> pdf(vec_size, 0);
			for (int i = 0; i < vec_size; i++) { 
				pdf[i] = 1 / (std_dev * sqrt(2 * PI))*exp(-pow((r_prop[i]-mean),2)/(2*pow(std_dev,2))); 
			}
			return pdf;
			break;
	}

	
}

