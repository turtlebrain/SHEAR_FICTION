#include "random_process.h"
#include "plastic_material.h"
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

#define PI 3.14159265358979

RandomProcess::RandomProcess()
{
	m_mean = 0;
	m_std_dev = 0;
	m_var = 0;
	m_distribution = eDistribution::NO_DISTRIBUTION;
}

RandomProcess::RandomProcess(double mean, eDistribution distribution)
{
	m_mean = mean;
	m_std_dev = 0;
	m_var = 0;
	m_distribution = distribution;
}

RandomProcess::RandomProcess(double mean, double stdDev, eDistribution distribution)
{
	m_mean = mean;
	m_std_dev = stdDev;
	m_var = stdDev * stdDev;
	m_distribution = distribution;
}

vector<double> RandomProcess::getInitialDistribution(const int vec_size) const
{
	vector<double> pdf(vec_size, 0);
	vector<double> r_prop(vec_size, 0);
	double r_min = m_mean - 5 * m_std_dev;
	double r_max = m_mean + 5 * m_std_dev;
	double r_delta = (r_max - r_min) / ((double)vec_size - 1);

	r_prop[0] = r_min;
	for (int i = 1; i < vec_size; i++) {
		r_prop[i] = r_prop[i - 1] + r_delta;
	}

	switch (m_distribution) {
		case eDistribution::NORMAL_DISTRIBUTION:
			vector<double> pdf(vec_size, 0);
			for (int i = 0; i < vec_size; i++) { 
				pdf[i] = 1 / (m_std_dev * sqrt(2 * PI))*exp(-pow((r_prop[i]-m_mean),2)/(2*pow(m_std_dev,2))); 
			}
			return pdf;
			break;
	}
}

