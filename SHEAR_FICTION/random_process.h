#pragma once
#include <iostream>
#include <vector>
#include "tensor_o2.h"
using namespace std;

enum class eDistribution {
	NO_DISTRIBUTION,
	NORMAL_DISTRIBUTION,
	LOGNORMAL_DISTRIBUTION,
};

class RandomProcess
{
public:
	RandomProcess();
	RandomProcess(double mean, eDistribution distribution = eDistribution::NO_DISTRIBUTION);
	RandomProcess(double mean, double stdDev, eDistribution dist);
	~RandomProcess() {}

public:
	vector<double> getInitialDistribution(const int vec_size) const;
	double getMean() const { return m_mean; }
	double getStdDev() const { return m_std_dev; }
	double getVar() const { return m_var; }
	eDistribution getDistribution() const { return m_distribution; }

protected:
	double m_mean;
	double m_std_dev;
	double m_var;
	eDistribution m_distribution;
};

class StochasticProcess :public RandomProcess
{
public:
	StochasticProcess() {}
	~StochasticProcess() {}

public:
	void SetTimeEvolutionPDF(vector<vector<double>> evolutionaryPDF) { m_3Dpdf = evolutionaryPDF; }

	vector<vector<double>> getTimeEvolutionPDF() { return m_3Dpdf; }
	vector<double> getInstantPDF(int t) { return m_3Dpdf[t]; }

public:
	TensorO2 m_deterministicTensor;
	double m_deterministicScalar = 0.0;

private:
	vector<vector<double>> m_3Dpdf;

};