#pragma once
#include <iostream>
#include "TensorOps.h"
#include <vector>
using namespace std;

class ProbabilisticMat
{
public:
	ProbabilisticMat() { }

	enum class eDistribution {
		NORMAL,
		LOGNORMAL,
	};

	vector<double> GetInitialDistribution(eDistribution dist, const double mean,const double std_dev, const int vec_size) const;

public:
	double mean;
	double std_dev;
	double var;
	eDistribution distribution;

};

class StochasticProcess :public ProbabilisticMat
{
public:
	StochasticProcess() {}

	void SetTimeEvolutionPDF(vector<vector<double>> EvolutionaryPDF) { m_3Dpdf = EvolutionaryPDF; }

	vector<vector<double>> GetTimeEvolutionPDF() { return m_3Dpdf; }
	vector<double> GetInstantPDF(int t) { return m_3Dpdf[t]; }
	double GetMean() { return mean; }

public:
	vector<double> DeterministicTensor;
	double DeterministicScalar;

private:
	vector<vector<double>> m_3Dpdf;

};