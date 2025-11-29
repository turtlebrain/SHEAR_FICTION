#pragma once
#include <iostream>
#include "Tensors.h"
#include "ProbabilisticMat.h"
#include <vector>
using namespace std;


class PlasticMat
{
public:
	PlasticMat() { }

	enum eStrengthType
	{
		VON_MISES,
		TRESCA,
		DRUCKER_PRAGER,
		MORH_COULOMB
	};

	void setphi(const double Phi) { m_phi.DeterministicScalar = Phi; }
	void setNvec(const vector<double>& N) { m_N_flow_vec.DeterministicTensor = N; }
	void setSigma(const vector<double>& sigma) { m_sigma.DeterministicTensor = sigma; }
	void setSigmaVolumetric(const vector<double>& sig_vol) { m_sigma_vol.DeterministicTensor = sig_vol; }
	void setSigmaDeviatoric(const vector<double>& sig_dev) { m_sigma_dev.DeterministicTensor = sig_dev; }
	void setTotalStrain(const vector<double>& epsilon) { m_eps_tot.DeterministicTensor = epsilon; }
	void setElasticStrain(const vector<double>& epsilon_e) { m_eps_e.DeterministicTensor = epsilon_e; }
	void setPlasticStrain(const double epsilon_p) { m_eps_p.DeterministicScalar = epsilon_p; }

	vector<double> getSigma() const { return m_sigma.DeterministicTensor; }
	vector<double> getSigmaVol() const { return m_sigma_vol.DeterministicTensor; }
	vector<double> getSigmaDev() const { return m_sigma_dev.DeterministicTensor; }
	vector<double> getTotalStrain() const { return m_eps_tot.DeterministicTensor; }
	vector<double> getElasticStrain() const { return m_eps_e.DeterministicTensor; }
	double getPlasticStrain() const { return m_eps_p.DeterministicScalar; }
	StochasticProcess GetStochasticStress() const { return m_sigma; }

private:
	StochasticProcess m_sigma;			//stress 		
	StochasticProcess m_sigma_vol;		//volumetric stress
	StochasticProcess m_sigma_dev;		//deviatoric stress
	StochasticProcess m_eps_tot;		//total strain
	StochasticProcess m_eps_e;			//elastic strain
	StochasticProcess m_eps_p;			//plastic strain
	StochasticProcess m_N_flow_vec;		//flow vector
	StochasticProcess m_phi;			//Yield function

private:
	eStrengthType m_strength_type;
	double m_mat_prop[9];

};


class UMAT_VonMisces :public PlasticMat
{
public:
	UMAT_VonMisces(double E, double v, double sigy, double H);

	enum class eHardeningType
	{
		LINEAR_ISOTROPIC
	};

	vector<vector<double>> getLinearElasticTensor(const double E, const double v) const;
	double getE() const { return m_young_modulus.mean; }
	double getv() const { return m_poisson_ratio.mean; }
	double getK() const { return m_bulk_modulus.mean; }
	double getG() const { return m_shear_modulus.mean; }
	double getH() const { return m_hardening_slope.mean; }

	double getSigmaY0() const { return m_sigma_y_init.mean; }

	void associative_flow_rule(const vector<double>& sigma);
	void hardening_law(eHardeningType H_type, double h, double sigma_y_0);
	double yield_function(const vector<double>& sigma,const double sigma_y);

	void EulerBackward(const vector<double>& delta_eps, const vector<double>& eps_tot_pstep, const vector<double>& eps_e_pstep, const vector<double>& sigma_pstep, const double eps_p_pstep);
	
	//-------------------------------------Probabilistic Constitutive Law----------------------------------------------/
public:
	vector<vector<double>> GetAdvectiveTensor();
	vector<vector<double>> GetDiffusiveTensor();

	void computeFPKE(vector<double> sig_0);

private:
	ProbabilisticMat m_young_modulus;
	ProbabilisticMat m_poisson_ratio;
	ProbabilisticMat m_shear_modulus;
	ProbabilisticMat m_bulk_modulus;
	ProbabilisticMat m_hardening_slope;
	ProbabilisticMat m_sigma_y_init;

};

