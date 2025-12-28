#pragma once
#include <iostream>
#include "random_process.h"
#include <vector>
#include "tensor_o2.h"
#include "tensor_o4.h"
using namespace std;

enum class eStrengthType
{
	VON_MISES,
	TRESCA,
	DRUCKER_PRAGER,
	MORH_COULOMB
};

enum class eHardeningType
{
	LINEAR_ISOTROPIC
};

class PlasticMat
{
public:
	PlasticMat() { }

public:
	void setPhi(const double phi) { m_phi.m_deterministicScalar = phi; }
	void setNvec(const TensorO2& n) { m_N_flow_vec.m_deterministicTensor = n; }
	void setSigma(const TensorO2& sigma) { m_sigma.m_deterministicTensor = sigma; }
	void setSigmaVolumetric(const TensorO2& sig_vol) { m_sigma_vol.m_deterministicTensor = sig_vol; }
	void setSigmaDeviatoric(const TensorO2& sig_dev) { m_sigma_dev.m_deterministicTensor = sig_dev; }
	void setTotalStrain(const TensorO2& epsilon) { m_eps_tot.m_deterministicTensor = epsilon; }
	void setElasticStrain(const TensorO2& epsilon_e) { m_eps_e.m_deterministicTensor = epsilon_e; }
	void setPlasticStrain(const double epsilon_p) { m_eps_p.m_deterministicScalar = epsilon_p; }

	TensorO2 getSigma() const { return m_sigma.m_deterministicTensor; }
	TensorO2 getSigmaVol() const { return m_sigma_vol.m_deterministicTensor; }
	TensorO2 getSigmaDev() const { return m_sigma_dev.m_deterministicTensor; }
	TensorO2 getTotalStrain() const { return m_eps_tot.m_deterministicTensor; }
	TensorO2 getElasticStrain() const { return m_eps_e.m_deterministicTensor; }
	double getPlasticStrain() const { return m_eps_p.m_deterministicScalar; }
	StochasticProcess getStochasticStress() const { return m_sigma; }

protected:
	StochasticProcess m_sigma;			//stress 		
	StochasticProcess m_sigma_vol;		//volumetric stress
	StochasticProcess m_sigma_dev;		//deviatoric stress
	StochasticProcess m_eps_tot;		//total strain
	StochasticProcess m_eps_e;			//elastic strain
	StochasticProcess m_eps_p;			//plastic strain
	StochasticProcess m_N_flow_vec;		//flow vector
	StochasticProcess m_phi;			//Yield function

protected:
	eStrengthType m_strength_type;
	double m_mat_prop[9];

};


class UMAT_VonMisces :public PlasticMat
{
public:
	UMAT_VonMisces(RandomProcess youngs_mod, RandomProcess poisson_ratio, RandomProcess sigy, RandomProcess h);

	TensorO4 getLinearElasticTensor(const double E, const double v) const;
	double getE() const { return m_young_modulus.getMean(); }
	double getv() const { return m_poisson_ratio.getMean(); }
	double getK() const { return m_bulk_modulus.getMean(); }
	double getG() const { return m_shear_modulus.getMean(); }
	double getH() const { return m_hardening_slope.getMean(); }

	double getSigmaY0() const { return m_sigma_y_init.getMean(); }

	void associative_flow_rule(const TensorO2& sigma);
	void hardening_law(eHardeningType h_type, double h, double sigma_y_0);
	double yield_function(const TensorO2& sigma, const double sigma_y);

	void EulerBackward(const TensorO2& delta_eps, const TensorO2& eps_tot_pstep, const TensorO2& eps_e_pstep, const TensorO2& sigma_pstep, const double eps_p_pstep);
	
	//-------------------------------------Probabilistic Constitutive Law----------------------------------------------/
public:
	TensorO4 GetAdvectiveTensor();
	TensorO4 GetDiffusiveTensor();

	void computeFPKE(TensorO2 sig_0);

protected:
	RandomProcess m_young_modulus;
	RandomProcess m_poisson_ratio;
	RandomProcess m_shear_modulus;
	RandomProcess m_bulk_modulus;
	RandomProcess m_hardening_slope;
	RandomProcess m_sigma_y_init;

};

