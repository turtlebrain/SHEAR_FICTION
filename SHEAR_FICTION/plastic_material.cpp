#include "plastic_material.h"
#include "tensor_operations.h"
#include <iostream>
#include <vector>
using namespace std;

UMAT_VonMisces::UMAT_VonMisces(RandomProcess E, RandomProcess v, RandomProcess sigy, RandomProcess H)
{
	m_young_modulus = E;
	m_poisson_ratio = v;
	m_sigma_y_init = sigy;
	m_hardening_slope = H;
	m_bulk_modulus = RandomProcess(E.getMean() / (3 * (1 - 2 * v.getMean())));
	m_shear_modulus = RandomProcess(E.getMean() / (2 * (1 + v.getMean())));
}

void UMAT_VonMisces::EulerBackward(const TensorO2& delta_eps, const TensorO2& eps_tot_pstep, const TensorO2& eps_e_pstep, const TensorO2& sigma_pstep, const double eps_p_pstep) 
{
	//These locals should probably be accessed from the class
	//Present values should be passed as arguments
	// Instantiation below can later become a SetMaterial function with input material model.
	//UMAT_VonMisces VM_Mat(E, v);
	TensorO2 eps_tot_fstep;											//Total strain at time t_n+1
	TensorO2 eps_e_fstep, eps_e_trial_fstep;						//elastic strain at time t_n+1 
	double eps_p_fstep, eps_p_trial_fstep;							//plastic strain at time t_n+1 
	TensorO2 p_fstep, s_fstep;										//volumetric and deviatoric stress at time t_n+1
	TensorO2 sigma_fstep;
 
	TensorO2 sigma_trial_fstep, p_trial_fstep, s_trial_fstep;
	TensorO2 eps_e_vol_trial, eps_e_dev_trial;

	double q_trial_fstep;											//trial Von-Mises effective stress
	double sigma_y;													//yield stress
	double sigma_y_init;											//Initial yield stress
	double H;														//Hardening slope
	double phi;														//yield function
	double delta_gamma = 0;											//plastic multiplier
	double resid_d;													//Residual derivative

	eps_p_trial_fstep = eps_p_pstep;
	TensorO2 eps_tot_pstep_copy = TensorO2(eps_tot_pstep);
	TensorO2 delta_eps_copy = TensorO2(delta_eps);
	for (int i = 0; i < 9; i++) {
		eps_tot_fstep.t2(i) = eps_tot_pstep_copy.t2(i) + delta_eps_copy.t2(i);
		if (i == 0) eps_e_trial_fstep.t2(i) = eps_tot_fstep.t2(i) - eps_p_trial_fstep;		
	}

	//Elastic trial state/Elastic predictor step
	//D_elastic = VM_Mat.getLinearElasticTensor(E, v);
	tie(eps_e_vol_trial, eps_e_dev_trial) = eps_e_trial_fstep.volDecomp();

	double K = getK();
	double G = getG();
	for (int i = 0; i < 9; i++) {
		//eps_e_trial_fstep[i] = eps_e_pstep[i] + delta_eps[i];
		p_trial_fstep.t2(i) =  K * eps_e_vol_trial.t2(i); 
		s_trial_fstep.t2(i) = 2 * G * eps_e_dev_trial.t2(i);
		sigma_trial_fstep.t2(i) = s_trial_fstep.t2(i) + p_trial_fstep.t2(i);
	}
	q_trial_fstep = sqrt((3.0 * (TensorOperations::contrac2X(s_trial_fstep, s_trial_fstep))/2.0));
	eps_p_trial_fstep = eps_p_pstep + delta_gamma;
	H = getH();
	sigma_y_init = getSigmaY0();
	hardening_law(eHardeningType::LINEAR_ISOTROPIC, H, sigma_y_init);	//H and sigma_y_0 are hardcoded for now
	sigma_y = sigma_y_init + H * eps_p_trial_fstep;
	phi = q_trial_fstep - sigma_y;//yield_function(sigma_trial_fstep, sigma_y);
	if (phi <= 0) {
		eps_e_fstep = eps_e_trial_fstep;
		eps_p_fstep = eps_p_trial_fstep;
		p_fstep = p_trial_fstep;
		s_fstep = s_trial_fstep;
		sigma_fstep = sigma_trial_fstep;
	}
	else {
		//Newton-Raphson algorithm
		//find hardening slope
		int iteration_counter = 0;
		double tol = 10e-5;
		double delta_gamma_new;
		delta_gamma = 0.0;


		while (iteration_counter <= 20 && abs(phi) >= tol) {

			phi = q_trial_fstep - 3.0 * G * delta_gamma - (sigma_y_init + (eps_p_pstep + delta_gamma) * H);

			resid_d = -3.0 * G - H;

			delta_gamma_new = delta_gamma - phi / resid_d;

			delta_gamma = delta_gamma_new;

			iteration_counter = iteration_counter + 1;
		}
		p_fstep = p_trial_fstep;
		eps_p_fstep = eps_p_pstep + delta_gamma;
		for (int i = 0; i < 9; i++) {
			s_fstep.t2(i) = (1.0 - (delta_gamma * 3.0 * G) / q_trial_fstep) * s_trial_fstep.t2(i);
			sigma_fstep.t2(i) = s_fstep.t2(i) + p_fstep.t2(i);
			eps_e_fstep.t2(i) = (1 / (2.0 * G)) * s_fstep.t2(i) +  eps_e_vol_trial.t2(i); // (1 / (2.0 * G)) * s_fstep[i] +  (1 / 3.0) * eps_e_vol_trial[i]
		}
	}

	setTotalStrain(eps_tot_fstep);
	setElasticStrain(eps_e_fstep);
	setPlasticStrain(eps_p_fstep);
	setSigmaVolumetric(p_fstep);
	setSigmaDeviatoric(s_fstep);
	setSigma(sigma_fstep);
}

double UMAT_VonMisces::yield_function(const TensorO2& sigma, const double sigma_y)
{
	TensorO2 s_ij, s_ji;
	TensorO2 sigmaCopy = TensorO2(sigma);
	double sigma_mean, J2, phi;
	int ind_ij;
	sigma_mean = (1 / 3.0) * (sigmaCopy.trace());
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			if (i == j) s_ij.t2(i,j) = sigmaCopy.t2(i,j) - sigma_mean;
			else s_ij.t2(i,j) = sigmaCopy.t2(i,j);
		}
	}
	s_ji = s_ij.transpose();
	J2 = (1 / 2.0) * (TensorOperations::contrac2X(s_ij, s_ji));
	
	phi = sqrt(3.0 * J2) - sigma_y;
	return phi;
}

void UMAT_VonMisces::associative_flow_rule(const TensorO2& sigma)
{
	TensorO2 s_ij, N_flow_vector;
	double sigma_mean;
	int ind_ij;
	TensorO2 sigmaCopy = TensorO2(sigma);
	sigma_mean = (1 / 3) * (sigmaCopy.trace());
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			if (i == j) s_ij.t2(i,j) = sigmaCopy.t2(i,j) - sigma_mean;
			else s_ij.t2(i,j) = sigmaCopy.t2(i,j);
		}
	}

	double mod_s = s_ij.normL2();
	for (int i = 0; i < 9; i++) N_flow_vector.t2(i) = (sqrt(3 / 2)) * s_ij.t2(i) / mod_s;

	setNvec(N_flow_vector);
}

TensorO4 UMAT_VonMisces::getLinearElasticTensor(const double E, const double v) const
{
	UnitTensorO4 I_d;
	TensorO4 D, I_x_I;
	TensorO2 I;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			if (i == j) I.t2(i,j) = 1;
		}
	}

	double G, K;	//Shear and bulk modulus
	double A = 1;
	G = E / (2 * (1 + v));
	K = E / (3 * (1 - 2 * v));

	I_d = UnitTensorO4(eIDTensor::UNIT_4_DEV);
	I_x_I = TensorOperations::dyadProduct(I, I);

	for (int ind_ij = 0; ind_ij < 9; ind_ij++) {
		for (int ind_kl = 0; ind_kl < 9; ind_kl++) {
			D.t4(ind_ij, ind_kl) = 2 * G * I_d.t4(ind_ij, ind_kl) + A * (K - (2 / 3) * G) * I_x_I.t4(ind_ij, ind_kl);
		}
	}

	return D;
}

void UMAT_VonMisces::hardening_law(eHardeningType H_type, double h, double sigma_y_0)
{
	switch (H_type) {
		case eHardeningType::LINEAR_ISOTROPIC: 
			m_hardening_slope = RandomProcess(h);
			m_sigma_y_init = RandomProcess(sigma_y_0);
			break;
		default:
			m_hardening_slope = RandomProcess(h);
			m_sigma_y_init = RandomProcess(sigma_y_0);
			break;
	}
}

void UMAT_VonMisces::computeFPKE(TensorO2 sig_0)
{
	vector<double> zeros(9, 0.0);
	vector<vector<double>> Ne_1_ij(9, zeros);		//Pre-yield Advective Tensor
	vector<vector<double>> Ne_2_ij(9, zeros);		//Pre-yield Diffusive Tensor
	vector<vector<double>> Np_1_ij(9, zeros);		//Post-yield Advective Tensor
	vector<vector<double>> Np_2_ij(9, zeros);		//Post-yield Diffusive Tensor

	vector<vector<double>> InitialStressPDF;
	for (int i = 0; i < 9; i++) {
		InitialStressPDF[i] = getStochasticStress().getInitialDistribution(100); // ProbabilisticMat::eDistribution::NORMAL, sig_0[i], 1e-5,
	}
		

	int ind_ij, ind_kl;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					ind_ij = 3 * (i - 1) + j;
					ind_kl = 3 * (k - 1) + l;
				}
			}
		}
	}
}