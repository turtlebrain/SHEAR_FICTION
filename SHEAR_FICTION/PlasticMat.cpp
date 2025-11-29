#include "PlasticMat.h"
#include <iostream>
#include "Tensors.h"
#include <vector>
using namespace std;

UMAT_VonMisces::UMAT_VonMisces(double E, double v, double sigy, double H)
{
	m_young_modulus.mean = E;
	m_poisson_ratio.mean = v;
	m_sigma_y_init.mean = sigy;
	m_hardening_slope.mean = H;
	m_bulk_modulus.mean = E / (3 * (1 - 2 * v));
	m_shear_modulus.mean = E / (2 * (1 + v));
}

void UMAT_VonMisces::EulerBackward(const vector<double>& delta_eps,const vector<double>& eps_tot_pstep,const vector<double>& eps_e_pstep, const vector<double>& sigma_pstep, const double eps_p_pstep) 
{
	//These locals should probably be accessed from the class
	//Present values should be passed as arguments
	TensorOps Tensor;
	// Instantiation below can later become a SetMaterial function with input material model.
	//UMAT_VonMisces VM_Mat(E, v);
	vector<double> eps_tot_fstep(9);									//Total strain at time t_n+1
	vector<double> eps_e_fstep(9), eps_e_trial_fstep(9);				//elastic strain at time t_n+1 
	double eps_p_fstep, eps_p_trial_fstep;								//plastic strain at time t_n+1 
	vector<double> p_fstep(9), s_fstep(9);								//volumetric and deviatoric stress at time t_n+1
	vector<double> sigma_fstep(9);
 
	vector<double> sigma_trial_fstep(9), p_trial_fstep(9), s_trial_fstep(9);;
	vector<double> eps_e_vol_trial(9), eps_e_dev_trial(9);

	double q_trial_fstep;											//trial Von-Mises effective stress
	double sigma_y;													//yield stress
	double sigma_y_init;											//Initial yield stress
	double H;														//Hardening slope
	double phi;														//yield function
	double delta_gamma = 0;											//plastic multiplier
	double resid_d;													//Residual derivative

	eps_p_trial_fstep = eps_p_pstep;
	for (int i = 0; i < 9; i++) {
		eps_tot_fstep[i] = eps_tot_pstep[i] + delta_eps[i];
		if (i == 0) eps_e_trial_fstep[i] = eps_tot_fstep[i] - eps_p_trial_fstep;
		
	}

	//Elastic trial state/Elastic predictor step
	//D_elastic = VM_Mat.getLinearElasticTensor(E, v);
	Tensor.tOps_volDecomp(eps_e_trial_fstep);
	eps_e_dev_trial = Tensor.m_dev_tensor;
	eps_e_vol_trial = Tensor.m_vol_tensor;
	double K = getK();
	double G = getG();
	for (int i = 0; i < 9; i++) {
		//eps_e_trial_fstep[i] = eps_e_pstep[i] + delta_eps[i];
		p_trial_fstep[i] =  K * eps_e_vol_trial[i]; 
		s_trial_fstep[i] = 2 * G * eps_e_dev_trial[i];
		sigma_trial_fstep[i] = s_trial_fstep[i] + p_trial_fstep[i];
	}
	q_trial_fstep = sqrt((3.0 * (Tensor.tOps_2XContrac(s_trial_fstep, s_trial_fstep))/2.0));
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
			s_fstep[i] = (1.0 - (delta_gamma * 3.0 * G) / q_trial_fstep) * s_trial_fstep[i];
			sigma_fstep[i] = s_fstep[i] + p_fstep[i];
			eps_e_fstep[i] = (1 / (2.0 * G)) * s_fstep[i] +  eps_e_vol_trial[i]; // (1 / (2.0 * G)) * s_fstep[i] +  (1 / 3.0) * eps_e_vol_trial[i]
		}
	}

	setTotalStrain(eps_tot_fstep);
	setElasticStrain(eps_e_fstep);
	setPlasticStrain(eps_p_fstep);
	setSigmaVolumetric(p_fstep);
	setSigmaDeviatoric(s_fstep);
	setSigma(sigma_fstep);

}

double UMAT_VonMisces::yield_function(const vector<double>& sigma, const double sigma_y)
{
	TensorOps Tensor;
	vector<double> s_ij(9), s_ji(9);
	double sigma_mean, J2, phi;
	int ind_ij;
	sigma_mean = (1 / 3.0) * (Tensor.tOps_trace(sigma));
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			ind_ij = 3 * (i - 1) + j;
			if (i == j) s_ij[ind_ij-1] = sigma[ind_ij-1] - sigma_mean;
			else s_ij[ind_ij-1] = sigma[ind_ij-1];
		}
	}
	s_ji = Tensor.tOps_transpose(s_ij);
	J2 = (1 / 2.0) * (Tensor.tOps_2XContrac(s_ij, s_ji));
	
	phi = sqrt(3.0 * J2) - sigma_y;
	return phi;
}

void UMAT_VonMisces::associative_flow_rule(const vector<double>& sigma)
{
	TensorOps Tensor;
	vector<double> s_ij(9), N_flow_vector(9);
	double sigma_mean;
	int ind_ij;

	sigma_mean = (1 / 3) * (Tensor.tOps_trace(sigma));
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			ind_ij = 3 * (i - 1) + j;
			if (i == j) s_ij[ind_ij] = sigma[ind_ij] - sigma_mean;
			else s_ij[ind_ij] = sigma[ind_ij];
		}
	}

	double mod_s = Tensor.tOps_L2norm(s_ij);
	for (int i = 0; i < 9; i++) N_flow_vector[i] = (sqrt(3 / 2)) * s_ij[i] / mod_s;

	setNvec(N_flow_vector);
}

vector<vector<double>> UMAT_VonMisces::getLinearElasticTensor(const double E, const double v) const
{
	vector<vector<double>> D, I_d, I_x_I;
	vector<double> I(9, 0);
	int ind_ij, ind_kl;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			ind_ij = 3 * (i - 1) + j;
			if (i == j) I[ind_ij] = 1;
		}
	}
	TensorOps Tensor;
	double G, K;	//Shear and bulk modulus
	double A = 1;

	G = E / (2 * (1 + v));
	K = E / (3 * (1 - 2 * v));

	I_d = Tensor.unit_4t(Tensor._Unit4Dev);
	I_x_I = Tensor.tOps_DyadProduct(I, I);

	for (ind_ij = 0; ind_ij < 9; ind_ij++) {
		for (ind_kl = 0; ind_kl < 9; ind_kl++) {
			D[ind_ij][ind_kl] = 2 * G * I_d[ind_ij][ind_kl] + A * (K - (2 / 3) * G) * I_x_I[ind_ij][ind_kl];
		}
	}

	return D;
}

void UMAT_VonMisces::hardening_law(eHardeningType H_type, double h, double sigma_y_0)
{
	switch (H_type) {
		case eHardeningType::LINEAR_ISOTROPIC: 
			m_hardening_slope.mean = h;
			m_sigma_y_init.mean = sigma_y_0;
			break;
		default:
			m_hardening_slope.mean = h;
			m_sigma_y_init.mean = sigma_y_0;
			break;
	}
}

void UMAT_VonMisces::computeFPKE(vector<double> sig_0)
{
	vector<double> zeros(9, 0.0);
	vector<vector<double>> Ne_1_ij(9, zeros);		//Pre-yield Advective Tensor
	vector<vector<double>> Ne_2_ij(9, zeros);		//Pre-yield Diffusive Tensor
	vector<vector<double>> Np_1_ij(9, zeros);		//Post-yield Advective Tensor
	vector<vector<double>> Np_2_ij(9, zeros);		//Post-yield Diffusive Tensor

	vector<vector<double>> InitialStressPDF;
	for (int i = 0; i < 9; i++) {
		InitialStressPDF[i] = GetStochasticStress().GetInitialDistribution(ProbabilisticMat::eDistribution::NORMAL, sig_0[i], 1e-5, 100);
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