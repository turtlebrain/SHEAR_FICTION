// SHEAR_FICTION.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "random_process.h"
#include "plastic_material.h"
#include "tensor_operations.h"
#include "tensor_o2.h"
#include "tensor_o4.h"
#include <vector>
#include <tuple>
using namespace std;

void print2DVec(vector<double> A);

int main()
{
	//------------------------------------------Tensor Operation Verification---------------------------------------------//
	
	TensorO2 s1(vector<double> { 50, 30, 20, 30, -20, -10, 20, -10, 10 });
	tuple<TensorO2, TensorO2> s1_vol_dev =  s1.volDecomp();
	TensorO2 s1vol = get<0>(s1_vol_dev);
	TensorO2 s1dev = get<1>(s1_vol_dev);
	print2DVec(s1vol.t2());
	cout << endl;
	print2DVec(s1dev.t2());
	cout << endl;
	//------------------------------------------ProbabilisticMat function Verification------------------------------------//
	
	RandomProcess StdNormalMat(0, 1, eDistribution::NORMAL_DISTRIBUTION);
	vector<double> pdf = StdNormalMat.getInitialDistribution(100);
	print2DVec(pdf);
	
	//------------------------------------------Von Mises Verification----------------------------------------------------//
	//Material Properties
	double Mat_Prop[9];
	Mat_Prop[0] = 220.5;			//Young's Modulus
	Mat_Prop[1] = 0.30;				//Poisson's ratio
	Mat_Prop[2] = 0.5;				//Initial yield stress
	Mat_Prop[3] = 0;				//Hardening slope
	UMAT_VonMisces Material_1(Mat_Prop[0], Mat_Prop[1], Mat_Prop[2], Mat_Prop[3]);

	// Initial state variables
	TensorO2 delta_eps;
	delta_eps.t2(0) = 0.00005;
	TensorO2 sigma_pstep, sig_dev_pstep, sig_vol_pstep;
	TensorO2 eps_tot_pstep, eps_e_pstep;
	double eps_p_pstep = 0;
	int counter = 0;

	Material_1.setSigma(sigma_pstep);
	Material_1.setTotalStrain(eps_tot_pstep);
	Material_1.setElasticStrain(eps_e_pstep);
	Material_1.setPlasticStrain(eps_p_pstep);
	//void setSigmaVolumetric(vector<double> sig_vol) { m_sigma_vol = sig_vol; }
	//void setSigmaDeviatoric(vector<double> sig_dev) { m_sigma_dev = sig_dev; }
	//void setElasticStrain(vector<double> epsilon_e) { m_eps_e = epsilon_e; }
	//void setPlasticStrain(double epsilon_p) { m_eps_p = epsilon_p; }


	while (counter < 100 && eps_e_pstep.t2(0) < 0.05) {
		Material_1.EulerBackward(delta_eps, eps_tot_pstep, eps_e_pstep, sigma_pstep, eps_p_pstep);

		eps_tot_pstep = Material_1.getTotalStrain();
		eps_p_pstep = Material_1.getPlasticStrain();
		eps_e_pstep = Material_1.getElasticStrain();
		sigma_pstep = Material_1.getSigma();
		sig_dev_pstep = Material_1.getSigmaDev();
		sig_vol_pstep = Material_1.getSigmaVol();

		cout << sigma_pstep.t2(0) << endl;
		//cout << eps_e_pstep[0] << endl;
		//cout << sig_vol_pstep[0] << endl;
		//cout << sig_dev_pstep[0] << endl;
		//cout << eps_e_pstep[0] << endl;
		//cout << sig_dev_pstep[0] << " + " << sig_vol_pstep[0] << " = " << sigma_pstep[0] << endl;
		counter = counter + 1;
	}
}

void print2DVec(vector<double> A) {
	int size_A = A.size();
	for (int i = 0; i < size_A; i++) cout << A[i] << " ";
}