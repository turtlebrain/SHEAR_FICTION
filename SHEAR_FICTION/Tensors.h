#pragma once
#include <vector>
using namespace std;

class TensorOps
{
public:
	TensorOps();
	~TensorOps();

	enum IDTensor
	{
		// Fourth-Order Unit Tensor
		_Unit4,
		_Unit4Sym,
		_Unit4AntiSym,
		_Unit4Vol,
		_Unit4Dev
	};

	vector<double> tOps_1XContrac(const vector<double>& A, const vector<double>& B);									//Single contraction of two 2nd order tensors -> A.B
	double tOps_2XContrac(const vector<double>& A, const vector<double>& B);											//Double contraction of two 2nd order tensors -> A:B	
	vector<vector<double>> tOps_2XContrac(const vector<vector<double>>& A, const vector<vector<double>>& B);			//Double contraction of a 4th order and 2nd order tensor
	vector<double> tOps_2XContrac(const vector<vector<double>>& A, const vector<double>& B);							//Double contraction of two 4th order tensors
	vector<double> tOps_2XContrac(const vector<double>& A, const vector<vector<double>>& B);							//Double contraction of 2nd order and 4th order tensor
	double tOps_trace(const vector<double>& A);																			//Trace Operation -> Aii = tr(A)
	vector<double> tOps_transpose(const vector<double>& A);																//Transpose Operation -> Aij= Aji
	double tOps_determinant(const vector<double>& A);																	//determinant Operation -> det(A)
	vector<double> tOps_inverse(const vector<double>& A);																//inverse Operation -> inv(A)
	vector<vector<double>> tOps_DyadProduct(const vector<double>& A, const vector<double>& B);							//Dyadic roduct of two 2nd order tensors -> AxB
	double tOps_L2norm(const vector<double>& A);																		//Euclidean norm (L2) ||A||
	void tOps_volDecomp(const vector<double>& A);																		//Volumetric-deviatoric tensor decomposition
	vector<vector<double>> unit_4t(IDTensor _U4);																		//Returns Fourth Oder Unit Tensor

public:
	vector<double> m_vol_tensor, m_dev_tensor;
};




