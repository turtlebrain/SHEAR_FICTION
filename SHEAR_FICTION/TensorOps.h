#pragma once
#include <vector>
using namespace std;
// 1st order tensor are represented by a double s_ 
// 2nd order tensors are represented by a 9x1 Array t_: {t11 t12 t13 t21 t22 t23 t31 t32 t33}
// 4th order tensors are represented by a 9x9 Arrary T_: {T1' T2' T3' T4' T5' T6' T7' T8' T9'} , T1...T9 are 9x1 Vectors
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

class TensorContainer
{
public:
	// There should be a RemoveT4 :https://stackoverflow.com/questions/39912/how-do-i-remove-an-item-from-a-stl-vector-with-a-certain-value
	void AddTensor(const vector<vector<double>>& T4);
	void AddTensor(const vector <double>& T2);
	void AddTensor(const double T1);
	vector<vector<double>> GetT4(int index) const;
	vector<double> GetT2(int index) const;
	double GetT1(int index) const;
	


protected:
	vector<vector<vector<double>>> T4_container;
	vector<vector<double>> T2_container;
	vector<double> T1_container;
};



