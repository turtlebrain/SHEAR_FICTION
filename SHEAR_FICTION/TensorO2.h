#pragma once
#include <vector>
#include "TensorO4.h"
#include <tuple>
using namespace std;
// 2nd order tensors are represented by a 9x1 Array t: {t11 t12 t13 t21 t22 t23 t31 t32 t33}
class TensorO2
{
public:
	TensorO2() { _m_t2.resize(9, 0.0); }
	TensorO2(const vector<double>& tensor2) : _m_t2(tensor2) {}
	TensorO2(const TensorO2& tensor2) : _m_t2(tensor2._m_t2) {}
	~TensorO2() {}

public:
	double& t2(size_t i, size_t j);
	double& t2(size_t ind0_ij);

public:
	// Tensor Operations
	TensorO2 contrac1X(TensorO2& tensor2);		//Single contraction of two 2nd order tensors -> A.B
	double contrac2X(TensorO2& tensor2);		//Double contraction of two 2nd order tensors -> A:B
	TensorO2 contrac2X(TensorO4& tensor4);		//Double contraction of 2nd order and 4th order tensor
	TensorO4 dyadProduct(TensorO2& tensor2);	//Dyadic roduct of two 2nd order tensors -> AxB
	double trace();								//Trace Operation -> Aii = tr(A)
	TensorO2 transpose();						//Transpose Operation -> Aij= Aji		
	double determinant();						//determinant Operation -> det(A)
	TensorO2 inverse();							//inverse Operation -> inv(A)
	tuple<TensorO2, TensorO2> volDecomp();		//Volumetric-deviatoric tensor decomposition

private:
	vector<double> _m_t2;
};