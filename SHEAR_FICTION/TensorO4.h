#pragma once
#include <vector>

using namespace std;
// 4th order tensors are represented by a 9x9 Array T: {T1' T2' T3' T4' T5' T6' T7' T8' T9'} , T1...T9 are 9x1 Vectors

enum class eIDTensor 
{
	// Fourth-Order Unit Tensor
	UNIT_4,
	UNIT_4_SYM,
	UNIT_4_ANTI_SYM,
	UNIT_4_VOL,
	UNIT_4_DEV
};

class TensorO4
{
public:
	TensorO4() {_m_t4.resize(9, vector<double>(9, 0.0));}
	TensorO4(const vector<vector<double>>& tensor4) : _m_t4(tensor4) {}
	TensorO4(const TensorO4& tensor4) : _m_t4(tensor4._m_t4) {}
	~TensorO4() {}

public:
	double& t4(size_t i, size_t j, size_t k, size_t l);
	double& t4(size_t ind0_ij, size_t ind0_kl);


private:
	vector<vector<double>> _m_t4;
};

class UnitTensorO4 : TensorO4
{
public:
	UnitTensorO4(eIDTensor u4 = eIDTensor::UNIT_4) : TensorO4(_unit4(u4)) {}
	~UnitTensorO4() {}

private:
	TensorO4 _unit4(eIDTensor u4);
};