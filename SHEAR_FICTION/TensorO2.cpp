#include "TensorO2.h"
#include <iostream>
using namespace std;

double& TensorO2::t2(size_t i, size_t j)
{
	// convert indices into a single index for 9x1 vector
	size_t ind_ij = 3 * (i - 1) + j;
	// convert index to zero-based index
	size_t ind0_ij = ind_ij - 1;
	return _m_t2[ind0_ij];
}

double& TensorO2::t2(size_t ind0_ij)
{
	return _m_t2[ind0_ij];
}

double TensorO2::trace()
{
	double C = 0.0;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			if (i == j) { C += this->t2(i,j); }
		}
	}
	return C;
}

TensorO2 TensorO2::transpose()
{
	TensorO2 A_trans;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			A_trans.t2(i, j) = this->t2(j, i);
		}
	}
	return A_trans;
}

double TensorO2::determinant()
{
	double det_A;
	det_A = _m_t2[0] * (_m_t2[4] * _m_t2[8] - _m_t2[5] * _m_t2[7]) - 
			_m_t2[1] * (_m_t2[3] * _m_t2[8] - _m_t2[5] * _m_t2[6]) + 
			_m_t2[2] * (_m_t2[3] * _m_t2[7] - _m_t2[4] * _m_t2[6]);
	return det_A;
}

TensorO2 TensorO2::inverse()
{
	double detA = this->determinant();
	TensorO2 Amin; TensorO2 Acof; TensorO2 invA;
	if (detA == 0) { cout << "Tensor is singular" << endl; exit(0); }

	Amin._m_t2[0] = (_m_t2[4] * _m_t2[8]) - (_m_t2[5] * _m_t2[7]);
	Amin._m_t2[1] = (_m_t2[3] * _m_t2[8]) - (_m_t2[5] * _m_t2[6]);
	Amin._m_t2[2] = (_m_t2[3] * _m_t2[8]) - (_m_t2[4] * _m_t2[6]);
	Amin._m_t2[3] = (_m_t2[1] * _m_t2[8]) - (_m_t2[2] * _m_t2[7]);
	Amin._m_t2[4] = (_m_t2[0] * _m_t2[8]) - (_m_t2[2] * _m_t2[6]);
	Amin._m_t2[5] = (_m_t2[0] * _m_t2[7]) - (_m_t2[1] * _m_t2[6]);
	Amin._m_t2[6] = (_m_t2[1] * _m_t2[5]) - (_m_t2[2] * _m_t2[4]);
	Amin._m_t2[7] = (_m_t2[0] * _m_t2[5]) - (_m_t2[2] * _m_t2[3]);
	Amin._m_t2[8] = (_m_t2[0] * _m_t2[4]) - (_m_t2[1] * _m_t2[3]);

	for (int i = 1; i <= 9; i++) {
		if ((i % 2) == 0) {
			Acof._m_t2[i - 1] = -Amin._m_t2[i - 1];
		}
		else {
			Acof._m_t2[i - 1] = Amin._m_t2[i - 1];
		}
	}

	TensorO2 Aadj = Acof.transpose();

	for (int i = 1; i <= 9; i++) {
		invA._m_t2[i - 1] = (1 / detA) * Aadj._m_t2[i - 1];
	}

	return invA;
}
