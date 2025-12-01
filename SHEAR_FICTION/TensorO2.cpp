#include "TensorO2.h"
#include <iostream>
#include <tuple>
using namespace std;

double& TensorO2::t2(size_t i, size_t j)
{
	// convert indices into a single index for 9x1 vector
	size_t ind_ij = 3 * (i - 1) + j;
	// convert index to zero-based index
	size_t ind0_ij = ind_ij - 1;
	return _m_t2[ind0_ij];
}

TensorO2 TensorO2::contrac1X(TensorO2& tensor2)
{
	TensorO2 C;
	for (size_t k = 1; k <= 3; k++) {
		for (size_t i = 1; i <= 3; i++) {
			for (size_t j = 1; j <= 3; j++) {
				C.t2(k, i) = this->t2(k, j) * tensor2.t2(j, i) + C.t2(k, i);
			}
		}
	}
	return C;
}

double TensorO2::contrac2X(TensorO2& tensor2)
{
	double C = 0;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			C += this->t2(i, j) * tensor2.t2(i, j);
		}
	}
	return C;
}

TensorO2 TensorO2::contrac2X(TensorO4& tensor4)
{
	TensorO2 C;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					C.t2(i, j) = this->t2(k,l) * tensor4.t4(k, l, i, j) + C.t2(i, j);
				}
			}
		}
	}
	return C;
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

	Amin._m_t2[0] = (this->_m_t2[4] * this->_m_t2[8]) - (this->_m_t2[5] * this->_m_t2[7]);
	Amin._m_t2[1] = (this->_m_t2[3] * this->_m_t2[8]) - (this->_m_t2[5] * this->_m_t2[6]);
	Amin._m_t2[2] = (this->_m_t2[3] * this->_m_t2[8]) - (this->_m_t2[4] * this->_m_t2[6]);
	Amin._m_t2[3] = (this->_m_t2[1] * this->_m_t2[8]) - (this->_m_t2[2] * this->_m_t2[7]);
	Amin._m_t2[4] = (this->_m_t2[0] * this->_m_t2[8]) - (this->_m_t2[2] * this->_m_t2[6]);
	Amin._m_t2[5] = (this->_m_t2[0] * this->_m_t2[7]) - (this->_m_t2[1] * this->_m_t2[6]);
	Amin._m_t2[6] = (this->_m_t2[1] * this->_m_t2[5]) - (this->_m_t2[2] * this->_m_t2[4]);
	Amin._m_t2[7] = (this->_m_t2[0] * this->_m_t2[5]) - (this->_m_t2[2] * this->_m_t2[3]);
	Amin._m_t2[8] = (this->_m_t2[0] * this->_m_t2[4]) - (this->_m_t2[1] * this->_m_t2[3]);

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

TensorO4 TensorO2::dyadProduct(TensorO2& tensor2)
{
	TensorO4 C;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					C.t4(i,j,k,l) = this->t2(i,j) * tensor2.t2(k,l);
				}
			}
		}
	}
	return C;
}

tuple<TensorO2, TensorO2> TensorO2::volDecomp()
{
	return { TensorO2(), TensorO2() };
}