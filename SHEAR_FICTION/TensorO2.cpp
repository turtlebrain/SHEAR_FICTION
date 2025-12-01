#include "TensorO2.h"
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

}