#include "TensorO4.h"

double& TensorO4::t4(size_t i, size_t j, size_t k, size_t l)
{
	// convert indices into a single index for 9x1 vector
	size_t ind_ij = 3 * (i - 1) + j;
	size_t ind_kl = 3 * (k - 1) + l;
	// convert index to zero-based index
	size_t ind0_ij = ind_ij - 1;
	size_t ind0_kl = ind_kl - 1;
	return _m_t4[ind0_ij][ind0_kl];
}