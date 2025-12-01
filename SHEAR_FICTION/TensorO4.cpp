#include "TensorO4.h"
#include "TensorO2.h"

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

double& TensorO4::t4(size_t ind0_ij, size_t ind0_kl)
{
	return _m_t4[ind0_ij][ind0_kl];
}

TensorO4 UnitTensorO4::_unit4(eIDTensor u4)
{
	TensorO2 kronekr;
	TensorO4 UT4;
	TensorO4 UT4Sym;
	TensorO4 C;
	TensorO4 C_temp;

	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			if (i == j) kronekr.t2(i,j) = 1.0;
			else kronekr.t2(i,j) = 0.0;
		}
	}

	for (int ind_ik = 0; ind_ik < 9; ind_ik++) {
		for (int ind_jl = 0; ind_jl < 9; ind_jl++) {
			if (ind_ik == ind_jl) { UT4.t4(ind_ik, ind_jl) = 1.0; }
			else { UT4.t4(ind_ik, ind_jl) = 0.0; }
		}
	}

	for (int i = 1; i <= 3; i++) {
		for (int k = 1; k <= 3; k++) {
			for (int j = 1; j <= 3; j++) {
				for (int l = 1; l <= 3; l++) {
					UT4Sym.t4(i, k, j, l) = (0.5) * (kronekr.t2(i,k) * kronekr.t2(j,l) + kronekr.t2(i,l) * kronekr.t2(j,k));
				}
			}
		}
	}

	TensorO4 dy_kronekr = kronekr.dyadProduct(kronekr);

	switch (u4) {
	case eIDTensor::UNIT_4:
		C = UT4;
		break;

	case eIDTensor::UNIT_4_SYM:
		C = UT4Sym;
		break;

	case eIDTensor::UNIT_4_ANTI_SYM:
		for (int ind_ik = 0; ind_ik < 9; ind_ik++) {
			for (int ind_jl = 0; ind_jl < 9; ind_jl++) {
				if (ind_ik == ind_jl) { C.t4(ind_ik, ind_jl) = 1.0; }
				else { C_temp.t4(ind_ik, ind_jl) = 0.0; }
			}
		}
		for (int i = 1; i <= 3; i++) {
			for (int k = 1; k <= 3; k++) {
				for (int j = 1; j <= 3; j++) {
					for (int l = 1; l <= 3; l++) {
						C.t4(i,k,j,l) = (0.5) * (kronekr.t2(i,k) * kronekr.t2(j,l) - kronekr.t2(i,l) * kronekr.t2(j,k));
					}
				}
			}
		}
		break;

	case eIDTensor::UNIT_4_VOL:
		for (int i = 1; i <= 3; i++) {
			for (int k = 1; k <= 3; k++) {
				for (int j = 1; j <= 3; j++) {
					for (int l = 1; l <= 3; l++) {
						C.t4(i,k,j,l) = (1 / 3.0) * dy_kronekr.t4(i,k,j,l);
					}
				}
			}
		}
		break;

	case eIDTensor::UNIT_4_DEV:
		for (int i = 1; i <= 3; i++) {
			for (int k = 1; k <= 3; k++) {
				for (int j = 1; j <= 3; j++) {
					for (int l = 1; l <= 3; l++) {
						C_temp.t4(i,k,j,l) = (1 / 3.0) * dy_kronekr.t4(i,k,j,l);
						C.t4(i,k,j,l) = UT4.t4(i,k,j,l) - C_temp.t4(i,k,j,l);
					}
				}
			}
		}
		break;
	}
	return C;
}