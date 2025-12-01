#include <iostream>
#include "Tensors.h"
#include <vector>
using namespace std;


TensorOps::TensorOps()
{

}

TensorOps::~TensorOps()
{

}


vector<double> TensorOps::tOps_1XContrac(const vector<double>& A, const vector<double>& B) {
	int ind_ji; int ind_ki; int ind_kj;
	vector<double> C(9);
	for (int i = 0; i < 9; i++) { C[i] = 0; }
	for (int k = 1; k <= 3; k++) {
		for (int i = 1; i <= 3; i++) {
			for (int j = 1; j <= 3; j++) {
				ind_ji = 3 * (j - 1) + i;
				ind_ki = 3 * (k - 1) + i;
				ind_kj = 3 * (k - 1) + j;
				C[ind_ki - 1] = A[ind_kj - 1] * B[ind_ji - 1] + C[ind_ki - 1];
			}
		}
	}
	return C;
}

double TensorOps::tOps_2XContrac(const vector<double>& A, const vector<double>& B) {
	double C = 0; int ind_ij;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			ind_ij = 3 * (i - 1) + j;
			C += A[ind_ij - 1] * B[ind_ij - 1];
		}
	}
	return C;
}

vector<vector<double>> TensorOps::tOps_2XContrac(const vector<vector<double>>& A, const vector<vector<double>>& B) {
	vector<double> temp(9, 0.0);
	vector<vector<double>> C(9, temp);

	for (int ij = 0; ij < 9; ij++) {
		for (int kl = 0; kl < 9; kl++) {
			for (int mn = 0; mn < 9; mn++) {
				C[ij][kl] = A[ij][mn] * B[mn][kl] + C[ij][kl];
			}
		}
	}
	return C;
}

vector<double> TensorOps::tOps_2XContrac(const vector<vector<double>>& A, const vector<double>& B) {
	int ind_ij; int ind_kl;
	vector<double> C(9);
	for (int i = 0; i < 9; i++) { C[i] = 0; }
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					ind_ij = 3 * (i - 1) + j;
					ind_kl = 3 * (k - 1) + l;
					C[ind_ij - 1] = A[ind_ij - 1][ind_kl - 1] * B[ind_kl - 1] + C[ind_ij - 1];
				}
			}
		}
	}
	return C;
}

vector<double> TensorOps::tOps_2XContrac(const vector<double>& A, const vector<vector<double>>& B) {
	vector<double> C(9);
	for (int i = 0; i < 9; i++) { C[i] = 0; }
	int ind_ij; int ind_kl;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					ind_ij = 3 * (i - 1) + j;
					ind_kl = 3 * (k - 1) + l;
					C[ind_ij - 1] = A[ind_kl - 1] * B[ind_kl - 1][ind_ij - 1] + C[ind_ij - 1];
				}
			}
		}
	}
	return C;
}

double TensorOps::tOps_trace(const vector<double>& A) {
	int ind_ij;
	double C = 0.0;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			ind_ij = 3 * (i - 1) + j;
			if (i == j) { C += A[ind_ij - 1]; }
		}
	}
	return C;
}

vector<double> TensorOps::tOps_transpose(const vector<double>& A) {
	int ind_ij;
	int ind_ji;
	vector<double> A_trans(9);

	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			ind_ji = 3 * (j - 1) + i;
			ind_ij = 3 * (i - 1) + j;
			A_trans[ind_ij - 1] = A[ind_ji - 1];
		}
	}
	return A_trans;
}

double TensorOps::tOps_determinant(const vector<double>& A) {
	double det_A;
	det_A = A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) + A[2] * (A[3] * A[7] - A[4] * A[6]);
	return det_A;
}


vector<double> TensorOps::tOps_inverse(const vector<double>& A) {
	double detA = tOps_determinant(A);
	vector<double> Amin(9); vector<double> Acof(9);
	vector<double> invA(9);
	if (detA == 0) { cout << "Tensor is singular" << endl; exit(0); }

	Amin[0] = A[4] * A[8] - A[5] * A[7];
	Amin[1] = A[3] * A[8] - A[5] * A[6];
	Amin[2] = A[3] * A[8] - A[4] * A[6];
	Amin[3] = A[1] * A[8] - A[2] * A[7];
	Amin[4] = A[0] * A[8] - A[2] * A[6];
	Amin[5] = A[0] * A[7] - A[1] * A[6];
	Amin[6] = A[1] * A[5] - A[2] * A[4];
	Amin[7] = A[0] * A[5] - A[2] * A[3];
	Amin[8] = A[0] * A[4] - A[1] * A[3];

	for (int i = 1; i <= 9; i++) {
		if ((i % 2) == 0) {
			Acof[i - 1] = -Amin[i - 1];
		}
		else {
			Acof[i - 1] = Amin[i - 1];
		}
	}

	vector<double> Aadj = tOps_transpose(Acof);

	for (int i = 1; i <= 9; i++) {
		invA[i - 1] = (1 / detA) * Aadj[i - 1];
	}

	return invA;
}

vector<vector<double>> TensorOps::tOps_DyadProduct(const vector<double>& A, const vector<double>& B) {

	vector<double> temp(9);
	vector<vector<double>> C(9, temp);

	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					int ind_ij = 3 * (i - 1) + j;
					int ind_kl = 3 * (k - 1) + l;
					C[ind_ij - 1][ind_kl - 1] = A[ind_ij - 1] * B[ind_kl - 1];
				}
			}
		}
	}
	return C;
}

double TensorOps::tOps_L2norm(const vector<double>& A) {
	double C = tOps_2XContrac(A, A);
	return sqrt(C);
}

void TensorOps::tOps_volDecomp(const vector<double>& A) {
	vector<vector<double>> U4vol, U4dev;
	U4vol = unit_4t(_Unit4Vol);
	U4dev = unit_4t(_Unit4Dev);
	m_vol_tensor = tOps_2XContrac(U4vol, A);
	m_dev_tensor = tOps_2XContrac(U4dev, A);
}

vector<vector<double>> TensorOps::unit_4t(IDTensor _U4) {
	vector<double> temp(9, 0.0);
	vector<double> kronekr(9);
	vector<vector<double>> UT4(9, temp);
	vector<vector<double>> UT4Sym(9, temp);
	vector<vector<double>> C(9, temp);
	vector<vector<double>> C_temp(9, temp);
	int ind_ij;  int ind_ik; int ind_jl; int ind_il; int ind_jk;

	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			ind_ij = 3 * (i - 1) + j;
			if (i == j) kronekr[ind_ij - 1] = 1.0;
			else kronekr[ind_ij - 1] = 0.0;
		}
	}

	for (int ind_ik = 0; ind_ik < 9; ind_ik++) {
		for (int ind_jl = 0; ind_jl < 9; ind_jl++) {
			if (ind_ik == ind_jl) { UT4[ind_ik][ind_jl] = 1.0; }
			else { UT4[ind_ik][ind_jl] = 0.0; }
		}
	}

	for (int i = 1; i <= 3; i++) {
		for (int k = 1; k <= 3; k++) {
			for (int j = 1; j <= 3; j++) {
				for (int l = 1; l <= 3; l++) {
					ind_ik = 3 * (i - 1) + k;
					ind_jl = 3 * (j - 1) + l;
					ind_il = 3 * (i - 1) + l;
					ind_jk = 3 * (j - 1) + k;
					UT4Sym[ind_ik - 1][ind_jl - 1] = (0.5) * (kronekr[ind_ik - 1] * kronekr[ind_jl - 1] + kronekr[ind_il - 1] * kronekr[ind_jk - 1]);
				}
			}
		}
	}

	vector<vector<double>> dy_kronekr = tOps_DyadProduct(kronekr, kronekr);

	switch (_U4) {
	case _Unit4:
		C = UT4;
		break;

	case _Unit4Sym:
		C = UT4Sym;
		break;

	case _Unit4AntiSym:
		for (int ind_ik = 0; ind_ik < 9; ind_ik++) {
			for (int ind_jl = 0; ind_jl < 9; ind_jl++) {
				if (ind_ik == ind_jl) { C[ind_ik][ind_jl] = 1.0; }
				else { C_temp[ind_ik][ind_jl] = 0.0; }
			}
		}
		for (int i = 1; i <= 3; i++) {
			for (int k = 1; k <= 3; k++) {
				for (int j = 1; j <= 3; j++) {
					for (int l = 1; l <= 3; l++) {
						ind_ik = 3 * (i - 1) + k;
						ind_jl = 3 * (j - 1) + l;
						ind_il = 3 * (i - 1) + l;
						ind_jk = 3 * (j - 1) + k;
						C[ind_ik - 1][ind_jl - 1] = (0.5) * (kronekr[ind_ik - 1] * kronekr[ind_jl - 1] - kronekr[ind_il - 1] * kronekr[ind_jk - 1]);
					}
				}
			}
		}
		break;

	case _Unit4Vol:
		for (int i = 1; i <= 3; i++) {
			for (int k = 1; k <= 3; k++) {
				for (int j = 1; j <= 3; j++) {
					for (int l = 1; l <= 3; l++) {
						ind_ik = 3 * (i - 1) + k;
						ind_jl = 3 * (j - 1) + l;
						C[ind_ik - 1][ind_jl - 1] = (1 / 3.0) * dy_kronekr[ind_ik - 1][ind_jl - 1];
					}
				}
			}
		}
		break;

	case _Unit4Dev:
		for (int i = 1; i <= 3; i++) {
			for (int k = 1; k <= 3; k++) {
				for (int j = 1; j <= 3; j++) {
					for (int l = 1; l <= 3; l++) {
						ind_ik = 3 * (i - 1) + k;
						ind_jl = 3 * (j - 1) + l;
						C_temp[ind_ik - 1][ind_jl - 1] = (1 / 3.0) * dy_kronekr[ind_ik - 1][ind_jl - 1];
						C[ind_ik - 1][ind_jl - 1] = UT4[ind_ik - 1][ind_jl - 1] - C_temp[ind_ik - 1][ind_jl - 1]; 
					}
				}
			}
		}
		break;
	}
	return C;
}

