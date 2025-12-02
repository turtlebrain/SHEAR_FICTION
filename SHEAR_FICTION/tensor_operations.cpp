#include "tensor_operations.h"

TensorO2 TensorOperations::contrac1X(const TensorO2& A, const TensorO2& B)
{
	TensorO2 C;
	TensorO2 Acopy = TensorO2(A);
	TensorO2 Bcopy = TensorO2(B);
	for (size_t k = 1; k <= 3; k++) {
		for (size_t i = 1; i <= 3; i++) {
			for (size_t j = 1; j <= 3; j++) {
				C.t2(k,i) = Acopy.t2(k,j) * Bcopy.t2(j,i) + C.t2(k,i);
			}
		}
	}
	return C;
}

double TensorOperations::contrac2X(const TensorO2& A, const TensorO2& B)
{
	double C = 0; 
	TensorO2 Acopy = TensorO2(A);
	TensorO2 Bcopy = TensorO2(B);
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			C += Acopy.t2(i,j) * Bcopy.t2(i,j);
		}
	}
	return C;
}

TensorO4 TensorOperations::contrac2X(const TensorO4& A, const TensorO4& B)
{
	TensorO4 C;
	TensorO4 Acopy = TensorO4(A);
	TensorO4 Bcopy = TensorO4(B);
	for (int ij = 0; ij < 9; ij++) {
		for (int kl = 0; kl < 9; kl++) {
			for (int mn = 0; mn < 9; mn++) {
				C.t4(ij,kl) = Acopy.t4(ij,mn) * Bcopy.t4(mn,kl) + C.t4(ij,kl);
			}
		}
	}
	return C;
}

TensorO2 TensorOperations::contrac2X(const TensorO4& A, const TensorO2& B)
{
	TensorO2 C;
	TensorO4 Acopy = TensorO4(A);
	TensorO2 Bcopy = TensorO2(B);
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					C.t2(i,j) = Acopy.t4(i,j,k,l) * Bcopy.t2(k,l) + C.t2(i,j);
				}
			}
		}
	}
	return C;
}

TensorO2 TensorOperations::contrac2X(const TensorO2& A, const TensorO4& B)
{
	TensorO2 C;
	TensorO2 Acopy = TensorO2(A);
	TensorO4 Bcopy = TensorO4(B);
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					C.t2(i,j) = Acopy.t2(k,l) * Bcopy.t4(k,l,i,j) + C.t2(i,j);
				}
			}
		}
	}
	return C;
}

TensorO4 TensorOperations::dyadProduct(const TensorO2& A, const TensorO2& B)
{
	TensorO4 C;
	TensorO2 Acopy = TensorO2(A);
	TensorO2 Bcopy = TensorO2(B);
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int k = 1; k <= 3; k++) {
				for (int l = 1; l <= 3; l++) {
					C.t4(i,j,k,l) = Acopy.t2(i,j) * Bcopy.t2(k,l);
				}
			}
		}
	}
	return C;
}