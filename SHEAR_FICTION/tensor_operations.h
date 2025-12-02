#pragma once
#include "TensorO4.h"
#include "TensorO2.h"
using namespace std;

class TensorOperations {
private:
	TensorOperations() {}

public:
	static TensorO2 contrac1X(const TensorO2& A, const TensorO2& B);			//Single contraction of two 2nd order tensors -> A.B
	static double contrac2X(const TensorO2& A, const TensorO2& B);				//Double contraction of two 2nd order tensors -> A:B	
	static TensorO4 contrac2X(const TensorO4& A, const TensorO4& B);			//Double contraction of two 4th order tensors
	static TensorO2 contrac2X(const TensorO4& A, const TensorO2& B);			//Double contraction of a 4th order and 2nd order tensor
	static TensorO2 contrac2X(const TensorO2& A, const TensorO4& B);			//Double contraction of 2nd order and 4th order tensor
	static TensorO4 dyadProduct(const TensorO2& A, const TensorO2& B);			//Dyadic roduct of two 2nd order tensors -> AxB
};