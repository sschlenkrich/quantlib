/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2016, Sebastian Schlenkrich

*/



#ifndef quantlib_templateauxilliaries_choleskyfactorisation_hpp
#define quantlib_templateauxilliaries_choleskyfactorisation_hpp

//#include <ql/types.hpp>
//#include <boost/function.hpp>

#include <stdlib.h>
#include <math.h>


namespace TemplateAuxilliaries {

    // This method is unsafe! Memory allocation need to be ensured by user
	// A(i,j) (and L) is stored row-wise as A[i*n+j]
	// L(i,j) is lower triangular matrix with L L^T = A
	template <typename Type> 
    inline void cholesky(Type *A, Type *L, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < (i+1); ++j) {
                Type s = 0;
                for (size_t k = 0; k < j; ++k) s += L[i * n + k] * L[j * n + k];
                if (i==j) {
                    if (A[i * n + i] < s) throw std::exception();
                    L[i * n + j] = sqrt(A[i * n + i] - s);
                }
                else L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
            }
        } 
    }

	template <typename Type>
	std::vector< std::vector< Type > > cholesky(const std::vector< std::vector< Type > >& A) {
		Type *arrayA = new Type[A.size()*A.size()];
		Type *arrayL = new Type[A.size()*A.size()];
		for (size_t i = 0; i < A.size(); ++i) {
			if (A.size()!=A[i].size()) throw std::exception();
			for (size_t j = 0; j < A[i].size(); ++j) arrayA[i*A[i].size() + j] = A[i][j];
		}
		cholesky(arrayA, arrayL, A.size());
		std::vector< std::vector< Type > > L(A.size(), std::vector< Type >(A.size(), 0.0));
		for (size_t i = 0; i < L.size(); ++i) {
			for (size_t j = 0; j <= i; ++j) L[i][j] = arrayL[i*L[i].size() + j];
		}
		delete arrayA;
		delete arrayL;
		return L;
	}

}

#endif  /* ifndef quantlib_templateauxilliaries_qrfactorisation_hpp */
