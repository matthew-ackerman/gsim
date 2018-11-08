#ifndef _LL_
#define _LL_

#include "phenotype.h"
#include "relatedness_data.h"
#include "pseudo.h"
#include "types.h"

#include <omp.h>

#include "Eigen/Core"
#include "Eigen/Dense"

#include <gsl/gsl_sf.h>

using namespace Eigen;

TYPE
get_ll(const MATRIX &, const Phenotype *, const Relatedness *, const int &, const int &);

TYPE
get_ll(const MATRIX &, const Phenotype &, const Relatedness &, const int &);

TYPE
get_ll(const MATRIX &, const Phenotype &, const Relatedness &, const pair &, const int &);

TYPE
get_ll_exact(const MATRIX &, const Phenotype &, const Relatedness &, const int &);

MATRIX 
makev(const size_t &);

MATRIX 
get_Sigma(const MATRIX &, const Relatedness &);

Matrix<TYPE, 6, 1>
get_vars(const MATRIX &, const Relatedness *, const int &);

MATRIX
get_var_norm(const Relatedness *, const int &);

Matrix<TYPE, 6, 1>
get_se(const MATRIX &, const MATRIX &, const Relatedness *, const int &);

Matrix<TYPE, 6, 1> 
sub5(const Matrix<TYPE, 5, 1> &, const Matrix<TYPE, 7, 1> &);

Matrix<TYPE, 5, 1> 
bus5(const Matrix<TYPE, 6, 1> &);

Matrix<TYPE, 6, 1> 
sub(const Matrix<TYPE, 6, 1> &);

Matrix<TYPE, 6, 1> 
bus(const Matrix<TYPE, 6, 1> &);

Matrix<TYPE, 7, 1> 
sub7(const Matrix<TYPE, 7, 1> &);

Matrix<TYPE, 7, 1> 
bus7(const Matrix<TYPE, 7, 1> &);

namespace ll
{
	void print_options(void);
}

#endif
