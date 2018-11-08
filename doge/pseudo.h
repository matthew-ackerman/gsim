#ifndef _PSEUDO_
#define _PSEUDO_

#include "phenotype.h"
#include "relatedness_data.h"
#include "types.h"

#include <omp.h>

#define EIGEN_NO_DEBUG

#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;

#define UINT uint32_t

struct pair
{
	MATRIX P;
	VECTOR s;
	bool pos_semi_def;
};

pair
get_pseudo(const MATRIX &);

pair
get_pseudo(const MATRIX &, const int &);

pair
get_pseudo_exact(const MATRIX &);

pair
get_inv(const MATRIX &);

namespace pseudo
{
 extern int drop_k;
}
#endif
