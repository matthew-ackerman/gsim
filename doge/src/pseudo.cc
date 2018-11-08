#include "pseudo.h"

namespace pseudo {
    int drop_k=0;
}

using namespace pseudo;

pair
get_inv(const MATRIX &S)
{
	pair ret=get_pseudo(S,0);
	return ret;
}

pair
get_pseudo(const MATRIX &S)
{
	pair ret=get_pseudo(S,1);
	return ret;
}

pair
get_pseudo_exact(const MATRIX &S)
{
	const int k=1;
	size_t n=S.rows();

	BDCSVD <MATRIX> svd(S, ComputeThinU | ComputeThinV);

	VECTOR singularValues_inv=svd.singularValues();

	for ( long i=0; i<n-1; ++i) 
	{
		singularValues_inv(i)=1.0/svd.singularValues()(i);
    	}

	for (long i=1;i<k+1; ++i)
		singularValues_inv(n-i)=0;

	pair ret;
	ret.pos_semi_def = true;

	for ( int i=0; i < n-k; ++i){ 
		for ( int j=0; j < n-k; ++j){
			if ( fabs(svd.matrixV()(i,j)-svd.matrixU()(i,j) ) > 0.0000001 ) 
			{
				ret.pos_semi_def = false;
				break;
			}
		}
	}

	ret.s = svd.singularValues();
	ret.P = (svd.matrixV()*singularValues_inv.asDiagonal()*svd.matrixU().transpose());
	return ret;
}

pair
get_pseudo(const MATRIX &S, const int &k)
{
	size_t n=S.rows();
	std::cerr << __FILE__ << ":" << __LINE__ << "getting pseudo.\n"; 

	BDCSVD<MATRIX> svd(S, ComputeThinU | ComputeThinV);

        Matrix <TYPE, 1, 1> lnD=MATRIX::Zero(1,1);
	VECTOR singularValues_inv=svd.singularValues();

	for ( long i=0; i<n-k; ++i) 
	{
		singularValues_inv(i)=1.0/svd.singularValues()(i);
		lnD(0,0)+=log(2*M_PI*svd.singularValues()(i) );
    	}

	for (long i=1;i<k+1; ++i)
		singularValues_inv(n-i)=0;

	for (long i=0;i<k; ++i)
		singularValues_inv(i)=0;

	std::cerr << svd.singularValues()(0) << std::endl;
	std::cerr << svd.singularValues()(n-1) << std::endl;
	std::cerr << svd.singularValues()(n-2) << std::endl;

	pair ret;
	ret.pos_semi_def = true;

	for ( int i=0; i < n-k; ++i){ 
		for ( int j=0; j < n-k; ++j){
			if ( fabs(svd.matrixV()(i,j)-svd.matrixU()(i,j) ) > 0.0000001 ) 
			{
				ret.pos_semi_def = false;
				break;
			}
		}
	}

	ret.s = svd.singularValues();
	ret.P = (svd.matrixV()*singularValues_inv.asDiagonal()*svd.matrixU().transpose());

	return ret;
}
