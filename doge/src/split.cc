#include "sketch.h"
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "type.h"

using namespace Eigen;

#define	t(A)	A.transpose()
#define s(A)	(A)(0,0)

//
MATRIX C
split(MATRIX A, MATRIX B)
{
	size_t n=A.rows();
	size_t a=1;
	Matrix C=MATRIX::Zero(n,n);
	for (int x=0; x<n; x+=a)
		for (int y=x; y<n; y+=a)
		{
			C.block(x,y,a,a)=A.block(x,0,n,a)*B.block(0,y,a,n);
			C.block(y,x,a,a)=A.block(y,0,n,a)*B.block(0,x,a,n);
			//C.block(y,x,a,a)=t(C.block(x,y,a,a));
		}
}

