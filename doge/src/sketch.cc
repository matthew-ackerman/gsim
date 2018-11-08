#include "sketch.h"
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"

using namespace Eigen;

#define s(a)	(a)(0,0)

void
sketch(void)
{
	MATRIX B, A=MATRIX::Random(3,3);
	A=A*t(A);
	MATRIX x1, x0=MATRIX::Random(3,1);
	double l0;

	for (int x=0; x<12; x++)
	{
 		x1=A*x0;
		std::cerr << x1/x1.norm() << ", " << s(t(x1)*x0)/s(t(x0)*x0) << std::endl;
		l0=s(t(x1)*x0)/s(t(x0)*x0);
		x0=x1/x1.norm();//(t(x1)*x1)(0,0);
	}
	B=A-(l0/x1.squaredNorm())*x1*t(x1);
	std::cerr << B;
	x0=MATRIX::Random(3,1);
	for (int x=0; x<12; x++)
	{
 		x1=B*x0;
		std::cerr << x1/x1.norm() << ", " << s(t(x1)*x0)/s(t(x0)*x0) << std::endl;
		l0=s(t(x1)*x0)/s(t(x0)*x0);
		x0=x1/x1.norm();//(t(x1)*x1)(0,0);
	}
	B=B-(l0/x1.squaredNorm())*x1*t(x1);

//	std::cerr << B << std::cerr;

	EigenSolver<MATRIX> es(A);
	std::cerr << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
	std::cerr << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl;

	EigenSolver<MATRIX> es2(B);

	std::cerr << "The eigenvalues of B are:" << std::endl << es2.eigenvalues() << std::endl;
	std::cerr << "The matrix of eigenvectors, B, is:" << std::endl << es2.eigenvectors() << std::endl;
}

