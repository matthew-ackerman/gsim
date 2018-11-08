#include "interface.h"
#include "map_file.h"
#include "state.h"
#include "sample_name.h"
#include "phenotype.h"
#include "relatedness_data.h"

#include <cstring>
#include <sstream>
#include <tuple>
#include <map>
#include <fstream>
#include <ctime>
#include <omp.h>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "pseudo.h"
#include "ll.h"
#include "sketch.h"

//#include "mpi.h"

using namespace Eigen;

#define UINT uint32_t
#define WORD 32

#define BLOCK 256

#define t(a)	a.transpose()
#define Av(a,b) VECTOR(a.array() * b.array() )
#define o(a,b) MATRIX( a.array() * b.array() )
#define d(a,b) MATRIX( a.array() / b.array() )
#define center(a,b) MATRIX( a.colwise() - b )
#define D(a,b) MATRIX( a.array().colwise() / b.array() )

//norm something?
#define NORM
//What is line_e?
#define LINE_E
#define MAX_TYPE	TYPE
//Use of Fastmath
#define FAST

std::string print_diag( Relatedness *r, const int &num_blocks)
{
	MATRIX MJM=r[0].Mt1*r[0].Mt1.transpose();
	MATRIX HJH=r[0].Ht1*r[0].Ht1.transpose();
	MATRIX MJH=(r[0].Mt1*r[0].Ht1.transpose()+r[0].Ht1*r[0].Mt1.transpose() )/2;

	return std::to_string(r[0].MtM.trace() )+", "+std::to_string(MJM.trace())+", "+std::to_string(r[0].HtH.trace())+", "+std::to_string(HJH.trace() )+", "+std::to_string(r[0].MtH.trace() )+", "+std::to_string(MJH.trace());
}

MATRIX qv(Matrix<TYPE, 6, 1> A, MATRIX B)
{
	Matrix<TYPE, 6, 1> ret;

	ret(1,0)=( A(1,0)*A(1,0) )*B(2,0);
	ret(3,0)=( A(3,0)*A(3,0) )*B(5,0);
	
	ret(0,0)=A(0,0)*B(0,0)-2*( A(0,0)*A(1,0)*A(1,0)/(A(0,0)+A(1,0)*A(1,0) ) )*B(1,0);
	ret(2,0)=A(2,0)*B(3,0)-2*( A(2,0)*A(3,0)*A(3,0)/(A(2,0)+A(3,0)*A(3,0) ) )*B(4,0);

	ret(4,0)=A(4,0)*B(6,0)-2*(A(4,0)*A(3,0)*A(1,0) )/(A(4,0)+A(3,0)*A(1,0) )*B(7,0)+(A(3,0)*A(1,0))*B(8,0);
	ret(5,0)=A(5,0)*B(9,0);
	return ret;
}

TYPE
q_norm(Matrix<TYPE, 6, 1> A)
{
	return fabs(A(0,0))
	+fabs(A(1,0))
	+fabs(A(2,0))
	+fabs(A(3,0))
	+A(4,0)
	+fabs(A(5,0));
}

struct ml
{
	VECTOR h;
	VECTOR se;
	MATRIX j;
};
//#define ∇ NB	 

MATRIX 
lm(Phenotype *p, Relatedness *r, const int &trait_n, const int &num_blocks)
{
	size_t *n=new size_t [num_blocks]; 
	size_t *sum=new size_t [num_blocks]; 
	size_t *sumsqr=new size_t [num_blocks]; 
	size_t nsqr=0;

	for (int z=0; z<num_blocks; z++)
	{
		n[z]=r[z].Mt1.rows();
		if(z>0)
			sum[z]=sum[z-1]+n[z-1];
		else 
			sum[z]=0;
		sumsqr[z]=nsqr;
		nsqr+=n[z]*n[z];
		//std::cerr << "lm:" << n[z] << std::endl;
	}

	size_t N=sum[num_blocks-1]+n[num_blocks-1];

//	std::cerr << "N: " << N << std::endl;
//	std::cerr << "n: " << n[0] << std::endl;
//	std::cerr << "n: " << n[1] << std::endl;

	MATRIX ret=MATRIX::Zero(6, 1);

	MATRIX x=MATRIX::Zero(N,2);
	MATRIX y=MATRIX::Zero(N,1);

	for (int z=0; z<num_blocks; z++)
	{
		x.block(sum[z],0,n[z],1)=r[z].Mt1;
		x.block(sum[z],1,n[z],1)=r[z].Ht1;
		y.block(sum[z],0,n[z],1)=p[z].value.block(0, trait_n,n[z],1);
	}
	
	MATRIX y2=y.array()-y.mean();

	std::cerr << "mean y: " << y2.mean() << std::endl;
	std::cerr << "mean x: " << x.mean() << std::endl;

	MATRIX V=( t(x)*x ).inverse()*t(x)*y2;

	std::cerr << "V: " << V << std::endl << std::endl;

	ret(1,0)=V(0,0);
	ret(3,0)=V(1,0);

	VECTOR cy=(y2-ret(1,0)*x.block(0,0,N,1)-ret(3,0)*x.block(0,1,N,1));

	MATRIX V2=( t(x)*x ).inverse()*t(x)*cy;
	
	std::cerr << "V2: " << V2 << std::endl << std::endl;

	cy=cy.array()-cy.mean();

	std::cerr << "Means: " << cy.mean() << ", " << y2.mean() << std::endl;

	MATRIX yyt=MATRIX::Zero(nsqr,1);

	for (int z=0; z<num_blocks; z++)
	{
		MATRIX temp=(cy.block(sum[z],0,n[z],1)*cy.block(sum[z],0,n[z],1).transpose() );
		yyt.block(sumsqr[z],0,n[z]*n[z],1)=Map<VECTOR>( temp.data(), n[z]*n[z]);
	}

	//std::cerr << yyt.block(nsqr-10,0,10,1) << std::endl;

	x=MATRIX::Zero(nsqr,4);

	for (int z=0; z<num_blocks; z++)
	{
		Eigen::Map<VECTOR> v1(r[z].MtM.data(), r[z].MtM.size());
		x.block(sumsqr[z],0,n[z]*n[z],1) = v1;
	}
	for (int z=0; z<num_blocks; z++)
	{
		Map<VECTOR> v2(r[z].HtH.data(), r[z].MtM.size());
		x.block(sumsqr[z],1,n[z]*n[z],1) = v2;
	}
	for (int z=0; z<num_blocks; z++)
	{
		Map<VECTOR> v3(r[z].MtH.data(), r[z].MtM.size());
		v3=v3*2.;
		x.block(sumsqr[z],2,n[z]*n[z],1) = v3;
	}
	for (int z=0; z<num_blocks; z++)
	{
		MATRIX v=makev(n[z]);
		Map<VECTOR> v4(v.data(), r[z].MtM.size());
		x.block(sumsqr[z],3,n[z]*n[z],1) = v4;
	}

	/*for (int z=0; z<nsqr; z++)
	{
	//	if (x(z,3)==1)
		{
			fprintf(stderr, "%f, %f, %f, %f, %f\n", yyt(z,0), x(z,0), x(z,1), x(z,2), x(z,3) );
		}
	}*/

	V=( t(x)*x ).inverse()*t(x)*yyt;

	ret(0,0)=V(0,0);//std::max(V(0,0), powf(10, -3) );
	ret(2,0)=V(1,0);//std::max(V(1,0), powf(10, -3) );
	ret(4,0)=V(2,0);
	ret(5,0)=V(3,0);//std::max(V(3,0), powf(10, -8) );

	ret=bus(ret);
	return ret;
}


struct pair2
{
	Matrix<TYPE, 5, 1> j;
	Matrix<TYPE, 5, 5> H;
};

struct trip
{
	MATRIX dS[5];	
	MATRIX dM[5];
	MATRIX ddS[5][5];	
};

MATRIX
get_drM(const MATRIX &h, const Relatedness &r, const size_t &x)
{
	switch(x)
	{
		case(1): 
			return -r.rMt1;
		case(3):
			return -r.rHt1;
		default:
			break;
	}
	MATRIX empty;
	return empty;
}

MATRIX
get_dlM(const MATRIX &h, const Relatedness &r, const size_t &x)
{
	switch(x)
	{
		case(1): 
			return -r.lMt1;
		case(3):
			return -r.lHt1;
		default:
			break;
	}
	MATRIX empty;
	return empty;
}

MATRIX
get_dM(const MATRIX &h, const Relatedness &r, const size_t &x)
{
	switch(x)
	{
		case(1): 
			return -r.Mt1;
		case(3):
			return -r.Ht1;
		default:
			break;
	}
	MATRIX empty;
	return empty;
}

/*
MATRIX
get_dS(const MATRIX &h, const Relatedness &r, const size_t &x)
{
	double a,b,b2,b3,b4,c,d,d2,d3,d4,e,f;
	size_t n=r.MtM.rows();
	a=h(0,0);
	b=h(1,0);
	c=h(2,0);
	d=h(3,0);
	e=h(4,0);
	f=h(5,0);
	b2=b*b;
	d2=d*d;
	b3=b2*b;
	d3=d2*d;
	b4=b3*b;
	d4=d3*d;

	MATRIX V=makev(n);
	MATRIX empty;

	switch(x)
	{
		case(0):
			return r.MtM*exp(a) + r.MJM*(-2*b2*exp(a)/(b2 + exp(a)) + 2*b2*exp(2*a)/pow((b2 + exp(a)),2) );
		case(1):
			return r.MJM*(4*b3*exp(a)/pow((b2 + exp(a)),2) - 4*b*exp(a)/(b2 + exp(a))) + r.MJH*(2*b*d2*e/pow((b*d + e),2) - 2*d*e/(b*d + e));
		case(2):
			return r.HtH*exp(c) + r.HJH*(-2*d2*exp(c)/(d2 + exp(c)) + 2*d2*exp(2*c)/pow((d2 + exp(c)),2) );
		case(3):
			return r.HJH*(4*d3*exp(c)/pow((d2 + exp(c)), 2) - 4*d*exp(c)/(d2 + exp(c))) + r.MJH*(2*b2*d*e/pow((b*d + e),2) - 2*b*e/(b*d + e));
		case(4):
			return r.MtH + r.MJH*(2*b*d*e/pow( (b*d + e), 2) - 2*b*d/(b*d + e));
		case(5):
			return V*exp(f);
		default:
			break;
	}
	return empty;
}*/


//To be called when make_outer() hasn't.
MATRIX
get_dS_small(const MATRIX &h, const Relatedness &r, const size_t &x)
{
	TYPE a,b,c,d,e,f;
	size_t n=r.MtM.rows();
	a=h(0,0);
	b=h(1,0);
	c=h(2,0);
	d=h(3,0);
	e=h(4,0);
	f=h(5,0);

	MATRIX empty;

	MATRIX v=makev(n);

	MATRIX MJM=r.Mt1*r.Mt1.transpose();
	MATRIX HJH=r.Ht1*r.Ht1.transpose();
	MATRIX MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;

	TYPE trMtM = r.tr(0,0);
	TYPE trMJM = r.tr(1,0);
	TYPE trHtH = r.tr(2,0);
	TYPE trHJH = r.tr(3,0);
	TYPE trMtH = r.tr(4,0);
	TYPE trMJH = r.tr(5,0);

switch(x)
{
	case(0):
		return MJM*(2*a*pow(b, 2)/pow(a + pow(b, 2), 2) - 2*pow(b, 2)/(a + pow(b, 2))) + r.MtM + v*(-2*a*pow(b, 2)*trMJM/pow(a + pow(b, 2), 2) + 2*pow(b, 2)*trMJM/(a + pow(b, 2)) - trMtM) ; // 0
	case(1):
		return MJH*(2*b*pow(d, 2)*e/pow(b*d + e, 2) - 2*d*e/(b*d + e)) + MJM*(4*a*pow(b, 3)/pow(a + pow(b, 2), 2) - 4*a*b/(a + pow(b, 2))) + v*(-4*a*pow(b, 3)*trMJM/pow(a + pow(b, 2), 2) + 4*a*b*trMJM/(a + pow(b, 2)) - 2*b*pow(d, 2)*e*trMJH/pow(b*d + e, 2) + 2*d*e*trMJH/(b*d + e)) ; // 1
	case(2):
		return HJH*(2*c*pow(d, 2)/pow(c + pow(d, 2), 2) - 2*pow(d, 2)/(c + pow(d, 2))) + r.HtH + v*(-2*c*pow(d, 2)*trHJH/pow(c + pow(d, 2), 2) + 2*pow(d, 2)*trHJH/(c + pow(d, 2)) - trHtH) ; // 2
	case(3):
		return HJH*(4*c*pow(d, 3)/pow(c + pow(d, 2), 2) - 4*c*d/(c + pow(d, 2))) + MJH*(2*pow(b, 2)*d*e/pow(b*d + e, 2) - 2*b*e/(b*d + e)) + v*(-2*pow(b, 2)*d*e*trMJH/pow(b*d + e, 2) + 2*b*e*trMJH/(b*d + e) - 4*c*pow(d, 3)*trHJH/pow(c + pow(d, 2), 2) + 4*c*d*trHJH/(c + pow(d, 2))) ; // 3
	case(4):
		return MJH*(2*b*d*e/pow(b*d + e, 2) - 2*b*d/(b*d + e)) + r.MtH + v*(-2*b*d*e*trMJH/pow(b*d + e, 2) + 2*b*d*trMJH/(b*d + e) - trMtH) ; // 4
	default:
		return empty;
}

	return empty;
}

MATRIX
get_ddS(const MATRIX &h, const Relatedness &r, const size_t &x, const size_t &y)
{
	TYPE a,b,b2,b3,b4,c,d,d2,d3,d4,e,e2,e3,f;
	size_t n=r.Mt1.rows();
	a=h(0,0);
	b=h(1,0);
	c=h(2,0);
	d=h(3,0);
	e=h(4,0);
	f=h(5,0);

	b2=b*b;
	d2=d*d;
	e2=e*e;
	b3=b2*b;
	d3=d2*d;
	e3=e2*e;
	b4=b3*b;
	d4=d3*d;

	MATRIX empty;

	MATRIX v=makev(n);

	MATRIX MJM=r.Mt1*r.Mt1.transpose();
	MATRIX HJH=r.Ht1*r.Ht1.transpose();
	MATRIX MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;

	TYPE trMtM = r.tr(0,0);
	TYPE trMJM = r.tr(1,0);
	TYPE trHtH = r.tr(2,0);
	TYPE trHJH = r.tr(3,0);
	TYPE trMtH = r.tr(4,0);
	TYPE trMJH = r.tr(5,0);

switch(x)
{
	case(0):
		switch(y)
		{
			case(0):
				return 4*pow(b, 2)*(MJM*(-a/(a + pow(b, 2)) + 1) + trMJM*v*(a/(a + pow(b, 2)) - 1))/pow(a + pow(b, 2), 2) ; // 0 0
			case(1):
				return 4*b*(MJM*(-2*a*pow(b, 2)/pow(a + pow(b, 2), 2) + a/(a + pow(b, 2)) + pow(b, 2)/(a + pow(b, 2)) - 1) + trMJM*v*(2*a*pow(b, 2)/pow(a + pow(b, 2), 2) - a/(a + pow(b, 2)) - pow(b, 2)/(a + pow(b, 2)) + 1))/(a + pow(b, 2)) ; // 0 1
			default:
				return empty;
		}
	case(1):
		switch(y)
		{
			case(1):
				return 4*MJH*(-b*pow(d, 3)*e/pow(b*d + e, 3) + pow(d, 2)*e/pow(b*d + e, 2)) + 4*MJM*(-4*a*pow(b, 4)/pow(a + pow(b, 2), 3) + 5*a*pow(b, 2)/pow(a + pow(b, 2), 2) - a/(a + pow(b, 2))) + 4*v*(4*a*pow(b, 4)*trMJM/pow(a + pow(b, 2), 3) - 5*a*pow(b, 2)*trMJM/pow(a + pow(b, 2), 2) + a*trMJM/(a + pow(b, 2)) + b*pow(d, 3)*e*trMJH/pow(b*d + e, 3) - pow(d, 2)*e*trMJH/pow(b*d + e, 2)) ; // 1 1
			case(3):
				return 2*e*(MJH*(-2*pow(b, 2)*pow(d, 2)/pow(b*d + e, 2) + 3*b*d/(b*d + e) - 1) + trMJH*v*(2*pow(b, 2)*pow(d, 2)/pow(b*d + e, 2) - 3*b*d/(b*d + e) + 1))/(b*d + e) ; // 1 3
			case(4):
				return 2*d*(MJH*(-2*b*d*e/pow(b*d + e, 2) + b*d/(b*d + e) + e/(b*d + e) - 1) + trMJH*v*(2*b*d*e/pow(b*d + e, 2) - b*d/(b*d + e) - e/(b*d + e) + 1))/(b*d + e) ; // 1 4
			default:
				return empty;
		}
	case(2):
		switch(y)
		{
			case(2):
				return 4*pow(d, 2)*(HJH*(-c/(c + pow(d, 2)) + 1) + trHJH*v*(c/(c + pow(d, 2)) - 1))/pow(c + pow(d, 2), 2) ; // 2 2
			case(3):
				return 4*d*(HJH*(-2*c*pow(d, 2)/pow(c + pow(d, 2), 2) + c/(c + pow(d, 2)) + pow(d, 2)/(c + pow(d, 2)) - 1) + trHJH*v*(2*c*pow(d, 2)/pow(c + pow(d, 2), 2) - c/(c + pow(d, 2)) - pow(d, 2)/(c + pow(d, 2)) + 1))/(c + pow(d, 2)) ; // 2 3
			default:
				return empty;
		}
	case(3):
		switch(y)
		{
			case(3):
				return 4*HJH*(-4*c*pow(d, 4)/pow(c + pow(d, 2), 3) + 5*c*pow(d, 2)/pow(c + pow(d, 2), 2) - c/(c + pow(d, 2))) + 4*MJH*(-pow(b, 3)*d*e/pow(b*d + e, 3) + pow(b, 2)*e/pow(b*d + e, 2)) + 4*v*(pow(b, 3)*d*e*trMJH/pow(b*d + e, 3) - pow(b, 2)*e*trMJH/pow(b*d + e, 2) + 4*c*pow(d, 4)*trHJH/pow(c + pow(d, 2), 3) - 5*c*pow(d, 2)*trHJH/pow(c + pow(d, 2), 2) + c*trHJH/(c + pow(d, 2))) ; // 3 3
			case(4):
				return 2*b*(MJH*(-2*b*d*e/pow(b*d + e, 2) + b*d/(b*d + e) + e/(b*d + e) - 1) + trMJH*v*(2*b*d*e/pow(b*d + e, 2) - b*d/(b*d + e) - e/(b*d + e) + 1))/(b*d + e) ; // 3 4
			default:
				return empty;
		}
	case(4):
		switch(y)
		{
			case(4):
				return 4*b*d*(MJH*(-e/(b*d + e) + 1) + trMJH*v*(e/(b*d + e) - 1))/pow(b*d + e, 2) ; // 4 4
			default:
				return empty;
		}
}


return empty;
}




Matrix<TYPE, 5, 1> 
nabla(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &trait_n)
#ifdef FAST
{
	clock_t t0, t1, t2, t3, t4;
	size_t n=r.Mt1.rows();
	t0=clock();
	MATRIX S=get_Sigma(h, r);
	t1=clock();
	pair a=get_pseudo(S);
	t2=clock();

	Matrix<TYPE, 5, 1> j=MATRIX::Zero(5,1);
	MATRIX cp=p.value.block(0,trait_n,n,1)-h(1,0)*r.Mt1-h(3,0)*r.Ht1;
	MATRIX ImSPcp=(MATRIX::Identity(n,n)-S*a.P)*cp;		//Orthogonal projector onto the kernel of S. . . 
	MATRIX Pc=a.P*cp;
	MATRIX tPc=t(Pc)*a.P;
	t3=clock();

	for (size_t x=0; x<5; x++)
	{
		MATRIX PdS=a.P*get_dS_small(h,r,x);
		TYPE dP2=(t(Pc)*PdS*ImSPcp)(0,0);   //4

		j(x,0)=(t(Pc)*get_dS_small(h,r,x)*Pc/2.)(0,0)+dP2-(2.*M_PI*PdS ).trace()/(4*M_PI) ; //5
		if (get_dM(h, r, x).cols()!=0 )
		{
			j(x,0)-= (t(get_dM(h, r, x))*Pc/2.)(0,0); //5
			j(x,0)-= (t(Pc)*(get_dM(h, r, x))/2.)(0,0); //5
		} 
	}
	t4=clock();
	//std::cerr << t1-t0 << ", " << t2-t1 << ", " << t3-t2 << ", " << t4-t3 << std::endl;
	return j*-1;
}
#else
{
	size_t n=r.MtM.rows();

	MATRIX S=get_Sigma(h, r);
	MATRIX P=get_pseudo(S).P;

	Matrix<TYPE, 5,1> j;

	MATRIX PP=P*t(P);	// PP(x,y) != P(x,y)*P(x,y)
	MATRIX SP=S*t(P);	// SP(x,y) != S(x,y)*P(x,y)


	MATRIX ImSP;		// S*t(p) is square. . . 
	ImSP=MATRIX::Identity(n,n)-SP; 		//Depends on where the diagonal is...

//	MATRIX ImPS=t(ImSP);					//Be careful with transposes too!

	VECTOR cp=p.value-h(1,0)*r.Mt1-h(3,0)*r.Ht1;

	for (size_t x=0; x<5; x++)
	{
		MATRIX dP=(-P*get_dS_small(h, r, x)*t(P) );	// dP  needs to be m x m
		MATRIX dP2=PP*get_dS_small(h, r, x)*ImSP;   // dP2 needs to be n x n 

		dP.noalias()+=dP2+t(dP2);

		j(x,0)=(-t(cp)*dP*(cp)/2.)(0,0);	// [X]
		j(x,0)-=(2.*M_PI*P*get_dS_small(h, r, x) ).trace()/(4*M_PI) ; // !=
		if (get_dlM(h, r, x).cols()!=0 )
		{
			j(x,0)-= (t(get_dM(h, r, x))*P*(cp)/2.)(0,0); // !=
			j(x,0)-= (t(cp)*P*(get_dM(h, r, x))/2.)(0,0); // !=
		} else {
			//
		}
	}
	return j*-1;
}
#endif

MATRIX
get_ddP(const MATRIX &P,
	const MATRIX &PP,
	const MATRIX &ImSP, 
	const MATRIX &D1P, 
	const MATRIX &D1ImSP, 
	const MATRIX &PPD1, 
	const MATRIX &dX2, 
	const MATRIX &dX2P, 
	const MATRIX &D2P, 
	const MATRIX &SdX2, 
	const MATRIX &D12)
{
	MATRIX dX2D1P=dX2*D1P;
	MATRIX mD2PmSdX2=-D2P-SdX2;
	//fprintf(stderr, "%d, dX2:%ldx%ld; D1ImSP:%ldx%ld; PPD1:%ldx%ld; mD2PmSdX2:%ldx%ld.\n", __LINE__, dX2P.rows(), dX2P.cols(), D1ImSP.rows(), D1ImSP.cols(), PPD1.rows(), PPD1.cols(), mD2PmSdX2.rows(), mD2PmSdX2.cols() );
	MATRIX RR = dX2P*D1ImSP;
	RR.noalias() += t(dX2P)*D1ImSP;
	RR.noalias() += PPD1*mD2PmSdX2;
	if (D12.cols() != 0)
	    RR.noalias() += PP*D12*ImSP;   // 

	MATRIX ret = ( -dX2D1P
              -t(dX2D1P)
		+RR
		+t(RR) );
	if(D12.cols() != 0)
              ret-=P*D12*P; //
	return ret;
}

Matrix<TYPE, 5,1> 
sub_nabla(const MATRIX &h, const VECTOR &lp, const VECTOR &rp,  const Relatedness &r, const bool &diagonal, const MATRIX &P, const MATRIX &S)
{
	size_t n=r.MtM.rows();

	Matrix<TYPE, 5,1> j;

	//let assume P is a n x m slice of P

	MATRIX PP=P*t(P);	// PP(x,y) != P(x,y)*P(x,y)
	MATRIX SP=S*t(P);	// SP(x,y) != S(x,y)*P(x,y)


	MATRIX ImSP;		// S*t(p) is square. . . 
	if (diagonal) ImSP=MATRIX::Identity(n,n)-SP; 		//Depends on where the diagonal is...
	else ImSP=-SP; 						//Depends on where the diagonal is...

//	MATRIX ImPS=t(ImSP);					//Be careful with transposes too!

	VECTOR lcp=lp-h(1,0)*r.lMt1-h(3,0)*r.lHt1;
	VECTOR rcp=rp-h(1,0)*r.rMt1-h(3,0)*r.rHt1;

	for (size_t x=0; x<5; x++)
	{

		MATRIX dP=(-P*get_dS_small(h, r, x)*t(P) );	// dP  needs to be m x m
		MATRIX dP2=PP*get_dS_small(h, r, x)*ImSP;   // dP2 needs to be n x n 

		dP.noalias()+=dP2+t(dP2);

		j(x,0)=(-t(lcp)*dP*(rcp)/2.)(0,0);	// [X]
		if (diagonal) j(x,0)-=(2.*M_PI*P*get_dS_small(h, r, x) ).trace()/(4*M_PI) ; // !=
		if (get_dlM(h, r, x).cols()!=0 )
		{
			j(x,0)-= (t(get_dlM(h, r, x))*P*(rcp)/2.)(0,0); // !=
			j(x,0)-= (t(lcp)*P*(get_drM(h, r, x))/2.)(0,0); // !=
		} else {
			//
		}
	}
	return j;//*-1;
}

Matrix<TYPE, 5, 1>
get_sub(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &step, const int &trait_n)
{
	size_t n=r.Mt1.rows();
	MATRIX S=get_Sigma(h, r);
	pair a=get_pseudo(S);
	
	Matrix<TYPE, 5, 1> ret=MATRIX::Zero(5,1);

	for (size_t x=0; x<n; x+=step)
	{
		int Nx;
		n-x > step ? Nx=step : Nx = n-x;
		for (size_t y=x; y<n; y+=step)
		{
			int Ny;
			n-y > step ? Ny=step : Ny = n-y;
			std::cerr << Nx << ", " << Ny << std::endl;
			if (x==y)
				ret  +=  sub_nabla(h, p.value.block(x,trait_n,Nx,1), p.value.block(y,trait_n,Ny,1), r, true,  a.P.block(0,y,n,Ny), S.block(0,y,n,Ny) );
			else
				ret += 2.*sub_nabla(h, p.value.block(x,trait_n,Nx,1), p.value.block(y,trait_n,Ny,1), r, false, a.P.block(0,y,n,Ny), S.block(0,y,n,Ny) );
		}
	}
	return ret;
}


Matrix<TYPE, 5, 5>
nabla2_forward(const Matrix<TYPE, 5, 1> &h, const Phenotype &p, const Relatedness &r, const Matrix<TYPE, 5, 1>k, const TYPE &delta, const int &trait_num)
{
	MATRIX H=MATRIX::Zero(5,5);

	for (size_t x=0; x<5; x++)
	{
		MATRIX dh=MATRIX::Zero(5,1);
		dh(x,0)=std::max(delta*fabs(h(x,0)), MAX_TYPE(0.000001) );
		Matrix<TYPE, 5, 1> k2=nabla(h+dh, p, r, trait_num);
		H.block(0,x,5,1)=(k2-k)/dh.norm();
	}

	return H;
}

Matrix<TYPE, 5, 5>
nabla2_center(const Matrix<TYPE, 5, 1> &h, const Phenotype &p, const Relatedness &r, const TYPE &delta, const int &trait_num)
{
	MATRIX H=MATRIX::Zero(5, 5);

	for (size_t x=0; x<5; x++)
	{
		MATRIX dh=MATRIX::Zero(5,1);
		//dh(x,0)=delta*fabs(h(x,0));
		dh(x,0)=std::max(delta*fabs(h(x,0)), MAX_TYPE(0.000001) );
		Matrix<TYPE, 5, 1> k1=nabla(h-dh, p, r, trait_num);
		Matrix<TYPE, 5, 1> k2=nabla(h+dh, p, r, trait_num);
		H.block(0,x,5,1)=(k2-k1)/(2.*dh.norm() );
	}

	return H;
}

pair2 
nablas(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &trait_num)
{
	size_t n=r.Mt1.rows();
	MATRIX S=get_Sigma(h, r);

//	std::cerr << "getting pseudo...";
	//THIS CANT JUST BE SUBDIVIDED
	pair a=get_pseudo(S);

//	std::cerr << "done in " << t1-s << std::endl;
//	std::cerr << get_ll(h,p,r,a) << std::endl;

	pair2 ret;

	MATRIX PP=a.P*a.P;
	MATRIX SP=S*a.P; 

	MATRIX ImSP=MATRIX::Identity(n,n)-S*a.P;
	MATRIX ImPS=t(ImSP);

	MATRIX cp=p.value.block(0,trait_num,n,1)-h(1,0)*r.Mt1-h(3,0)*r.Ht1;

	MATRIX dP[5];
	MATRIX dPP[5], DdP[5], PPdS[5], dSP[5], dSImSP[5];


	for (size_t x=0; x<5; x++)
	{
		MATRIX dSx=get_dS_small(h, r, x);
		dP[x]=(-a.P*dSx*a.P); //2
		MATRIX dP2=PP*dSx*ImSP;   //4
		dP[x].noalias() += dP2+t(dP2);
		dPP[x]=dP[x]*a.P; //i.e. dX2P=dPP[y]
		DdP[x]=S*dP[x]; //i.e. SdX2=DdP[y]
		PPdS[x]=PP*dSx; //i.e. PPD1=PPdS[x]
		dSP[x]=dSx*a.P; //i.e. D1P=dSP[x] and D2P=dSP[y]
		dSImSP[x]=dSx*ImSP;	//i.e. D1ISP=dSImSP[x]
	}

	for (size_t x=0; x<5; x++)
	{
		MATRIX dMx=get_dM(h,r,x);
		for (size_t y=x; y<5; y++)
		{

			MATRIX ddSxy=get_ddS(h,r,x,y);
			MATRIX dSy=get_dS_small(h,r,y);
			MATRIX dMy=get_dM(h,r,y);

			MATRIX ddPxy=get_ddP(a.P, PP, ImSP, dSP[x], dSImSP[x], PPdS[x], dP[y], dPP[y], dSP[y], DdP[y], ddSxy); 
			{
			TYPE d=(dP[x]*dSy).trace()/2.;
			if (ddSxy.cols()!=0)
				d+=(a.P*ddSxy).trace()/2.;
			ret.H(x,y)=d+( t(cp)*ddPxy*(cp)/2.)(0,0);
			}

			if (dMx.cols() != 0 )
			{
				TYPE d=(t(dMx)*dP[y]*(cp) )(0,0);
				ret.H(x,y) += d;//+t(x);
			}
			if (dMy.cols() != 0 )
			{
				TYPE d =(t(dMy)*dP[x]*(cp) )(0,0);
				ret.H(x,y) += d;
			}
			if (dMx.cols() != 0 &&  dMy.cols() !=0 )
			{
				TYPE d = (t(dMx)*a.P*(dMy) )(0,0);
				ret.H(x,y) += d;
			}
			if (x!=y) ret.H(y,x)=ret.H(x,y);
		}
		{
			MATRIX dSx=get_dS_small(h,r,x);
			ret.j(x,0)=(-t(cp)*dP[x]*(cp)/2.)(0,0)-(2.*M_PI*a.P*dSx ).trace()/(4*M_PI); //5
			if (dMx.cols()!=0 )
			{
				ret.j(x,0)-= (t(dMx)*a.P*(cp)/2.)(0,0); //5
				ret.j(x,0)-= (t(cp)*a.P*(dMx)/2.)(0,0); //5
			} else {
//				std::cerr << "dM[" << x << "] empty.\n" << std::endl;
			}
		}
	}
	return ret;
}

void
print_options(void)
{
	std::cerr << __FILE__;
#ifdef __FAST_MATH__
	std::cerr << ", FAST MATH";
#endif
	std::cerr << "\n";
}

int main (int argc, char **argv){

//	sketch();
//	exit(0);

	bool hessian=false, block=false;
	int trait_num;
	int num_blocks=1;
	TYPE loci=0;
	print_options();
	ll::print_options();	

	std::string rel_name="", pheno_name="", names_name="";

	Environment env;
	env.set_name("mle");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Get mle fit. Please direct questions to matthew.s.ackerman@gmail.com");
	env.optional_arg('i',"input",  rel_name,      "please .", ".");
	env.optional_arg('p',"pheno",  pheno_name,      "please .", ".");
	env.optional_arg('t',"trait",  trait_num,      "please .", ".");
	env.optional_arg('k',"dropk",  pseudo::drop_k,   "please .", ".");
	env.optional_arg('l',"loci",  loci,   "please .", ".");
	env.flag(	'b',"block",  	&block,   	&flag_set,	"please .", ".");
	env.flag(	'h',"hessian", &hessian,	&flag_set, 	"an error occurred while displaying the version message", "use hessian");		//DONE


//	env.positional_arg('n',"names",  names_name,      "please provide a filename.", "names of individuals in the sample.");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);
	
	if (!hessian) 
	{
		std::cerr << "Maximizing w/o hessian. Results may not converge.\n";
	} else {
		std::cerr << "Maximizing using hessian. Computation may proceed slowly.\n";
	}

	initParallel();
	setNbThreads(4);
	if (block) num_blocks=2;
	std::cerr << "using " << nbThreads( ) << " threads.\n";


	Relatedness *r=new Relatedness[num_blocks];

	std::cerr << "reading Relatedness..\n";
	for (int z=0; z<num_blocks; z++)
	{
		Flat_file <Relatedness> rel_file;
		std::cerr << "opening " << z << std::endl;
		rel_file.open(rel_name.c_str(), READ);
		r[z]=rel_file.read_header();
		rel_file.read(r[z]);
		rel_file.close();
		r[z].make_trace();
	}

	if (loci==0)
		loci=r[0].mean_sites();
	
	std::cerr << "diagonals: " << print_diag(r, num_blocks) << std::endl;

	std::cerr << "mean sites: " << loci << std::endl;

	Phenotype *p=new Phenotype[num_blocks];

	std::cerr << "reading Phenotypes..\n";
	for (int z=0; z<num_blocks; z++)
	{
		Flat_file <Phenotype> pheno_file;
		std::cerr << "opening " << z << std::endl;
		pheno_file.open(pheno_name.c_str(), READ);
		p[z]=pheno_file.read_header();
		while(pheno_file.table_is_open() )
		{
			pheno_file.read(p[z]);
		}
		pheno_file.close(); 
	}

	std::cerr << __LINE__ << "Files read " << std::endl;

	if (num_blocks==2)
	{
		int N=p[0].value.rows();
		Eigen::Matrix<int,Eigen::Dynamic,1>  drops=Matrix<int, Eigen::Dynamic, 1>::Zero(N/2,1);
		for (int z=0; z<N/2; z++)
			drops(z,0)=z;
		r[0].drop_missing(drops);
		p[0].drop_missing(drops);
		for (int z=0; z<N/2; z++)
			drops(z,0)=z+N/2;
		r[1].drop_missing(drops);
		p[1].drop_missing(drops);
	}

	Matrix<TYPE, 5, 1> h=MATRIX::Zero(5,1);

/*  10961
0.0235893
  9082.31
0.0349932
   8663.2
0.0369896 */

	for (int x=10; x<30; x++)
	{
		h(0,0)=0.01*x;
		h(1,0)=10000.;
		h(2,0)=0.01*x;
		h(3,0)=10000.;
		h(4,0)=0;

		get_ll(h, p[0], r[0], 0);
	}

	h=lm(p, r, trait_num, num_blocks);


	delete [] r;
	delete [] p;

	return 0;
}
