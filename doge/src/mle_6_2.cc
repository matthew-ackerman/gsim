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
#include "types.h"
#include "matrix.h"

//#include "mpi.h"

using namespace Eigen;

#define UINT uint32_t
#define WORD 32

#define BLOCK 256

//These should just agree about type

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
	VectorXd h;
	VectorXd se;
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

	MATRIX ret=MATRIX::Zero(6,1);

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
		yyt.block(sumsqr[z],0,n[z]*n[z],1) = Map<VECTOR>( temp.data(), n[z]*n[z]);
	}

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

	std::cerr << "ret: " << ret << std::endl;

	ret=bus(ret);
	return ret;
}


struct pair2
{
	Matrix<TYPE, 6, 1> j;
	Matrix<TYPE, 6, 6> H;
};

struct trip
{
	MATRIX dS[6];	
	MATRIX dM[6];
	MATRIX ddS[6][6];	
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
			return r.MJM*(4*b3*exp(a)/pow((b2 + exp(a)),2) - 4*b*exp(a)/(b2 + exp(a))) + r.MJH*(2*b*d2*e/pow((D1),2) - 2*d*e/(D1));
		case(2):
			return r.HtH*exp(c) + r.HJH*(-2*d2*exp(c)/(d2 + exp(c)) + 2*d2*exp(2*c)/pow((d2 + exp(c)),2) );
		case(3):
			return r.HJH*(4*d3*exp(c)/pow((d2 + exp(c)), 2) - 4*d*exp(c)/(d2 + exp(c))) + r.MJH*(2*b2*d*e/pow((D1),2) - 2*b*e/(D1));
		case(4):
			return r.MtH + r.MJH*(2*b*d*e/pow( (D1), 2) - 2*b*d/(D1));
		case(5):
			return V*exp(f);
		default:
			break;
	}
	return empty;
}*/

#define LINE_E

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

	TYPE D1=b*d+e;

switch(x)
{
	case(0):
		return MJM*(2*a*pow(b, 2)/pow(a + pow(b, 2), 2) - 2*pow(b, 2)/(a + pow(b, 2) ) ) + r.MtM ; // 0
	case(1):
		if (D1!=0)
			return MJH*(2*b*pow(d, 2)*e/pow(D1, 2) - 2*d*e/D1 ) + MJM*(4*a*pow(b, 3)/pow(a + pow(b, 2), 2) - 4*a*b/(a + pow(b, 2) ) ) ; // 1
		return MJM*(4*a*pow(b, 3)/pow(a + pow(b, 2), 2) - 4*a*b/(a + pow(b, 2) ) ) ; // 1
	case(2):
		return HJH*(2*c*pow(d, 2)/pow(c + pow(d, 2), 2) - 2*pow(d, 2)/(c + pow(d, 2) ) ) + r.HtH ; // 2
	case(3):
		if (D1!=0)
			return HJH*(4*c*pow(d, 3)/pow(c + pow(d, 2), 2) - 4*c*d/(c + pow(d, 2) ) ) + MJH*(2*pow(b, 2)*d*e/pow(D1, 2) - 2*b*e/D1 ) ; // 3
		return HJH*(4*c*pow(d, 3)/pow(c + pow(d, 2), 2) - 4*c*d/(c + pow(d, 2) ) ); // 3
	case(4):
		if (D1!=0)
			return MJH*(2*b*d*e/pow(D1, 2) - 2*b*d/(D1) ) + r.MtH ; // 4
		return r.MtH ; // 4
	case(5):
		return v; // 5
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

	MATRIX MJM=r.Mt1*r.Mt1.transpose();
	MATRIX HJH=r.Ht1*r.Ht1.transpose();
	MATRIX MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;

	TYPE D1=b*d+e;

switch(x)
{
	case(0):
		switch(y)
		{
			case(0):
				return 4*MJM*pow(b, 2)*(-a/(a + pow(b, 2)) + 1)/pow(a + pow(b, 2), 2) ; // 0 0
			case(1):
				return 4*MJM*b*(-2*a*pow(b, 2)/pow(a + pow(b, 2), 2) + a/(a + pow(b, 2)) + pow(b, 2)/(a + pow(b, 2)) - 1)/(a + pow(b, 2)) ; // 0 1
			default:
				return empty;
		}
	case(1):
		switch(y)
		{
			case(1):
				if(D1!=0)
					return 4*MJH*(-b*pow(d, 3)*e/pow(D1, 3) + pow(d, 2)*e/pow(D1, 2) ) + 4*MJM*(-4*a*pow(b, 4)/pow(a + pow(b, 2), 3) + 5*a*pow(b, 2)/pow(a + pow(b, 2), 2) - a/(a + pow(b, 2))) ; // 1 1
				return 4*MJM*(-4*a*pow(b, 4)/pow(a + pow(b, 2), 3) + 5*a*pow(b, 2)/pow(a + pow(b, 2), 2) - a/(a + pow(b, 2))) ; // 1 1
			case(3):
				if(D1!=0)
					return 2*MJH*e*(-2*pow(b, 2)*pow(d, 2)/pow(D1, 2) + 3*b*d/(D1) - 1)/(D1) ; // 1 3
				return empty; // 1 3
			case(4):
				if(D1!=0)
					return 2*MJH*d*(-2*b*d*e/pow(D1, 2) + b*d/(D1) + e/(D1) - 1)/(D1) ; // 1 4
				return empty; // 1 4
			default:
				return empty;
		}
	case(2):
		switch(y)
		{
			case(2):
				return 4*HJH*pow(d, 2)*(-c/(c + pow(d, 2)) + 1)/pow(c + pow(d, 2), 2) ; // 2 2
			case(3):
				return 4*HJH*d*(-2*c*pow(d, 2)/pow(c + pow(d, 2), 2) + c/(c + pow(d, 2)) + pow(d, 2)/(c + pow(d, 2)) - 1)/(c + pow(d, 2)) ; // 2 3
			default:
				return empty;
		}
	case(3):
		switch(y)
		{
			case(3):
				if (D1!=0)
					return 4*HJH*(-4*c*pow(d, 4)/pow(c + pow(d, 2), 3) + 5*c*pow(d, 2)/pow(c + pow(d, 2), 2) - c/(c + pow(d, 2))) + 4*MJH*(-pow(b, 3)*d*e/pow(D1, 3) + pow(b, 2)*e/pow(D1, 2)) ; // 3 3
				return 4*HJH*(-4*c*pow(d, 4)/pow(c + pow(d, 2), 3) + 5*c*pow(d, 2)/pow(c + pow(d, 2), 2) - c/(c + pow(d, 2))); // 3 3
			case(4):
				if (D1!=0)
					return 2*MJH*b*(-2*b*d*e/pow(D1, 2) + b*d/(D1) + e/(D1) - 1)/(D1) ; // 3 4
				return empty;
			default:
				return empty;
		}
	case(4):
		switch(y)
		{
			case(4):
				if (D1!=0)
					return 4*MJH*b*d*(-e/(D1) + 1)/pow(D1, 2) ; // 4 4
				return empty;
			default:
				return empty;
		}
	case(5):
		switch(y)
		{
			default:
				return empty;
		}
}

return empty;
}



#define FAST

Matrix<TYPE, 6, 1> 
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

	Matrix<TYPE, 6, 1> j=MATRIX::Zero(6,1);
	MATRIX cp=p.value.block(0,trait_n,n,1)-h(1,0)*r.Mt1-h(3,0)*r.Ht1;
	MATRIX ImSPcp=(MATRIX::Identity(n,n)-S*a.P)*cp;		//Orthogonal projector onto the kernel of S. . . 
	MATRIX Pc=a.P*cp;
	MATRIX tPc=t(Pc)*a.P;
	t3=clock();

	for (size_t x=0; x<6; x++)
	{
		MATRIX PdS=a.P*get_dS_small(h,r,x);
		TYPE dP2=(t(Pc)*PdS*ImSPcp)(0,0);   //4

		j(x,0)=(t(Pc)*get_dS_small(h,r,x)*Pc/2.)(0,0)+dP2;
		j(x,0)-=(2.*M_PI*PdS ).trace()/(4*M_PI) ; //5 //include the thing?
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

	Matrix<TYPE, 6,1> j;

	MATRIX PP=P*t(P);	// PP(x,y) != P(x,y)*P(x,y)
	MATRIX SP=S*t(P);	// SP(x,y) != S(x,y)*P(x,y)


	MATRIX ImSP;		// S*t(p) is square. . . 
	ImSP=MATRIX::Identity(n,n)-SP; 		//Depends on where the diagonal is...

//	MATRIX ImPS=t(ImSP);					//Be careful with transposes too!

	VectorXd cp=p.value-h(1,0)*r.Mt1-h(3,0)*r.Ht1;

	for (size_t x=0; x<6; x++)
	{
		MATRIX dP=(-P*get_dS_small(h, r, x)*t(P) );	// dP  needs to be m x m
		MATRIX dP2=PP*get_dS_small(h, r, x)*ImSP;   // dP2 needs to be n x n 

		dP.noalias()+=dP2+t(dP2);

		j(x,0)=(-t(cp)*dP*(cp)/2.)(0,0);	// [X]

		j(x,0)-=(2.*M_PI*P*get_dS_small(h, r, x) ).trace()/(4*M_PI) ; // include the thing?

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

Matrix<TYPE, 6,1> 
sub_nabla(const MATRIX &h, const VectorXd &lp, const VectorXd &rp,  const Relatedness &r, const bool &diagonal, const MATRIX &P, const MATRIX &S)
{
	size_t n=r.MtM.rows();

	Matrix<TYPE, 6,1> j;

	//let assume P is a n x m slice of P

	MATRIX PP=P*t(P);	// PP(x,y) != P(x,y)*P(x,y)
	MATRIX SP=S*t(P);	// SP(x,y) != S(x,y)*P(x,y)


	MATRIX ImSP;		// S*t(p) is square. . . 
	if (diagonal) ImSP=MATRIX::Identity(n,n)-SP; 		//Depends on where the diagonal is...
	else ImSP=-SP; 						//Depends on where the diagonal is...

//	MATRIX ImPS=t(ImSP);					//Be careful with transposes too!

	VectorXd lcp=lp-h(1,0)*r.lMt1-h(3,0)*r.lHt1;
	VectorXd rcp=rp-h(1,0)*r.rMt1-h(3,0)*r.rHt1;

	for (size_t x=0; x<6; x++)
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

Matrix<TYPE, 6, 1>
get_sub(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &step, const int &trait_n)
{
	size_t n=r.Mt1.rows();
	MATRIX S=get_Sigma(h, r);
	pair a=get_pseudo(S);
	
	Matrix<TYPE, 6,1> ret=MATRIX::Zero(6,1);

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

#define MAX_TYPE	double

Matrix<TYPE, 6, 6>
nabla2_forward(const Matrix<TYPE, 6, 1> &h, const Phenotype &p, const Relatedness &r, const Matrix<TYPE, 6, 1>k, const TYPE &delta, const int &trait_num)
{
	MATRIX H=MATRIX::Zero(6,6);

	for (size_t x=0; x<6; x++)
	{
		MATRIX dh=MATRIX::Zero(6,1);
		dh(x,0)=std::max(delta*fabs(h(x,0)), MAX_TYPE(0.000001) );
		Matrix<TYPE, 6, 1> k2=nabla(h+dh, p, r, trait_num);
		H.block(0,x,6,1)=(k2-k)/dh.norm();
	}

	return H;
}

Matrix<TYPE, 6, 6>
nabla2_center(const Matrix<TYPE, 6, 1> &h, const Phenotype &p, const Relatedness &r, const TYPE &delta, const int &trait_num)
{
	MATRIX H=MATRIX::Zero(6,6);

	for (size_t x=0; x<6; x++)
	{
		MATRIX dh=MATRIX::Zero(6,1);
		dh(x,0)=std::max(delta*fabs(h(x,0)), MAX_TYPE(0.000001) );
		Matrix<TYPE, 6, 1> k1=nabla(h-dh, p, r, trait_num);
		Matrix<TYPE, 6, 1> k2=nabla(h+dh, p, r, trait_num);
		H.block(0,x,6,1)=(k2-k1)/(2.*dh.norm() );
	}

	return H;
}

Matrix<TYPE, 6, 1>
nabla_center(const Matrix<TYPE, 6, 1> &h, const Phenotype &p, const Relatedness &r, const TYPE &delta, const int &trait_num)
{
	MATRIX j=MATRIX::Zero(6, 1);

	for (size_t x=0; x<6; x++)
	{
		MATRIX dh=MATRIX::Zero(6,1);
		dh(x,0)=std::max(delta*fabs(h(x,0)), MAX_TYPE(0.000001) );
		TYPE k1=get_ll(h-dh, p, r, trait_num);
		TYPE k2=get_ll(h+dh, p, r, trait_num);
		j(x,0)=(k2-k1)/(2.*dh.norm() );
	}
	return j;
}


pair2 
nablas(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &trait_num)
{
	size_t n=r.Mt1.rows();
	MATRIX S=get_Sigma(h, r);

//	std::cerr << "getting pseudo...";
	pair a=get_pseudo(S);
//	pair a=get_inv(S);

//	std::cerr << "done in " << t1-s << std::endl;
//	std::cerr << get_ll(h,p,r,a) << std::endl;

	pair2 ret;

	MATRIX PP=a.P*a.P;
	MATRIX SP=S*a.P; 

	MATRIX ImSP=MATRIX::Identity(n,n)-S*a.P;
	MATRIX ImPS=t(ImSP);

	MATRIX cp=p.value.block(0,trait_num,n,1)-h(1,0)*r.Mt1-h(3,0)*r.Ht1;

	MATRIX dP[6];
	MATRIX dPP[6], DdP[6], PPdS[6], dSP[6], dSImSP[6];


	for (size_t x=0; x<6; x++)
	{
		MATRIX dSx=get_dS_small(h, r, x);
		std::cerr << x << "=" << dSx.trace() << std::endl;
		dP[x]=(-a.P*dSx*a.P); //2
		MATRIX dP2=PP*dSx*ImSP;   //4
		dP[x].noalias() += dP2+t(dP2);
		dPP[x]=dP[x]*a.P; //i.e. dX2P=dPP[y]
		DdP[x]=S*dP[x]; //i.e. SdX2=DdP[y]
		PPdS[x]=PP*dSx; //i.e. PPD1=PPdS[x]
		dSP[x]=dSx*a.P; //i.e. D1P=dSP[x] and D2P=dSP[y]
		dSImSP[x]=dSx*ImSP;	//i.e. D1ISP=dSImSP[x]
	}

	for (size_t x=0; x<6; x++)
	{
		MATRIX dMx=get_dM(h,r,x);
		for (size_t y=x; y<6; y++)
		{

			MATRIX ddSxy=get_ddS(h,r,x,y);
			MATRIX dSy=get_dS_small(h,r,y);
			MATRIX dMy=get_dM(h,r,y);

			MATRIX ddPxy=get_ddP(a.P, PP, ImSP, dSP[x], dSImSP[x], PPdS[x], dP[y], dPP[y], dSP[y], DdP[y], ddSxy); 
			{
			double d=(dP[x]*dSy).trace()/2.;
			if (ddSxy.cols()!=0)
				d+=(a.P*ddSxy).trace()/2.;
				ret.H(x,y)=d+( t(cp)*ddPxy*(cp)/2.)(0,0);
			}

			if (dMx.cols() != 0 )
			{
				double d=(t(dMx)*dP[y]*(cp) )(0,0);
				ret.H(x,y) += d;//+t(x);
			}
			if (dMy.cols() != 0 )
			{
				double d =(t(dMy)*dP[x]*(cp) )(0,0);
				ret.H(x,y) += d;
			}
			if (dMx.cols() != 0 &&  dMy.cols() !=0 )
			{
				double d = (t(dMx)*a.P*(dMy) )(0,0);
				ret.H(x,y) += d;
			}
			if (x!=y) ret.H(y,x)=ret.H(x,y);
		}
		{
			MATRIX dSx=get_dS_small(h,r,x);
			ret.j(x,0)=(-t(cp)*dP[x]*(cp)/2.)(0,0)-(2.*M_PI*a.P*dSx ).trace()/(4*M_PI); //5
			if (dMx.cols()!=0 )
			{
				ret.j(x,0) -= (t(dMx)*a.P*(cp)/2.)(0,0); //5
				ret.j(x,0) -= (t(cp)*a.P*(dMx)/2.)(0,0); //5
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
	env.optional_arg('b',"block",  num_blocks, "please .", ".");
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
	Phenotype *p;
	Matrix<TYPE, 6, 1> norm=MATRIX::Zero(6,1);
	MATRIX K=MATRIX::Zero(6,1);
	MATRIX H_n;
	Matrix<TYPE, 6, 1> h=MATRIX::Zero(6,1);
	MATRIX iH=MATRIX::Zero(6,6);
	MATRIX j=MATRIX::Zero(6,1);

//	loop to here
//	for (int e=0; e<150; e++)
//	{

	std::cerr << "reading Relatedness..\n";
	for (int z=0; z<num_blocks; z++)
	{
		Flat_file <Relatedness> rel_file;
		std::cerr << "opening " << z << std::endl;
		rel_file.open(rel_name.c_str(), READ);
		r[z]=rel_file.read_header();
		rel_file.read(r[z]);
		rel_file.close();
	}

	if (loci==0)
		loci=r[0].mean_sites();
	
	p=new Phenotype[num_blocks];

	for (int z=0; z<num_blocks; z++)
	{
		Flat_file <Phenotype> pheno_file;
		pheno_file.open(pheno_name.c_str(), READ);
		p[z]=pheno_file.read_header();
		while(pheno_file.table_is_open() )
		{
			pheno_file.read(p[z]);
		}
		pheno_file.close(); 
	}

	if (num_blocks!=1)
	{
		srand(time(0));

		int N=p[0].value.rows();
		int U=double(N)/(double)(num_blocks) ;

		std::vector<int> list;
		for (int z=0; z<N; z++) list.push_back(z);
		std::random_shuffle (list.begin(), list.end() );
		Eigen::Matrix<int,Eigen::Dynamic,1>  drops=Matrix<int, Eigen::Dynamic, 1>::Zero(N-U, 1);


		for (int x=0; x<num_blocks; x++)
		{
			for (int z=0; z<x*U; z++)
				drops(z,0)=list[z];
			for (int z=(x+1)*U; z<N; z++)
				drops(z-U,0)=list[z];
			r[x].drop_missing(drops);
			p[x].drop_missing(drops);
		}
	}

	if (hessian)
	{
		for(int z=0; z<num_blocks; z++)
			r[z].make_outer_prod();
	}

/*
	{
		size_t N=0;

		for (int z=0; z<num_blocks; z++)
			N+=p[z].value.rows();

		MATRIX y=MATRIX::Zero(N, 1);

		size_t n=0;

		for (int z=0; z<num_blocks; z++)
		{
			y.block(n,0,p[z].value.rows(),1)=p[z].value.block(0,trait_num,p[z].value.rows(),1);
			n+=p[z].value.rows();
		}

		n=0;

		for (int z=0; z<num_blocks; z++)
		{
			p[z].value.block(0,trait_num,p[z].value.rows(),1)=y.block(n,0,p[z].value.rows(),1).array()-y.mean();
			n+=p[z].value.rows();
		}

		y=y.array()-y.mean();

		double var=TYPE( (y.adjoint() * y)(0,0) ) / TYPE(N - 1.);

#define NORM

#ifdef NORM
		for (int z=0; z<num_blocks;z++)
			p[z].value.block(0,trait_num,p[z].value.rows(),1).array()*=(1./sqrt(var));

		norm(0,0)=var;
		norm(1,0)=var;
		norm(2,0)=var;
		norm(3,0)=var;
		norm(4,0)=var;
		norm(5,0)=1.;
#else 
		norm(0,0)=1.;
		norm(1,0)=1.;
		norm(2,0)=1.;
		norm(3,0)=1.;
		norm(4,0)=1.;
		norm(5,0)=1.;
#endif
	}

	K(0,0)=0.579007;
	K(1,0)=317.897;
	K(2,0)=0.473052;
	K(3,0)= 0.886349;
	K(4,0)=-1.46819;
	K(5,0)=0.721093;

	std::cerr << "norm" << t(norm) << std::endl;
	
	h=d(K, norm);

	iH=MATRIX::Zero(6,6);
	j=MATRIX::Zero(6,1);


	for (int z=0; z<num_blocks; z++)
	{
		pair2 q=nablas(h, p[z], r[z], trait_num);	
		iH+=q.H;
		j+=q.j;
	}

	//std::cout << "H:" << iH << std::endl;
	std::cout << "j:" << t(j) << std::endl;
	
	H_n=iH.inverse();

	}

	delete [] p;

	return 0;
*/

	MATRIX S;

	MATRIX var_norm=MATRIX::Zero(9,1);

	{
		size_t N=0;

		for (int z=0; z<num_blocks; z++)
			N+=p[z].value.rows();

		MATRIX y=MATRIX::Zero(N, 1);

		size_t n=0;

		for (int z=0; z<num_blocks; z++)
		{
			y.block(n,0,p[z].value.rows(),1)=p[z].value.block(0,trait_num,p[z].value.rows(),1);
			n+=p[z].value.rows();
		}

		n=0;

		for (int z=0; z<num_blocks; z++)
		{
			p[z].value.block(0,trait_num,p[z].value.rows(),1)=y.block(n,0,p[z].value.rows(),1).array()-y.mean();
			n+=p[z].value.rows();
		}

		y=y.array()-y.mean();

		double var=TYPE( (y.adjoint() * y)(0,0) ) / TYPE(N - 1.);

#define NORM

#ifdef NORM
		for (int z=0; z<num_blocks;z++)
			p[z].value.block(0,trait_num,p[z].value.rows(),1).array()*=(1./sqrt(var));

		std::cerr << "use this var: " << var << std::endl;

		norm(0,0)=var;
		norm(1,0)=var;
		norm(2,0)=var;
		norm(3,0)=var;
		norm(4,0)=var;
		norm(5,0)=1.;
#else 
		norm(0,0)=1.;
		norm(1,0)=1.;
		norm(2,0)=1.;
		norm(3,0)=1.;
		norm(4,0)=1.;
		norm(5,0)=1.;
#endif
		var_norm=get_var_norm(r, num_blocks);
	}

//h: 2.03823 -887.418  -7.7534 -288.935 -3.03126  1.81161

	double sh0=0.00005;
	double sh1=1.;
	double sh2=0.00005;
	double sh3=1.;
	double sh4=0.00005;
	double sh5=0.1;
	
	MATRIX h_t;

	h=lm(p, r, trait_num, num_blocks);
	std::cerr << "h is: " << t(h) << std::endl;

	h(0,0)=0.000015;
	h(1,0)=0;
	h(2,0)=0.000015;
	h(3,0)=0;
	h(4,0)=0.000000001;
	h(5,0)=0.5;

//	std::cerr << t(nabla(h, p[0], r[0], trait_num) ) << std::endl;
//	std::cerr << t(nabla_center(h, p[0], r[0], 0.01, trait_num) ) << std::endl;

//	return 0;

	/*
	//h_t=h;
	for (int j=0; j<1; j++){
	for (int x0=2; x0<5; x0++){
	for (int x1=0; x1<5; x1++){
	for (int x2=0; x2<5; x2++){
	for (int x3=0; x3<5; x3++){
	for (int x4=0; x4<5; x4++){
	for (int x5=0; x5<5; x5++){
		h_t=h;
		h_t(0,0)=h(0,0)+sh0*double(x0-2);
		h_t(1,0)=h(1,0)+sh1*double(x1-2);
		h_t(2,0)=h(2,0)+sh2*double(x2-2);
		h_t(3,0)=h(3,0)+sh3*double(x3-2);
		h_t(4,0)=h(4,0)+sh4*double(x4-2);
		h_t(5,0)=h(5,0)+sh5*double(x5-2);
		std::cerr << "OUT:" << t(h_t) << '\t' << get_ll(h_t, p[0], r[0], trait_num) << std::endl; 
	}
	}
	}
	}
	}
	return 0;
	}
	}*/
	


	std::cerr << "lm: " << t( o(sub(h), norm) ) << std::endl;

	std::cout << "Un-normed : " << t(sub(h) ) << std::endl;
	std::cout << "LL: " << get_ll(h, p[0], r[0], trait_num) << std::endl;
	std::cerr << "vars: " << t(qv(h, var_norm) ) << std::endl;
	std::cerr << "Sigma trace: " << get_Sigma(h, r[0]).trace() << std::endl;

	std::cerr << "Should be 0: " << p[0].value.sum() << std::endl;
	TYPE this_var=0;
	for (int x=0; x<p[0].value.rows(); x++)
	{
		this_var+=powf(p[0].value(x,0), 2);
	}

	std::cerr << "var p: " << this_var << std::endl;
	std::cerr << "Normed: " <<  qv(h, var_norm).sum() << std::endl;

	std::cout << "LL: " << get_ll(h, p, r, trait_num, num_blocks) << std::endl;
	std::cerr << "h is: " << t(h) << std::endl;
	std::cerr << "var_norm is: " << t(var_norm) << std::endl;
	std::cerr << "vars: " << t(qv(h, var_norm) ) << std::endl;

	{
	std::cerr << "for h:" << t( o(h, norm) ) << std::endl;
	}

	MATRIX f_n, f_m, f_m2, h_m, h_n, h_o, df, dh, H_m;

	h_m=h;

	clock_t t0, t1;
	t0=clock();

	if (!hessian)
	{

	f_m=MATRIX::Zero(6,1);

	for (int x=0; x<num_blocks; x++)
	{
		f_m+=nabla(h_m, p[x], r[x], trait_num);
	}

	{

		H_m=MATRIX::Zero(6,6);

		for (int x=0; x<num_blocks; x++)
			H_m+=nabla2_forward(h, p[x], r[x], f_m, 0.05, trait_num);
		H_m=H_m.inverse();	
		std::cerr << "done: " << h_m.rows() << ", " << H_m.rows() << ", " << f_m.rows() << std::endl;
		h_n=h_m-H_m*f_m;
		std::cerr << f_m.squaredNorm() << "\n"; 
		std::cerr << __LINE__ << "f_m:" << t(o(f_m,norm)) << std::endl;
		std::cerr << __LINE__ << "h_m:" << t(o(h_m,norm)) << std::endl;
		std::cerr << __LINE__ << "h_n:" << t(o(h_n,norm)) << std::endl;
	}
	}

	for (size_t x=0; x<150; x++)
	{

	if (!hessian)
	{
		
		f_n=MATRIX::Zero(6,1);
		for (int z=0; z<num_blocks; z++)
			f_n+=nabla(h_n, p[z], r[z], trait_num);

		df=f_n-f_m;
		dh=h_n-h_m;
		
		H_n=H_m+(dh-H_m*df)*(t(dh)*H_m*df).inverse()*t(dh)*H_m;

		VectorXd step=(H_n*f_n);

		std::cerr << "step size (var):" << qv(step,var_norm).norm() << " : " << t(qv(step,var_norm) ) << std::endl;
		std::cerr << "f(theta):" << f_n << std::endl;

		if ( qv(step,var_norm).norm() < pow(10,-5) ) break;
		//step/=(1.+10.*pow(q_norm(qv(h_n-step,var_norm) )-1.,2) );
		std::cerr << pow(q_norm(qv(h_n-step,var_norm) )-1.,2) << std::endl;
		std::cerr << q_norm(qv(h_n-step,var_norm) )-1. << std::endl;
		std::cerr << qv(h_n-step,var_norm) << std::endl;
		std::cerr << "step size 2: " << qv(step,var_norm).norm() << " : " << t(qv(step,var_norm) ) << std::endl;

		h_o=h_n-step;

		f_m=f_n;
		h_m=h_n;
		h_n=h_o;
		H_m=H_n;

		t1=clock();
		std::cerr << x << ", " << t1-t0 << ", " << f_n.squaredNorm() << "\nh: " << t( o(sub(h_o), norm) ) << std::endl << std::endl;
//		std::cout << "LL:" << get_ll(h_n, p[0], r[0], trait_num) << std::endl;
		std::cout << "LL: ~" << get_ll(h_n, p, r, trait_num, num_blocks) << std::endl;

		std::cerr << "vars: " << t(qv(h_n, var_norm) ) << "="  << qv(h_n, var_norm).sum() << std::endl;

	} else {

		std::cerr << __LINE__ << std::endl;

		MATRIX iH=MATRIX::Zero(6,6);
		MATRIX j=MATRIX::Zero(6,1);

		for (int z=0; z<num_blocks; z++)
		{
			std::cerr << __LINE__ << std::endl;
			pair2 q=nablas(h, p[z], r[z], trait_num);	
			for (int z=0; z<num_blocks; z++)
			{
			}
			std::cerr << "H rows:" << q.H.rows() << std::endl;
			std::cerr << "H cols:" << q.H.cols() << std::endl;
			iH+=q.H;
			j+=q.j;
		}

		std::cerr << "H:" << iH << std::endl;
		std::cout << "j:" << j << std::endl;
	
		H_n=iH;
		iH=iH.inverse();

		//std::cerr << __LINE__ << H_n << std::endl << std::endl;
		//std::cerr << __LINE__ << iH << std::endl << std::endl;
		//std::cerr << __LINE__ << nabla2_forward(h, p[0], r[0], j, 0.001, trait_num) << std::endl;

		VECTOR D=iH.diagonal();
		VECTOR step=(iH*j);
		std::cerr << "Inverse:" << iH << std::endl;
		std::cerr << "vars: " << t(qv(h, var_norm) ) << "="  << qv(h, var_norm).sum() << std::endl;
		std::cerr << "step size (before):" << step.norm() << std::endl;
		if ( step.norm() < pow(10,-12) ) break;
		step/=std::min(20., (1.+5.*pow(q_norm(qv(h-step,var_norm) )-1.,2) ) );
		std::cerr << "step size (after):" << step.norm() << std::endl;
		//step/=log(pow(step.norm(),2)+exp(1) );
		h+=(step);

		t1=clock();
		std::cerr << x << ", " << t1-t0 << ", " << j.squaredNorm() << std::endl << "j:" << t(j) << std::endl << "h:" << t( o(sub(h), norm) ) << std::endl << std::endl;	
		std::cout << "LL: ~" << get_ll(h, p, r, trait_num, num_blocks) << std::endl;

		}
	}

	if(!hessian)
	{
		H_n=MATRIX::Zero(6,6);
		for (int x=0; x<num_blocks; x++)
			H_n=nabla2_center(h_n,p[x],r[x], 0.005, trait_num);
		H_n=H_n.inverse();
	}
	else
	{
		h_n=h;
		H_n=H_n.inverse();
	}


	MATRIX v0=qv(h_n-VectorXd(ArrayXd(H_n.diagonal()).abs().sqrt())*2, var_norm);
	MATRIX v1=qv(h_n, var_norm);
	MATRIX v2=qv(h_n+VectorXd(ArrayXd(H_n.diagonal()).abs().sqrt())*2, var_norm);

	std::cout << var_norm << std::endl;
	std::cout << norm << std::endl;
	std::cout << "text, :, , a, A, d, D, r, e" << std::endl;
	std::cout << "vars 0 : " << t(v0) << std::endl;
	std::cout << "vars 1 : " << t(v1) << std::endl;
	std::cout << "vars 2 : " << t(v2) << std::endl;
	std::cout << "h0 : " << t(o(sub(h_n-VectorXd(ArrayXd(H_n.diagonal()).abs().sqrt())*2 ), norm) ) << std::endl;
	std::cout << "h1 : " << t(o(sub(h_n), norm) ) << std::endl;
 	std::cout << "h2 : " << t(o(sub(h_n+VectorXd(ArrayXd(H_n.diagonal()).abs().sqrt())*2 ), norm) ) << std::endl;
	std::cout << "Un-normed : " << t(sub(h_n) ) << std::endl;

	delete [] r;
	delete [] p;

	return 0;
}
