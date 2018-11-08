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

using namespace Eigen;

#define UINT uint32_t
#define WORD 32

#define BLOCK 256

#define t(a)	a.transpose()
#define Av(a,b) VectorXf(a.array() * b.array() )
#define Am(a,b) MATRIX( a.array() * b.array() )
#define center(a,b) MATRIX( a.colwise() - b )
#define D(a,b) MATRIX( a.array().colwise() / b.array() )

/*
Source	Variance	SE
Vp(G)	835.201822	98.811209
V(e)	547.330978	88.294858
Vp	1382.532800	38.509201
V(G)/Vp	0.604110	0.065580
logL	-11210.978
logL0	-11260.607
LRT	99.258
df	1
Pval	0.0000e+00
n	2733
*/

struct ml
{
	VectorXf h;
	VectorXf se;
	MATRIX j;
};
//#define ∇ NB	 

MATRIX 
lm(Phenotype &p, Relatedness &r, const int &trait_n)
{
	size_t n=r.Mt1.rows();
	std::cerr << n << std::endl;
	size_t nsqr=n*n;

	MATRIX ret=MATRIX::Zero(6,1);
	MATRIX x=MATRIX::Zero(n,2);

	x.block(0,0,n,1)=r.Mt1;
	x.block(0,1,n,1)=r.Ht1;

	MATRIX V=( t(x)*x ).inverse()*t(x)*p.value.block(0,trait_n,n,1);

	ret(1,0)=V(0,0);
	ret(3,0)=V(1,0);

	VectorXf cy=(p.value.block(0,trait_n,n,1)-ret(1,0)*r.Mt1-ret(3,0)*r.Ht1);
	//MATRIX cy=(p.value-ret(1,0)*r.Mt1-ret(3,0)*r.Ht1);
	MATRIX Myyt=(cy*cy.transpose() );
	Map<RowVectorXf> yyt(Myyt.data(), Myyt.size() );
	MATRIX v=makev(n);

	x=MATRIX::Zero(nsqr,4);
	Eigen::Map<RowVectorXf> v1(r.MtM.data(), r.MtM.size());
	x.block(0,0,nsqr,1) = v1.transpose();
	Map<RowVectorXf> v2(r.HtH.data(), r.MtM.size());
	x.block(0,1,nsqr,1) = v2.transpose();
	Map<RowVectorXf> v3(r.MtH.data(), r.MtM.size());
	x.block(0,2,nsqr,1) = v3.transpose();
	Map<RowVectorXf> v4(v.data(), r.MtM.size());
	x.block(0,3,nsqr,1) = v4.transpose();
	V=( t(x)*x ).inverse()*t(x)*yyt.transpose();
	ret(0,0)=std::max(V(0,0), powf(10, -3) );
	ret(2,0)=std::max(V(1,0), powf(10, -3) );
	ret(4,0)=V(2,0);
	ret(5,0)=std::max(V(3,0), powf(10, -3) );

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


trip 
get_partials(const MATRIX &h, const Relatedness &r)
{
	double a,b,b2,b3,b4,c,d,d2,d3,d4,e,f;
	size_t n=r.Mt1.rows();
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

	trip q;

	q.dM[1]=-r.Mt1;
	q.dM[3]=-r.Ht1;

	MATRIX V=makev(n);

	q.dS[0] = r.MtM*exp(a) + r.MJM*(-2*b2*exp(a)/(b2 + exp(a)) + 2*b2*exp(2*a)/pow((b2 + exp(a)),2) );
	q.dS[1] = r.MJM*(4*b3*exp(a)/pow((b2 + exp(a)),2) - 4*b*exp(a)/(b2 + exp(a))) + r.MJH*(2*b*d2*e/pow((b*d + e),2) - 2*d*e/(b*d + e));
	q.dS[2] = r.HtH*exp(c) + r.HJH*(-2*d2*exp(c)/(d2 + exp(c)) + 2*d2*exp(2*c)/pow((d2 + exp(c)),2) );
	q.dS[3] = r.HJH*(4*d3*exp(c)/pow((d2 + exp(c)), 2) - 4*d*exp(c)/(d2 + exp(c))) + r.MJH*(2*b2*d*e/pow((b*d + e),2) - 2*b*e/(b*d + e));
	q.dS[4] = r.MtH + r.MJH*(2*b*d*e/pow( (b*d + e), 2) - 2*b*d/(b*d + e));
	q.dS[5] = V*exp(f);

	q.ddS[0][0] = (r.MtM + r.MJM*(-2*b2/(b2 + exp(a)) + 6*b2*exp(a)/pow((b2 + exp(a)),2) - 4*b2*exp(2*a)/pow((b2 + exp(a)),3) ))*exp(a);
	q.ddS[0][1] = 4*r.MJM*b*(b2/(b2 + exp(a)) - 2*b2*exp(a)/pow((b2 + exp(a)),2) - 1 + exp(a)/(b2 + exp(a)))*exp(a)/(b2 + exp(a));
	q.ddS[1][1] = 4*r.MJM*(-4*b4*exp(a)/pow((b2 + exp(a)),3) + 5*b2*exp(a)/pow((b2 + exp(a)),2) - exp(a)/(b2 + exp(a))) + 4*r.MJH*(-b*d3*e/pow((b*d + e),3) + d2*e/pow((b*d + e),2));
	q.ddS[1][3] = 2*r.MJH*e*(-2*b2*d2/pow((b*d + e),2) + 3*b*d/(b*d + e) - 1)/(b*d + e);
	q.ddS[1][4] = 2*r.MJH*d*(-2*b*d*e/pow((b*d + e),2) + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e);
	q.ddS[2][2] = (r.HtH + r.HJH*(-2*d2/(d2 + exp(c)) + 6*d2*exp(c)/pow((d2 + exp(c)),2) - 4*d2*exp(2*c)/pow((d2 + exp(c)),3)))*exp(c);
	q.ddS[2][3] = 4*r.HJH*d*(d2/(d2 + exp(c)) - 2*d2*exp(c)/pow((d2 + exp(c)),2) - 1 + exp(c)/(d2 + exp(c)))*exp(c)/(d2 + exp(c));
	q.ddS[3][3] = 4*r.HJH*(-4*d4*exp(c)/pow((d2 + exp(c)),3) + 5*d2*exp(c)/pow((d2 + exp(c)),2) - exp(c)/(d2 + exp(c))) + 4*r.MJH*(-b3*d*e/pow((b*d + e),3) + b2*e/pow((b*d + e),2));
	q.ddS[3][4] = 2*r.MJH*b*(-2*b*d*e/pow((b*d + e),2) + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e);
	q.ddS[4][4] = 4*r.MJH*b*d*(-e/(b*d + e) + 1)/pow((b*d + e),2);
	q.ddS[5][5] = V*exp(f);

	return q;
}

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
}

//To be called when make_outer() hasn't.
MATRIX
get_dS_small(const MATRIX &h, const Relatedness &r, const size_t &x)
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
	MATRIX MJM, MJH, HJH;
	switch(x)
	{
		case(0):
			MJM=r.Mt1*r.Mt1.transpose();
			return r.MtM*exp(a) + MJM*(-2*b2*exp(a)/(b2 + exp(a)) + 2*b2*exp(2*a)/pow((b2 + exp(a)),2) );
		case(1):
			MJM=r.Mt1*r.Mt1.transpose();
			MJH=r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose();
			return MJM*(4*b3*exp(a)/pow((b2 + exp(a)),2) - 4*b*exp(a)/(b2 + exp(a))) + MJH*(2*b*d2*e/pow((b*d + e),2) - 2*d*e/(b*d + e));
		case(2):
			HJH=r.Ht1*r.Ht1.transpose();
			return r.HtH*exp(c) + HJH*(-2*d2*exp(c)/(d2 + exp(c)) + 2*d2*exp(2*c)/pow((d2 + exp(c)),2) );
		case(3):
			HJH=r.Ht1*r.Ht1.transpose();
			MJH=r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose();
			return HJH*(4*d3*exp(c)/pow((d2 + exp(c)), 2) - 4*d*exp(c)/(d2 + exp(c))) + MJH*(2*b2*d*e/pow((b*d + e),2) - 2*b*e/(b*d + e));
		case(4):
			MJH=r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose();
			return r.MtH + r.MJH*(2*b*d*e/pow( (b*d + e), 2) - 2*b*d/(b*d + e));
		case(5):
			return V*exp(f);
		default:
			break;
	}
	return empty;
}

MATRIX
get_ddS(const MATRIX &h, const Relatedness &r, const size_t &x, const size_t &y)
{
	double a,b,b2,b3,b4,c,d,d2,d3,d4,e,f;
	size_t n=r.Mt1.rows();
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
			switch (y)
			{
				case(0):
					return (r.MtM + r.MJM*(-2*b2/(b2 + exp(a)) + 6*b2*exp(a)/pow((b2 + exp(a)),2) - 4*b2*exp(2*a)/pow((b2 + exp(a)),3) ))*exp(a);
				case(1):
					return 4*r.MJM*b*(b2/(b2 + exp(a)) - 2*b2*exp(a)/pow((b2 + exp(a)),2) - 1 + exp(a)/(b2 + exp(a)))*exp(a)/(b2 + exp(a));
				default:
					return empty;
			}
			return empty;
		case(1):
			switch (y)
			{
				case(1):
					return 4*r.MJM*(-4*b4*exp(a)/pow((b2 + exp(a)),3) + 5*b2*exp(a)/pow((b2 + exp(a)),2) - exp(a)/(b2 + exp(a))) + 4*r.MJH*(-b*d3*e/pow((b*d + e),3) + d2*e/pow((b*d + e),2));
				case(3):
					return 2*r.MJH*e*(-2*b2*d2/pow((b*d + e),2) + 3*b*d/(b*d + e) - 1)/(b*d + e);
				case(4):
					return 2*r.MJH*d*(-2*b*d*e/pow((b*d + e),2) + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e);
				default:
					return empty;
			}
		case(2):
			switch (y)
			{
				case(2):
					return (r.HtH + r.HJH*(-2*d2/(d2 + exp(c)) + 6*d2*exp(c)/pow((d2 + exp(c)),2) - 4*d2*exp(2*c)/pow((d2 + exp(c)),3)))*exp(c);
				case(3):
					return 4*r.HJH*d*(d2/(d2 + exp(c)) - 2*d2*exp(c)/pow((d2 + exp(c)),2) - 1 + exp(c)/(d2 + exp(c)))*exp(c)/(d2 + exp(c));
				default:
					return empty;
			}
		case(3):
			switch (y)
			{
				case(3):
					return 4*r.HJH*(-4*d4*exp(c)/pow((d2 + exp(c)),3) + 5*d2*exp(c)/pow((d2 + exp(c)),2) - exp(c)/(d2 + exp(c))) + 4*r.MJH*(-b3*d*e/pow((b*d + e),3) + b2*e/pow((b*d + e),2));
				case(4):
					return 2*r.MJH*b*(-2*b*d*e/pow((b*d + e),2) + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e);
			}
		case(4):
			if(y==4) return 4*r.MJH*b*d*(-e/(b*d + e) + 1)/pow((b*d + e),2);
			return empty;
		case(5):
			if(y==5) return  V*exp(f);
			return empty;
		default:
			break;
	}
	return empty;
}

Matrix<TYPE, 6, 1> 
nabla(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &trait_n)
/*{
	size_t n=r.Mt1.rows();
	MATRIX S=get_Sigma(h, r);
	pair a=get_pseudo(S);

	Matrix<TYPE, 6, 1> j=MATRIX::Zero(6,1);
	MATRIX cp=p.value.block(0,trait_n,n,1)-h(1,0)*r.Mt1-h(3,0)*r.Ht1;
	MATRIX ImSP=MATRIX::Identity(n,n)-S*a.P;		//Orthogonal projector onto the kernel of S. . . 
	MATRIX Pc=a.P*cp;
	MATRIX tPcP=t(Pc)*a.P;
	MATRIX ImSPcp=ImSP*cp;

	for (size_t x=0; x<6; x++)
	{
		TYPE dP2=(tPcP*get_dS(h, r, x)*ImSPcp)(0,0);   //4

		j(x,0)=(t(Pc)*get_dS(h,r,x)*Pc/2.)(0,0)+dP2-(2.*M_PI*a.P*get_dS(h, r, x) ).trace()/(4*M_PI) ; //5
		if (get_dM(h, r, x).cols()!=0 )
		{
			j(x,0)-= (t(get_dM(h, r, x))*Pc/2.)(0,0); //5
			j(x,0)-= (t(Pc)*(get_dM(h, r, x))/2.)(0,0); //5
		} else {
			//
		}
	}
	return j*-1;
}*/
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

	VectorXf cp=p.value-h(1,0)*r.Mt1-h(3,0)*r.Ht1;

	for (size_t x=0; x<6; x++)
	{
		MATRIX dP=(-P*get_dS(h, r, x)*t(P) );	// dP  needs to be m x m
		MATRIX dP2=PP*get_dS(h, r, x)*ImSP;   // dP2 needs to be n x n 

		dP.noalias()+=dP2+t(dP2);

		j(x,0)=(-t(cp)*dP*(cp)/2.)(0,0);	// [X]
		j(x,0)-=(2.*M_PI*P*get_dS(h, r, x) ).trace()/(4*M_PI) ; // !=
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
sub_nabla(const MATRIX &h, const VectorXf &lp, const VectorXf &rp,  const Relatedness &r, const bool &diagonal, const MATRIX &P, const MATRIX &S)
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

	VectorXf lcp=lp-h(1,0)*r.lMt1-h(3,0)*r.lHt1;
	VectorXf rcp=rp-h(1,0)*r.rMt1-h(3,0)*r.rHt1;

	for (size_t x=0; x<6; x++)
	{
		std::cerr << "HI:" << __LINE__ << std::endl;

		MATRIX dP=(-P*get_dS(h, r, x)*t(P) );	// dP  needs to be m x m
		MATRIX dP2=PP*get_dS(h, r, x)*ImSP;   // dP2 needs to be n x n 

		dP.noalias()+=dP2+t(dP2);

		j(x,0)=(-t(lcp)*dP*(rcp)/2.)(0,0);	// [X]
		if (diagonal) j(x,0)-=(2.*M_PI*P*get_dS(h, r, x) ).trace()/(4*M_PI) ; // !=
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

Matrix<TYPE, 6, 6>
nabla2_forward(const Matrix<TYPE, 6, 1> &h, const Phenotype &p, const Relatedness &r, const Matrix<TYPE, 6, 1>k, const TYPE &delta, const int &trait_num)
{
	MATRIX H=MATRIX::Zero(6,6);

	for (size_t x=0; x<6; x++)
	{
		MATRIX dh=MATRIX::Zero(6,1);
		dh(x,0)=delta*fabs(h(x,0));
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
		dh(x,0)=delta*fabs(h(x,0));
		Matrix<TYPE, 6, 1> k1=nabla(h-dh, p, r, trait_num);
		Matrix<TYPE, 6, 1> k2=nabla(h+dh, p, r, trait_num);
		H.block(0,x,6,1)=(k2-k1)/(2.*dh.norm() );
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

	MATRIX dP[6];
	MATRIX dPP[6], DdP[6], PPdS[6], dSP[6], dSImSP[6];


	for (size_t x=0; x<6; x++)
	{
		MATRIX dSx=get_dS(h, r, x);
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
			MATRIX dSy=get_dS(h,r,y);
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
			MATRIX dSx=get_dS(h,r,x);
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

	int trait_num;
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


//	env.positional_arg('n',"names",  names_name,      "please provide a filename.", "names of individuals in the sample.");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);

	//initParallel();
	setNbThreads(4);
	std::cerr << "using " << nbThreads( ) << " threads.\n";

	Flat_file <Relatedness> rel_file;
	Relatedness r;

	if (rel_name=="")
	{
		rel_file.open(READ);
	}
	else
	{
		rel_file.open(rel_name.c_str(), READ);
	}
	r=rel_file.read_header();
	rel_file.read(r);
	rel_file.close();

	Flat_file <Phenotype> pheno_file;
	Phenotype p;

	if (pheno_name=="")
	{
		pheno_file.open(READ);
	}
	else
	{
		pheno_file.open(pheno_name.c_str(), READ);
	}
	p=pheno_file.read_header();
	while(pheno_file.table_is_open() )
	{
		pheno_file.read(p);
	}
	pheno_file.close(); 

	std::cerr << __LINE__ << "Files read " << std::endl;
	


	Eigen::Matrix<int,Eigen::Dynamic,1>  drops=p.get_missing(trait_num);

	std::cerr << "got \n" << std::endl;

	r.drop_missing(drops);
	p.drop_missing(drops);

	int N(p.value.rows() );
	std::cerr << N << std::endl;

	std::cerr << p.value.block(0,trait_num, N, 1) << std::endl;

	//MATRIX V = makev(N);

	//std::cerr << p.value.block(0,0,N,1) << std::endl;
	//MATRIX P = p.value.rowwise()- p.value.colwise().mean();

	//MATRIX var=(P.adjoint()*P)/double( N - 1. );
	//std::cerr << "Var: " << var << std::endl;

	r.make_outer_prod();

	MATRIX S=MATRIX::Zero(N,N);

	Matrix<TYPE, 6, 1> h=MATRIX::Zero(6,1);

	h=lm(p,r, trait_num);

/*
	h(0, 0) =  2.06062;
	h(1, 0) =  1085.79;
	h(2, 0) =  2.57764;
	h(3, 0) = -613.632;
	h(4, 0) =-0.687815;
	h(5, 0) =  64.8181;
*/
	std::cerr << h << std::endl;

	MATRIX f_n, f_m, f_m2, h_m, h_n, h_o, df, dh, H_n, H_m;

	h_m=h;
/*


	h_n=bus(h_n);

	H_n=nabla2_center(h_n,p,r,0.005, 0).inverse();

	MATRIX v0=get_se(h_n, h_n-VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2, r);
	MATRIX v1=get_vars(h_n, r);
	MATRIX v2=get_se(h_n, h_n+VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2, r);

	std::cerr << "vars 0 : " << v0 << std::endl;
	std::cerr << "vars 1 : " << v1 << std::endl;
	std::cerr << "vars 2 : " << v2 << std::endl;*/

#define BROYDEN

	clock_t t0, t1;
	t0=clock();

#ifdef 	BROYDEN
	f_m=nabla(h_m, p, r, trait_num);
	std::cerr << f_m.squaredNorm() << "\n"; 
	{
		std::cerr << "calling nabla2.\n"; 
		H_m=nabla2_forward(h,p,r,f_m, 0.005, trait_num).inverse();
		std::cerr << "done.\n"; 
		h_n=h_m-H_m*f_m;
	}
#endif


	for (size_t x=0; x<50; x++)
	{

#ifdef 	BROYDEN
		
		f_n=nabla(h_n, p, r, trait_num);

		if (f_n.squaredNorm() > f_m.squaredNorm()*2 )
		{
			std::cerr << "whoops! re-doing step. " << f_m.squaredNorm() << " < " << f_n.squaredNorm() << "\n"; 
//			h_m=
//			f_m=
			H_m=nabla2_forward(h_m,p,r,f_m, 0.005, trait_num).inverse();
			std::cerr << "done.\n"; 
			h_n=h_m-H_m*f_m;
			f_n=nabla(h_n, p, r, trait_num);
		} 

		df=f_n-f_m;
		dh=h_n-h_m;
		
		//H_n=H_m+(dh-H_m*df)/df.squaredNorm()*t(df);
		H_n=H_m+(dh-H_m*df)*(t(dh)*H_m*df).inverse()*t(dh)*H_m;

		h_o=h_n-H_n*f_n;

		f_m=f_n;
		h_m=h_n;
		h_n=h_o;
		H_m=H_n;

		t1=clock();
		std::cerr << x << ", " << t1-t0 << ", " << f_n.squaredNorm() << "\nh: " << sub(h_n) << "\nse1: " << sub(h_n-VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2 ) << "\nse2:" << sub(h_n+VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2 ) << std::endl;
		MATRIX v=get_vars(h_n, r);
		std::cerr << "vars: " << v;

		if (f_n.squaredNorm() < 0.000005 ) break;
#else
		pair2 q=nablas(h, p, r, trait_num);	

		MATRIX iH=q.H.inverse();
		VectorXf D=iH.diagonal();
		h+=iH*q.j;

		t1=clock();
		std::cerr << x << ", " << t1-t0 << ", " << q.j.squaredNorm() << std::endl;
	
		if (q.j.squaredNorm() < 0.000005 ) break;
#endif
	}
	H_n=nabla2_center(h_n,p,r,0.005, trait_num).inverse();
	std::cerr << "done: " << t1-t0 << ", " << f_n.squaredNorm() << "\nh: " << sub(h_n) << "\nse1: " << sub(h_n-VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2 ) << "\nse2:" << sub(h_n+VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2 ) << std::endl;

	MATRIX v0=get_se(h_n, h_n-VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2, r);
	MATRIX v1=get_vars(h_n, r);
	MATRIX v2=get_se(h_n, h_n+VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2, r);

	std::cerr << "vars 0 : " << v0 << std::endl;
	std::cerr << "vars 1 : " << v1 << std::endl;
	std::cerr << "vars 2 : " << v2 << std::endl;
	return 0;
}
