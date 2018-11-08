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
#define Av(a,b) VectorXf(a.array() * b.array() )
#define o(a,b) MATRIX( a.array() * b.array() )
#define d(a,b) MATRIX( a.array() / b.array() )
#define center(a,b) MATRIX( a.colwise() - b )
#define D(a,b) MATRIX( a.array().colwise() / b.array() )

/*
//GLOBALS
MJM
HJH
MJH
float trMtM=r.MtM.trace();
float trMJM=MJM.trace();
float trHtH=r.HtH.trace();
float trHJH=HJH.trace();
float trMtH=r.MtH.trace();
float trMJH=MJH.trace();
float trE=
*/
double trMtM, trMJM, trHtH, trHJH, trMtH, trMJH, trE;

void set_trace(Relatedness *r, const int &num_blocks)
{
	trMtM=0;
	trMJM=0;
	trHtH=0;
	trHJH=0;
	trMtH=0;
	trMJH=0;
	trE=0;

	size_t n=r[0].Mt1.rows();

	for (int x=0; x<num_blocks; x++)
	{
		MATRIX MJM=r[x].Mt1*r[x].Mt1.transpose();
		MATRIX HJH=r[x].Ht1*r[x].Ht1.transpose();
		MATRIX MJH=(r[x].Mt1*r[x].Ht1.transpose()+r[x].Ht1*r[x].Mt1.transpose() )/2;

		trMtM+=double(r[x].MtM.trace() );
		trMJM+=double(MJM.trace() );
		trHtH+=double(r[x].HtH.trace() );
		trHJH+=double(HJH.trace() );
		trMtH+=double(r[x].MtH.trace() );
		trMJH+=double(r[x].MJH.trace() );
		trE+=double(n);
	}
}

std::string print_diag( Relatedness *r, const int &num_blocks)
{
	MATRIX MJM=r[0].Mt1*r[0].Mt1.transpose();
	MATRIX HJH=r[0].Ht1*r[0].Ht1.transpose();
	MATRIX MJH=(r[0].Mt1*r[0].Ht1.transpose()+r[0].Ht1*r[0].Mt1.transpose() )/2;

	return std::to_string(r[0].MtM.trace() )+", "+std::to_string(MJM.trace())+", "+std::to_string(r[0].HtH.trace())+", "+std::to_string(HJH.trace() )+", "+std::to_string(r[0].MtH.trace() )+", "+std::to_string(MJH.trace());
}


//calculates variance components, of which there are 6.
MATRIX qv(Matrix<float, 7, 1> A, MATRIX B)
{
	Matrix<float, 6, 1> ret;

	ret(1,0)=( A(1,0)*A(1,0) )*B(2,0);
	ret(3,0)=( A(3,0)*A(3,0) )*B(5,0);
	
	ret(0,0)=A(0,0)*B(0,0)-2*( A(0,0)*A(1,0)*A(1,0)/(A(0,0)+A(1,0)*A(1,0) ) )*B(1,0);
	ret(2,0)=A(2,0)*B(3,0)-2*( A(2,0)*A(3,0)*A(3,0)/(A(2,0)+A(3,0)*A(3,0) ) )*B(4,0);

	ret(4,0)=A(4,0)*B(6,0)-2*(A(4,0)*A(3,0)*A(1,0) )/(A(4,0)+A(3,0)*A(1,0) )*B(7,0)+(A(3,0)*A(1,0))*B(8,0);
	ret(5,0)=A(5,0)*B(9,0);
	return ret;
}

//Heuristic function to help keep step sizes small.
float
q_norm(Matrix<float, 7, 1> A)
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
	VectorXf h;
	VectorXf se;
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

	MATRIX ret=MATRIX::Zero(7,1);

	MATRIX x=MATRIX::Zero(N,2);
	MATRIX y=MATRIX::Zero(N,1);

	for (int z=0; z<num_blocks; z++)
	{
		x.block(sum[z],0,n[z],1)=r[z].Mt1;
		x.block(sum[z],1,n[z],1)=r[z].Ht1;
		y.block(sum[z],0,n[z],1)=p[z].value.block(0, trait_n, n[z], 1);
	}
	
	MATRIX y2=y.array()-y.mean();

	std::cerr << "mean y: " << y2.mean() << std::endl;
	std::cerr << "mean x: " << x.mean() << std::endl;

	MATRIX V=( t(x)*x ).inverse()*t(x)*y2;

	std::cerr << "V: " << V << std::endl << std::endl;

	ret(1,0)=V(0,0);
	ret(3,0)=V(1,0);

	VectorXf cy=(y2-ret(1,0)*x.block(0,0,N,1)-ret(3,0)*x.block(0,1,N,1));

	MATRIX V2=( t(x)*x ).inverse()*t(x)*cy;
	
	std::cerr << "V2: " << V2 << std::endl << std::endl;

	cy=cy.array()-cy.mean();

	std::cerr << "Means: " << cy.mean() << ", " << y2.mean() << std::endl;

	MATRIX yyt=MATRIX::Zero(nsqr,1);

	for (int z=0; z<num_blocks; z++)
	{
		MATRIX temp=(cy.block(sum[z],0,n[z],1)*cy.block(sum[z],0,n[z],1).transpose() );
		yyt.block(sumsqr[z],0,n[z]*n[z],1)=Map<VectorXf>( temp.data(), n[z]*n[z]);
	}

	//std::cerr << yyt.block(nsqr-10,0,10,1) << std::endl;

	x=MATRIX::Zero(nsqr,4);

	for (int z=0; z<num_blocks; z++)
	{
		Eigen::Map<VectorXf> v1(r[z].MtM.data(), r[z].MtM.size());
		x.block(sumsqr[z],0,n[z]*n[z],1) = v1;
	}
	for (int z=0; z<num_blocks; z++)
	{
		Map<VectorXf> v2(r[z].HtH.data(), r[z].MtM.size());
		x.block(sumsqr[z],1,n[z]*n[z],1) = v2;
	}
	for (int z=0; z<num_blocks; z++)
	{
		Map<VectorXf> v3(r[z].MtH.data(), r[z].MtM.size());
		v3=v3*2.;
		x.block(sumsqr[z],2,n[z]*n[z],1) = v3;
	}
	for (int z=0; z<num_blocks; z++)
	{
		MATRIX v=makev(n[z]);
		Map<VectorXf> v4(v.data(), r[z].MtM.size());
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
	Matrix<float, 7, 1> j;
	Matrix<float, 7, 7> H;
};

struct trip
{
	MATRIX dS[7];	
	MATRIX dM[7];
	MATRIX ddS[7][7];	
};

/*
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
}*/

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

#define LINE_E

//To be called when make_outer() hasn't.
MATRIX
get_dS_small(const MATRIX &h, const Relatedness &r, const size_t &x)
{
	double a,b,c,d,e,f,l;
	size_t n=r.MtM.rows();
	a=h(0,0);
	b=h(1,0);
	c=h(2,0);
	d=h(3,0);
	e=h(4,0);
	f=h(5,0);

	MATRIX empty;

	MATRIX v=makev(n);


switch(x)
{
	case(0):
	{
		MJM=r.Mt1*r.Mt1.transpose();
		return MJM*(2*a*pow(b, 2)/pow(a + pow(b, 2), 2) - 2*pow(b, 2)/(a + pow(b, 2))) - l*(-2*a*pow(b, 2)*trMJM/pow(a + pow(b, 2), 2) + 2*pow(b, 2)*trMJM/(a + pow(b, 2)) - trMtM) + r.MtM ; // 0
	}
	case(1):
	{
		MJM=r.Mt1*r.Mt1.transpose();
		MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;
		return MJH*(2*b*pow(d, 2)*e/pow(b*d + e, 2) - 2*d*e/(b*d + e)) + MJM*(4*a*pow(b, 3)/pow(a + pow(b, 2), 2) - 4*a*b/(a + pow(b, 2))) - l*(-4*a*pow(b, 3)*trMJM/pow(a + pow(b, 2), 2) + 4*a*b*trMJM/(a + pow(b, 2)) - 2*b*pow(d, 2)*e*trMJH/pow(b*d + e, 2) + 2*d*e*trMJH/(b*d + e)) ; // 1
	}
	case(2):
	{
		HJH=r.Ht1*r.Ht1.transpose();
		return HJH*(2*c*pow(d, 2)/pow(c + pow(d, 2), 2) - 2*pow(d, 2)/(c + pow(d, 2))) - l*(-2*c*pow(d, 2)*trHJH/pow(c + pow(d, 2), 2) + 2*pow(d, 2)*trHJH/(c + pow(d, 2)) - trHtH) + r.HtH ; // 2
	}
	case(3):
	{
		HJH=r.Ht1*r.Ht1.transpose();
		MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;
		return HJH*(4*c*pow(d, 3)/pow(c + pow(d, 2), 2) - 4*c*d/(c + pow(d, 2))) + MJH*(2*pow(b, 2)*d*e/pow(b*d + e, 2) - 2*b*e/(b*d + e)) - l*(-2*pow(b, 2)*d*e*trMJH/pow(b*d + e, 2) + 2*b*e*trMJH/(b*d + e) - 4*c*pow(d, 3)*trHJH/pow(c + pow(d, 2), 2) + 4*c*d*trHJH/(c + pow(d, 2))) ; // 3
	}
	case(4):
	{
		MJH=(r.Mt1*r.Ht1.transpose()+r.Ht1*r.Mt1.transpose() )/2;
		return MJH*(2*b*d*e/pow(b*d + e, 2) - 2*b*d/(b*d + e)) - l*(-2*b*d*e*trMJH/pow(b*d + e, 2) + 2*b*d*trMJH/(b*d + e) - trMtH) + r.MtH ; // 4
	}
	case(5):
	{
		return l*trE + v ; // 5
	}
	case(6):
	{
		return -2*a*pow(b, 2)*trMJM/(a + pow(b, 2)) + a*trMtM - 2*b*d*e*trMJH/(b*d + e) - 2*c*pow(d, 2)*trHJH/(c + pow(d, 2)) + c*trHtH + e*trMtH + f*trE - 1.0 ; // 6
	}
	default:
		return empty;
}
	return empty;
}

MATRIX
get_ddS(const MATRIX &h, const Relatedness &r, const size_t &x, const size_t &y)
{
	double a,b,b2,b3,b4,c,d,d2,d3,d4,e,e2,e3,f;
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
				return 4*MJH*(-b*pow(d, 3)*e/pow(b*d + e, 3) + pow(d, 2)*e/pow(b*d + e, 2)) + 4*MJM*(-4*a*pow(b, 4)/pow(a + pow(b, 2), 3) + 5*a*pow(b, 2)/pow(a + pow(b, 2), 2) - a/(a + pow(b, 2))) ; // 1 1
			case(3):
				return 2*MJH*e*(-2*pow(b, 2)*pow(d, 2)/pow(b*d + e, 2) + 3*b*d/(b*d + e) - 1)/(b*d + e) ; // 1 3
			case(4):
				return 2*MJH*d*(-2*b*d*e/pow(b*d + e, 2) + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e) ; // 1 4
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
				return 4*HJH*(-4*c*pow(d, 4)/pow(c + pow(d, 2), 3) + 5*c*pow(d, 2)/pow(c + pow(d, 2), 2) - c/(c + pow(d, 2))) + 4*MJH*(-pow(b, 3)*d*e/pow(b*d + e, 3) + pow(b, 2)*e/pow(b*d + e, 2)) ; // 3 3
			case(4):
				return 2*MJH*b*(-2*b*d*e/pow(b*d + e, 2) + b*d/(b*d + e) + e/(b*d + e) - 1)/(b*d + e) ; // 3 4
			default:
				return empty;
		}
	case(4):
		switch(y)
		{
			case(4):
				return 4*MJH*b*d*(-e/(b*d + e) + 1)/pow(b*d + e, 2) ; // 4 4
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

Matrix<float, 7, 1> 
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

	Matrix<float, 7, 1> j=MATRIX::Zero(7,1);
	MATRIX cp=p.value.block(0,trait_n,n,1)-h(1,0)*r.Mt1-h(3,0)*r.Ht1;
	MATRIX ImSPcp=(MATRIX::Identity(n,n)-S*a.P)*cp;		//Orthogonal projector onto the kernel of S. . . 
	MATRIX Pc=a.P*cp;
	MATRIX tPc=t(Pc)*a.P;
	t3=clock();

	for (size_t x=0; x<7; x++)
	{
		MATRIX PdS=a.P*get_dS_small(h,r,x);
		float dP2=(t(Pc)*PdS*ImSPcp)(0,0);   //4

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

	Matrix<float, 7,1> j;

	MATRIX PP=P*t(P);	// PP(x,y) != P(x,y)*P(x,y)
	MATRIX SP=S*t(P);	// SP(x,y) != S(x,y)*P(x,y)


	MATRIX ImSP;		// S*t(p) is square. . . 
	ImSP=MATRIX::Identity(n,n)-SP; 		//Depends on where the diagonal is...

//	MATRIX ImPS=t(ImSP);					//Be careful with transposes too!

	VectorXf cp=p.value-h(1,0)*r.Mt1-h(3,0)*r.Ht1;

	for (size_t x=0; x<7; x++)
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

Matrix<float, 7,1> 
sub_nabla(const MATRIX &h, const VectorXf &lp, const VectorXf &rp,  const Relatedness &r, const bool &diagonal, const MATRIX &P, const MATRIX &S)
{
	size_t n=r.MtM.rows();

	Matrix<float, 7,1> j;

	//let assume P is a n x m slice of P

	MATRIX PP=P*t(P);	// PP(x,y) != P(x,y)*P(x,y)
	MATRIX SP=S*t(P);	// SP(x,y) != S(x,y)*P(x,y)


	MATRIX ImSP;		// S*t(p) is square. . . 
	if (diagonal) ImSP=MATRIX::Identity(n,n)-SP; 		//Depends on where the diagonal is...
	else ImSP=-SP; 						//Depends on where the diagonal is...

//	MATRIX ImPS=t(ImSP);					//Be careful with transposes too!

	VectorXf lcp=lp-h(1,0)*r.lMt1-h(3,0)*r.lHt1;
	VectorXf rcp=rp-h(1,0)*r.rMt1-h(3,0)*r.rHt1;

	for (size_t x=0; x<7; x++)
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

Matrix<float, 7, 1>
get_sub(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &step, const int &trait_n)
{
	size_t n=r.Mt1.rows();
	MATRIX S=get_Sigma(h, r);
	pair a=get_pseudo(S);
	
	Matrix<float, 7, 1> ret=MATRIX::Zero(7,1);

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

Matrix<float, 7, 7>
nabla2_forward(const Matrix<float, 7, 1> &h, const Phenotype &p, const Relatedness &r, const Matrix<float, 7, 1>k, const float &delta, const int &trait_num)
{
	MATRIX H=MATRIX::Zero(7, 7);

	for (size_t x=0; x<7; x++)
	{
		MATRIX dh=MATRIX::Zero(7,1);
		dh(x,0)=std::max(delta*fabs(h(x,0)), MAX_TYPE(0.000001) );
		Matrix<float, 7, 1> k2=nabla(h+dh, p, r, trait_num);
		H.block(0,x,7,1)=(k2-k)/dh.norm();
	}

	return H;
}

Matrix<float, 7, 7>
nabla2_center(const Matrix<float, 7, 1> &h, const Phenotype &p, const Relatedness &r, const float &delta, const int &trait_num)
{
	MATRIX H=MATRIX::Zero(7, 7);

	for (size_t x=0; x<7; x++)
	{
		MATRIX dh=MATRIX::Zero(7,1);
		//dh(x,0)=delta*fabs(h(x,0));
		dh(x,0)=std::max(delta*fabs(h(x,0)), MAX_TYPE(0.000001) );
		Matrix<float, 7, 1> k1=nabla(h-dh, p, r, trait_num);
		Matrix<float, 7, 1> k2=nabla(h+dh, p, r, trait_num);
		H.block(0,x,7,1)=(k2-k1)/(2.*dh.norm() );
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

	MATRIX dP[7];
	MATRIX dPP[7], DdP[7], PPdS[7], dSP[7], dSImSP[7];


	for (size_t x=0; x<7; x++)
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

	for (size_t x=0; x<7; x++)
	{
		MATRIX dMx=get_dM(h,r,x);
		for (size_t y=x; y<7; y++)
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
	float loci=0;
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
	//std::cerr << __LINE__ << p[0].value << std::endl;


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

	/*
	for (int x=0; x<num_blocks; x++)
	{
		Eigen::Matrix<int,Eigen::Dynamic,1>  drops=p[x].get_missing(trait_num);
		r[x].drop_missing(drops);
		p[x].drop_missing(drops);

		std::vector<int> pre_drop;
		for (int z=0; z<r[x].Ht1.rows(); z++)
		{
			if (r[x].Ht1(z,0) < -1.5 || r[x].Ht1(z,0) > 2 || r[x].Mt1(z,0) < -2 || r[x].Mt1(z,0) > 2 ) pre_drop.push_back(z);
		}

		drops=Eigen::Matrix<int, Eigen::Dynamic, 1>::Zero(pre_drop.size(), 1);

		for (int z=0; z<pre_drop.size(); z++)
		{
			drops(z,0)=pre_drop[z];
		}

		r[x].drop_missing(drops);
		p[x].drop_missing(drops);
	}*/

	if (hessian)
	{
		for(int z=0; z<num_blocks; z++)
			r[z].make_outer_prod();
	}

	MATRIX S;

	Matrix<float, 7, 1> h=MATRIX::Zero(7,1);
	Matrix<float, 7, 1> norm=MATRIX::Zero(7,1);
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

		double var=float( (y.adjoint() * y)(0,0) ) / float(N - 1.);

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
		var_norm=get_var_norm(r, num_blocks);
	}

	h=lm(p, r, trait_num, num_blocks);

	std::cerr << "lm: " << t( o(sub7(h), norm) ) << std::endl;

	MATRIX K=MATRIX::Zero(7,1);

/*	K(0,0) = 1.;
	K(1,0) = -0.25;
	K(2,0) = 1.;
	K(3,0) = -0.25;
	K(4,0) = -0.25;
	K(5,0) = 0.5;

	h=d(K, norm);*/

//	std::cout << get_Sigma( sub(h), r[0]);
//	return 0;

	std::cerr << "lm: " << t( o(sub7(h), norm) ) << std::endl;

	std::cout << "Un-normed : " << t(sub7(h) ) << std::endl;
	std::cout << "LL: " << get_ll(h, p[0], r[0], trait_num) << std::endl;
	std::cerr << "vars: " << t(qv(h, var_norm) ) << std::endl;

/*	K(0,0) = 1.014;
	K(1,0) =-0.237;
	K(2,0) = 0.969;
	K(3,0) =-0.238;
	K(4,0) =-0.096;
	K(5,0) = 0.5;
	h=d(K, norm);*/

	std::cout << "LL: " << get_ll(h, p[0], r[0], trait_num) << std::endl;
	std::cerr << "vars: " << t(qv(h, var_norm) ) << std::endl;

	{
	std::cerr << "vars: " << t(qv(h, var_norm) ) << std::endl;
	std::cerr << "for h:" << t( o(h, norm) ) << std::endl;
	}

	MATRIX f_n, f_m, f_m2, h_m, h_n, h_o, df, dh, H_n, H_m;

	h_m=h;

	clock_t t0, t1;
	t0=clock();

	if (!hessian)
	{

	f_m=MATRIX::Zero(7,1);

	for (int x=0; x<num_blocks; x++)
	{
		f_m+=nabla(h_m, p[x], r[x], trait_num);
	}

	{

		H_m=MATRIX::Zero(7,7);

		for (int x=0; x<num_blocks; x++)
			H_m+=nabla2_forward(h,p[x],r[x],f_m, 0.05, trait_num);
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
		
		f_n=MATRIX::Zero(7,1);
		for (int z=0; z<num_blocks; z++)
			f_n+=nabla(h_n, p[z], r[z], trait_num);

		df=f_n-f_m;
		dh=h_n-h_m;
		
		H_n=H_m+(dh-H_m*df)*(t(dh)*H_m*df).inverse()*t(dh)*H_m;

		VectorXf step=(H_n*f_n);

		std::cerr << "step size (var):" << qv(step,var_norm).norm() << " : " << t(qv(step,var_norm) ) << std::endl;
		if ( qv(step,var_norm).norm() < pow(10,-5) ) break;
		step/=(1.+10.*pow(q_norm(qv(h_n-step,var_norm) )-1.,2) );
		std::cerr << pow(q_norm(qv(h_n-step,var_norm) )-1.,2) << std::endl;
		std::cerr << q_norm(qv(h_n-step,var_norm) )-1. << std::endl;
		std::cerr << qv(h_n-step,var_norm) << std::endl;
		std::cerr << "step size:" << qv(step,var_norm).norm() << " : " << t(qv(step,var_norm) ) << std::endl;

		h_o=h_n-step;

		f_m=f_n;
		h_m=h_n;
		h_n=h_o;
		H_m=H_n;

		t1=clock();
		std::cerr << x << ", " << t1-t0 << ", " << f_n.squaredNorm() << "\nh: " << t( o(sub7(h_o), norm) ) << std::endl << std::endl;
		std::cout << "LL:" << get_ll(h_n, p[0], r[0], trait_num) << std::endl;

		std::cerr << "vars: " << t(qv(h_n, var_norm) ) << "="  << qv(h_n, var_norm).sum() << std::endl;


		} else {

		std::cerr << __LINE__ << std::endl;

		MATRIX iH=MATRIX::Zero(7,7);
		MATRIX j=MATRIX::Zero(7,1);

		for (int z=0; z<num_blocks; z++)
		{
			std::cerr << __LINE__ << std::endl;
			pair2 q=nablas(h, p[z], r[z], trait_num);	
			std::cerr << q.H.rows() << std::endl;
			std::cerr << q.H.cols() << std::endl;
			iH+=q.H;
			j+=q.j;
		}
	
		H_n=iH;
		iH=iH.inverse();

		//std::cerr << __LINE__ << H_n << std::endl << std::endl;
		//std::cerr << __LINE__ << iH << std::endl << std::endl;
		//std::cerr << __LINE__ << nabla2_forward(h, p[0], r[0], j, 0.001, trait_num) << std::endl;

		VectorXf D=iH.diagonal();
		VectorXf step=(iH*j);
		std::cerr << "step size:" << step.norm() << std::endl;
		if ( step.norm() < pow(10,-5) ) break;
		step/=log(pow(step.norm(),2)+exp(1) );
		h+=(step);

		t1=clock();
		std::cerr << x << ", " << t1-t0 << ", " << j.squaredNorm() << std::endl << "j:" << t(j) << std::endl << "h:" << t( o(sub7(h), norm) ) << std::endl << std::endl;	
		std::cout << "LL:" << get_ll(h, p[0], r[0], trait_num) << std::endl;

		}
	}

	if(!hessian)
	{
		H_n=MATRIX::Zero(7,7);
		for (int x=0; x<num_blocks; x++)
			H_n=nabla2_center(h_n,p[x],r[x], 0.005, trait_num);
		H_n=H_n.inverse();
	}
	else
		h_n=h;	


	MATRIX v0=qv(h_n-VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2, var_norm);
	MATRIX v1=qv(h_n, var_norm);
	MATRIX v2=qv(h_n+VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2,  var_norm);

	std::cout << var_norm << std::endl;
	std::cout << "text, :, , a, A, d, D, r, e" << std::endl;
	std::cout << "vars 0 : " << t(v0) << std::endl;
	std::cout << "vars 1 : " << t(v1) << std::endl;
	std::cout << "vars 2 : " << t(v2) << std::endl;
	std::cout << "h0 : " << t(o(sub7(h_n-VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2 ), norm) ) << std::endl;
	std::cout << "h1 : " << t(o(sub7(h_n), norm) ) << std::endl;
 	std::cout << "h2 : " << t(o(sub7(h_n+VectorXf(ArrayXf(H_n.diagonal()).abs().sqrt())*2 ), norm) ) << std::endl;
	std::cout << "Un-normed : " << t(sub7(h_n) ) << std::endl;

	delete [] r;
	delete [] p;

	return 0;
}
