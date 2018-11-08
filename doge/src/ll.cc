#include "ll.h"

#define t(a) a.transpose()

void
ll::print_options(void)
{
	std::cerr << __FILE__;
#ifdef __FAST_MATH__
	std::cerr << ", FAST MATH";
#endif
	std::cerr << std::endl;
}


MATRIX makev(const size_t &n)
{
	MATRIX I=MATRIX::Identity(n,n); 
	MATRIX J=MATRIX::Constant(n,n, -1./(TYPE(n)-1.) );
	return TYPE(n)/(TYPE(n)-1)*I+J;
}

TYPE
get_ll(const MATRIX &h, const Phenotype *p, const Relatedness *r, const int &trait_num, const int &num_blocks)
{
	TYPE ret=0;
	for (int x=0; x<num_blocks; x++)
	{
		size_t k=r[x].Mt1.rows()-1;
		MATRIX S=get_Sigma(h, r[x]);
		pair a=get_pseudo_exact(S);
		if(!a.pos_semi_def)
		{
			std::cerr << S << std::endl;
			exit(0);
			return std::numeric_limits<double>::quiet_NaN();
		}
		ret+=get_ll(h, p[x],r[x], a, trait_num);
	}
	return ret;
}

TYPE
get_ll(const MATRIX &h, const Phenotype &p, const Relatedness &r, const int &trait_num)
{
	size_t k=r.Mt1.rows()-1;
	MATRIX S=get_Sigma(h, r);
	pair a=get_pseudo_exact(S);
	if(!a.pos_semi_def)
		return std::numeric_limits<double>::quiet_NaN();
//	std::cerr << "Sigma:" << S << std::endl;
	return get_ll(h,p,r,a,trait_num);
}

TYPE
get_ll(const MATRIX &h, const Phenotype &p, const Relatedness &r, const pair &a, const int &trait_num)
{
	size_t k=r.Mt1.rows()-1;
	size_t n2=(k+1)/2.;
	size_t n=(k+1);
	double lnD=0;

	for (size_t x=0; x<k; x++)
	{
		lnD+=log(a.s(x,0));
	}
	lnD+=log(2*M_PI)*k;
	//std::cerr << p.value.block(0,trait_num,k+1,1) << std::endl;
	MATRIX cp=p.value.block(0,trait_num,k+1,1)-h(1,0)*r.Mt1-h(3,0)*r.Ht1;
	//std::cerr << a.P.block(0,0,10,10) << std::endl;

//	std::cerr << "part2: " << lnD/2. << std::endl;
//	std::cerr << "part1: " << (cp.transpose()*a.P*cp)(0,0)/2. << std::endl;

//	std::cerr << "A B" << std::endl;
//	std::cerr << "C D" << std::endl;
//	std::cerr << "A: " << (cp.block(0,0,n2,1).transpose()*a.P.block(0,0,n2,n2)*cp.block(0,0,n2,1) )(0,0)/2. << std::endl;
//	std::cerr << "B: " << (cp.block(n2,0,n2,1).transpose()*a.P.block(n2,0,n2,n2)*cp.block(0,0,n2,1) )(0,0)/2. << std::endl;
//	std::cerr << "C: " << (cp.block(0,0,n2,1).transpose()*a.P.block(0,n2,n2,n2)*cp.block(n2,0,n2,1) )(0,0)/2. << std::endl;
//	std::cerr << "D: " << (cp.block(n2,0,n2,1).transpose()*a.P.block(n2,n2,n2,n2)*cp.block(n2,0,n2,1) )(0,0)/2. << std::endl;

	std::cerr << "Condition: " << a.s(0,0)/a.s(k-1,0) << ", " << a.s(k,0) << std::endl;
//	std::cerr << "Sigma^(+):" <<  a.P << std::endl;
//	std::cerr << "t(x)*Sigma^(+) x:" <<  cp.transpose()*a.P*cp/2. << std::endl;

	return lnD/2.+ (cp.transpose()*a.P*cp)(0,0) /2.;
//	return (cp.transpose()*a.P*cp)(0,0) /2.;
}

MATRIX 
get_Sigma(const MATRIX &h, const Relatedness &r)
{
	double a,b,b2,c,d,d2,e,f;

	MATRIX subh;

	if ( h.rows()==5 )
	{
		subh=sub5(h, r.tr);
	} else if ( h.rows()==6 ) {
		subh=sub(h);
	} else {
		subh=sub7(h);
	}
		
	a=subh(0,0);
	b=subh(1,0);
	c=subh(2,0);
	d=subh(3,0);
	e=subh(4,0);
	f=subh(5,0);
	b2=b*b;
	d2=d*d;
	size_t n=r.Mt1.rows();
	MATRIX V=makev(n);
	std::cerr << "H" << a << ", " << b << ", " << c << ", " << d << ", " << e << ", " << f << std::endl;
	MATRIX ret=a*r.MtM+c*r.HtH+e*r.MtH+f*V;
	if (a+b2 > 0)
	{
		MATRIX MJM=2*r.Mt1*t(r.Mt1);
		ret-=(a*b2)/(a+b2)*MJM;
	}
	if (c+d2 > 0)
	{
		MATRIX HJH=2*r.Ht1*t(r.Ht1);
		ret-=(c*d2)/(c+d2)*HJH;
	}
	if ( (e+b*d) != 0 )
	{
		MATRIX MJH=r.Mt1*t(r.Ht1)+r.Ht1*t(r.Mt1);
		//ret-=2*(e*b*d)/(e+b*d)*MJH;
		ret-=(e*b*d)/(e+b*d)*MJH;
	}
	return ret;
}

MATRIX
get_var_norm(const Relatedness *r, const int &num_blocks)
{
	MATRIX ret=MATRIX::Zero(10,1);
	TYPE N=0;
	for (int x=0; x<num_blocks; x++)
	{
		size_t n=r[x].Mt1.rows();
		N+=n;
		MATRIX V=makev(n);

		ret(0,0)+=(r[x].MtM).trace();
		ret(1,0)+=(r[x].Mt1*t(r[x].Mt1) ).trace();
		ret(2,0)+=(t(r[x].Mt1)*V*r[x].Mt1)(0,0);

		ret(3,0)+=(r[x].HtH).trace();
		ret(4,0)+=(r[x].Ht1*t(r[x].Ht1) ).trace();
		ret(5,0)+=(t(r[x].Ht1)*V*r[x].Ht1)(0,0);

		ret(6,0)+=2.*(r[x].MtH).trace();
		ret(7,0)+=(r[x].Mt1*t(r[x].Ht1)+r[x].Ht1*t(r[x].Mt1) ).trace();
		ret(8,0)+=( t(r[x].Ht1)*V*r[x].Mt1+t(r[x].Mt1)*V*r[x].Ht1)(0,0);

		ret(9,0)+=(V).trace();
	}
	ret.array()/=N;
	return ret;
}

Matrix<TYPE, 6, 1>
get_vars(const MATRIX &h, const Relatedness *r, const int &num_blocks)
{
	double a,b,b2,c,d,d2,e,f;
	MATRIX subh=sub(h);
	a=subh(0,0);
	b=subh(1,0);
	c=subh(2,0);
	d=subh(3,0);
	e=subh(4,0);
	f=subh(5,0);
	b2=b*b;
	d2=d*d;
	Matrix<TYPE, 6, 1> ret=Matrix<TYPE, 6, 1>::Zero(6,1);
	TYPE N=0;
	for (int x=0; x<num_blocks; x++)
	{
		size_t n=r[x].Mt1.rows();
		N+=n;
		MATRIX V=makev(n);
		ret(0,0)+=(a*r[x].MtM).trace();
		ret(1,0)-=2*((a*b2)/(a+b2)*r[x].Mt1*t(r[x].Mt1) ).trace();
		ret(2,0)+=(c*r[x].HtH).trace();
		ret(3,0)-=2*((c*d2)/(c+d2)*r[x].Ht1*t(r[x].Ht1) ).trace();
		ret(4,0)+=2*( (e*r[x].MtH).trace()-(e*b*d)/(e+b*d)*(r[x].Mt1*t(r[x].Ht1)+r[x].Ht1*t(r[x].Mt1) ).trace() );
		ret(5,0)+=(f*V).trace();
	
	}
	ret.array()/=N;
	return ret;
}

Matrix<TYPE, 6, 1>
get_se(const MATRIX &h, const MATRIX &se, const Relatedness *r, const int &num_blocks)
{
	double a,b,b2,c,d,d2,e,f;
	double a_,b_,b2_,c_,d_,d2_,e_,f_;
	MATRIX subh=sub(h);

	a=subh(0,0);
	b=subh(1,0);
	c=subh(2,0);
	d=subh(3,0);
	e=subh(4,0);
	f=subh(5,0);

	b2=b*b;
	d2=d*d;
	
	subh=sub(se);

	a_=subh(0,0);
	b_=subh(1,0);
	c_=subh(2,0);
	d_=subh(3,0);
	e_=subh(4,0);
	f_=subh(5,0);

	b2_=b_*b_;
	d2_=d_*d_;
	
	Matrix<TYPE, 6, 1> ret=Matrix<TYPE, 6, 1>::Zero(6,1);
	Matrix<TYPE, 6, 1> den=Matrix<TYPE, 6, 1>::Zero(6,1);

	for (int x=0; x<num_blocks; x++)
	{
		size_t n=r[x].Mt1.rows();

		MATRIX V=makev(n);

		ret(0,0)=(a*r[x].MtM).trace();
		den(0,0)=(a_*r[x].MtM).trace();
		ret(1,0)=((b*r[x].Mt1).transpose()*V*b*r[x].Mt1);
		den(1,0)=((b_*r[x].Mt1).transpose()*V*b_*r[x].Mt1);
		ret(2,0)=(c*r[x].HtH).trace();
		den(2,0)=(c_*r[x].HtH).trace();
		ret(3,0)=((d*r[x].Ht1).transpose()*V*d*r[x].Ht1);
		den(3,0)=((d_*r[x].Ht1).transpose()*V*d_*r[x].Ht1);
		ret(4,0)=(e*r[x].MtH).trace();
		den(4,0)=(e_*r[x].MtH).trace();
		ret(5,0)=(f*V).trace();
		den(5,0)=(f_*V).trace();

		if (a+b2 > 0)
		{
			MATRIX MJM=r[x].Mt1*t(r[x].Mt1);
			ret(0,0)-=((a*b2)/(a+b2)*MJM).trace();
		}
		if (a_+b2_ > 0)
		{
			MATRIX MJM=r[x].Mt1*t(r[x].Mt1);
			den(0,0)-=((a_*b2)/(a_+b2_)*MJM).trace();
		}
		if (c+d2 > 0)
		{
			MATRIX HJH=r[x].Ht1*t(r[x].Ht1);
			ret(0,2)-=((c*d2)/(c+d2)*HJH).trace();
		}
		if (c_+d2_ > 0)
		{
			MATRIX HJH=r[x].Ht1*t(r[x].Ht1);
			den(0,2)-=((c_*d2)/(c_+d2_)*HJH).trace();
		}
		if ( (e+b*d) != 0 )
		{
			MATRIX MJH=r[x].Mt1*t(r[x].Ht1)+r[x].Ht1*t(r[x].Mt1);
			ret(0,4)-=((e*b*d)/(e+b*d)*MJH).trace();
		}
		if ( (e_+b_*d_) != 0 )
		{
			MATRIX MJH=r[x].Mt1*t(r[x].Ht1)+r[x].Ht1*t(r[x].Mt1);
			den(0,4)-=((e_*b_*d_)/(e_+b_*d_)*MJH).trace();
		}
	}

	//den/=ret.sum();

	return den;
}


Matrix<TYPE, 7, 1> sub7(const Matrix<TYPE, 7, 1> &a)
{
	Matrix<TYPE, 7, 1> b;
	b <<    a(0,0),
		a(1,0),
		a(2,0),
		a(3,0),
		a(4,0),
		a(5,0),
		a(6,0);
	return b;
}

Matrix<TYPE, 7, 1> bus7(const Matrix<TYPE, 7, 1> &a)
{
	Matrix<TYPE, 7, 1> b;
	b <<    a(0,0),
		a(1,0),
		a(2,0),
		a(3,0),
		a(4,0),
		a(5,0),
		a(6,0);
	return b;
}

Matrix<TYPE, 6, 1> sub(const Matrix<TYPE, 6, 1> &a)
{
	Matrix<TYPE, 6, 1> b;
	b <<    a(0,0),
		a(1,0),
		a(2,0),
		a(3,0),
		a(4,0),
		a(5,0);
	return b; 
}

Matrix<TYPE, 6, 1> bus(const Matrix<TYPE, 6, 1> &a)
{
	Matrix<TYPE, 6, 1> b;
	b << 	a(0,0),
		a(1,0),
		a(2,0),
		a(3,0),
		a(4,0),
		a(5,0);
	return b; 
}

Matrix<TYPE, 6, 1> sub5(const Matrix<TYPE, 5, 1> &a, const Matrix<TYPE, 7, 1> &b )
{
	std::cerr << b << std::endl;
	Matrix<TYPE, 6, 1> c;
	c <<   a(0,0),
		a(1,0),
		a(2,0),
		a(3,0),
		a(4,0),
		(b(6,0)-a(0,0)*b(0,0)-a(2,0)*b(2,0)-a(4,0)*b(4,0)
		 +2*(a(0,0)*a(1,0)*a(1,0) )/(a(0,0)+a(1,0)*a(1,0) )*b(1,0)
		 +2*(a(2,0)*a(3,0)*a(3,0) )/(a(2,0)+a(3,0)*a(3,0) )*b(3,0)
		 +2*(a(4,0)*a(1,0)*a(3,0) )/(a(4,0)+a(1,0)*a(3,0) )*b(5,0) )/b(6,0);
	std::cerr << c << std::endl;
	return c; 
}

Matrix<TYPE, 5, 1> bus5(const Matrix<TYPE, 6, 1> &a)
{
	Matrix<TYPE, 5, 1> b;
	b << 	a(0,0),
		a(1,0),
		a(2,0),
		a(3,0),
		a(4,0);
	return b; 
}

