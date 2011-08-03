/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/nlbit.h"
#include "mg/PPRep.h"
#include "mg/LBRepEndC.h"
#include "mg/LBRep.h"
#include "mg/LinearEquation.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGLBRep.cpp
//
// Implement MGLBRep class.
//******Smoothing function fo Shoenberg and Reinch variation.
//This smoothing function keeps boundary condition, start and end point
//and 1st derivative at start and end point.

//Q-transpose is the tridiagonal matrix of order n*n with general row
//{1/dxim1 -(1/dxim1+1/dxi) 1/dxi} for i=0,...,n-1,
//D is the diagonal matrix of dy[i] for i=0,...,n-1,
//constructQTD will construct dx, QtDDQ, and QtY of below.
void constructQTD(
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	std::vector<double>& dx,//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2, will be output
	MGBPointSeq& QtDDQ,	//three bands of Q-transpose*D*D*Q at and above the diagonal
					//will be output.
					//QtDDQ(i,0)=(i,i),
					//QtDDQ(i,1)=(i,i+1) and (i+1,i)
					//QtDDQ(i,2)=(i,i+2) and (i+2,i) for i=0,...,n-1
	MGBPointSeq& QtY)	//Q-transpose*y will be output
		//QtY(.,0) and (.,n-1) will be modified to take the end conditions into account.
{
	size_t n=tau.length();assert(n>=3);
	size_t nm1=n-1, nm2=n-2;
	size_t sd=y.sdim();

	size_t i,j;
	dx.resize(nm1);
	dx[0]=tau[1]-tau[0];

	MGBPointSeq QtD(n,3);//QtD(i,0) is (i,i-1), (i,1) is (i,i) , and (i,2) is (i,i+1).
	QtD(0,0)=0.; QtD(0,1)=-dy[0]/dx[0]; QtD(0,2)=dy[1]/dx[0];
	for(i=1; i<nm1; i++){
		double dxim1=dx[i-1];
		double dxi=tau[i+1]-tau[i];
		dx[i]=dxi;
		QtD(i,0)=dy[i-1]/dxim1;
		QtD(i,1)=-dy[i]/dxim1-dy[i]/dxi;
		QtD(i,2)=dy[i+1]/dxi;
	}
	QtD(nm1,0)=dy[nm2]/dx[nm2]; QtD(nm1,1)=-dy[nm1]/dx[nm2]; QtD(nm1,2)=0.;
	//cout<<endl<<"*****dx="<<dx; cout<<"****QtD="<<QtD<<endl;

	QtDDQ.resize(n,3);//QtDDQ(i,0)=(i,i),
					//QtDDQ(i,1)=(i,i+1) and (i+1,i)
					//QtDDQ(i,2)=(i,i+2) and (i+2,i)
	for(i=0; i<n; i++){
		double di0=QtD(i,0), di1=QtD(i,1), di2=QtD(i,2);
		QtDDQ(i,0)=di0*di0+di1*di1+di2*di2;
	}
	for(i=1; i<n; i++){
		size_t im1=i-1;
		QtDDQ(im1,1)=QtD(im1,1)*QtD(i,0)+QtD(im1,2)*QtD(i,1);
	}
	QtDDQ(nm1,1)=0.;
	for(i=2; i<n; i++){
		QtDDQ(i-2,2)=QtD(i-2,2)*QtD(i,0);
	}
	QtDDQ(nm2,2)=0.;
	QtDDQ(nm1,2)=0.;

	QtY.resize(n,sd);
	for(j=0; j<sd; j++){
		double prev=(y(1,j)-y(0,j))/dx[0];
		QtY(0,j)=prev;
		for(i=1; i<nm1; i++){
			double diff=(y(i+1,j)-y(i,j))/dx[i];
			QtY(i,j)=diff-prev;
			prev=diff;
		}
		QtY(nm1,j)=-prev;
	}
}

//Q-transpose is the tridiagonal matrix of order (n-2)*n with general row
//{1/dxim1 -(1/dxim1+1/dxi) 1/dxi} for i=1,...,n-2,
//D is the diagonal matrix of dy[i] for i=0,...,n-1,
//constructQTD will construct dx, QtDDQ, and QtY of below for free end condition SR.
void constructQTDFreeE(
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2, will be output
	MGBPointSeq& QtDDQ,	//three bands of Q-transpose*D*D*Q at and above the diagonal
					//will be output.
					//QtDDQ(i,0)=(i,i),
					//QtDDQ(i,1)=(i,i+1) and (i+1,i)
					//QtDDQ(i,2)=(i,i+2) and (i+2,i) for i=0,...,n-1
	MGBPointSeq& QtY)	//Q-transpose*y will be output
		//QtY(.,0) and (.,n-1) will be modified to take the end conditions into account.
{
	size_t n=tau.length();assert(n>=3);
	size_t nm1=n-1, nm2=n-2, nm3=n-3;
	size_t sd=y.sdim();

	size_t i,j;

	MGBPointSeq QtD(nm2,3);//QtD(i,0) is (i,i-1), (i,1) is (i,i) , and (i,2) is (i,i+1).
	dx.resize(nm1);
	dx[0]=tau[1]-tau[0];
	for(i=1; i<=nm2; i++){
		size_t im1=i-1, ip1=i+1;
		double dxim1=dx[im1];
		double dxi=tau[ip1]-tau[i];
		dx[i]=dxi;
		QtD(im1,0)=dy[im1]/dxim1;
		QtD(im1,1)=-dy[i]/dxim1-dy[i]/dxi;
		QtD(im1,2)=dy[ip1]/dxi;
	}

	QtDDQ.resize(nm2,3);//QtDDQ(i,0)=(i,i),
					//QtDDQ(i,1)=(i,i+1) and (i+1,i)
					//QtDDQ(i,2)=(i,i+2) and (i+2,i)
	for(i=0; i<nm2; i++){
		double di0=QtD(i,0), di1=QtD(i,1), di2=QtD(i,2);
		QtDDQ(i,0)=di0*di0+di1*di1+di2*di2;
	}
	for(i=1; i<nm2; i++){
		size_t im1=i-1;
		QtDDQ(im1,1)=QtD(im1,1)*QtD(i,0)+QtD(im1,2)*QtD(i,1);
	}
	QtDDQ(nm3,1)=0.;
	for(i=2; i<nm2; i++){
		QtDDQ(i-2,2)=QtD(i-2,2)*QtD(i,0);
	}
	if(n>3)  QtDDQ(n-4,2)=0.;
	QtDDQ(nm3,2)=0.;
	//cout<<"****QtDDQ="<<QtDDQ<<endl;

	QtY.resize(nm2,sd);
	for(j=0; j<sd; j++){
		double prev=(y(1,j)-y(0,j))/dx[0];
		for(i=1; i<=nm2; i++){
			double diff=(y(i+1,j)-y(i,j))/dx[i];
			QtY(i-1,j)=diff-prev;
			prev=diff;
		}
	}
	//cout<<"****QtY="<<QtY<<endl;
}

//Q-transpose is the tridiagonal matrix of order n*n with general row
//{1/dxim1 -(1/dxim1+1/dxi) 1/dxi} for i=0,...,n-1,
//D is the diagonal matrix of dy[i] for i=0,...,n-1,
//constructQTD will construct dx, QtDDQ, and QtY of below.
void constructQTD(
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2, will be output
	MGBPointSeq& QtDDQ,	//three bands of Q-transpose*D*D*Q at and above the diagonal
					//will be output.
					//QtDDQ(i,0)=(i,i),
					//QtDDQ(i,1)=(i,i+1) and (i+1,i)
					//QtDDQ(i,2)=(i,i+2) and (i+2,i) for i=0,...,n-1
	MGBPointSeq& QtY	//Q-transpose*y will be output
		//QtY(.,0) and (.,n-1) will be modified to take the end conditions into account.
){
	constructQTD(tau,y,dy,dx,QtDDQ,QtY);

	//Set the end coditions in QtY.
	const MGVector& deriS=begin.first();
	const MGVector& deriE=end.first();
	size_t sd=y.sdim(), nm1=tau.length()-1;
	for(size_t j=0; j<sd; j++){
		QtY(0,j)-=deriS[j]; QtY(nm1,j)+=deriE[j];
	}
}

//Given p, dx, dy, and U,  constructQU() will construct QU.
//Consequently sfp, and max_dev_at_p will be obtained.
void constructQU(
	double p,	//value between 0. and 1.
	const std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	const MGBPointSeq& U,//input U
	MGBPointSeq& QU,//Q*U will be output.
	double& sfp,	//sum of all the coordinate deviation squares at point y(.,.).
	double& max_dev_at_p// square of maximum deviation will be output.
){
//construct Q*U in QU and compute sfp.
	size_t sd=U.sdim(), n=U.length();
	size_t nm1=n-1;

	QU.resize(n,sd);
	max_dev_at_p=sfp=0.;
	MGVector prev(sd); prev.clear();
	MGVector qui;
	double dyi, qui2, max_devi;
	double six1mp=6.*(1.-p);
	size_t i;
	for(i=0; i<nm1; i++){
		MGVector diffi=U(i+1)-U(i); diffi/=dx[i];
		qui=diffi-prev;
		QU.store_at(i,qui);
		dyi=dy[i];
		qui*=dyi*six1mp;
		qui2=qui%qui;
		sfp+=qui2;
		if(dyi>1.) max_devi=qui2;
		else max_devi=qui2*dyi*dyi;
		if(max_devi>max_dev_at_p)
			max_dev_at_p=max_devi;
		prev=diffi;
	}
	qui=(-prev);
	QU.store_at(nm1,qui);
	dyi=dy[i];
	qui*=dyi*six1mp;
	qui2=qui%qui;
	sfp+=qui2;
	max_devi=qui2*dyi*dyi;
	if(max_devi>max_dev_at_p)
		max_dev_at_p=max_devi;
}

//Given p, dx, QtDDQ, QtY, and dy, solve() will solve the linear equation
// (6*(1-p)*QtDDQ+p*R)*U=QtY, where R is the tridiagonal matrix whose
// general row is dx[i-1] 2(dx[i-1]+dx[i]) dx[i] for i=0, ..., n-1.
//Consequently, U, QU, sfp, and max_dev_at_p will be obtained.
void solveFreeE(
	double p,	//value between 0. and 1. If this is not the case,
				//p will be set into the range.
	const std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2
	const MGBPointSeq& QtDDQ,
		//three bands of Q-transpose*D*D*Q at and above the diagonal are input.
		//QtDDQ(i,0)=(i,i),
		//QtDDQ(i,1)=(i,i+1) and (i+1,i)1
		//QtDDQ(i,2)=(i,i+2) and (i+2,i) for i=0,...,n-1
	const MGBPointSeq& QtY,	//Q-transpose*y are input.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	MGBPointSeq& U,	//solved U will be output.
	MGBPointSeq& QU,//Q*U will be output.
	double& sfp,	//sum of all the coordinate deviation squares at point y(.,.).
	double& max_dev_at_p// square of maximum deviation will be output.
){
	if(p<0.) p=0.;if(p>1.) p=1.;
	int i, n=QtDDQ.length()+2, sd=QtY.sdim();
	assert(n>2);
	int nm1=n-1, nm2=n-2;

	//construct W=6*(1-p)*Qt*D*D*Q + pR and solve the linear equation W*U=QtY.
	double six1mp=6.*(1.-p);
	double twop=2.*p;
	MGBPointSeq W(nm2,3);//W will be 6*(1-p)*Qt*D*D*Q + pR;
	for(i=0; i<nm2; i++){
		W(i,0)=six1mp*QtDDQ(i,0)+twop*(dx[i]+dx[i+1]);
		W(i,1)=six1mp*QtDDQ(i,1)+p*dx[i+1];
		W(i,2)=six1mp*QtDDQ(i,2);
	}
	solveSymetricTridiagonal(W,QtY,U);//cout<<"out of solveSymetricTridiagonal6="<<endl<<U<<endl;
	U.reshape(n,1); U.set_length(n);//U.length() was n-2.
	for(int j=0; j<sd; j++){
		U(0,j)=0.; U(nm1,j)=0.;//set free end condition.
	}

//construct Q*U in QU and compute sfp.
	constructQU(p,dx,dy,U,QU,sfp,max_dev_at_p);
}

//Given p, dx, QtDDQ, QtY, and dy, solve() will solve the linear equation
// (6*(1-p)*QtDDQ+p*R)*U=QtY, where R is the tridiagonal matrix whose
// general row is dx[i-1] 2(dx[i-1]+dx[i]) dx[i] for i=0, ..., n-1.
//Consequently, U, QU, sfp, and max_dev_at_p will be obtained.
void solve(
	double p,	//value between 0. and 1. If this is not the case,
				//p will be set into the range.
	const std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2
	const MGBPointSeq& QtDDQ,
		//three bands of Q-transpose*D*D*Q at and above the diagonal are input.
		//QtDDQ(i,0)=(i,i),
		//QtDDQ(i,1)=(i,i+1) and (i+1,i)1
		//QtDDQ(i,2)=(i,i+2) and (i+2,i) for i=0,...,n-1
	const MGBPointSeq& QtY,	//Q-transpose*y are input.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	MGBPointSeq& U,	//solved U will be output.
	MGBPointSeq& QU,//Q*U will be output.
	double& sfp,	//sum of all the coordinate deviation squares at point y(.,.).
	double& max_dev_at_p// square of maximum deviation will be output.
){
	if(p<0.) p=0.;if(p>1.) p=1.;
	int i, n=QtDDQ.length(), sd=QtY.sdim();
	assert(n>2);
	int nm1=n-1, nm2=n-2;

	//construct W=6*(1-p)*Qt*D*D*Q + pR and solve the linear equation W*U=QtY.
	double six1mp=6.*(1.-p);
	double twop=2.*p;
	MGBPointSeq W(n,3);//W will be 6*(1-p)*Qt*D*D*Q + pR;
	W(0,0)=six1mp*QtDDQ(0,0)+twop*dx[0];
	W(0,1)=six1mp*QtDDQ(0,1)+p*dx[0];
	W(0,2)=six1mp*QtDDQ(0,2);
	for(i=1; i<nm1; i++){
		W(i,0)=six1mp*QtDDQ(i,0)+twop*(dx[i-1]+dx[i]);
		W(i,1)=six1mp*QtDDQ(i,1)+p*dx[i];
		W(i,2)=six1mp*QtDDQ(i,2);
	}
	W(nm1,0)=six1mp*QtDDQ(nm1,0)+twop*(dx[nm2]);
	W(nm1,1)=0.;
	W(nm1,2)=0.;
	solveSymetricTridiagonal(W,QtY,U);//cout<<"out of solveSymetricTridiagonal6="<<endl<<U<<endl;

//construct Q*U in QU and compute sfp.
	constructQU(p,dx,dy,U,QU,sfp,max_dev_at_p);
}

//A dedicated class for compute_smoothed_p, provides a functional object for mgNlbit,
//to solve the equation sfp_of_solve(p) -max_sfp=0.
class MGsfpDiff{
	//all of the following member varialbles are the same as compute_smoothed_p.
	//See compute_smoothed_p.
	bool m_sum;
	double m_dev;
	const std::vector<double>& m_dx;
	const MGBPointSeq& m_QtDDQ;
	const MGBPointSeq& m_QtY;
	const double* m_dy;//Weights  at tau[i] for i=0,..., tau.length()-1.
	MGBPointSeq& m_U;
	MGBPointSeq& m_QU;
	bool m_freeEnd;
public:
	MGsfpDiff(
		bool dev_is_sum,	//dev is upper bound of sum if true, is max if false.
		double dev,
		const std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2
		const MGBPointSeq& QtDDQ,
		const MGBPointSeq& QtY,
		const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
		MGBPointSeq& U,
		MGBPointSeq& QU,
		bool freeEnd=false
	):m_sum(dev_is_sum),m_dev(dev*1.01),m_dx(dx),m_QtDDQ(QtDDQ),
	m_QtY(QtY),m_dy(dy),m_U(U),m_QU(QU),m_freeEnd(freeEnd){;};

	double operator()(double p)const{
		double sfp,max_dev_at_p;
		if(m_freeEnd)
			solveFreeE(p,m_dx,m_QtDDQ,m_QtY,m_dy,m_U,m_QU,sfp,max_dev_at_p);
		else
			solve(p,m_dx,m_QtDDQ,m_QtY,m_dy,m_U,m_QU,sfp,max_dev_at_p);
		if(m_sum) return m_dev-sfp;
		return m_dev-max_dev_at_p;
	};

	double dev()const{return m_dev;};

};

#define zero 1.e-7
double compute_smoothed_p(
	const MGLBRepEndC& begin,//Begin end condition
	const MGLBRepEndC& end,	//End end conditoion
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	double deviation,//if dev_is_sum is true,
		//deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		//if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		//dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		//dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		//for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum,
	std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2, will be output
	MGBPointSeq& U,	//df2(tau[i])/(6.*p) will be returned, where df2 is 2nd derivative
					//and p is the output of the compute_smoothed_p.
	MGBPointSeq& QU	//Q*U will be output where Q is the tridiagonal matrix of order n
					//havine the general row [1/dx[i-1], -(1/dx[i-1]+1/dx[i]), dx[i]]
					//for i=0,...,n-1
){
	MGBPointSeq QtDDQ,QtY;
	constructQTD(begin,end,tau,y,dy,dx,QtDDQ,QtY);
	//cout<<dx<<endl<<"***QtDDQ="<<QtDDQ<<endl<<"QtY="<<QtY<<endl;

    double p=zero;
	double sfp,max_dev_at_p;
	solve(p,dx,QtDDQ,QtY,dy,U,QU,sfp,max_dev_at_p);//get U, QU, and sfp.

	MGsfpDiff sfpdiff(dev_is_sum,deviation,dx,QtDDQ,QtY,dy,U,QU);
	double dev=sfpdiff.dev();
    if(dev_is_sum){
		if(sfp<=dev) return p;
	}else{
		if(max_dev_at_p<=dev) return p;
	}

	int ier;
	const double dpmin = .002;
	return mgNlbit(sfpdiff, zero,1., dev*dpmin, 20, ier);
}

//compute_smoothed_p for free end condition SR.
double compute_smoothed_p(
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	double deviation,//if dev_is_sum is true,
		//deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		//if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		//dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		//dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		//for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum,
	std::vector<double>& dx,	//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2, will be output
	MGBPointSeq& U,	//df2(tau[i])/(6.*p) will be returned, where df2 is 2nd derivative
					//and p is the output of the compute_smoothed_p.
	MGBPointSeq& QU	//Q*U will be output where Q is the tridiagonal matrix of order n
					//havine the general row [1/dx[i-1], -(1/dx[i-1]+1/dx[i]), dx[i]]
					//for i=0,...,n-1
){
	MGBPointSeq QtDDQ,QtY;
	constructQTDFreeE(tau,y,dy,dx,QtDDQ,QtY);
	//cout<<dx<<endl<<"***QtDDQ="<<QtDDQ<<endl<<"QtY="<<QtY<<endl;

    double p=zero;
	double sfp,max_dev_at_p;
	solveFreeE(p,dx,QtDDQ,QtY,dy,U,QU,sfp,max_dev_at_p);//get U, QU, and sfp.
	//cout<<"U="<<endl<<U<<endl; cout<<"QU="<<endl<<QU<<endl;

	MGsfpDiff sfpdiff(dev_is_sum,deviation,dx,QtDDQ,QtY,dy,U,QU,true);
	double dev=sfpdiff.dev();
    if(dev_is_sum){
		if(sfp<=dev) return p;
	}else{
		if(max_dev_at_p<=dev) return p;
	}

	int ier;
	const double dpmin = .002;
	return mgNlbit(sfpdiff, zero,1., dev*dpmin, 20, ier);
}

//Build the MGLBRep from the output of compute_smoothed_p.
void build_LBRep(
	double p,
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const std::vector<double>& dx,//dx[i]=tau[i+1]-tau[i]] for i=0,.., n-2, is input
	const double* dy,		// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	const MGBPointSeq& U,	//df2(tau[i])/(6.*p) is input, where df2 is 2nd derivative
					//and p is the output of the compute_smoothed_p.
	const MGBPointSeq& QU,	//Q*U is input where Q is the tridiagonal matrix of order n
					//havine the general row [1/dx[i-1], -(1/dx[i-1]+1/dx[i]), dx[i]]
					//for i=0,...,n-1
	MGLBRep& lbrep	//LBRep will be output.
){
	const int k=4;//order is 4.
	int i,j, n=tau.length();
	int nm1=n-1;
	int sd=U.sdim();
    double six1mp=(1.-p)*6.;
    double sixp=6.*p;
	MGPPRep pp(k,n,sd);
	for(i=0; i<n; i++)
		pp.break_point(i)=tau[i];
    for(j=0; j<sd; ++j){
		for(i=0; i<=nm1; ++i){
		    double dyi = dy[i];
			pp(0,i,j)=y(i,j)-six1mp*dyi*dyi*QU(i,j);//postional data.
			pp(2,i,j)=sixp*U(i,j);// Computing 2nd power/2.
		}
		for(i=0; i<nm1; ++i){
			double dxi=dx[i];
			pp(3,i,j)=(pp(2,i+1,j)-pp(2,i,j))/dxi;
			pp(1,i,j)=(pp(0,i+1,j)-pp(0,i,j))/dxi -(pp(2,i,j)+pp(3,i,j)/3.*dxi)/2.*dxi;
		}
    }

	int brep_dim=n+2;
	MGKnotVector t(k,brep_dim);//Order is 4.
	for(i=0; i<k; i++) {
		t(i)=pp.break_point(0);
		t(i+brep_dim)=pp.break_point(nm1);
	}
	j=1; i=k;
	while(i<brep_dim)
		t(i++)=pp.break_point(j++);
	MGLBRep lbpp(pp,t);
	lbrep.knot_vector()=lbpp.knot_vector();
	lbrep.line_bcoef()=lbpp.line_bcoef();
}

//Build line B-Rep by Schoenberg and Reinsch smoothing function, given
//1st derivatives on the start and end points, data points (tau,y),
//weights dy at data points, and a mean deviation deviation.
//If dy[i] gets larger, deviation at tau(i) gets larger.
//n can be any number greater than or equal to 2.
void MGLBRep::buildSRSmoothedLB_of_1stDeriv(
	const MGLBRepEndC& begin,//Begin end condition
	const MGLBRepEndC& end,	//End end conditoion.
		//begin.cond() and end.cond() must be MGENDC_1D or MGENDC_12D.
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const double* dy,// dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	double deviation,//if dev_is_sum is true,
		//deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		//if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		//dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		//dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		//for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum
){
	if(tau.length()<=2 || deviation<=MGTolerance::wc_zero()*.1){
		int error;
		*this=MGLBRep(begin,end,tau,y,error);
		return;
	}

	std::vector<double> dx;
	MGBPointSeq U, QU;
	double p=compute_smoothed_p(begin,end,tau,y,dy,deviation,dev_is_sum,dx,U,QU);

// CORRECT VALUE OF P HAS BEEN FOUND.
// COMPUTE POL.COEFFICIENTS FROM  Q*U AND Build B-Rep.
	build_LBRep(p,tau,y,dx,dy,U,QU,*this);
}

//Build line B-Rep by Schoenberg and Reinsch smoothing function, supposing the end
//conditions are free end conditions, given
//data points (tau,y), weights dy at data points, and a deviation.
//If dy[i] gets larger, deviation at tau(i) gets larger.
//n can be any number greater than or equal to 2.
//***End conditions are free end condition.***
void MGLBRep::buildSRSmoothedLB_of_FreeEnd(
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& y,	//Data point ordinates.
	const double* dy,//dy[i] is the weights  at tau[i] for i=0,..., tau.length()-1.
	double deviation,//if dev_is_sum is true,
		//deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		//if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		//dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		//dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		//for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum
){
	if(tau.length()<=2 || deviation<=MGTolerance::wc_zero()*.1){
		*this=MGLBRep(tau,y);
		return;
	}

	std::vector<double> dx;
	MGBPointSeq U, QU;
	double p=compute_smoothed_p(tau,y,dy,deviation,dev_is_sum,dx,U,QU);

// CORRECT VALUE OF P HAS BEEN FOUND.
// COMPUTE POL.COEFFICIENTS FROM  Q*U AND Build B-Rep.
	build_LBRep(p,tau,y,dx,dy,U,QU,*this);
}

