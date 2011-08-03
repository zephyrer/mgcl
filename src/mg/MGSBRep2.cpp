/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Vector.h"
#include "mg/Position.h"
#include "mg/KnotArray.h"
#include "mg/KnotVector.h"
#include "mg/BPointSeq.h"
#include "mg/SPointSeq.h"
#include "mg/LBRep.h"
#include "mg/Surface.h"
#include "mg/SBRep.h"
#include "mg/Plane.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/Bsgsmt.h"
#include "cskernel/Bluprt.h"
#include "cskernel/Blumor.h"
#include "cskernel/Blgl2a.h"
#include "cskernel/Bludkt.h"
#include "cskernel/Bluakt.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSBRep2.cpp
//
// Implements Surface B-Representation class MGSBRep.

//<< Constructor >>
								  
//**** 1. Approximation Constructor ****

MGSBRep::MGSBRep(			//BLUDKT
	const MGSBRep& old_brep,	//Original B-Rep.
	const MGKnotVector& tu,		//knot vector of u-direction
	const MGKnotVector& tv,		//knot vector of v-direction
	int &error)					//Error flag.
//Approximate an original B-Rep by a new knot configuration.
//The new knot config must be inside the range of the original B-Rep
//parameter. However, new knots may be coarse or fine.
:MGSurface(old_brep)
,m_uknot(tu), m_vknot(tv)
{
	update_mark();
	assert(tu.order()==old_brep.order_u()
		&& tv.order()==old_brep.order_v());

	size_t lenu=tu.bdim(), lenv=tv.bdim(), ncd=old_brep.sdim();
	unsigned ku=tu.order(), kv=tv.order();
	unsigned kmax=ku; if(kmax<kv) kmax=kv;
	size_t nu, nv, sizeu, sizev;
	old_brep.m_surface_bcoef.length(nu, nv);
	old_brep.m_surface_bcoef.capacity(sizeu, sizev);
	size_t m;if(kmax<=5) m=9;else m=2*kmax-1;
	size_t n2=lenu; if(n2<lenv) n2=lenv;
	double* work1=new double[n2];
	double* work2=new double[n2*m];

	MGSPointSeq surf1(lenv, lenu, ncd);
	MGSPointSeq surf2(lenu, lenv, ncd);
	MGBPointSeq temp(n2,ncd);

	size_t irc=sizeu*sizev;
	for(size_t j=0; j<nv; j++){
		bludkt_(ku,nu,old_brep.knot_data_u(),old_brep.coef_data(0,j,0),
				irc,ncd,lenu,tu.data(),n2,
				work1,work2,&temp(0,0),&error);
		if(error!=1) break;
		for(size_t i=0; i<lenu; i++) for(size_t k=0; k<ncd; k++) surf1(j,i,k)=temp(i,k);
	}
	if(error==1){
		irc=lenu*lenv;
		for(size_t i=0; i<lenu; i++){
			bludkt_(kv,nv,old_brep.knot_data_v(),surf1.data(0,i,0),
					irc,ncd,lenv,tv.data(),n2,
					work1,work2,&temp(0,0),&error);
			if(error!=1) break;
			for(size_t j=0; j<lenv; j++) for(size_t k=0; k<ncd; k++) surf2(i,j,k)=temp(j,k);
		}
	}
							   
	delete[] work1; delete[] work2;
	if(error==1){
		error=0;		//Return of BLUDKT:error=1 means normal.
		m_surface_bcoef=surf2;
	}
	else{
		*this=old_brep;
	}
}

//**** 2.Conversion Constructor.****

MGSBRep::MGSBRep(		//BLUAKT
	const MGSBRep& old_brep,	//Original B-Rep.
	const MGKnotArray& uknots,	//Knots to add for u-direction
	const MGKnotArray& vknots)	//Knots to add for v-direction.
//Gets new B-Rep by adding knots to an original B-Rep.
:MGSurface(old_brep)
{
	update_mark();
	size_t ncd=old_brep.sdim();
	unsigned ku=old_brep.order_u(), kv=old_brep.order_v();
	size_t nu, nv, sizeu, sizev;
	old_brep.m_surface_bcoef.length(nu, nv);
	old_brep.m_surface_bcoef.capacity(sizeu, sizev);

	unsigned kmax=ku; if(kmax<kv) kmax=kv;
	double* work=new double[kmax*kmax];

	const int nadu=uknots.length();
	double* tadu=new double[nadu];
	int* mltu=new int[nadu];
	unsigned mltu_total=0;
	for(int i=0; i<nadu; i++){
		tadu[i]=uknots[i].value();
		mltu[i]=uknots[i].multiplicity();
		mltu_total+=mltu[i];
	}
	size_t lenu=old_brep.bdim_u()+mltu_total;

	const int nadv=vknots.length();
	double* tadv=new double[nadv];
	int* mltv=new int[nadv];
	unsigned mltv_total=0;
	for(int j=0; j<nadv; j++){
		tadv[j]=vknots[j].value();
		mltv[j]=vknots[j].multiplicity();
		mltv_total+=mltv[j];
	}
	size_t lenv=old_brep.bdim_v()+mltv_total;
 
	MGSPointSeq surf1(lenv, lenu, ncd);
	MGSPointSeq surf2(lenu, lenv, ncd);
	size_t lmax=lenu; if(lmax<lenv) lmax=lenv;
	MGBPointSeq temp(lmax,ncd);
	MGKnotVector tu(ku,lenu), tv(kv,lenv);

	int lenunew,lenvnew;
	size_t irc=sizeu*sizev;
	for(size_t j=0; j<nv; j++){
		bluakt_(ku, nu, old_brep.knot_data_u(), old_brep.coef_data(0,j,0),
			irc, ncd, nadu, tadu, mltu, lmax, work,
			&lenunew, &tu(0), &temp(0,0));
		for(int i=0; i<lenunew; i++) for(size_t k=0; k<ncd; k++) surf1(j,i,k)=temp(i,k);
	}
	tu.set_bdim(lenunew);

	irc=lenu*lenv;
	for(int i=0; i<lenunew; i++){
		bluakt_(kv, nv, old_brep.knot_data_v(), surf1.data(0,i,0),
			irc, ncd, nadv, tadv, mltv, lmax, work,
			&lenvnew, &tv(0), &temp(0,0));
		for(int j=0; j<lenvnew; j++) for(size_t k=0; k<ncd; k++) surf2(i,j,k)=temp(j,k);
	}
	tv.set_bdim(lenvnew);
	surf2.set_length(lenunew,lenvnew);
							   
	delete[] tadu; delete[] mltu; delete[] work;
	delete[] tadv; delete[] mltv;

	m_uknot=tu; m_vknot=tv;
	m_surface_bcoef=surf2;
}

//**** 3.Approximation Constructor.****

MGSBRep::MGSBRep(			//Least square approximation.
	const MGNDDArray& utau,	//Data point abscissa of u-direction
	const MGNDDArray& vtau,	//Data point abscissa of v-direction
	const MGSPointSeq& points,	//Point seq data
	const MGKnotVector& tu,	//knot vector of u-direction
	const MGKnotVector& tv,	//knot vector of v-direction
	const MGSPointSeq& weight)	//Weight for each point(space dimension 1).
			//weight(i,j,0) is the weight for points(i,j,.).
:MGSurface(),
m_surface_bcoef(tu.bdim(),tv.bdim(),points.sdim()),
m_uknot(tu), m_vknot(tv)
{
	const size_t nutau=utau.length(),nvtau=vtau.length();
	const size_t ip1=points.capacity_u(),iw=weight.capacity_u();
	const size_t nu=tu.bdim(),nv=tv.bdim();
	size_t nmax=nu; if(nmax<nv) nmax=nv;
	const unsigned ku=tu.order(),kv=tv.order();
	unsigned kmax=ku; if(kmax<kv) kmax=kv;
	const unsigned sdim=points.sdim();
	double* work=new double[2*nmax];
	double* q=new double[kmax*nmax];
	double* swork=new double[nu*nvtau];		//swork[nu][nvtau].

	//Generate weight of u B-coef.
	size_t nubynvtau=nu*nvtau;
	double* vweight=new double[nubynvtau];	//Actually vweight[nu][nvtau].
	for(size_t i=0; i<nubynvtau; i++) vweight[i]=0.0;
	size_t id;
	for(size_t j=0; j<nvtau; j++){
		for(size_t i=0; i<nutau; i++)	{
			id=m_uknot.eval_coef(utau(i),q);
			for(size_t m=0; m<ku; m++) vweight[j+nvtau*(id+m)]+=q[m]*weight(i,j,0);
		}
	}

	int kw=1;
	for(size_t m=0; m<sdim; m++){
		blgl2a_(nutau,utau.data(),points.data(0,0,m),weight.data(),
			ip1,iw,nvtau,kw,tu.data(),nu,ku,nvtau,work,q,swork);

		blgl2a_(nvtau,vtau.data(),swork,vweight,
			nvtau,nvtau,nu,kw,tv.data(),nv,kv,nu,work,q,
			&m_surface_bcoef(0,0,m));
	}

	delete[] work; delete[] q; delete[] swork; delete[] vweight;
}

MGSBRep::MGSBRep(		//Shoenberg and Reinch's smoothing function
						//approximation.
	const MGNDDArray& utau,	//Data point abscissa of u-direction
	const MGNDDArray& vtau,	//Data point abscissa of v-direction
	const MGSPointSeq& points,	//Point seq data
	const MGSPointSeq& delp,	//Error estimate for each point
		//(space dimension 1).
		// delp(i,j,0) is for points(i,j,.) at (utau(i),vtau(j)).
	double deviation)		//Mean deviation of each point.
//If delp(i,j,0) becomes bigger, deviation at (utau(i),vtau(j)) becomes bigger.
:MGSurface(),
m_surface_bcoef(utau.length()+2,vtau.length()+2,points.sdim()),
m_uknot(4,utau.length()+2), m_vknot(4,vtau.length()+2)
{
	size_t nu=utau.length(),nv=vtau.length();
	size_t npmax=nu;if(npmax<nv) npmax=nv;
	size_t nup2=nu+2, nvp2=nv+2;
	size_t nup2bynv=nup2*nv;
	size_t ip1,ip2; points.capacity(ip1,ip2);
	size_t ip12=ip1*ip2;
	size_t ncd=points.sdim();
	double* v=new double[npmax*7];
	double* a4=new double[npmax*4];
	double* a3=new double[npmax*ncd];
	double* pwork=new double[(npmax+2)*ncd];	//pwork[ncd][npmax+2]
	double* swork=new double[nup2bynv*ncd];		//swork[ncd][nv][nup2]

	int lud;
	for(size_t j=0; j<nv; j++){	//Process of u-direction.
		bsgsmt_(ncd,nu,utau.data(),points.data(0,j,0),delp.data(0,j,0),
		deviation,ip12,nv,nup2,v,a4,a3,pwork,&lud,&m_uknot(0),
		swork+j);
	}

	//Generate delp of u B-coef.
	double* delpu=new double[nup2bynv];	//Actually delpu[nup2][nv].
	for(size_t i=0; i<nup2bynv; i++) delpu[i]=0.0;
	size_t id;
	double coef[4];
	for(size_t j=0; j<nv; j++){
		for(size_t i=0; i<nu; i++)	{
			id=m_uknot.eval_coef(utau(i),coef);
			for(size_t m=0; m<4; m++) delpu[j+nv*(id+m)]+=coef[m]*delp(i,j,0);
		}
	}
							//Process of v-direction.
	int nvi,lvd;
	for(size_t i=0; i<nup2; i++){
		nvi=nv*i;
		bsgsmt_(ncd,nv,vtau.data(),swork+nvi,delpu+nvi,
		deviation,nup2bynv,nup2,nvp2,v,a4,a3,pwork,&lvd,&m_vknot(0),
		&m_surface_bcoef(i,0,0));
	}
	delete[] v;	delete[] a4; delete[] a3; delete[] pwork;
	delete[] swork; delete[] delpu;
}

//Gets new B-Rep by connecting two B-Rep to one.
MGSBRep::MGSBRep(
	const MGSBRep& brep1,	//B-Rep 1.
	int which1,				//which perimeter of brep1.
	int continuity,			//continuity.
	const MGSBRep& brep2,	//B-Rep 2.
	int which2,				//which perimeter of brep2.
	int opposite			// Input if parameter direction of which2
					// is the same as which1 along common edge.
					// If opposite is true, the direction is opposite.
):MGSurface(brep1){
	update_mark();
	assert(continuity>=0);
	assert(which1>=0 && which1<4 && which2>=0 && which2<4);
	assert(brep1.sdim()<=3 && brep2.sdim()<=3);

	size_t nu,nv;
	MGSBRep sb2(brep2);
	if(opposite){
		if(which2==0 || which2==2) sb2.negate(1); //Reverse u-direction.
		else                       sb2.negate(0); //Reverse v-direction.
	}
	if(which1==0 || which1==2){
		if(which2==1 || which2==3){
			sb2.exchange_uv(); which2++; if(which2>3) which2=0;
		}
	}else{
		if(which2==0 || which2==2){
			sb2.exchange_uv(); which2--; if(which2<0) which2=3;
		}
	}

	int lwhich;
	MGLBRep lb, lb1, lb2;
	size_t sd1=brep1.sdim(),sd2=sb2.sdim();
	if(which1==0 || which1==2){
		m_uknot=brep1.knot_vector_u();
		size_t nv1=brep1.bdim_v(), nv2=sb2.bdim_v();
		MGBPointSeq bp1(nv1, sd1); MGBPointSeq bp2(nv2, sd2);
		for(size_t j=0; j<nv1; j++) for(size_t m=0; m<sd1; m++) bp1(j,m)=brep1.coef(0,j,m);
		for(size_t j=0; j<nv2; j++) for(size_t m=0; m<sd2; m++) bp2(j,m)=sb2.coef(0,j,m);
		lb1=MGLBRep(brep1.knot_vector_v(), bp1);
		lb2=MGLBRep(sb2.knot_vector_v(), bp2);
		if(which1==0){
			if(which2==0) lwhich=0;
			else          lwhich=1;
		}else{
			if(which2==0) lwhich=2;
			else          lwhich=3;
		}
		lb=MGLBRep(lb1,continuity,lwhich,lb2); 
		nu=m_uknot.bdim(); nv=lb.bdim();
		size_t ncd=lb.sdim();
		m_surface_bcoef=MGSPointSeq(nu,nv,ncd);
		m_vknot=lb.knot_vector();
		for(size_t j=0; j<nv; j++)
			for(size_t m=0; m<ncd; m++) m_surface_bcoef(0,j,m)=lb.coef(j,m);
		for(size_t i=1; i<nu; i++){
			for(size_t j=0; j<nv1; j++) 
				for(size_t m=0; m<sd1; m++) bp1(j,m)=brep1.coef(i,j,m);
			for(size_t j=0; j<nv2; j++) 
				for(size_t m=0; m<sd2; m++) bp2(j,m)=sb2.coef(i,j,m);
			lb1=MGLBRep(brep1.knot_vector_v(), bp1);
			lb2=MGLBRep(sb2.knot_vector_v(), bp2);
			lb=MGLBRep(lb1,continuity,lwhich,lb2); 
			for(size_t j=0; j<nv; j++)
				for(size_t m=0; m<ncd; m++) m_surface_bcoef(i,j,m)=lb.coef(j,m);
		}
	}else{
		m_vknot=brep1.knot_vector_v();
		size_t i;
		size_t nu1=brep1.bdim_u(), nu2=sb2.bdim_u();
		MGBPointSeq bp1(nu1, sd1); MGBPointSeq bp2(nu2, sd2);
		for(i=0; i<nu1; i++)
			for(size_t m=0; m<sd1; m++)
				bp1(i,m)=brep1.coef(i,0,m);
		for(i=0; i<nu2; i++)
			for(size_t m=0; m<sd2; m++)
				bp2(i,m)=sb2.coef(i,0,m);
		lb1=MGLBRep(brep1.knot_vector_u(), bp1);
		lb2=MGLBRep(sb2.knot_vector_u(), bp2);
		if(which1==3){
			if(which2==3) lwhich=0;
			else          lwhich=1;
		}else{
			if(which2==3) lwhich=2;
			else          lwhich=3;
		}
		lb=MGLBRep(lb1,continuity,lwhich,lb2); 
		nu=lb.bdim(); nv=m_vknot.bdim();
		size_t ncd=lb.sdim();
		m_surface_bcoef=MGSPointSeq(nu,nv,ncd);
		m_uknot=lb.knot_vector();
		for(i=0; i<nu; i++)
			for(size_t m=0; m<ncd; m++)
				m_surface_bcoef(i,0,m)=lb.coef(i,m);
		for(size_t j=1; j<nv; j++){
			for(i=0; i<nu1; i++) 
				for(size_t m=0; m<sd1; m++)
					bp1(i,m)=brep1.coef(i,j,m);
			for(i=0; i<nu2; i++) 
				for(size_t m=0; m<sd2; m++)
					bp2(i,m)=sb2.coef(i,j,m);
			lb1=MGLBRep(brep1.knot_vector_u(), bp1);
			lb2=MGLBRep(sb2.knot_vector_u(), bp2);
			lb=MGLBRep(lb1,continuity,lwhich,lb2); 
			for(i=0; i<nu; i++)
				for(size_t m=0; m<ncd; m++)
					m_surface_bcoef(i,j,m)=lb.coef(i,m);
		}
	}
}

// Gets new B-Rep by computing a part of the original. New one is exactly
// the same as the original except that it is partial.
MGSBRep::MGSBRep(
	const MGBox& uvrange,		//u and v parameter range.
	const MGSBRep& old_brep,	//Original B-Rep.
	int multiple)	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
//If multiple==true(!=0), knot_u(i)=t1 and knot_u(n+i)=t2 for i=0,..., k-1
//will be guaranteed. Here, n=bdim_u(), k=order_u(),
//t1=uvrange(0).low_point(), and t2=uvrange(0).high_point().
//About knot_v(j), the same.
// Both u-range and v-range must be inside the range of old_brep.
:MGSurface(old_brep){
	update_mark();
	MGBox b2=old_brep.param_range();
	MGBox b3=b2&uvrange;;
	if(b3 == b2){ *this=old_brep; return; }

	size_t ncd=old_brep.sdim();
	unsigned ku=old_brep.order_u(), kv=old_brep.order_v();
	unsigned kmax=ku; if(kmax<kv) kmax=kv;
	size_t nu, nv, sizeu, sizev;
	old_brep.m_surface_bcoef.length(nu, nv);
	old_brep.m_surface_bcoef.capacity(sizeu, sizev);
	size_t nmax=nu; if(nmax<nv) nmax=nv;
	double* work=new double[kmax*kmax];

	MGSPointSeq surf1(nv, nu, ncd);
	MGSPointSeq surf2(nu, nv, ncd);
	MGBPointSeq temp(nmax,ncd);
	MGKnotVector tu(ku, nu), tv(kv,nv);

	int nunew, nvnew;
	double t1,t2;
	size_t irc=sizeu*sizev;
	t1=b3(0).low_point(); t2=b3(0).high_point();
	for(size_t j=0; j<nv; j++){
		bluprt_(ku, nu, old_brep.knot_data_u(), old_brep.coef_data(0,j,0),
			irc,ncd,t1,t2,nmax,work,&nunew,&tu(0), &temp(0,0),multiple);
		for(int i=0; i<nunew; i++) for(size_t k=0; k<ncd; k++) surf1(j,i,k)=temp(i,k);
	}
	tu.set_bdim(nunew); tu.reshape(nunew+ku);

	irc=nu*nv;
	t1=b3(1).low_point(); t2=b3(1).high_point();
	for(int i=0; i<nunew; i++){
		bluprt_(kv,nv,old_brep.knot_data_v(),surf1.data(0,i,0),
			irc,ncd,t1,t2,nmax,work,&nvnew,&tv(0),&temp(0,0),multiple);
		for(int j=0; j<nvnew; j++) for(size_t k=0; k<ncd; k++) surf2(i,j,k)=temp(j,k);
	}
	tv.set_bdim(nvnew); tv.reshape(nvnew+kv);
	surf2.set_length(nunew, nvnew); surf2.reshape(nunew, nvnew);
							   
	delete[] work;

	m_uknot=tu; m_vknot=tv;
	m_surface_bcoef=surf2;
}

// Construct a Surface B-Rep by changing space dimension and
//ordering of coordinate.
MGSBRep::MGSBRep(
	size_t dim,				// New space dimension.
	const MGSBRep& sbrep,	// Original Surface B-rep.
	size_t start1, 			// Destination order of new Surface.
	size_t start2) 			// Source order of original Surface.
:MGSurface(sbrep)
,m_uknot(sbrep.knot_vector_u()), m_vknot(sbrep.knot_vector_v())
,m_surface_bcoef(dim,sbrep.surface_bcoef(),start1,start2){
	update_mark();
}

MGSBRep::MGSBRep(
		const MGPlane& plane,	// Original Plane
		const MGBox& prange)	// parameter range of new Surface.
// Construct a Surface B-Rep (order = 2 ) from Plane and Parameter ranges.
:MGSurface(plane),m_surface_bcoef(2,2,plane.sdim())
,m_uknot(2,2)
,m_vknot(2,2){
	update_mark();
// m_surface_bcoef‚ÉControl Polygon‚Ìƒf[ƒ^‚ð‘ã“ü
	double u[2] = {prange(0).low_point(), prange(0).high_point()};
	double v[2] = {prange(1).low_point(), prange(1).high_point()};		
	
	for(int i = 0; i < 2; i++){
		m_uknot(2*i)=m_uknot(2*i+1) = u[i];
		m_vknot(2*i)=m_vknot(2*i+1) = v[i];
		m_surface_bcoef.store_at(0, i, plane.eval(u[0], v[i]));
		m_surface_bcoef.store_at(1, i, plane.eval(u[1], v[i]));
	}
}

//Destructor
//MGSBRep::~MGSBRep(){;}

// Exchange parameter u and v.
MGSurface& MGSBRep::exchange_uv(){
	size_t nu=bdim_u(), nv=bdim_v(), dim=sdim();
	// Exchange surface coefficients.
	MGSPointSeq temp(nv,nu,dim);
	for(size_t k=0; k<dim; k++)
		for(size_t j=0; j<nv; j++)
			for(size_t i=0; i<nu; i++) temp(j,i,k)=coef(i,j,k);
	return *this=MGSBRep(temp, m_vknot, m_uknot);
}

//Modify the original Surface by extrapolating the specified perimeter.
//The extrapolation is C2 continuous if the order >=4.
//The extrapolation is done so that extrapolating length is "length"
//at the position of the parameter value "param" of the perimeter.
MGSBRep& MGSBRep::extend(
	int perimeter,	//perimeter number of the Surface.
					// =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	// parameter value of above perimeter.
	double length,	//chord length to extend at the parameter param of the perimeter.
	double dk      //Coefficient of how curvature should vary at
//    extrapolation start point. When dk=0, curvature keeps same, i.e.
//    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
//    i.e. dK/dS=-K/length at extrapolation start point.
//    (S=parameter of arc length, K=Curvature at start point)
//    That is, when dk reaches to 1 from 0, curve changes to flat.
){
	assert(sdim()<=3);
	assert(perimeter>=0 && perimeter<4);
	
	const size_t ncd=surface_bcoef().sdim();
	int at_start=1;//starting perimeter
	size_t nu, nv;
	size_t order; size_t n,m; MGKnotVector* t;
	if(perimeter==1 || perimeter==3){	// Extrapolate to u-direction
		order=order_u();
		n=bdim_u();
		t=&(knot_vector_u());
		if(perimeter==1) at_start=0;//ending perimeter
		nu=n+1; m=nv=bdim_v();
	}else{
		// Extrapolate to v-direction
		order=order_v();
		n=bdim_v();
		t=&(knot_vector_v());
		if(perimeter==2) at_start=0;//ending perimeter
		nv=n+1; m=nu=bdim_u();
	}
	//(nu,nv) are new surface B-Rep dimensions of u and v.
	//(order,n,t) is line B-rep to extrapolate.
	//m is the number of line B-reps to extrapolate.

	double* dcoef=new double[order];
	size_t id;
	if(at_start)
		id=t->eval_coef(t->param_s(),dcoef,1);	//To get 1st derivative.
	else
		id=t->eval_coef(t->param_e(),dcoef,1);	//To get 1st derivative.

	MGSPointSeq surf(nu,nv,ncd);
	MGLBRep lbtemp(n,order,ncd);
	MGKnotVector& t1=lbtemp.knot_vector();
	MGBPointSeq& coeftemp=lbtemp.line_bcoef();

	MGPosition uv=perimeter_uv(perimeter,param);//Surface parameter value of param.
	size_t ndu=0,ndv=0;
	if(perimeter==0 || perimeter==2) ndv=1;
	else                             ndu=1;
	double slen=length/(eval(uv,ndu,ndv)).len();

	int nnew,cid; double firstd_len,data,dlen;
	for(size_t i=0; i<m; i++){
		if(perimeter==0 || perimeter==2){
			for(size_t j=0; j<n; j++)
				for(size_t k=0; k<ncd; k++) coeftemp(j,k)=coef(i,j,k);
		}else{
			for(size_t j=0; j<n; j++)
				for(size_t k=0; k<ncd; k++) coeftemp(j,k)=coef(j,i,k);
		}
		coeftemp.set_length(n);

	//Compute first derivative length at the end of the extrapolating line.
		firstd_len=0.;
		for(size_t k=0; k<ncd; k++){
			data=0.;
			for(size_t j=0; j<order; j++){
				cid=id+j; data=data+coeftemp(cid,k)*dcoef[j];
			}
			firstd_len+=data*data;
		}
		firstd_len=sqrt(firstd_len);

		t1=*t; dlen=firstd_len*slen;
		//std::cout<<"before:"<<lbtemp<<std::endl;///////
		lbtemp.extend(at_start,dlen,dk);
		//std::cout<<"after:"<<lbtemp<<std::endl;///////
		nnew=lbtemp.bdim();
		if(perimeter==0 || perimeter==2){
			for(int j=0; j<nnew; j++)
				for(size_t k=0; k<ncd; k++) surf(i,j,k)=coeftemp(j,k);
		}else{
			for(int j=0; j<nnew; j++)
				for(size_t k=0; k<ncd; k++) surf(j,i,k)=coeftemp(j,k);
		}
	}

	*t=t1;
	surf.set_length(nu,nv);
	surface_bcoef()=surf;

	delete[] dcoef;
	update_mark();
	return *this;
}

MGSBRep& MGSBRep::move(		//BLUMOR
	int move_kind_u,		//Indicates how to move Surface for u direction.
	int move_kind_v,		//Indicates how to move Surface for v direction.
	const MGPosition& move_point_param, //indicates object point to move
							// by the (u,v) parameter value.
	const MGPosition& to_point,	//destination point of the abve source point.
	const MGPosition fix_point[2])
//Modify the original Surface by moving move_point to to_point. fix_point can be
//applied according to move_kind.
{
	int i,id;

	double u=move_point_param(0), v=move_point_param(1);
	if(u<param_s_u() || u>param_e_u() || v<param_s_v() || v>param_e_v())
		return *this;

	MGVector delta=to_point-eval(move_point_param);
	//delta is a vector from original to destination point.

	int ku=order_u(), kv=order_v();
	size_t nu=bdim_u(), nv=bdim_v();
	double* ucoef=new double[ku]; double* vcoef=new double[kv];
	double* uratio=new double[nu]; double* vratio=new double[nv];
	int ui,uj,vi,vj,umovek, vmovek;
	int uki=m_uknot.eval_coef(u,ucoef);
	int vki=m_vknot.eval_coef(v,vcoef);
	//ucoef[i] and vcoef[j] are coefficients that are multiplied to
	//B-Coef's. uki and vki are the ids of first B-Coef's. 0<=i<ku, 0<=j<kv. 
	double us,ue,vs,ve;			//For fixed point parameter values.
	int uis,uie,vis,vie;	//For start and end id's of B-Coef's that
								//should be modified.
	double usum,vsum;

	//Compute u-direction ratio in uratio.
	switch(move_kind_u){
	case (1): 
		ui=ku; uj=nu; umovek=1; break;
	case (2):
		us=fix_point[0](0); us=m_uknot.range(us);
		ui=m_uknot.locate(us)+1;
		if(us<=u){umovek=2; uj=nu;}
		else     {umovek=3; uj=ui; ui=1;}
		break;
	case (3):
		us=fix_point[0](0); us=m_uknot.range(us);
		ue=fix_point[1](0); ue=m_uknot.range(ue);
		if(us>ue){ double uu=us; us=ue; ue=uu;}
		ui=m_uknot.locate(m_uknot.range(us))+1;
		uj=m_uknot.locate(m_uknot.range(ue))+1;
		if(ue<=u)    {umovek=2; ui=uj; uj=nu;}
		else if(u<us){umovek=3; uj=ui; ui=1;}
		else         {umovek=1;}
		break;
	default:
		umovek=1; uj=ui=(uki+ku); break;
	}
	blumor_(umovek,&ui,&uj,u,ku,nu,m_uknot.data(),&uis,&uie,uratio);
	uis-=1; uie-=1; 
	id=0; usum=0;
	for(i=uki; i<uki+ku; i++){
		if(i>=uis && i<=uie) usum=usum+ucoef[id]*uratio[i];
		id+=1;
	}
	if(MGMZero(usum)) goto end_p;

	//Compute v-direction ratio in vratio.
	switch(move_kind_v){
	case (1): 
		vi=kv; vj=nv; vmovek=1; break;
	case (2):
		vs=fix_point[0](1); vs=m_vknot.range(vs);
		vi=m_vknot.locate(vs)+1;
		if(vs<=v){vmovek=2; vj=nv;}
		else     {vmovek=3; vj=vi; vi=1;}
		break;
	case (3):
		vs=fix_point[0](1); vs=m_vknot.range(vs);
		ve=fix_point[1](1); ve=m_vknot.range(ve);
		if(vs>ve){ double vv=vs; vs=ve; ve=vv;}
		vi=m_vknot.locate(m_vknot.range(vs))+1;
		vj=m_vknot.locate(m_vknot.range(ve))+1;
		if(ve<=v)    {vmovek=2; vi=vj; vj=nv;}
		else if(v<vs){vmovek=3; vj=vi; vi=1;}
		else         {vmovek=1;}
		break;
	default:
		vmovek=1; vj=vi=(vki+kv); break;
	}
	blumor_(vmovek,&vi,&vj,v,kv,nv,m_vknot.data(),&vis,&vie,vratio);
	vis-=1; vie-=1;
	id=0; vsum=0.;
	for(int j=vki; j<vki+kv; j++){
		if(j>=vis && j<=vie) vsum=vsum+vcoef[id]*vratio[j];
		id+=1;
	}
	if(MGMZero(vsum)) goto end_p;

	//New B-Coefficients.
	double ratio1,ratio2;

	for(i=uis; i<=uie; i++){
		ratio1=uratio[i]/usum;
		for(int j=vis; j<=vie; j++){
			ratio2=ratio1*vratio[j]/vsum;
			for(size_t m=0; m<sdim(); m++)
				m_surface_bcoef(i,j,m)+=delta.ref(m)*ratio2;
		}
	}

end_p:
	delete[] ucoef; delete[] vcoef; delete[] uratio; delete[] vratio;
	update_mark();
	return *this;
}

void MGSBRep::negate(	//Change direction of the surface.
	int is_u)				// Negate along u-direction if is_u is ture,
							// else along v-direction.
{
	if(is_u){		//u-direction.
		m_uknot.reverse(); m_surface_bcoef.reverse(is_u);
	}else{
		m_vknot.reverse(); m_surface_bcoef.reverse(is_u);
	}
}

//Obtain parameter value if this surface is negated by "negate()".
// Negate along u-direction if is_u is ture,
// else along v-direction.
MGPosition MGSBRep::negate_param(const MGPosition& uv, int is_u)const{
	double u=uv(0), v=uv(1);
	if(is_u){		//u-direction.
		u=m_uknot.reverse_param(u);
	}else{
		v=m_vknot.reverse_param(v);
	}
	return MGPosition(u,v);
}

//Shrink this surface to the part limitted by the parameter range of uvbx.
//New parameter range uvbx2 is so determined that uvbx2 is the smallest
//box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
//the values of u or v knots of the surface knotvector.
//uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
//That is u range is from us to ue , and so on.
void MGSBRep::shrink_to_knot(
	const MGBox& uvbx,
	int multiple	//Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
){
	MGBox uvb2=box_param();
	uvb2&=uvbx;

	MGInterval& uspan=uvb2[0];
	double u0=uspan.low_point(), u1=uspan.high_point();
	int idu0=m_uknot.locate(u0), idu1=m_uknot.locate(u1)+1;
	uspan.set_low(m_uknot[idu0]);
	uspan.set_high(m_uknot[idu1]);

	MGInterval& vspan=uvb2[1];
	double v0=vspan.low_point(), v1=vspan.high_point();
	int idv0=m_vknot.locate(v0), idv1=m_vknot.locate(v1)+1;
	vspan.set_low(m_vknot[idv0]);
	vspan.set_high(m_vknot[idv1]);

	MGSBRep sbnew(uvb2,*this,multiple);
	m_surface_bcoef=sbnew.m_surface_bcoef;
	m_uknot=sbnew.m_uknot;
	m_vknot=sbnew.m_vknot;
}
