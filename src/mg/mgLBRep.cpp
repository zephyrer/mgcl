/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/KnotArray.h"
#include "mg/CParam_list.h"
#include "mg/BPointSeq.h"
#include "mg/OscuCircle.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"
#include "mg/PPRep.h"
#include "mg/LBRepEndC.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/SurfCurve.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/Bluprt.h"
#include "cskernel/Blumix.h"
#include "cskernel/bkdtpg.h"
#include "cskernel/Bludkt.h"
#include "cskernel/Bluakt.h"
#include "cskernel/Blgl2a.h"
#include "cskernel/Blgint.h"
#include "cskernel/Bkdnp.h"
#include "cskernel/bkdtkt.h"
#include "cskernel/Blctpb.h"
#include "cskernel/blg4sc.h"
#include "cskernel/blg4sp2.h"
#include "cskernel/blgcs.h"
#include "cskernel/mgblgsq.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// MGLBRep.cpp
//
// Implement MGLBRep class.
//
// This file contains all of the constructors of MGLBRep.

//<< Constructor >>

//Copy constructor.
//MGLBRep::MGLBRep(const MGLBRep& lb2)
//:MGCurve(lb2), m_knot_vector(lb2.m_knot_vector)	,m_line_bcoef(lb2.m_line_bcoef)
//{;}

//Construct Line B-Representation, providing all the data.
MGLBRep::MGLBRep(
	const MGKnotVector& t,		//Knot Vector.
	const MGBPointSeq& bcoef)	//Line B-Coef.
:MGCurve(), m_knot_vector(t),m_line_bcoef(bcoef){
	assert(t.bdim()==bcoef.length());
}

	//**** 1. Interpolation Constructor ****

MGLBRep::MGLBRep(
	const MGBPointSeq& points,	//Point seq data
	int& error,					//Error flag.
	unsigned order,				// Order
	int circular)				//Circular flag
// Construct Line B-rep by intepolation from Point data only.
//If circular is true, start and end points are connected smoothly.
//If circular is true, order will be always 4, and input order is neglected.
:MGCurve(), m_knot_vector(order,points.length()+2)
	,m_line_bcoef(points.length()+2,points.sdim())
{
	assert(order>0 && points.length()>=2);

	MGBPointSeq points_temp(points);//Copy, since BLG4SC, BLGSPL may change
									//input data.
	size_t nv=points_temp.length();	//Number of input points.
	const size_t iv=points_temp.capacity();
	const unsigned ncd=points_temp.sdim();
	const size_t irc=m_line_bcoef.capacity();
	int n;

	if(circular){
		if(order!=4) m_knot_vector.change_order(4);
		double* work1=new double[irc];
		double* work2=new double[irc*9];
		blg4sc_(nv, &points_temp(0,0), iv, ncd, irc, work1,work2,
				&n, &m_knot_vector(0), &m_line_bcoef(0,0), &error);
		delete[] work1;
		delete[] work2;
	}else{
		if(order==2){	//Generate polyline B-Rep. of order 2.
			m_knot_vector=MGKnotVector(2,nv,0.,1.);
			m_line_bcoef=points;
			error=0;
			return;
		}else{
			if(nv<order){
				order=nv;
				m_knot_vector.size_change(order,nv);
			}
			double* work=new double[nv*order*2];
			n=nv;
			//   GENERATE DATA POINTS IN WORK(I,1) 0<=I<=N-1 . 
			error=2;
			double* val=&points_temp(0,0);
			bkdtpg_(val, n, ncd, iv, work);
			//   DISCARD TOO NEAR POINTS 
			int nnew = n;
		    bkdnp_(&nnew,work,val,iv,ncd,1,MGTolerance::max_knot_ratio());
			if(nnew>=int(order)){
				//   GENERATE KNOT VECTOR IN T FROM DATA POINTS. 
				double* t=&m_knot_vector(0);
				bkdtkt_(work, nnew, order, t);
				//   GENERATE B-REP. 
				error=blgint_(work,val,t,order,nnew,ncd,iv,irc,work+n,&m_line_bcoef(0,0));
				n = nnew;
			}
			delete[] work;
		}
	}

	if(error==1) {
		error=0;		//Return of BLG4SC,blgint_ error=1 means normal.
		m_knot_vector.set_bdim(n);
		m_line_bcoef.set_length(n);
	}
}

// Construct Line B-rep of a specified order, given data point abscissa and
//the ordinates.
MGLBRep::MGLBRep(
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& points,	//Point seq data
	size_t order,				//order
	double ratio)			//Maximum of data point ratio of pre and after spans.
			// Let d(i)=tau[i]-tau[i-1], then if d(i)/d(i-1)>ratio or
			// d(i-1)/d(i)>ratio, either tau[i] or tau[i-1] will be removed.
			// This is done to prevent control polygon computation error.
			// When ratio<0. no data point removal will be done.
:MGCurve(),m_line_bcoef(points.length(),points.sdim()){
	int IMLT=1;
	MGNDDArray tau2(tau);
	MGBPointSeq points2(points);
	size_t n=points2.length(), pointSize=points2.capacity(), pointDim=points2.sdim();
	if(ratio>0.){
		bkdnp_((int*)&n,tau2.data(),points2.data(),pointSize,pointDim,IMLT,ratio);
		tau2.set_length(n); points2.set_length(n);
	}
	int error;
	size_t k=order;//Order
	if(k>n) k=n;
	m_knot_vector=MGKnotVector(tau2,k);
	double* work=new double[n*(2*k-1)];
	const size_t irc=m_line_bcoef.capacity();
	error=blgint_(tau2.data(), points2.data(), m_knot_vector.data(), k, n, pointDim,
			pointSize, irc, work, &m_line_bcoef(0,0));

	delete[] work;
	if(error==1){
		error=0;		//Return of BLGINT error=1 means normal.
		m_line_bcoef.set_length(n);
	}
}

MGLBRep::MGLBRep(
	const MGNDDArray& tau,	//Data point abscissa
	const MGBPointSeq& points,	//Point seq data
	const MGKnotVector& t,	//knot vector
	int &error)				//Error flag.
// Construct Line B-rep of any order number by interpolation 
//from data point and knot vector.
:MGCurve(), m_knot_vector(t)
	,m_line_bcoef(points.length(),points.sdim())
{
	const unsigned k=t.order();
	const size_t n=points.length();
	const unsigned ncd=points.sdim();
	const size_t iv=points.capacity();
	const size_t irc=m_line_bcoef.capacity();
	double* work=new double[n*(2*k-1)];

	error=blgint_(tau.data(), points.data(), t.data(), k, n, ncd,
			iv, irc, work, &m_line_bcoef(0,0));

	delete[] work;
	if(error==1){
		error=0;		//Return of BLGINT error=1 means normal.
		m_line_bcoef.set_length(n);
	}
}

void build_mgblgsp_parameter(
	size_t& order,				//Order of the target B-Spline.
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& value,	//Data point ordinate
	MGENDCOND& beginc,
	MGENDCOND& endc,
	MGNDDArray& tau_new,	//Data point abscissa
	MGBPointSeq& value_new	//Data point ordinate.
							//value_new.size() wil be value.length().
){
	int i,j;
	int ie=1, is=0;
	const int n=value.length();
	const int nm1=n-1;
	int nnew=n;
	beginc=begin.cond(); endc=end.cond();
	if(beginc==MGENDC_1D || beginc==MGENDC_2D) {nnew+=1;is=1;}
	else if(beginc==MGENDC_12D) {nnew+=2; is=2;}
	if(endc==MGENDC_1D || endc==MGENDC_2D){nnew+=1; ie=2;}
	else if(endc==MGENDC_12D) {nnew+=2; ie=3;}
	if(int(order)>nnew) order=nnew;
	const int ncd=value.sdim();

	value_new=value;
	value_new.reshape(nnew,is);
	int nnewm1=nnew-1;
	for(j=0; j<ncd; j++) {
		value_new(0,j)=value(0,j);
		value_new(nnewm1,j)=value(nm1,j);
	}
	tau_new=tau;
	tau_new.reshape(nnew,is);
	for(int k=0; k<is; k++) tau_new(k)=tau(0);
	for(i=n+is;i<nnew; i++) tau_new(i)=tau(nm1);

	if(beginc==MGENDC_1D || beginc==MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(1,j)=(begin.first()).ref(j);
	if(beginc==MGENDC_2D || beginc==MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(is,j)=(begin.second()).ref(j);
	if(endc==MGENDC_1D || endc==MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(nnew-2,j)=(end.first()).ref(j);
	if(endc==MGENDC_2D || endc==MGENDC_12D)
		for(j=0; j<ncd; j++) value_new(nnew-ie,j)=(end.second()).ref(j);
	value_new.set_length(nnew); tau_new.set_length(nnew);
}

// Construct Line B-rep of order 4 by interpolation from Point data
//and end condition.
// Inner point may include derivative inf.
MGLBRep::MGLBRep(
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& value,	//Data point ordinate
	int &error)					//Error flag.
:MGCurve(){
	*this=MGLBRep(4,begin,end,tau,value,error);
}

// Construct Line B-rep of input order by interpolation from Point data
//and end condition.
// Inner point may include derivative inf.
MGLBRep::MGLBRep(
	size_t order,				//Order of the target B-Spline.
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& value,	//Data point ordinate
	int &error)					//Error flag.
:MGCurve(){
	MGENDCOND beginc, endc;
	MGBPointSeq valuen;	//Data point ordinate
	MGNDDArray taun;	//Data point abscissa
	build_mgblgsp_parameter(order,begin,end,tau,value,beginc,endc,taun,valuen);
	int nnew=valuen.length();
	int ncd=valuen.sdim();

	m_line_bcoef.resize(nnew,ncd);
	m_knot_vector.size_change(order,nnew);
	int irc=nnew;
	int iv=nnew;
	double* work=new double[nnew*(order*2+1)];
	mgblgsq((int)order,beginc,endc,taun.data(),valuen.data(),iv,nnew,ncd,
			irc,work,&m_knot_vector(0),&m_line_bcoef(0,0),&error);

	delete[] work;
	if(error==1){
		error=0;		//Return of mgblgsq error=1 means normal.
		m_knot_vector.set_bdim(nnew);
		m_line_bcoef.set_length(nnew);
	}
}

//Construct Line B-rep of any order by interpolation from Point data
//with end condition and the knot vector for the B-rep to construct.
// (tau(i), value(i,.)) for 0<=i<=n(the length of value).
//For the start and end point, tau does not have multiplicity. However,
//if tau has multiplicity at inner point, this means 1st derivative data is
//provided for the associated value, i.e.:
//If tau has multiplicity 2 as tau(i)=tau(i+1), value(i,.) is 1st derivative
//at tau(i) and value(i+1,.) is positional data at tau(i)(=tau(i+1)).
//If tau has multiplicity 3 as tau(i)=tau(i+1)=tau(i+2),
//value(i,.) is 1st derivative at tau(i)- ,
//value(i+1,.) is positional data at tau(i)(=tau(i+1)), 
//value(i+2,.) is 1st derivative at tau(i)+.
//Maximum multiplicity allowed is 3.
//Data points for the start and end point will be added according to begin and end,
//and tau does not have to have the data point multiplicity at the start and end.
MGLBRep::MGLBRep(
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& value,	//Data point ordinate
	const MGKnotVector& t,		//knot vector.
	int &error)				//Error flag.
:MGCurve(), m_knot_vector(t){
	MGENDCOND beginc, endc;
	MGBPointSeq valuen;	//Data point ordinate
	MGNDDArray taun;	//Data point abscissa
	size_t k=t.order();
	build_mgblgsp_parameter(k,begin,end,tau,value,beginc,endc,taun,valuen);
	int nnew=valuen.length();
	int ncd=valuen.sdim();

	m_line_bcoef.resize(nnew,ncd);
	double* work=new double[nnew*(k*2+1)];
	error = 2;
	for(int i=0; i<ncd; ++i){
		blg4sp2_((int)k,&error,beginc,endc,taun.data(),valuen.data()+i*nnew,nnew,nnew,1,
			m_knot_vector.data(),1,work,work+nnew,work+2*nnew,&m_line_bcoef(0,i));
			if (error!=1) break;
	}

	delete[] work;
	if(error==1) error=0;//Return of mgblgsq error=1 means normal.
}

// Construct Line B-rep of order 4 from point and point-kind followed by
//osculating circle data.
MGLBRep::MGLBRep(
	const MGLBRepEndC& begin,	//Begin end condition
	const MGLBRepEndC& end,		//End end conditoion
	const MGBPointSeq& points,	//Point seq data
	const int*  point_kind,		// Point kind of above point.
	const MGOscuCircle& circle,	//Provides osculating circle data.
	int &error)					//Error flag.
:MGCurve(), m_knot_vector(4)
	,m_line_bcoef(0,points.sdim())
{
	assert(points.sdim()==2 || points.sdim()==3);
	size_t j;

	MGENDCOND ibc[2]; ibc[0]=begin.cond(), ibc[1]=end.cond();
	if(ibc[0]==MGENDC_12D) ibc[0]=MGENDC_1D;
	if(ibc[1]==MGENDC_12D) ibc[1]=MGENDC_1D;
	const unsigned ncd=points.sdim();
	const size_t nv=points.length();

	double begin_deriv[3], end_deriv[3];
	if(ibc[0]==MGENDC_1D)
		for(j=0; j<ncd; j++) begin_deriv[j]=(begin.first()).ref(j);
	else if(ibc[0]==MGENDC_2D)
		for(j=0; j<ncd; j++) begin_deriv[j]=(begin.second()).ref(j);
	if(ibc[1]==MGENDC_1D)
		for(j=0; j<ncd; j++) end_deriv[j]=(end.first()).ref(j);
	else if(ibc[1]==MGENDC_2D)
		for(j=0; j<ncd; j++) end_deriv[j]=(end.second()).ref(j);

	size_t ncir=circle.length();	//Construct index and radius data
	int* idk=new int[ncir];	//from circle input.
	double* rcir=new double[ncir];
	for(size_t i=0; i<ncir; i++){
		idk[i]=circle(i).index()+1;
		rcir[i]=circle(i).radius();
	}

	int maxn=nv+2+7*ncir+nv;
	m_line_bcoef.reshape(maxn);
	m_knot_vector.reshape(maxn+4);

	const size_t iv=points.capacity();
	const size_t irc=maxn; 
	double* work1=new double[maxn];
	int m1=5*irc+105; int m2=9*maxn;
	int m=(m1>m2)?m1:m2;
	double* work2=new double[m];

	int n;
	blgcs_((const int*)ibc, begin_deriv, end_deriv, ncd, nv, point_kind,
			points.data(), ncir, idk, rcir, iv, irc, work1,work2,
			&n, &m_knot_vector(0), &m_line_bcoef(0,0), &error);
							   
	delete[] work1; delete[] work2; delete[] idk; delete[] rcir;
	if(error==1){
		error=0;		//Return of BLGCS error=1 means normal.
		m_knot_vector.set_bdim(n);
		m_knot_vector.reshape(n+4);
		m_line_bcoef.set_length(n);
		m_line_bcoef.reshape(n);
	}
}

//Construct MGLBRep that interpolates data points (tau(i), rlb.eval(tau(i)))
// for 0<=i<tau.length().
MGLBRep::MGLBRep(const MGRLBRep& rlb, const MGNDDArray& tau)
:MGCurve(rlb){
	int error;
	size_t n=tau.length();
	size_t sd=rlb.sdim();
	MGBPointSeq point(n,sd);
	for(size_t i=0; i<n; i++) point.store_at(i,rlb.eval(tau(i)));
	MGLBRepEndC sderi(MGENDC_1D, rlb.eval(tau(0),1));
	MGLBRepEndC ederi(MGENDC_1D, rlb.eval(tau(n-1),1));
	*this=MGLBRep(sderi,ederi,tau,point,error);
}

//**** 2. Approximation Constructor ****

//Construct Curve B-Rep ep.
//This is an approximation, and the tolerance is MGTolerance::line_zero().
MGLBRep::MGLBRep(
	const MGCurve& crv,	//Original Curve.
	size_t order)//Order. When order=0 is input, and crv was a MGLBRep,
		//the original order will be used. Otherwise(order=0 and crv was not an MGLBRep)
		//order will be set to 4.
:MGCurve(crv){
	const MGLBRep *brep = dynamic_cast<const MGLBRep*>(&crv);
	if(brep)
		if(!order || brep->order() == order){
			*this = *brep;
			return;
		}
	const MGTrimmedCurve *tcrv = dynamic_cast<const MGTrimmedCurve*>(&crv);
	if(tcrv){
		const MGLBRep *bcrv = dynamic_cast<const MGLBRep*>(tcrv->curve());
		if(bcrv) if(!order || bcrv->order()==order){
			if(tcrv->m_sameRange){*this = *bcrv; return;}
			*this=MGLBRep(tcrv->m_range.low_point(), tcrv->m_range.high_point(), *bcrv);
			return;
		}
	}

	update_mark();
	size_t knew=order;
	if(knew==0)
		knew=4;
	const MGSurfCurve *scrv = dynamic_cast<const MGSurfCurve*>(&crv);
	if(scrv){
		const MGLBRep *bcrv = dynamic_cast<const MGLBRep*>(scrv->m_curve.curve());
		if(bcrv&&!order) knew=bcrv->order();
	}

	double ts=crv.param_s(), te=crv.param_e();
	const MGKnotVector& t=crv.knot_vector();
	size_t k=t.locate(ts)+1, n=t.locate(te)+1;
		//k!=t.order() or n!=t.bdim() of t occurs when crv is a trimmedcurve.
	size_t index, multi_found;
	size_t km1=t.order()-1, start=k;
	do{	//Approximation by dividing to parts of continuity>=C0.
		multi_found=t.locate_multi(start,km1,index);//Locate C0 continuity point.
		if(start==k && index==n){
			crv.approximate_as_LBRep(*this,knew,start-1,index);	//Use original.
			if(ts>param_s() || te<param_e()){
				MGLBRep lb(ts, te, *this); *this=lb;
			}
			return;
		}else{
			if(start==k){//For the 1st span.
				crv.approximate_as_LBRep(*this,knew,start-1,index);//1st approximation.
				if(ts>param_s()){
					MGLBRep lb(ts, param_e(), *this); *this=lb;
				}
			}else{//For the span from the 2nd.
				MGLBRep lbt;
				crv.approximate_as_LBRep(lbt,knew,start-1,index);//from the 2nd approximation.
				if(te<lbt.param_e()){
					MGLBRep lbt2(lbt.param_s(), te, lbt); lbt=lbt2;
				}
				//double ratio; int which;
				//int cont=continuity(lbt,which,ratio);
				//cout<<(*this)<<lbt;
				connect(0,2,lbt);
				//cout<<(*this);
			}
			start=index+multi_found;
		}
	}while(index<n);
}

// Construct 3D B-Rep by mixing two 2D B-Rep.
//The two 2D B-Rep's directions and start and end points must be the same.
MGLBRep::MGLBRep(
	unsigned coordinate1,	//Missing oordinate kind of the brep1 
							// 0:x, 1:y, 2:z.
	const MGLBRep& brep1,	//Original 2D B-Rep1. Coordinates are
							//(y,z), (z,x), (x,y) according to coordinate1. 
	unsigned coordinate2,	//Missing coordinate kind of the brep2.
							// 0:x, 1:y, 2:z, and 3:girth rep.
	const MGLBRep& brep2)	//Original 2D B-Rep2. Coordinates are
		//(y,z), (z,x), (x,y) and (t, g2) according to coordinate2.
		//t is parameter of brep1 and g2 is x, y, or z according to coordinate1.
//Second 2D B-Rep can be girth representaion. That is,
//Let brep1 is f(t)=(f1(t),f2(t)) and brep2 is g(s)=(g1(s),g2(s)), where
//f1,f2 are two coordinates, g1 is parameter t of f(t), and 
//g2(s) is the missing coordinate of f(t). Given s value,
// ( f1(g1(s)), f2(g1(s)), g2(s) ) is a 3D space point.
:MGCurve(brep1){
	update_mark();
	assert(brep1.sdim()==2 && brep2.sdim()==2);
	assert(coordinate1<=2 && coordinate2<=3 && coordinate1!=coordinate2);

	int i;
	size_t kcod1=coordinate1+1;
	const unsigned k1=brep1.order();
	const size_t n1=brep1.m_line_bcoef.length();
	const size_t irc1=brep1.m_line_bcoef.capacity();

	size_t kcod2=coordinate2+1;
	const unsigned k2=brep2.order();
	const size_t n2=brep2.m_line_bcoef.length();
	const size_t irc2=brep2.m_line_bcoef.capacity();

	size_t k=k1; if(k<k2) k=k2; k=4*k*k+3*k;
	size_t irc=n1+n2; if(irc<11) irc=11;
	if(irc*4<k) irc=int(k/4)+1;

	double error=MGTolerance::line_zero();
	int* kseq=new int[irc];
	double* wk2=new double[irc*4];
	int n; MGNDDArray tau(irc); MGBPointSeq bp(irc,3);
	int iflag; MGVector vec_s(3),vec_e(3);

	//Compute initial data by BLUMIX.
	blumix_(error,kcod1,k1,n1,brep1.knot_data(),brep1.coef_data(),irc1,
		kcod2,k2,n2,brep2.knot_data(),brep2.coef_data(),irc2,
		irc,kseq,wk2,&vec_s(0),&vec_e(0),&n,&tau(0),&bp(0,0),&iflag);
	tau.set_length(n); bp.set_length(n);
	delete[] kseq;
	delete[] wk2;
	MGLBRepEndC endc_s(MGENDC_1D,vec_s);
	MGLBRepEndC endc_e(MGENDC_1D,vec_e);

	//Check tolerance.
	size_t cod11,cod12,cod21;
	cod11=kcod1;if(cod11>=3) cod11=0;
	cod12=cod11+1;if(cod12>=3) cod12=0;
	unsigned com21=1;
	if(coordinate2<3){
		cod21=kcod2;if(cod21>=3) cod21=0;
		if(cod21==coordinate1)com21=0;
	}
	//cod11,12,21 are 3D coordinate id of brep1 and 2.
	//com21 is id of 2nd line coordinate that is missing in 1st line.

	MGVector F,G(3); MGPosition b1,b2,p1,p2; 
	int added;
	MGLBRep line;
	double tmid,s1,s2,diff;
	for(size_t loop=0; loop<3; loop++){
		//Maximum loop counter is 3(3 times as much as output of BLUMIX)
		added=0;
		n=bp.length();
		line=MGLBRep(endc_s,endc_e,tau,bp,iflag); if(iflag) return;
		for(i=1; i<n-1; i++){
			while((tau(i)==tau(i+1)) && i<n-1) i+=1;
			tmid=(tau(i)+tau(i+1))/2.;
			F=line.eval(tmid);
			b1=MGPosition(2,F,0,cod11);
			brep1.on(b1,s1); p1=brep1.eval(s1);
			diff=(p1-b1).len();
			if(diff>error){			//brep1 differs from line, add point.
				G=brep2.closest_mix(coordinate2,coordinate1,s1,p1,F);
				bp.insert_at(i+1,G); tau.add_data(tmid);
				i+=1; n+=1; added+=1;
			} else{
				if(coordinate2<3) b2=MGPosition(2,F,0,cod21);
				else			  b2=MGPosition(s1,F(coordinate1));
				brep2.on(b2,s2); p2=brep2.eval(s2);
				diff=(p2-b2).len();
				if(diff>error){		//brep2 differs from line, add point.
					G.set(cod11)=p1.ref(0); G.set(cod12)=p1.ref(1);
					G.set(coordinate1)=p2.ref(com21);
					bp.insert_at(i+1,G); tau.add_data(tmid);
					i+=1; n+=1; added+=1;
				}
			}
		}
		if(!added) break;
	}
	if(added){
		n=bp.length();
		line=MGLBRep(endc_s,endc_e,tau,bp,iflag);
		if(iflag) return;
	}
	(*this)=line;
}

//Function for BLUMIX constructor. Given 3D point F,
//compute correct point of F that is closest to F.
MGPosition MGLBRep::closest_mix(
	unsigned coordinate2,	//Missing coordinate kind of this(say 2nd line).
							//coordinate2 can be 3.
	unsigned coordinate1,	//Missing coordinate kind of P(say, of 1st line).
							//coordinate1 cannot be 3.
	double tau,				//Parameter value of P, of 1st line,
							//used only when coordinate2=3.
	const MGPosition& P,	//Point of 1st line(correct coordinates).
	const MGPosition& F		//3D point, used to coose closest point to F
							//when more than one point are found.
	) const
{
	MGCParam_list list; unsigned kcod;

	unsigned com2=0, com21=1;
	if(coordinate2<3){
		kcod=coordinate2+1; if(kcod>=3) kcod=0;
		if(kcod==coordinate1){com2=1; com21=0;}
	}
	unsigned com1=com2+1; if(com1>=2) com1=0;
	//com2 is id of 2nd line coordinate that is common to P(1st line).
	//com21 is id of 2nd line coordinate that is missing in P(1st line).
	//com1 is id of P coordinate that is common to 2nd line.

	if(coordinate2==3) list=isect_1D(tau);
	else{
		double data=P.ref(com1); list=isect_1D(data,com2);
		double tol=MGTolerance::wc_zero();
		if(!list.entries()) list=isect_1D(data+tol,com2); 
		if(!list.entries()) list=isect_1D(data-tol,com2);
	}

	kcod=coordinate1+1; if(kcod>=3) kcod=0;
	MGPosition Q(3,P,kcod);
	size_t n=list.entries();
	if(!n){
		Q(coordinate1)=F(coordinate1);
		return Q;
	}
	MGPosition A=eval(list.removeFirst()); Q(coordinate1)=A(com21);
	if(n==1) return Q;

	double dist=(Q-F).len(), dist2;
	MGPosition Q2(Q);
	for(size_t i=1; i<n; i++){
		A=eval(list.removeFirst()); Q2(coordinate1)=A(com21);
		dist2=(Q2-F).len();
		if(dist2<dist){ dist=dist2; Q=Q2;}
	}
	return Q;
}	

// Construct Line B-rep of any order number by least square approximation
//from Point data with approximation weights and knot vector of B-Rep.
MGLBRep::MGLBRep(
	const MGNDDArray& tau,	//Data point abscissa.
	const MGBPointSeq& points,	//Data Point ordinates.
	const double* weight,	//Weights for each points 
	const MGKnotVector& t)	//knot vector
	:MGCurve(), m_knot_vector(t)
	,m_line_bcoef(t.bdim(),points.sdim())
{
	const size_t ntau=tau.length();
	const size_t ig=points.capacity();
	const size_t iw=1; const size_t irc=1; const int kw=1;
	const size_t n=t.bdim();
	unsigned k=t.order();
	const unsigned sdim=points.sdim();
	const unsigned m=1;
	double* work=new double[2*n];
	double* q=new double[k*n];

	for(size_t i=0; i<sdim; i++)
		blgl2a_(ntau,tau.data(),points.data(0,i),weight,ig,iw,irc,
				kw,t.data(),n,k,m, work, q, &m_line_bcoef(0,i));
	m_line_bcoef.set_length(n);

	delete[] work; delete[] q;
}

//Approximate an original B-Rep by a new knot configuration.
//The new knot config must be inside the range of the original B-Rep
//parameter. However new knots may be coarse or fine.
MGLBRep::MGLBRep(
	const MGLBRep& old_brep,//Original B-Rep.
	const MGKnotVector& t,	//knot vector
	int &error)				//Error flag.
:MGCurve(old_brep), m_knot_vector(t)
	,m_line_bcoef(t.bdim(),old_brep.sdim())
{
	update_mark();
	const unsigned k=t.order();
	const size_t n1=old_brep.m_line_bcoef.length();
	const size_t irc1=old_brep.m_line_bcoef.capacity();
	const unsigned ncd=m_line_bcoef.sdim();
	const size_t n2=t.bdim();
	const size_t irc2=n2;
	size_t m;if(k<=5) m=9;else m=2*k-1;
	double* work1=new double[n2];
	double* work2=new double[n2*m];

	bludkt_(k,n1,old_brep.knot_data(),old_brep.coef_data(),
			irc1,ncd,n2, m_knot_vector.data(),irc2,
			work1, work2, &m_line_bcoef(0,0), &error);
							   
	delete[] work1; delete[] work2;
	if(error==1){
		error=0;		//Return of BLUDKT:error=1 means normal.
		m_knot_vector.set_bdim(n2);
		m_line_bcoef.set_length(n2);
	}
}

//Gets new B-Rep by a new knots.
MGLBRep::MGLBRep(
	const MGCurve& old_curve,	//Original curve.
	const MGKnotVector& t	//knot vector
):MGCurve(old_curve),m_knot_vector(t){
	update_mark();
	MGNDDArray tau; tau.update_from_knot(t);
	MGBPointSeq bp;
	old_curve.eval_line(tau,bp);
	int error;
	m_line_bcoef=MGLBRep(tau,bp,t,error).line_bcoef();
}

//**** 3.Conversion Constructor.****

MGLBRep::MGLBRep(
	const MGPPRep& pprep)	//PP-rep
//Convert PP-Rep to B-rep.
:MGCurve(){
	size_t i,j;

	const unsigned ncd=pprep.sdim();
	const unsigned k=pprep.order();	//order.
	const size_t l=pprep.nbreak()-1;//Number of spans.
	size_t n=l+k-1;					//Initial B-Rep dimension.
	MGPPRep pp(pprep); pp.normalize();
	//Generate initial knot vector.
	MGKnotVector t(k,n);
	for(i=0; i<k; i++) {
		t(i)=pp.break_point(0);
		t(i+n)=pp.break_point(l);
	}
	j=1; i=k;
	while(i<n) t(i++)=pp.break_point(j++);
	//Add multiple knot, taking continuity into account.
	MGVector v1,v2; size_t imk,imkp1; double brk;
	for(i=k; i<n; i++){
		imk=i-k; imkp1=imk+1;
		brk=pp.break_point(imkp1);
		v1=pp.eval_i(imk,brk,0);
		v2=pp.coef(0,imkp1);
		j=0;
		if(v1==v2){
			for(j=1; j<=k-2; j++){
				v1=pp.eval_i(imk,brk,j);
				v2=pp.coef(j,imkp1);
				if(v1!=v2) break;
			}
		}
		int mult=k-1-j;	 //// mult=k-j
		if(mult>0)		 //// mult>1
			t=MGKnotVector(t,MGKnotArray(brk,mult));
	}

	*this=MGLBRep(pp,t);
}

MGLBRep::MGLBRep(
	const MGLBRep& old_brep,	//Original B-Rep.
	const MGKnotArray& knots)	//Knots to add.
//Gets new B-Rep by adding knots to an original B-Rep.
:MGCurve(old_brep)
	,m_knot_vector(old_brep.order(),old_brep.bdim())
	,m_line_bcoef(old_brep.bdim(),old_brep.sdim())
{
	update_mark();
	const unsigned k=old_brep.order();
	const size_t n1=old_brep.bdim();
	const size_t irc1=old_brep.m_line_bcoef.capacity();
	const unsigned ncd=m_line_bcoef.sdim();

	const size_t nad=knots.length();
	double* tad=new double[nad];
	int* mlt=new int[nad];
	unsigned mlt_total=0;
	const MGKnotVector& t=old_brep.knot_vector();
	size_t km1=k-1;
	for(size_t i=0; i<nad; i++){
		double tadi=tad[i]=knots[i].value();
		int j=t.locate(tadi);
		size_t mlt1=0,mlt2=knots[i].multiplicity();
		if(tadi==t[j]){
			size_t j1=j;
			while(j1>=km1 && t[j1-1]==tadi){j1--; mlt1++;};
		}
		mlt2+=mlt1;
		if(mlt2>=k) mlt2=km1;
		mlt[i]=mlt2;
		mlt_total+=mlt[i];
	}
	const size_t irc2=n1+mlt_total;
	m_knot_vector.reshape(irc2+k);
	m_line_bcoef.reshape(irc2);
	double* work=new double[k*k];

	int n2;
	bluakt_( k, n1, old_brep.knot_data(), old_brep.coef_data(),
			irc1, ncd, nad, tad, mlt, irc2, work,
			&n2, &m_knot_vector(0), &m_line_bcoef(0,0));
	m_knot_vector.set_bdim(n2);
	m_line_bcoef.set_length(n2);

	delete[] tad; delete[] mlt, delete[] work;
}

//Construct LBRep by connecting brep1 and brep2 to make one B-Representation.
//brep1's parameter range will not be changed, instead brep2's range
//will be so modified that brep2 has the same 1st derivative magnitude
//as the original brep1's at the connecting point
//(start or end point of brep1).
//continuity and which can be obtained using the fucntion continuity().
MGLBRep::MGLBRep(
	const MGLBRep& brep1,	//B-Rep 1.
	int continuity,			//continuity.
	int which,				//which point of brep1 to which of brep2.
	const MGLBRep& brep2)	//B-Rep 2.
:MGCurve(brep1){
	assert(0<=which && which<=3);
	*this=brep1;
	connect(continuity,which,brep2);
}

//Gets new B-Rep by computing a part of the original. New one is exactly the same
//as the original except that it is partial.
MGLBRep::MGLBRep(
	double t1, double t2,	//New parameter range. t1 must be less than t2.
	const MGLBRep& old_brep,//Original B-Rep.
	int multiple)   //Indicates if start and end knot multiplicities
					//are necessary. =0:unnecessary, !=0:necessary.
//If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1 will be
//guaranteed.
//Both t1 and t2 must be inside te range of old_brep.
:MGCurve(old_brep)
	,m_knot_vector(old_brep.order(),old_brep.bdim())
	,m_line_bcoef(old_brep.bdim(),old_brep.sdim())
{
	update_mark();
	int i1=old_brep.m_knot_vector.locate(t1), i2=old_brep.m_knot_vector.locate(t2,1);
	double error=old_brep.m_knot_vector.param_error();
	double ti1=old_brep.m_knot_vector[i1];
	if((t1-ti1)<=error) t1=ti1;
//	else{
//		double ti1p1=old_brep.m_knot_vector[i1+1];
//		if((ti1p1-t1)<=error) t1=ti1p1;
//	}
//	double ti2=old_brep.m_knot_vector[i2];
//	if((t2-ti2)<=error) t2=ti2;
//	else{
		double ti2p1=old_brep.m_knot_vector[i2+1];
		if((ti2p1-t2)<=error) t2=ti2p1;
//	}
	const unsigned k=old_brep.order();
	const size_t n1=old_brep.bdim();
	const size_t irc1=old_brep.m_line_bcoef.capacity();
	const size_t ncd=old_brep.sdim();
	const size_t irc2=old_brep.bdim();
	double* work=new double[k*k];
	int n2;

	bluprt_(k, n1, old_brep.knot_data(), old_brep.coef_data(),
			irc1, ncd, t1, t2, irc2, work, &n2, &knot(0), &coef(0,0)
			,multiple);
	m_knot_vector.set_bdim(n2);
	m_line_bcoef.set_length(n2);

	delete[] work;//cout<<(*this);///////////
}

// Construct a Line B-Rep by changing space dimension and order of coordinate.
MGLBRep::MGLBRep(
		size_t dim,					// New space dimension.
		const MGLBRep& old_brep,	// Original Line B-rep.
		size_t start1, 				// Destination order of new line.
		size_t start2)	 			// Source order of original line.
:MGCurve(old_brep)
	,m_knot_vector(old_brep.knot_vector())
	,m_line_bcoef(dim,old_brep.m_line_bcoef,start1,start2)
{
	update_mark();
}

//	MGLBRep(const MGLBRep&);  //Copy constructor.
//  We can use default copy constructor.

//Destructor
//	~MGLBRep();	//We can use default destructor.

//This constructor constructs B-Rep, converting from PP-Rep.
//Knot Vector is input. Each knot of the knot vector is break point of
//pprep. The continuities at all the break points must be C(k-2) where
//k is the order of pprep.
MGLBRep::MGLBRep(
	const MGPPRep& pp,	//PP-rep
	const MGKnotVector t)	//Knot Vector
:MGCurve(),m_knot_vector(t)
,m_line_bcoef(t.bdim(),pp.sdim()){
	assert(t.order()==pp.order());

	const unsigned ncd=pp.sdim();
	const unsigned k=pp.order();	//order.
	const size_t ipc=pp.nbreak();
	const size_t l=ipc-1;//Number of spans.
	size_t n=t.bdim();
	const size_t ipcw=n+k-2;
	double* pcwork=new double[ipcw*k*ncd];
	const size_t irc=m_line_bcoef.capacity();
	blctpb_(pp.break_point_data(),pp.coef_data(),l,k,ncd,n,
			t.data(),ipc,ipcw,irc, pcwork, &m_line_bcoef(0,0));
	delete[] pcwork;
}
