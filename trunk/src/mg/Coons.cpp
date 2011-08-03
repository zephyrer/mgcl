/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/LBRep.h"
#include "mg/Coons.h"
#include "mg/SPointSeq.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Constructor.
MGCoons::MGCoons(
	MGPvector<MGCurve>& perimeters,
	MGPvector<MGCurve>& derivatives
):m_perimeters(perimeters), m_derivatives(derivatives){
	eval_corner();
}

MGCoons::MGCoons(
	MGPvector<MGLBRep>& perimeters,
	MGPvector<MGLBRep>& derivatives
){
	m_perimeters.resize(4);
	m_derivatives.resize(4);
	for(size_t i=0; i<4; i++){
		m_perimeters.reset(i,perimeters.release(i));
		m_derivatives.reset(i,derivatives.release(i));
	}
	eval_corner();
}

//get the derivative data d2f/((dv)(du)) at corner i.
const MGVector& MGCoons::d2fdvu(size_t i)const{
	switch(i){
		case 0: return m_df2dvu00;
		case 1: return m_df2dvu10;
		case 2: return m_df2dvu11;
		default: return m_df2dvu01;
	}
}

void MGCoons::eval_corner(){
	m_f00=(m_perimeters[0]->eval(0.)+m_perimeters[3]->eval(0.))*.5;
	m_f01=(m_perimeters[3]->eval(1.)+m_perimeters[2]->eval(0.))*.5;
	m_f10=(m_perimeters[0]->eval(1.)+m_perimeters[1]->eval(0.))*.5;
	m_f11=(m_perimeters[2]->eval(1.)+m_perimeters[1]->eval(1.))*.5;

	m_dfdu00=m_derivatives[3]->eval(0.);
	m_dfdu01=m_derivatives[3]->eval(1.);
	m_dfdu10=m_derivatives[1]->eval(0.);
	m_dfdu11=m_derivatives[1]->eval(1.);

	m_dfdv00=m_derivatives[0]->eval(0.);
	m_dfdv01=m_derivatives[2]->eval(0.);
	m_dfdv10=m_derivatives[0]->eval(1.);
	m_dfdv11=m_derivatives[2]->eval(1.);

	m_df2duv00=m_derivatives[3]->eval(0.,1);
	m_df2dvu00=m_derivatives[0]->eval(0.,1);
	m_df2duv01=m_derivatives[3]->eval(1.,1);
	m_df2dvu01=m_derivatives[2]->eval(0.,1);

	m_df2duv10=m_derivatives[1]->eval(0.,1);
	m_df2dvu10=m_derivatives[0]->eval(1.,1);
	m_df2duv11=m_derivatives[1]->eval(1.,1);
	m_df2dvu11=m_derivatives[2]->eval(1.,1);
}

//Copy constructor.
//MGCoons(const MGCoons& rhs);

//////////Destructor////////

//~MGCoons();
									
////////// Operator Overload //////////

//MGCoons& operator=(const MGCoons&); //Assignment

//Friend Function

//Debug Function
std::ostream& operator<< (std::ostream& out, const MGCoons& coons){
	out<<"MGCoons::"<<&coons<<", m_perimeters[i] and m_derivatives[i]="<<std::endl;
	for(size_t i=0; i<3; i++){
		out<<i<<":: "<<coons.m_perimeters[i]<<coons.m_derivatives[i];
	}
	return out;
}

////////// Member Function //////////

MGVector MGCoons::eval_inner(
	double u, double v	// Parameter value of the surface.
						// must be 0<=u,v<=1.
)const{
	if(MGMZero(v)) return m_perimeters[0]->eval(u);
	if(MGMZero(u)) return m_perimeters[3]->eval(v);
	if(MGMZero(1.-u))return m_perimeters[1]->eval(v);
	if(MGMZero(1.-v)) return m_perimeters[2]->eval(u);

	double h0u=mgHermite0(u),h1u=mgHermite1(u),h2u=mgHermite2(u),h3u=mgHermite3(u);
	double h0v=mgHermite0(v),h1v=mgHermite1(v),h2v=mgHermite2(v),h3v=mgHermite3(v);

	MGVector f0v=m_perimeters[3]->eval(v);
	MGVector f1v=m_perimeters[1]->eval(v);
	MGVector dfdu0v=m_derivatives[3]->eval(v);
	MGVector dfdu1v=m_derivatives[1]->eval(v);
	MGVector fa=(h0u*f0v+h2u*dfdu0v)+(h1u*f1v+h3u*dfdu1v);

	MGVector fu0=m_perimeters[0]->eval(u);
	MGVector fu1=m_perimeters[2]->eval(u);
	MGVector dfdvu0=m_derivatives[0]->eval(u);
	MGVector dfdvu1=m_derivatives[2]->eval(u);
	MGVector fb=(h0v*fu0+h2v*dfdvu0)+(h1v*fu1+h3v*dfdvu1);

	MGVector df2dudv00=(u*m_df2duv00+v*m_df2dvu00)/(u+v);
	MGVector df2dudv01=(-u*m_df2duv01+(v-1.)*m_df2dvu01)/(v-u-1.);
	MGVector df2dudv10=((1.-u)*m_df2duv10+v*m_df2dvu10)/(1.-u+v);
	MGVector df2dudv11=((u-1.)*m_df2duv11+(v-1.)*m_df2dvu11)/(u+v-2.);
	MGVector m0=(h0u*m_f00 +h2u*m_dfdu00)  +(h1u*m_f10 +h3u*m_dfdu10);
	MGVector m1=(h0u*m_f01 +h2u*m_dfdu01)  +(h1u*m_f11 +h3u*m_dfdu11);
	MGVector m2=(h0u*m_dfdv00+h2u*df2dudv00)+(h1u*m_dfdv10+h3u*df2dudv10);
	MGVector m3=(h0u*m_dfdv01+h2u*df2dudv01)+(h1u*m_dfdv11+h3u*df2dudv11);

	MGVector fc=(m0*h0v+m2*h2v)+(m1*h1v+m3*h3v);

	return fa+(fb-fc);
}

//Evaluate surface data.
//Currently ndu=ndv=0 is assumed.
MGVector MGCoons::eval(
	double u, double v	// Parameter value of the surface.
						// must be 0<=u,v<=1.
	, size_t ndu		// Order of derivative along u.
	, size_t ndv		// Order of derivative along v.
)const{
	return eval_inner(u,v);
}

void MGCoons::eval(
	const MGNDDArray&	utau,		//u方向のデータポイント
	const MGNDDArray&	vtau,		//v方向のデータポイント
	MGSPointSeq&		spoint)const//evaluated data will be output to spoint.
{
	size_t nu=utau.length(), nv=vtau.length();
	for(size_t i=0; i<nu; i++){
		double u=utau[i];
		for(size_t j=0; j<nv; j++){
			spoint.store_at(i,j,eval(u,vtau[j]));
		}
	}
}

size_t MGCoons::sdim()const{return m_perimeters[0]->sdim();}

double mgHermite0(double t){
	double s=(1.-t);
	return s*s*(2.*t+1.);
}
double mgHermite1(double t){
	return t*t*(-2.*t+3.);
}
double mgHermite2(double t){
	double s=(1.-t);
	return s*s*t;
}
double mgHermite3(double t){
	return t*t*(t-1.);
}
