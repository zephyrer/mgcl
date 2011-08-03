/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/SBRepVecTP.h"
#include "mg/LBRep.h"
#include "mg/Surface.h"
#include "mg/SBRep.h"
#include "mg/Interval.h"
#include "mg/Position.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implements Tangent Plane Line B-Representation.
// A vector version of MGSBRepTP

////////Constructor////////

//Copy Constructor
//all of the ownership of m_TP of tp2 will be
//transfered to this object.
MGSBRepVecTP::MGSBRepVecTP(const MGSBRepVecTP& vectp2){
	MGSBRepVecTP* rhsp=const_cast<MGSBRepVecTP*>(&vectp2);
	for(size_t i=0; i<4; i++){
		m_TP[i].push_back(rhsp->m_TP[i]);
		size_t i2=i*2;
		m_to_SE[i2]=rhsp->m_to_SE[i2];
		m_to_SE[i2+1]=rhsp->m_to_SE[i2+1];
		m_prange[i]=rhsp->m_prange[i];
	}
}

MGSBRepVecTP::MGSBRepVecTP(const MGSBRepTP& tp2){
	for(size_t i=0; i<4; i++){
		if(!(tp2.specified(i))) continue;
		m_TP[i].push_back(new MGLBRep(tp2.TP(i)));
		m_prange[i]=MGInterval(tp2.TP(i).param_range());
		m_to_SE[i*2]=m_to_SE[i*2+1]=true;
	}
}

//Assignment.
//all of the ownership of m_TP of tp2 will be
//transfered to this object.
MGSBRepVecTP& MGSBRepVecTP::operator=(const MGSBRepVecTP& vectp2){
	MGSBRepVecTP* rhsp=const_cast<MGSBRepVecTP*>(&vectp2);
	for(size_t i=0; i<4; i++){
		m_TP[i].clear();
		if(vectp2.specified(i)) m_TP[i].push_back(rhsp->m_TP[i]);
		size_t i2=i*2;
		m_to_SE[i2]=rhsp->m_to_SE[i2];
		m_to_SE[i2+1]=rhsp->m_to_SE[i2+1];
		m_prange[i]=rhsp->m_prange[i];
	}
	return *this;
}

//////////////// Member Function ///////////////////////

//change the parameter range to [t0,t1] if t0<t1,
//                           to [t1,t0] if t1<t0.
//when t1<t0, the direction will be reversed.
void MGSBRepVecTP::change_range(
	bool along_u,	//the objective range is along u parameter or v.
	double t0, double t1
){
	double lennew=t1-t0;
	size_t j=0; if(!along_u) j=1;
	for(size_t m=0; m<2; m++, j+=2){
		double lenold=m_prange[j].length().value();
		double ostart=m_prange[j][0].value();
		double ratio=lennew/lenold;
		MGPvector<MGLBRep>& tpvecm=m_TP[j];
		size_t i,n=tpvecm.size();
		for(i=0; i<n; i++){
			MGLBRep& lb=*(tpvecm[i]);
			double t20=t0+(lb.param_s()-ostart)*ratio;
			double t21=t0+(lb.param_e()-ostart)*ratio;
			lb.change_range(t20,t21);
		}
		if(t0>t1) tpvecm.reverse_sequence();
	}
	m_prange[j].change_range(t0,t1);
}

void MGSBRepVecTP::change_range(
	bool along_u,	//the objective range is along u parameter or v.
	const MGInterval& prange
){
	change_range(along_u,prange[0].value(), prange[1].value());
}

//evaluate TP at the perimeter i's parameter t.
//Function's return value is:
//true if t was inside the parameter range of a tangent plane of m_TP[i].
//false if t was outside the parameter range of the m_TP[i] for all i.
bool MGSBRepVecTP::eval(
	size_t i,		//perimeter numeber.
	double t,		//parameter vaule of perimeter i.
	MGVector& normal//evaluated normal will be returned.
)const{
	if(!(m_prange[i].includes(t))) return false;
	size_t n=m_TP[i].size();
	for(size_t j=0; j<n; j++){
		MGLBRep& tp=*(m_TP[i][j]);
		if(tp.in_range(t)){
			normal=tp.eval(t);
			return true;
		}
	}
	return false;
}

// Compute the maximum (absolute) cos value of between vector deris[i](t) 
// and vector this->TP(i)(t) for i=0,1,2,3, where t is a common
// parameter of the data point obtained from deris[i]'s knot vector.
//Function's return value is the max out of cosmax[.].
double MGSBRepVecTP::get_perimeters_max_cos(
	const MGPvector<MGLBRep>& deris,
	double taumax[4],
	double cosmax[4]
)const{
	assert(deris.size() == 4);
	MGVector N(3), T(3);
	MGNDDArray tau;
	double max=0.;
	for(int i = 0; i < 4; i++){
		if(!specified(i)){
			taumax[i] = deris[i]->param_s();
			cosmax[i] = 0.;
			continue;
		}

		tau.update_from_knot(deris[i]->knot_vector());
		double taus=tau[0];
		double cmi=0.;
		double tmi = taus;
		if(eval(i,taus,T)){
			N = deris[i]->eval(taus);
			cmi = fabs(N.cangle(T));
			tmi = taus;
		}

		int ntau=tau.length();
		for(int j = 1; j < ntau; j++){
			double tauj=tau[j];
			double tmid = (tau[j-1]+tauj)*.5;
			double cm;
			if(eval(i,tmid,T)){
				N = deris[i]->eval(tmid);
				cm = fabs(N.cangle(T));
				if(cmi<cm){	cmi = cm; tmi = tmid;}
			}

			if(!eval(i,tauj,T)) continue;
			N = deris[i]->eval(tauj);
			cm = fabs(N.cangle(T));
			if(cmi<cm){cmi = cm; tmi = tauj;}
		}
		taumax[i] = tmi;
		cosmax[i] = cmi;
		if(max<cmi) max=cmi;
	}
	return max;
}

// Compute the maximum (absolute) sin value of between vector srf.normal(uv(t))
// and vector this->TP(i)(t) for i=0,1,2,3, where perim[i] is 
// the same as srf.perimeter_curve(i), and t is a common parameter
// of deris[i] and TP(i).
double MGSBRepVecTP::get_perimeters_max_sin(
	const MGSurface& srf,
	double         taumax[4],
	double         sinmax[4],
	bool*          evalf	//indicates perimeters to evalate if evalf!=null
			//When evalf[i] is true, perimeter i is evaluated for 0<=i<=3.
)const{
	MGVector N(3), T(3);
	MGPosition uv(2);
	MGNDDArray tau;
	double max=0.;
	for(int i = 0; i < 4; i++){
		if(!specified(i) || (evalf && !evalf[i])){
			taumax[i] = (i % 2 == 0) ? srf.param_s_u() : srf.param_s_v();
			sinmax[i] = 0.;
			continue;
		}

		int id = i%2;
		switch(i){
		case 0:
			uv(1) = srf.param_s_v();
			tau.update_from_knot(srf.knot_vector_u());
			break;
		case 1:
			uv(0) = srf.param_e_u();
			tau.update_from_knot(srf.knot_vector_v());
			break;
		case 2:
			uv(1) = srf.param_e_v();
			tau.update_from_knot(srf.knot_vector_u());
			break;
		case 3:
			uv(0) = srf.param_s_u();
			tau.update_from_knot(srf.knot_vector_v());
			break;
		};
		//std::cout<<tau<<std::endl;
		double taus=tau[0];
		uv(id) = taus;

		double cmi=0.;
		double tmi = taus;
		if(eval(i,taus,T)){
			N = srf.normal(uv);
			cmi = fabs(N.sangle(T));
			tmi = uv(id);
		}

		int ntau=tau.length();
		for(int j = 1; j < ntau; j++){
			double tauj=tau[j];
			double tmid = (tau[j-1]+tauj)*.5;
			uv(id) = tmid;
			if(eval(i,tmid,T)){
				N = srf.normal(uv);
				double cm = fabs(N.sangle(T));
				if(cmi<cm){cmi = cm; tmi = tmid;}
			}

			uv(id) = tauj;
			if(eval(i,tauj,T)){
				N = srf.normal(uv);
				double cm = fabs(N.sangle(T));
				if(cmi<cm){	cmi = cm; tmi = tauj;}
			}
		}
		taumax[i] = tmi;
		sinmax[i] = cmi;
		if(max<cmi) max=cmi;
	}
	return max;
}

//Set i-th perimeter's TP.
//vectp[i] must be newed objects, and all of the ownership will be transferer to
//this instance.
void MGSBRepVecTP::set_TP(
	size_t i,					//perimeter numeber.
	MGPvector<MGLBRep>& vectp,
	const MGInterval& prange		//Whole perimeter's parameter range.
){
	assert(i<4);
	MGPvector<MGLBRep>& tpi=m_TP[i];
	tpi.clear();
	tpi.push_back(vectp);
	m_prange[i]=prange;
	size_t n=tpi.size();
	if(!n) return;

	double plen=prange.length().value();
	m_to_SE[2*i]=MGREqual_base(prange[0].value(),tpi[0]->param_s(),plen);
	m_to_SE[2*i+1]=MGREqual_base(prange[1].value(),tpi[n-1]->param_e(),plen);
}

//Set i-th perimeter's TP(std::vector version).
//vectp[i] must be newed objects, and all of the ownership will be transferer to
//this instance.
void MGSBRepVecTP::set_TP(
	size_t i,					//perimeter number.
	std::vector<MGLBRep*>& vectp,
	const MGInterval& prange		//Whole perimeter's parameter range.
){
	assert(i<4);
	MGPvector<MGLBRep>& tpi=m_TP[i];
	tpi.clear();

	size_t n=vectp.size();
	for(size_t j=0; j<n; j++) tpi.push_back(vectp[j]);
	vectp.clear();
	m_prange[i]=prange;
	if(!n) return;

	double plen=prange.length().value();
	m_to_SE[2*i]=MGREqual_base(prange[0].value(),tpi[0]->param_s(),plen);
	m_to_SE[2*i+1]=MGREqual_base(prange[1].value(),tpi[n-1]->param_e(),plen);
}

//Set i-th perimeter's TP as a null, as an unspecified one.
void MGSBRepVecTP::set_TP_null(size_t i){
	assert(i<4);
	m_TP[i].clear();
}

//Debug Function.
std::ostream& operator<<(
	std::ostream& ostrm, const MGSBRepVecTP& vectp
){
	ostrm<<"MGSBRepVecTP::"<<&vectp;
	for(size_t i=0; i<4; i++){
		ostrm<<std::endl;
		ostrm<<"Perimeter "<<i<<":";
		ostrm<<",m_prange="<<vectp.m_prange[i];
		ostrm<<",m_to_SE=["; if(vectp.m_to_SE[2*i])ostrm<<"1"; else ostrm<<"0";
		ostrm<<",";if(vectp.m_to_SE[2*i+1])ostrm<<"1"; else ostrm<<"0";
		ostrm<<"]"<<std::endl;
		//ostrm<<vectp.m_TP[i];		
	}
	return ostrm;
}
