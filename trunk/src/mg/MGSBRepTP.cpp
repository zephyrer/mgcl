/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/SBRepTP.h"
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

//Member Data
//	MGLBRep* m_TP[4];	//Tangent Plane will be stored.
			//Tangent plane m_TP[i] is a line b-representation of
			//(unit)normal vector along the i-th perimeter.
			//Parameter range of the TP is the same as u or v parameter
			//range of the corresponding surface representation.

//Friend Function

////////Constructor////////

//Default Constructor, will be set as no TPs' are specified.
MGSBRepTP::MGSBRepTP()
{m_TP[0]=m_TP[1]=m_TP[2]=m_TP[3]=0;}

//Copy Constructor
MGSBRepTP::MGSBRepTP(const MGSBRepTP& tp)
{
	m_TP[0]=m_TP[1]=m_TP[2]=m_TP[3]=0;
	for(size_t i=0; i<4; i++)
		if(tp.m_TP[i]) m_TP[i]=new MGLBRep(*(tp.m_TP[i]));
}

MGSBRepTP::MGSBRepTP(const MGSurface& brep)
{
	double param_s_u = brep.param_s_u(), param_e_u = brep.param_e_u(),
		   param_s_v = brep.param_s_v(), param_e_v = brep.param_e_v();

	MGTolerance::push();
	size_t i;
	for(size_t iperim = 0; iperim < 4; iperim++){
		MGTolerance::set_line_zero(MGTolerance::angle_zero() * 0.25);
		MGCurve* crv = brep.perimeter_curve(iperim);
		const MGKnotVector& tempKnotVector = crv->knot_vector();
		MGKnotVector knotVector(tempKnotVector);	//指定オーダーのノットベクトルに作り替える
		for(i = tempKnotVector.order() - 1; i < tempKnotVector.bdim(); i++){
			double	tpara = 0.0,					//テンポラリ
				spara = tempKnotVector(i),		//スパンの始点
				epara = tempKnotVector(i + 1);	//スパンの終点
			if(epara - spara < crv->param_error())continue;	//マルチノットのときの処理(RLBRepのみ)
			//1スパンの分割数を決定する
			MGInterval interval(spara, epara);
			int ndiv = crv->calc_div_num(interval);			//オフセット曲線の分割数
			double shortspan = (epara - spara) / ndiv;
			tpara = spara;
			for (int j = 0; j < ndiv; j++){knotVector.add_data(tpara); tpara += shortspan;}
		}

		//制御点を生成する
		MGNDDArray dataPoint;
		dataPoint.update_from_knot(knotVector);
		size_t len = dataPoint.length();
		MGBPointSeq bp1(len, brep.sdim());
		for(i = 0; i < len; i++){
			MGPosition pos;
			switch (iperim){
			case 0:
				pos = brep.unit_normal(dataPoint(i), param_s_v);
				break;
			case 1:
				pos = brep.unit_normal(param_e_u, dataPoint(i));
				break;
			case 2:
				pos = brep.unit_normal(dataPoint(i), param_e_v);
				break;
			case 3:
				pos = brep.unit_normal(param_s_u, dataPoint(i));
				break;
			}
			bp1.store_at(i, pos);
		}
		int error = 0;
		//精度十分の曲線を生成する
		MGLBRep* tangentPL=new MGLBRep(dataPoint, bp1, knotVector, error);

		//余分なKnotを削除
		// 2.0* sinθ/2 ≒ θ とする。
		MGTolerance::set_line_zero(MGTolerance::angle_zero());
		tangentPL->remove_knot();
		m_TP[iperim] = tangentPL;
		delete crv;
	}
	MGTolerance::pop();
}

//Destructor
MGSBRepTP::~MGSBRepTP(){
	for(size_t i=0; i<4; i++) if(m_TP[i]) delete m_TP[i];
}

//Assignment.
MGSBRepTP& MGSBRepTP::operator=(const MGSBRepTP& tp){
	for(size_t i=0; i<4; i++){
		if(m_TP[i]){delete m_TP[i]; m_TP[i]=0;}
		if(tp.m_TP[i]) m_TP[i]=new MGLBRep(*(tp.m_TP[i]));
	}
	return *this;
}

//Member Function

//Set i-th perimeter's TP.
void MGSBRepTP::set_TP(size_t i, const MGLBRep& tp)
{
	assert(i<4);
	if(m_TP[i]) delete m_TP[i];
	m_TP[i]=new MGLBRep(tp);
}

//Set i-th perimeter's TP(auto_ptr version).
void MGSBRepTP::set_TP(size_t i, std::auto_ptr<MGLBRep>& tp){
	assert(i<4);
	if(m_TP[i]) delete m_TP[i];
	m_TP[i]=tp.release();
}

//Set i-th perimeter's TP as a null, as an unspecified one.
void MGSBRepTP::set_TP_null(size_t i)
{
	assert(i<4);
	if(m_TP[i]){delete m_TP[i]; m_TP[i]=0;}
}

// Compute the maximum (absolute) cos value of between vector deris[i](t) 
// and vector this->TP(i)(t) for i=0,1,2,3, where t is a common
// parameter of the data point obtained from deris[i]'s knot vector.
//Function's return value is the max out of cosmax[.].
double MGSBRepTP::get_perimeters_max_cos(
	const MGPvector<MGLBRep>& deris,
	double                    taumax[4],
	double                    cosmax[4]
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
		N = deris[i]->eval(taus);
		T = TP(i).eval(taus);
		double cmi = fabs(N.cangle(T));
		double tmi = taus;

		int ntau=tau.length();
		for(int j = 1; j < ntau; j++){
			double tauj=tau[j];
			double tmid = (tau[j-1]+tauj)*.5;
			N = deris[i]->eval(tmid);
			T = TP(i).eval(tmid);
			double cm = fabs(N.cangle(T));
			if(cmi < cm){
				cmi = cm; tmi = tmid;
			}

			N = deris[i]->eval(tauj);
			T = TP(i).eval(tauj);
			cm = fabs(N.cangle(T));
			if(cmi < cm){
				cmi = cm; tmi = tauj;
			}
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
double MGSBRepTP::get_perimeters_max_sin(
	const MGSurface& srf,
	double         taumax[4],
	double         sinmax[4],
	bool*          eval	//indicates perimeters to evalate if eval!=null
			//When eval[i] is true, perimeter i is evaluated for 0<=i<=3.
)const{
	MGVector N(3), T(3);
	MGPosition uv(2);
	MGNDDArray tau;
	double max=0.;
	for(int i = 0; i < 4; i++){
		if(!specified(i) || (eval && !eval[i])){
			taumax[i] = (i % 2 == 0) ? srf.param_s_u() : srf.param_s_v();
			sinmax[i] = 0.;
			continue;
		}

		const int id = (i % 2 == 0) ? 0 : 1;
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

		N = srf.normal(uv);
		T = TP(i).eval(taus);
		double cmi = fabs(N.sangle(T));
		double tmi = uv(id);

		int ntau=tau.length();
		for(int j = 1; j < ntau; j++){
			double tauj=tau[j];
			double tmid = (tau[j-1]+tauj)*.5;
			uv(id) = tmid;
			N = srf.normal(uv);
			T = TP(i).eval(tmid);
			double cm = fabs(N.sangle(T));
			if(cmi < cm){
				cmi = cm; tmi = tmid;
			}

			uv(id) = tauj;
			N = srf.normal(uv);
			T = TP(i).eval(tauj);
			cm = fabs(N.sangle(T));
			if(cmi < cm){
				cmi = cm; tmi = tauj;
			}
		}
		taumax[i] = tmi;
		sinmax[i] = cmi;
		if(max<cmi) max=cmi;
	}
	return max;
}

//Compute maximun abs(cons(theta)), where theta=angle of TP(i) and  corresponding 
//perimeter[i]'s start and end points' tangent vector.
double MGSBRepTP::max_cos(
	const MGCurve*	perimeter[4]//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
)const{
	size_t cid02[2]={3,1};
	size_t cid31[2]={0,2};
	double t[8]={perimeter[3]->param_s(), perimeter[1]->param_s(),
					perimeter[0]->param_e(), perimeter[2]->param_e(),
					perimeter[3]->param_e(), perimeter[1]->param_e(),
					perimeter[0]->param_s(), perimeter[2]->param_s()};

	double cosmax=0.;
	for(size_t i=0; i<4; i++){//for all the perimeters.

	size_t* cid;
	if(i%2) cid=cid31; else cid=cid02;
	if(specified(i)){
		const MGLBRep& ncrv=TP(i);
		double s[2]={ncrv.param_s(), ncrv.param_e()};
		for(size_t m=0; m<2; m++){//for start and end
			const MGCurve& peri=*(perimeter[cid[m]]);
			MGVector tan=peri.eval(t[i*2+m],1);
			MGVector N=ncrv.eval(s[m]);
			double cos=fabs(N.cangle(tan));
			if(cos>cosmax) cosmax=cos;
		}
	}

	}
	return cosmax;
}
double MGSBRepTP::max_cos(
	const MGPvector<MGLBRep>& perimeters//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
)const{
	const MGCurve*	peris[4]=
		{perimeters[0], perimeters[1], perimeters[2], perimeters[3]};
	return max_cos(peris);
}
