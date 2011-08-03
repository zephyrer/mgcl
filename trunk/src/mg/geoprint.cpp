/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/CCisect.h"
#include "mg/CCisect_list.h"
#include "mg/CompositeCurve.h"
#include "mg/CParam_list.h"
#include "mg/CSisect.h"
#include "mg/CSisect_list.h"
#include "mg/Ellipse.h"
#include "mg/EReal.h"
#include "mg/Interval.h"
#include "mg/Knot.h"
#include "mg/KnotArray.h"
#include "mg/KnotVector.h"
#include "mg/LBRep.h"
#include "mg/LBRepEndC.h"
#include "mg/Matrix.h"
#include "mg/NDDArray.h"
#include "mg/Object.h"
#include "mg/OscuCircle.h"
#include "mg/OscuCircleData.h"
#include "mg/Plane.h"
#include "mg/Point.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/PPRep.h"
#include "mg/RLBRep.h"
#include "mg/RSBRep.h"
#include "mg/SBRep.h"
#include "mg/BSumSurf.h"
#include "mg/SBRepTP.h"
#include "mg/SBRepEndC.h"
#include "mg/SPointSeq.h"
#include "mg/SSisect.h"
#include "mg/SSisect_list.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/Transf.h"
#include "mg/TrimmedCurve.h"
#include "mg/FSurface.h"
#include "mg/MGStl.h"
#include "mg/Tolerance.h"
#include "mgGL/Appearance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include <iomanip>
using namespace std;

//////**********Object output*********//////////

//////////// MGObject ////////////
// Output function.
ostream& MGObject::out(std::ostream& ostrm) const{
	if(m_appearance)
		ostrm<<std::endl<<","<<(*m_appearance);
	//else
	//	ostrm<<",m_appearance="<<m_appearance;
	return ostrm;
}

//////////// MGGeometry ////////////
ostream& MGGeometry::out(ostream& ostrm) const{
	MGObject::out(ostrm);
	if(m_box) ostrm<<",m_box="<<(*m_box);
	return ostrm;
}

// Output virtual function.
std::ostream& MGCurve::out(std::ostream& ostrm) const{
	MGGeometry::out(ostrm);
	return ostrm;
}

// Output virtual function.
std::ostream& MGSurface::out(std::ostream& ostrm) const{
	MGGeometry::out(ostrm);
	return ostrm;
}

//////////// MGCompositeCurve ////////////
ostream& MGCompositeCurve::out(ostream& ostrm) const{
	ostrm<<"MGCompositeCurve::"<<this;
	ostrm<<",number_of_curves="<<number_of_curves()<<endl;
	deque<MGCurve*>::const_iterator
		i=m_composite.begin(), ie=m_composite.end();
	size_t j=0;
	for(; i!=ie; i++,j++)
		ostrm<<"Curve"<<j<<"::"<<(**i);
	MGCurve::out(ostrm);
	return ostrm;
}

//////////// MGEllipse ////////////
ostream& MGEllipse::out(ostream& ostrm) const{
//	ostrm.setf(ios::scientific,ios::floatfield);
//	ostrm.precision(10);
	double t0=m_prange[0], t1=m_prange[1];
	ostrm<<"MGEllipse::"<<this;
	ostrm<<",m_circle="<<(int)m_circle
		<<",m_r="<<m_r
		<<",m_center="<<m_center<<endl
		<<",m_normal="<<m_normal
		<<",m_m="<<m_m
		<<",m_n="<<m_n<<endl
		<<",m_prange=("<<t0<<","<<t1<<")"
		<<" ,start="<<eval(t0)
		<<" ,end="<<eval(t1);
	if(m_gprange)
		ostrm<<endl<<",m_gprange=("<<m_gprange[0]<<","<<m_gprange[1]<<")";
	MGCurve::out(ostrm);
	return ostrm;
}

//////////// MGLBRep ////////////
ostream& MGLBRep::out(ostream& ostrm) const{
	ostrm<<"MGLBRep::"<<this;
	ostrm<<",Space Dimension="<<sdim()<<endl;
	ostrm<<","<<knot_vector();
	ostrm<<","<<line_bcoef();
	MGCurve::out(ostrm);
	return ostrm;
}

//////////// MGPlane ////////////
ostream& MGPlane::out(ostream& ostrm) const{
//	o.setf(ios::scientific,ios::floatfield);
//	o.precision(10);
	ostrm<<"MGPlane::"<<this;
	ostrm<<",m_root_point="<<m_root_point
		<<",m_normal="<<m_normal<<",m_d="<<m_d
		<<",m_uderiv="<<m_uderiv<<",m_vderiv="<<m_vderiv<<endl;
	MGSurface::out(ostrm);
	return ostrm;
}

//////////// MGPoint ////////////
ostream& MGPoint::out(ostream& ostrm)const{
	ostrm<<"MGPoint::"<<this;
	size_t dim=sdim();
	ostrm <<",Q(";
	for(size_t i=0; i<dim; i++){
		ostrm<<m_point(i);
		if(i<dim-1) ostrm<<",";
	}
	ostrm<<")";
	MGGeometry::out(ostrm);
	return ostrm;
}

//////////// MGRLBRep ////////////
ostream& MGRLBRep::out(ostream& ostrm) const{
	ostrm<<"MGRLBRep(non homogeneous form)::"<<this;
	ostrm<<",Space Dimension="<<sdim()<<endl;
	ostrm<<","<<knot_vector();
	ostrm<<","<<non_homogeneous_bcoef();
	MGCurve::out(ostrm);
	return ostrm;
}

//////////// MGRSBRep ////////////
ostream& MGRSBRep::out(ostream& ostrm) const{
	ostrm<<"MGRSBRep(non homogeneous form)::"<<this;
	ostrm<<",Space Dimension="<<sdim()<<endl;
	ostrm<<","<<knot_vector_u();
	ostrm<<","<<knot_vector_v();
	ostrm<<","<<non_homogeneous_bcoef();
	MGSurface::out(ostrm);
	return ostrm;
}

//////////// MGSBRep ////////////
ostream& MGSBRep::out(ostream& ostrm) const{
	ostrm<<"MGSBRep::"<<this;
	ostrm<<",m_uknot="<<m_uknot;
	ostrm<<",m_vknot="<<m_vknot;
	ostrm<<",m_surface_bcoef="<<m_surface_bcoef;
	MGSurface::out(ostrm);
	return ostrm;
}

//////////// MGStraight ////////////
ostream& MGStraight::out(ostream& ostrm) const{
	ostrm<<"MGStraight::"<<this;
	ostrm<<",m_root_point="<<m_root_point
		 <<",m_direction="<<m_direction
		 <<",m_sparam="<<m_sparam<<",m_endparam="<<m_endparam<<endl;
	if(m_sparam.finite())
		ostrm<<", start="<<eval(m_sparam.value());
	if(m_endparam.finite())
		ostrm<<", end="<<eval(m_endparam.value());
	ostrm<<endl;
	MGCurve::out(ostrm);
	return ostrm;
}

//////////// MGSurfCurve ////////////
ostream& MGSurfCurve::out(ostream& ostrm) const{
	ostrm<<"MGSurfCurve::"<<this;
	ostrm<<",start="<<start_point()<<",end="<<end_point();
	ostrm<<",m_surface="<<m_surface
		<<",m_curve="<<m_curve << endl;
	MGCurve::out(ostrm);
	return ostrm;
}

//////////// MGTrimmedCurve ////////////
ostream& MGTrimmedCurve::out(ostream& ostrm) const{
	ostrm<<"MGTrimmedCurve::"<<this;
	ostrm<<",m_curve="<<m_curve;
	if(m_curve)
		ostrm<<",m_sameRange="<<m_sameRange
		<<",m_range="<<m_range
		<<",start="<<start_point()<<",end="<<end_point()
		<<",*m_curve="<<endl<<(*m_curve);
	MGCurve::out(ostrm);
	return ostrm;
}

//Debug Function
// Output virtual function.
ostream& MGBSumSurf::out(std::ostream& ostrm) const{
	ostrm<<"MGBSumSurf::"<<this;
	ostrm<<",m_g1="<<m_g1;if(m_g1) ostrm<<","<<(*m_g1);
	ostrm<<",m_g2="<<m_g1;if(m_g2) ostrm<<","<<(*m_g2);
	ostrm<<",m_g12="<<m_g1;if(m_g12) ostrm<<","<<(*m_g12);
	MGSurface::out(ostrm);
	return ostrm;
}

ostream& MGStl::out(std::ostream& ostrm) const{
	ostrm<<"MGStl::"<<this;

	// ボックス枠の座標を出力
	ostrm<<",m_box="<<m_box;

	//Print out triangle data.
	ostrm<<std::endl;
	int nIndices = m_indices.size();
	int nTri=nIndices/3;
	assert(nTri*3==nIndices);
	ostrm<<"Number of triangles="<<nTri<<endl;
	for(int i=0, j=0; j<nTri; i+=3, j++){
		int i0=m_indices[i], i1=m_indices[i+1], i2=m_indices[i+2];
		const MGPosition& P0=m_vecPos[i0];
		const MGPosition& P1=m_vecPos[i1];
		const MGPosition& P2=m_vecPos[i2];
		ostrm<<j<<": Normal="<<m_vecNormlTriang[j]<<endl;
		ostrm<<"0("<<P0[0]<<","<<P0[1]<<","<<P0[2]<<"), 1(";
		ostrm<<P1[0]<<","<<P1[1]<<","<<P1[2]<<"), 2(";
		ostrm<<P2[0]<<","<<P2[1]<<","<<P2[2]<<")"<<endl;
	}
	MGObject::out(ostrm);
	return ostrm;
}

////////************Non Object output************///////////////

//////////// Box ////////////
ostream& operator<< (ostream& out, const MGBox& box){
	size_t dim=box.sdim();
	out <<"MGBox(";
	for(size_t i=0; i<dim; i++){
		out <<box.m_range[i];
		if(i<dim-1) out<<",";
	}
	out<<")";

	return out;
}

//////////// BPointSeq ////////////
ostream& operator<<(ostream& out, const MGBPointSeq& bps){
	int space_dim=bps.sdim();
	int length=bps.length();
	out<<"BPointSeq:: m_capacity="<<bps.capacity()<<", m_sdim="<<space_dim;
	out<<", m_length="<<length;
	int m;
	for(int i=0; i<length ; i++){
		m=space_dim*i;
		if((m-int(m/12)*12)==0) out<<endl;
		out<<i<<"(";
		for(int j=0; j<(space_dim-1) ; j++) out<<bps(i,j)<<",";
		out<<bps(i,space_dim-1)<<") ";
	}
	out<<endl;
	return out;
}

//////////// MGCCisect_list ////////////
ostream& operator<< (ostream& out, const MGCCisect_list& lst){
	out<<"MGCCisect_list::m_curve1="<<(lst.m_curve1);
	out<<",m_curve2="<<(lst.m_curve2);
	size_t n=lst.entries();
	out<<",number of isect="<<n<<endl;
	list<MGCCisect>::const_iterator itr; size_t i=0;
	for(itr=lst.begin(); itr!=lst.end(); itr++)
		out<<i++<<":"<<(*itr);
	return out;
}

//////////// MGCParam_list ////////////
ostream& operator<< (ostream& out, const MGCParam_list& lst){
	out<<"MGCParam_list::m_curve="<<(lst.m_curve);
	size_t n=lst.entries();
	out<<", number of param="<<n<<endl;
	MGCParam_list::const_Citerator itr; size_t i=0;
	for(itr=lst.begin(); itr!=lst.end(); itr++)
		out<<i++<<":"<<(*itr)<<" ";
	out<<endl;
	return out;
}

//////////// MGCSisect_list ////////////
ostream& operator<< (ostream& out, const MGCSisect_list& lst){
	out<<"MGCSisect_list::m_curve="<<(lst.m_curve);
	out<<", m_surface="<<(lst.m_surface);
	size_t n=lst.entries();
	out<<", number of isect="<<n<<endl;
	MGCSisect_list::const_CSiterator itr; size_t i=0;
	for(itr=lst.begin(); itr!=lst.end(); itr++)
		out<<i++<<":"<<(*itr);
	return out;
}

//////////// MGDefault ////////////
ostream& operator<<(ostream& o, const MGDefault& def){
	o<<"mgNULL_VEC="<<mgNULL_VEC;
	o<<", mgZERO_VEC="<<mgZERO_VEC;
	o<<", mgZERO_VEC_2D="<<mgZERO_VEC_2D<<","<<endl;
	o<<"mgNULL_Pos="<<mgNULL_Pos;
	o<<", mgORIGIN="<<mgORIGIN;
	o<<", mgORIGIN_2D="<<mgORIGIN_2D<<","<<endl;
	o<<"mgX_UVEC="<<mgX_UVEC;
	o<<", mgY_UVEC="<<mgY_UVEC;
	o<<", mgZ_UVEC="<<mgZ_UVEC<<","<<endl;
	o<<"mgX_UVEC_2D="<<mgX_UVEC_2D;
    o<< ", mgY_UVEC_2D="<<mgY_UVEC_2D<<","<<endl;
	o<<"mgEMP_INTERV="<<mgEMP_INTERV;
	o<<", mgNULL_BOX="<<mgNULL_BOX<<","<<endl;
	o<<"mgEMP_BOX="<<mgEMP_BOX<<","<<endl;
	o<<"mgEMP_BOX_2D="<<mgEMP_BOX_2D<<","<<endl;
	o<<"mgNULL_MATR="<<mgNULL_MATR<<","<<endl;
	o<<"mgUNIT_MATR="<<mgUNIT_MATR<<","<<endl;
	o<<"mgUNIT_MATR_2D="<<mgUNIT_MATR_2D<<","<<endl;
	o<<"mgNULL_TRANSF="<<mgNULL_TRANSF<<","<<endl;
	o<<"mgID_TRANSF="<<mgID_TRANSF<<","<<endl;
	o<<"mgID_TRANSF_2D="<<mgID_TRANSF_2D<<endl;
	o<<"mgNULL_KNOT_VECTOR="<<mgNULL_KNOT_VECTOR;
	o<<"mgINFINITE_KNOT_VECTOR="<<mgINFINITE_KNOT_VECTOR;
	return o;
}

//////////// MGEReal ////////////
ostream& operator<< (ostream& out, const MGEReal& er){
	if(er.finite()) out<<er.value();
	else if(er.minus_infinite()) out<<"MGINFINITE_MINUS";
	else out<<"MGINFINITE_PLUS";
	return out;
}

//////////// MGInterval ////////////
ostream& operator<< (ostream &o, const MGInterval& i2){
//		o.setf( ios::scientific, ios::floatfield );
//        o.precision( 10 );
	if(i2.empty()) o<<"I"<<"[MGINTERVAL_EMPTY]";
	else           o<<"I"<<"["<<i2.m_low<<","<<i2.m_high<<"]";
	return o;
}

//////////// MGKnot ////////////
ostream& operator<<(ostream& out, const MGKnot& knot){
	out<< "(" << knot.value() <<" " << knot.multiplicity() <<") ";
	return out;
}

//////////// MGKnotArray ////////////
ostream& operator<<(ostream& out, const MGKnotArray& klist){
	size_t n=klist.length();
	out<<"MGKnotArray::length="<<n<<endl;
	for(size_t i=0; i<n; i++) out<<klist[i]<<" ";
	out<<endl;
	return out;
}

//////////// MGKnotVector ////////////
ostream & operator << (ostream & out, const MGKnotVector & t){
	unsigned order=t.order(); size_t bdim=t.bdim();
	out<<"MGKnotVector::order="<< order<<",bdim=" <<bdim;
	out<<",current=" << t.m_current <<endl;
	int n=order+bdim; int nm1=n-1; int mult=1;
	if(n<=0) return out;
	double tprev=t(0)-(t(nm1)-t(0));
	for(int i=0; i<n; i++){
		if(t(i)==tprev) mult+=1;
		else if(i){
			if((i-int(i/12)*12)==0) out<<endl;
			out<<mult<<")"<<t(i-1)<<" ";
			mult=1;
		}
		tprev=t(i);
		if(i==nm1){
			if((i-int(i/12)*12)==0) out<<endl;
			out<<mult<<")";
			out<<tprev<<" ";
		}
	}
	out<<endl;
	return out;
}

//////////// MGLBRepEndC ////////////
ostream& operator<< (ostream& out, const MGLBRepEndC& endc){
	out<<"MGLBRepEndC::m_cond="<<endc.cond();	
	if(endc.cond()==MGENDC_1D || endc.cond()==MGENDC_12D)
		out<<", 1st deriv="<<endc.first();
	if(endc.cond()==MGENDC_2D || endc.cond()==MGENDC_12D)
		out<<", 2nd deriv="<<endc.second();
	out<<endl;
	return out;
}

//////////// MGMatrix ////////////
ostream& operator<< (ostream& out, const MGMatrix& mat){
	size_t dim=mat.sdim();
    out<<endl<<"MGMatrix(";
	for(size_t i=0; i<dim; i++){
		if(i>0) out<<"         ";
		for(size_t j=0; j<dim; j++){
			out<<mat(i,j); if(j<dim-1) out<<",";
		}
		if(i<dim-1) out<<endl;
	}
	out<<")";
	return out;
}

//////////// MGNDDArray ////////////
ostream& operator<<(ostream & out, const MGNDDArray& ndda){
	int n=ndda.length();
	out<<"MGNDDArray::capacity="<<ndda.capacity()<<",length="<<n<<endl;
	for(int i=0; i<n; i++) out<<i<<":"<<ndda(i)<<" ";
	out<<endl;
	return out;
}

//////////// MGOscuCircle ////////////
ostream& operator<<(ostream& out, const MGOscuCircle& circle){
	size_t n=circle.length();
	out<<"MGOscuCircle::length="<<n<<endl;
	for(size_t i=0; i<n; i++) out<<circle(i)<<" ";
	out<<endl;
	return out;
}

//////////// MGOscuCircleData ////////////
ostream& operator<<(ostream& out, const MGOscuCircleData& circle){
	out<< "(" << circle.index() <<" " << circle.radius() <<") ";
	return out;
}

//////////// MGPosition ////////////
ostream& operator<< (ostream& out, const MGPosition& p){
	size_t dim=p.sdim();
	out <<"P(";
	for(size_t i=0; i<dim; i++){
		out<<p.m_element(i);
		if(i<dim-1) out<<",";
	}
	out<<")";
	return out;
}

//////////// MGPosition_list ////////////
ostream& operator<< (ostream& out, const MGPosition_list& lst){
	size_t n=lst.size();
	out<<"MGPosition_list::number of position="<<n<<endl;
	MGPosition_list::const_iterator i;
	size_t j=0;
	for(i=lst.begin(); i!=lst.end(); i++) out<<j++<<":"<<(*i)<<" ";
	out<<endl;
	return out;
}

//////////// MGPPRep ////////////
ostream& operator<< (ostream& out, const MGPPRep& pprep)
//	friend istream& operator>> (istream&, MGPPRep& );
{
	unsigned order=pprep.order();
	int nbreak=pprep.nbreak();
	size_t sdim=pprep.sdim();
	out<<"MGPPRep::"<<&pprep<<"order=" <<order <<", nbreak=" <<nbreak;
	out<<", sdim=" <<sdim <<endl;
	if(!nbreak) return out;

	out<<"break_id break_point:: sdim:pp-coef(form order 0 to order-1)" <<endl;
	int j;
	for(j=0; j<nbreak-1; j++){
		out<<j<<" "<<pprep.break_point(j)<<"::";
		for(size_t k=0; k<sdim; k++){
			out<<" "<<k <<":";
			for(unsigned i=0; i<order ; i++){
				out<<pprep.coef(i,j,k)<<" ";
			}
		}
		out<<endl;
	}
	out<<j<<" "<<pprep.break_point(j)<<"::"<<endl;
	return out;
}

//////////// MGSBRepEndC ////////////
ostream& operator<< (ostream& out, const MGSBRepEndC& endc){
	out<<"MGSBRepEndC::m_cond={"
		<<endc.m_cond[0]<<","<<endc.m_cond[1]<<","
		<<endc.m_cond[2]<<","<<endc.m_cond[3]<<"}"
		<<endl;	
	for(size_t i=0; i<4; i++){
		out<<" m_11d["<<i<<"]="<<endc.m_11d[i];
		out<<" m_12d["<<i<<"]="<<endc.m_12d[i];
		out<<" m_21d["<<i<<"]="<<endc.m_21d[i];
		out<<" m_22d["<<i<<"]="<<endc.m_22d[i];
		if(endc.cond(i)==MGENDC_1D || endc.cond(i)==MGENDC_12D)
			out<<" 1st deriv("<<i<<")="<<endc.first(i);
		if(endc.cond(i)==MGENDC_2D || endc.cond(i)==MGENDC_12D)
			out<<" 2nd deriv("<<i<<")="<<endc.second(i);
	}
	return out;
}

//////////// MGSBRepTP ////////////
ostream& operator<< (ostream& out, const MGSBRepTP& tp)
{
	out<<"MGSBRepTP::"<<endl;
	for(size_t i=0; i<4; i++)
		if(tp.specified(i)) out<<" TP"<<i<<":"<<tp.TP(i);
	return out;
}

//////////// MGSPointSeq ////////////
ostream& operator<<(ostream& out, const MGSPointSeq& sp){
	size_t w=8;
	int space_dim=sp.sdim();
	int lenu=sp.length_u(), lenv=sp.length_v();
	out<<"MGSPointSeq::m_capacityu="<<sp.capacity_u()<<", m_capacityv="<<sp.capacity_v()
	   <<", m_sdim="<<space_dim;
	out<<", m_lengthu="<<lenu<<", m_lengthv="<<lenv;
	int m;

	out <<setiosflags( ios::fixed )<<setprecision(3);
	for(int j=0; j<lenv ; j++){
		for(int i=0; i<lenu ; i++){
			m=space_dim*i;
			if((m-int(m/12)*12)==0) out<<endl<<i<<","<<j<<"(";
			else out<<i<<"(";
			for(int k=0; k<(space_dim-1) ; k++){
				out<<setw(w);
				out<<sp(i,j,k)<<",";
			}
			out<<setw(w);
			out<<sp(i,j,space_dim-1)<<") ";
		}
	}
	out<<setprecision(6)<<resetiosflags(ios::fixed)<<endl;
	return out;
}

//////////// MGSSisect_list ////////////
ostream& operator<< (ostream& out, const MGSSisect_list& lst){
	out<<"MGSSisect_list::m_surface1="<<(lst.m_surface1);
	out<<",surface2="<<(lst.m_surface2);
	size_t n=lst.entries(); size_t j=0;
	out<<",number of isects="<<n<<endl;
	MGSSisect_list::const_SSiterator i;
	for(i=lst.begin(); i!=lst.end(); i++) out<<j++<<":"<<(*i)<<endl;
	return out;
}
//////////// MGTolerance ////////////
ostream & operator << (ostream & o, const MGTolerance & tl){
	o << "MGTolerance::m_mach_zero="<<tl.m_mach_zero
	  << ", m_wc_zero="<<tl.m_wc_zero
		<<", m_wc_zero_sqr="<<tl.m_wc_zero_sqr<<endl;
	o <<", m_rc_zero="<<tl.m_rc_zero<<", m_rc_zero_sqr="<<tl.m_rc_zero_sqr
		<<", m_angle_zero="<<tl.m_angle_zero<<endl;
	o <<", m_line_zero="<<tl.m_line_zero<<", m_max_knot_ratio="
		<<tl.m_max_knot_ratio<<", m_count="<<tl.m_count<<endl;
	return o;
}

//////////// MGTransf ////////////
ostream& operator<< (ostream& o, const MGTransf& trn1){
    o<<"MGTransf::"<<trn1.m_affine<<", "
		<<trn1.m_translation;
	return o;
}
//////////// MGUnit_vector ////////////
ostream& operator<< (ostream& out, const MGUnit_vector& unit) {
	size_t dim=unit.sdim();
	out <<"U(";
	for(size_t i=0; i<dim; i++){
		out<<unit.m_element[i];
		if(i<dim-1) out<<",";
	}
	out<<")";
	return out;
}

//////////// MGVector ////////////
ostream& operator<< (ostream& out, const MGVector& vec){
	size_t dim=vec.sdim();
	out <<"V(";
	for(size_t i=0; i<dim; i++){
		out<<vec.m_element[i];
		if(i<dim-1) out<<",";
	}
	out<<")";
	return out;
}

// MGINTERVAL_TYPEを標準出力に出力する。
ostream & operator << (ostream & out, MGINTERVAL_TYPE rel){
	const char* name[5]={
     "MGINTERVAL_EMPTY",             // Interval は空集合
     "MGINTERVAL_FINITE",            // 上方、下方とも有限
     "MGINTERVAL_FINITE_ABOVE",      // 上方有限、下方無限
     "MGINTERVAL_FINITE_BELOW",      // 上方無限、下方有限
     "MGINTERVAL_INFINITE"};         // 上方、下方ともに無限 
	 int i=(int)rel;
	 if(i>=5) i=4; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGPSRELATIONを標準出力に出力する。
ostream & operator << (ostream & out, MGPSRELATION rel){
	const char* name[6]={
     "MGPSREL_UNKNOWN",              // 不明
     "MGPSREL_TORSION",              // ねじれ
     "MGPSREL_ISECT",               // 交差
     "MGPSREL_PARALLEL",             // 平行
     "MGPSREL_COIN",                 // 包含
     "MGPSREL_VIRTUAL_ISECT"};      // 指定直線の延長上の点での交差
	int i=(int)rel;
	if(i>=6) i=5; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGCURVE_TYPEを標準出力に出力する。
ostream & operator << (ostream & out, MGCURVE_TYPE rel){
	const char* name[7]={
	"MGCURVE_UNKNOWN",		// 不明
	"MGCURVE_STRAIGHT",		// 直線
	"MGCURVE_ELLIPSE",		// 楕円
	"MGCURVE_SPLINE",		//B-Spline(LBRep). スプライン
	"MGCURVE_RSPLINE",		//B-Spline(Rational)
	"MGCURVE_USER1",		//Auxiliary curve type 1
	"MGCURVE_USER2"};		//Auxiliary curve type 2
	int i=(int)rel;
	if(i>=7 || i<0) i=0;
	out<<name[i];
	return out;
}

// MGELLIPSE_TYPEを標準出力に出力する。
ostream & operator << (ostream & out, MGELLIPSE_TYPE rel){
	const char* name[3]={
	 "MGELLIPSE_EMPTY"		// 空
	,"MGELLIPSE_SEGMENT"	// セグメント
	,"MGELLIPSE_CLOSED"};	// 閉じた楕円
	int i=(int)rel;
	if(i>=3) i=2; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGSTRAIGHT_TYPEを標準出力に出力する。
ostream & operator << (ostream & out, MGSTRAIGHT_TYPE rel){
	const char* name[4]={
	 "MGSTRAIGHT_EMPTY" 		// 空
	,"MGSTRAIGHT_SEGMENT"		// 線分
	,"MGSTRAIGHT_HALF_LIMIT"	// 半直線
	,"MGSTRAIGHT_UNLIMIT"};		// 無限直線
	 int i=(int)rel;
	 if(i>=4) i=3; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGSURFACE_TYPEを標準出力に出力する。
ostream & operator << (ostream & out, MGSURFACE_TYPE rel){
	const char* name[10]={
	"MGSURFACE_UNKNOWN",	// 不明
	"MGSURFACE_PLANE",		// 平面
	"MGSURFACE_CONE",		// 円錐
	"MGSURFACE_SPHERE",		// 球面
	"MGSURFACE_TORUS",		// 円環面
	"MGSURFACE_SPLINE"		//Free form surface 自由曲面                                 
							//(Tensor product surface of LBRep)
	"MGSURFACE_RSPLINE",	//Free form surface                             
							//(Tensor product surface of Rational B-Spline) 
	"MGSURFACE_CYLINDER",	//Cylinder surface(A special case of MGSURFACE_CONE)
	"MGSURFACE_USER1",		//Auxiliary surface type 1
	"MGSURFACE_USER2"		//Auxiliary surface type 2
	};
	int i=(int)rel;
	if(i>=10 || i<0) i=0;
	out<<name[i];
	return out;
}

// MGCCRELATIONを標準出力に出力する。
ostream & operator << (ostream & out, MGCCRELATION rel){
	const char* name[5]={
	"MGCCREL_UNKNOWN",	// 未知（調べていなくて、わからない、以下のいずれか）
	"MGCCREL_ISECT",	// 交差（直交、接する、一致以外）
	"MGCCREL_NORMAL",	// 直交
	"MGCCREL_TANGENT",	// 接する
	"MGCCREL_COIN"};	// 一致する
	 int i=(int)rel;
	 if(i>=5) i=4; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGCSRELATIONを標準出力に出力する。
ostream & operator << (ostream & out, MGCSRELATION rel){
	const char* name[6]={
	"MGCSREL_UNKNOWN",		// 未知
	"MGCSREL_IN",			// 曲面の内側からの交点
	"MGCSREL_OUT",			// 曲面の外側からの交点
	"MGCSREL_IN_TAN",		// 曲面の内側から接している
	"MGCSREL_OUT_TAN",		// 曲面の外側から接している
	"MGCSREL_COIN"};		// 曲線が曲面に含まれている
	 int i=(int)rel;
	 if(i>=6) i=5; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGSSRELATIONを標準出力に出力する。
ostream & operator << (ostream & out, MGSSRELATION rel){
	const char* name[4]={
	"MGSSREL_UNKNOWN",	// 未知（調べていなくて、わからない、以下のいずれか）
	"MGSSREL_ISECT",	// 交差（接する、一致以外）
	"MGSSREL_TANGENT",	// 接する
	"MGSSREL_COIN"};	// 一致する
	 int i=(int)rel;
	 if(i>=4) i=3; if(i<0) i=0;
	out<<name[i];
	return out;
}

// MGENDCONDを標準出力に出力する。
ostream & operator << (ostream & out, MGENDCOND rel){
	const char* name[5]={
	"MGENDC_UNKNOWN",	// 未知
	"MGENDC_1D",		// 1st deravative provided.
	"MGENDC_2D",		// 2nd deravative provided.
	"MGENDC_NO",		// no end cond(only positional data)
	"MGENDC_12D"};		// both 1st and 2nd deravatives provided.
	 int i=(int)rel;
	 if(i>=5) i=4; if(i<0) i=0;
	out<<name[i];
	return out;
}
