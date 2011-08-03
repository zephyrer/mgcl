/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mg/KnotVector.h"
#include "mg/FSurface.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "Tl/TLInputParam.h"
#include "Tl/TLparameter.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//mgTLparameter is a proprietry class for Face tessellation.
//mgTLparameter represents one intesection of Loop and u=const(or v=const)
//line in the parameter space of a face.
//When intersection with u=const, m_t is v value, or vise versa.

////////Constructor/////////////
mgTLparameter::mgTLparameter(
	const MGFace& f,
	mgTLPoints* tlpoints,
	double crvTol,
	double surfTol,
	double max_ratio,
	double max_edge_len)
:m_face(&f), m_surface(f.surface()),
m_tess_crvError(crvTol), m_tess_srfError(surfTol), m_max_ratio(max_ratio),
m_max_edge_len(max_edge_len),m_points(tlpoints),m_texture(false){
	double wczero=MGTolerance::wc_zero();
	if(crvTol<wczero) crvTol=wczero;
	if(surfTol<wczero) surfTol=wczero;

	const MGSurface& srf=*m_surface;
	m_surface_error=srf.parameter_error();
	m_surface_error*=m_surface_error;
	const MGBox& uvbox=f.box_param();
	double u0=uvbox[0].low_point(), u1=uvbox[0].high_point();
	double v0=uvbox[1].low_point(), v1=uvbox[1].high_point();
	m_surf_prange[0]=u0; m_surf_prange[1]=u1;
	m_surf_prange[2]=v0; m_surf_prange[3]=v1;
	double rczero=MGTolerance::rc_zero();
	m_puerror=rczero*(u1-u0);
	m_pverror=rczero*(v1-v0);

	const double minimum_ratio_span=.001;
	//Subdivision beyond the 1/1000 will never be done.
	m_uError=srf.knot_vector_u().param_span()*minimum_ratio_span;
	m_vError=srf.knot_vector_v().param_span()*minimum_ratio_span;

	if(m_tess_crvError<wczero) m_tess_crvError=wczero;
	if(m_tess_srfError<wczero) m_tess_srfError=wczero;
	if(m_max_ratio<1.5) m_max_ratio=1.5;
		//1.5 is an approximation of sqrt(2). If ratio is less than sqrt(2), 
		//the subdivision will not improve the ratio squareness.
	m_max_ratio_sqr=m_max_ratio*m_max_ratio;
	m_error=m_tess_crvError;
	if(m_tess_crvError<m_tess_srfError) m_error=m_tess_srfError;
	m_error*=5.; m_error*=m_error;

	size_t nedgPerim=edgPerimLen();
	m_edgPerim=new std::vector<bool>[nedgPerim];	
		//Array of number of perimeter loops.
		//If outer loop exists, the array size is one.
		//m_edgPerim[j][i] indicates if j-th loop's i-th edge is a perimeter
		//edge or not. That is if on a perimeter or not. 
		//If true, it means the edge is a perimeter edge.

	size_t nploop=m_face->number_of_perimeter_boundaries();
	size_t i;
	if(nploop){
		//When perimeter boudaries exist.
		for(size_t j=0; j<nploop; j++){
			const MGLoop& lp=*(f.loop(j));
			size_t nedge=lp.number_of_pcells();
			m_edgPerim[j].resize(nedge);
			for(i=0; i<nedge; i++) m_edgPerim[j][i]=false;
		}
	}else if(!f.no_outer_boundaries()){
		//When one outer boudary exists.
		const MGLoop& olp=*(f.loop(size_t(0)));
		MGComplex::const_pcellItr ei=olp.pcell_begin(), ee=olp.pcell_end();
		size_t nedge=olp.number_of_pcells();
		m_edgPerim[0].resize(nedge);
		for(i=0; i<nedge; ei++, i++){
			const MGEdge* edg=edge_from_iterator(ei);
			if(edg->on_surface_perimeter()) m_edgPerim[0][i]=true;
			else m_edgPerim[0][i]=false;
		}
	}
	mgTLisectsList& uList=get_uList();
	mgTLisectsList& vList=get_vList();

	uList.push_front(u0); uList.push_back(u1);
	uList.set_middle((u0+u1)*0.5);
	vList.push_front(v0); vList.push_back(v1);
	vList.set_middle((v0+v1)*0.5);
	//Actual data of uList and vList are initialized in mgTLRects constructor
	//(by mgTLRects::init()).
}

mgTLparameter::mgTLparameter(
	const MGSurface& srf,
	mgTLPoints* tlpoints,
	double crvTol,
	double surfTol,
	double max_ratio,
	double max_edge_len)
:m_face(0), m_surface(&srf), m_edgPerim(0),
m_tess_crvError(crvTol), m_tess_srfError(surfTol), m_max_ratio(max_ratio),
m_max_edge_len(max_edge_len),m_points(tlpoints),m_texture(false){
	double wczero=MGTolerance::wc_zero();
	if(crvTol<wczero) crvTol=wczero;
	if(surfTol<wczero) surfTol=wczero;

	m_surface_error=srf.parameter_error();
	m_surface_error*=m_surface_error;
	double u0=srf.param_s_u(), u1=srf.param_e_u();
	double v0=srf.param_s_v(), v1=srf.param_e_v();
	m_surf_prange[0]=u0; m_surf_prange[1]=u1;
	m_surf_prange[2]=v0; m_surf_prange[3]=v1;
	double rczero=MGTolerance::rc_zero();
	m_puerror=rczero*(u1-u0);
	m_pverror=rczero*(v1-v0);

	const double minimum_ratio_span=.001;
	//Subdivision beyond the 1/1000 will never be done.
	m_uError=srf.knot_vector_u().param_span()*minimum_ratio_span;
	m_vError=srf.knot_vector_v().param_span()*minimum_ratio_span;

	if(m_tess_crvError<wczero) m_tess_crvError=wczero;
	if(m_tess_srfError<wczero) m_tess_srfError=wczero;
	if(m_max_ratio<1.5) m_max_ratio=1.5;
		//1.5 is an approximation of sqrt(2). If ratio is less than sqrt(2), 
		//the subdivision will not improve the ratio squareness.
	m_max_ratio_sqr=m_max_ratio*m_max_ratio;
	m_error=m_tess_crvError;
	if(m_tess_crvError<m_tess_srfError) m_error=m_tess_srfError;
	m_error*=5.; m_error*=m_error;

	mgTLisectsList& uList=get_uList();
	mgTLisectsList& vList=get_vList();

	uList.push_front(u0); uList.push_back(u1);
	uList.set_middle((u0+u1)*0.5);
	vList.push_front(v0); vList.push_back(v1);
	vList.set_middle((v0+v1)*0.5);
	//Actual data of uList and vList are initialized in mgTLRects constructor
	//(by mgTLRects::init()).
}
mgTLparameter::mgTLparameter(
	const MGFSurface& obj,
	mgTLPoints* tlpoints,
	double crvTol,
	double surfTol,
	double max_ratio,
	double max_edge_len,
	bool texture,
	double tex_surfTol,
	double tex_angleTol,
	double tex_max_edge_len
):m_edgPerim(0){
	const MGFace* face=obj.get_face_pointer();
	if(face)
		*this=mgTLparameter(*face,tlpoints,crvTol,surfTol,max_ratio,max_edge_len);
	else{
		const MGSurface* surf=obj.get_surface_pointer();
		*this=mgTLparameter(*surf,tlpoints,crvTol,surfTol,max_ratio,max_edge_len);
	}
	m_texture=texture;
	if(texture){
		m_tex_surfTol=tex_surfTol;
		if(m_tex_surfTol<m_surface_error)
			m_tex_surfTol=m_surface_error;

		double azero=MGTolerance::angle_zero();
		m_tex_angleTol=tex_angleTol;
		if(m_tex_angleTol<azero)
			m_tex_angleTol=azero;
		m_tex_max_edge_len=tex_max_edge_len;
		if(max_edge_len<0. || max_edge_len>tex_max_edge_len)
			m_max_edge_len=tex_max_edge_len;
	}
}

mgTLparameter::mgTLparameter(
	const MGFSurface& obj,	//テセレーションするフェイス
							//Must be MGFace or MGSurface.
	mgTLPoints* tlpoints,
	const mgTLInputParam& param//parameter for the tessellation.
):m_edgPerim(0){
	*this=mgTLparameter(obj,tlpoints,
		param.crvTol(),param.surfTol(),param.max_ratio(),param.max_edge_len(),
		param.texture(),param.tex_surfTol(),param.tex_angleTol(),param.tex_max_edge_len());
}

//Copy constructor.
mgTLparameter::mgTLparameter(const mgTLparameter& param2)
:m_face(param2.m_face),m_surface(param2.m_surface),
 m_surface_error(param2.m_surface_error),m_max_edge_len(param2.m_max_edge_len),
m_tess_crvError(param2.m_tess_crvError), m_tess_srfError(param2.m_tess_srfError),
m_error(param2.m_error), m_uError(param2.m_uError), m_vError(param2.m_vError),
m_edgPerim(0),m_max_ratio(param2.m_max_ratio),
m_uList(param2.m_uList), m_vList(param2.m_vList),
m_puerror(param2.m_puerror), m_pverror(param2.m_pverror),m_points(param2.m_points),
m_texture(param2.m_texture), m_tex_surfTol(param2.m_tex_surfTol),
m_tex_angleTol(param2.m_tex_angleTol),m_tex_max_edge_len(param2.m_tex_max_edge_len){
	m_max_ratio_sqr=m_max_ratio*m_max_ratio;
	for(size_t i=0; i<4; i++) m_surf_prange[i]=param2.m_surf_prange[i];

	size_t nedgPerim=edgPerimLen();
	if(nedgPerim){
		m_edgPerim=new std::vector<bool>[nedgPerim];
		for(size_t j=0; j<nedgPerim; j++){
			size_t nedge=param2.m_edgPerim[j].size();
			m_edgPerim[j].resize(nedge);
			for(size_t i=0; i<nedge; i++)
				m_edgPerim[j][i]=param2.m_edgPerim[j][i];
		}
	}
}

//////////Destructor/////////////
mgTLparameter::~mgTLparameter(){
	delete[] m_edgPerim;
}

///////Operator oveload///////
mgTLparameter& mgTLparameter::operator=(const mgTLparameter& param2){
	m_face=param2.m_face;
	m_surface=param2.m_surface;
	m_points=param2.m_points;

	delete[] m_edgPerim;
	size_t nedgPerim=edgPerimLen();
	if(nedgPerim){
		m_edgPerim=new std::vector<bool>[nedgPerim];
		for(size_t j=0; j<nedgPerim; j++){
			size_t nedge=param2.m_edgPerim[j].size();
			m_edgPerim[j].resize(nedge);
			for(size_t i=0; i<nedge; i++) m_edgPerim[j][i]=param2.m_edgPerim[j][i];
		}
	}else m_edgPerim=0;

	m_uList=param2.m_uList;
	m_vList=param2.m_vList;

	m_puerror=param2.m_puerror;
	m_pverror=param2.m_pverror;
	m_surface_error=param2.m_surface_error;
	for(size_t i=0; i<4; i++) m_surf_prange[i]=param2.m_surf_prange[i];

	m_tess_crvError=param2.m_tess_crvError;
	m_tess_srfError=param2.m_tess_srfError;

	m_uError=param2.m_uError;
	m_vError=param2.m_vError;
	m_error=param2.m_error;

	m_max_ratio=param2.m_max_ratio;
	m_max_ratio_sqr=m_max_ratio*m_max_ratio;
	m_max_edge_len=param2.m_max_edge_len;

	m_texture=param2.m_texture;
	m_tex_surfTol=param2.m_tex_surfTol;
	m_tex_angleTol=param2.m_tex_angleTol;
	m_tex_max_edge_len=param2.m_tex_max_edge_len;

	return *this;
}

size_t mgTLparameter::edgPerimLen()const{
	if(!is_face()) return 0;
	size_t nedgPerim=1;
	size_t nploop=m_face->number_of_perimeter_boundaries();
	if(!m_face->no_outer_boundaries()) if(nploop) nedgPerim=nploop;
	return nedgPerim;
}

void mgTLparameter::get_texture_parameter(
	double& tex_surfTol, double& tex_angleTol, double& tex_max_edge_len
)const{
	tex_surfTol=m_tex_surfTol;
	tex_angleTol=m_tex_angleTol;
	tex_max_edge_len=m_tex_max_edge_len;
}

ostream& operator<< (ostream& out, const mgTLparameter& para){
	out<<"TLparam::m_face="<<para.m_face<<",m_surface="<<para.m_surface
		<<",m_tess_crvError="<<para.m_tess_crvError;
	out<<",m_tess_srfError="<<para.m_tess_srfError
		<<",m_max_ratio="<<para.m_max_ratio<<endl;
	out<<",m_surface_error="<<para.m_surface_error
		<<",m_uError="<<para.m_uError<<",m_vError="<<para.m_vError
		<<",m_error="<<para.m_error<<endl;
	out<<",m_puerror="<<para.m_puerror<<",m_pverror="<<para.m_pverror;
	out<<", m_points="<<para.m_points;
	out<<", m_max_edge_len="<<para.m_max_edge_len<<endl;
	out<<", m_surf_prange=[";
	size_t i;
	for(i=0; i<4; i++){
		out<<para.m_surf_prange[i];
		if(i!=3) out<<",";
	}
	out<<"]"<<endl;
	size_t ne=para.edgPerimLen();
	for(i=0; i<ne; i++) {
		out<<"m_edgPerim["<<i<<"]"<<"::";
		size_t m=para.m_edgPerim[i].size();
		for(size_t j=0; j<m; j++) out<<para.m_edgPerim[i][j]<<" ";
		out<<endl;
	}
	out<<"****** m_uList::"<<para.m_uList;
	out<<"****** m_vList::"<<para.m_vList;

	return out;
}
