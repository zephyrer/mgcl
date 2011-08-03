/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/

#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/HHisect_vector.h"
#include "tl/TLTriangles.h"
#include "tl/TLPoints.h"
#include "tl/TLInputParam.h"
#include "tl/TLparameter.h"
#include "tl/TLRect.h"
#include "tl/TLRects.h"
#include "tl/TLData.h"
#include "tl/TLDataVector.h"
#include "tl/TLTexPoints.h"
#include "tl/TLTexPlane.h"
#include "tl/TLFans.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////// A class that contains one Face or Surface's tessellated data //////////

///////////Constructor//////////
mgTLData::mgTLData(const mgTLData& tlData)//copy constructor.
:m_tlparam(tlData.m_tlparam), m_rects(tlData.m_rects),
m_triangles(tlData.m_triangles), m_tlpoints(tlData.m_tlpoints),
m_tex_coords(tlData.m_tex_coords),m_textured(NOT_TEXTURED),
m_tex_distorted(tlData.m_tex_distorted){
	mgTLData* tld=const_cast<mgTLData*>(&tlData);
	tld->m_triangles=0;
	tld->m_rects=0;
	tld->m_tlpoints=0;
	tld->m_tlparam=0;
	tld->m_tex_coords=0;
	tld->m_tex_distorted=0;
}

mgTLData::mgTLData(
	const MGFSurface& obj,	//テセレーションするフェイス
							//Must be MGFace or MGSurface.
	const mgTLInputParam& param//parameter for the tessellation.
):m_rects(0),m_triangles(new mgTLTriangles()),m_tlpoints(new mgTLPoints()),
m_tlparam(0),m_tex_coords(0),m_textured(NOT_TEXTURED),m_tex_distorted(0){
	m_tlparam=new mgTLparameter(obj,m_tlpoints,param);
	m_rects=new mgTLRects(*m_tlparam,param.minimum_tri()/2);//Plannerな矩形のクラスを作成
	//std::cout<<*m_rects;
	triangulate(param.fanKind());
	//std::cout<<*m_rects;
}

mgTLData::mgTLData(
	const MGFSurface& obj,	//テセレーションするフェイス
	double crvTol,			//バウンダリのトレランス
	double surfTol,			//平面とみなすトレランス
	double max_ratio,		//最大アスペクト比
	MGCL::fan_kind fk,
		//fk=SINGLE_TRIANGLE:   1 triangle/FAN
		//fk=MULTIPLE_TRIANGLES: as many triangles as possible/FAN
	size_t minimum_tri,		//Specify minimum number of triangles.
	double max_edge_len
):m_rects(0),m_triangles(new mgTLTriangles()),m_tlpoints(new mgTLPoints()),
m_tlparam(0),m_tex_coords(0),m_textured(NOT_TEXTURED),m_tex_distorted(0){
	m_tlparam=new mgTLparameter(obj,m_tlpoints,crvTol,surfTol,max_ratio,max_edge_len);
	m_rects=new mgTLRects(*m_tlparam,minimum_tri/2);//Plannerな矩形のクラスを作成
	//std::cout<<(*m_tlpoints)<<std::endl;
	//std::cout<<(*m_rects)<<std::endl;
	triangulate(fk);
}

//////////// Destructor.////////
mgTLData::~mgTLData(){
	if(m_rects) delete m_rects;
	if(m_triangles) delete m_triangles;
	if(m_tlpoints) delete m_tlpoints;
	if(m_tlparam) delete m_tlparam;
	if(m_tex_coords) delete m_tex_coords;
	if(m_tex_distorted) delete m_tex_distorted;
}

// Assignment
mgTLData& mgTLData::operator=(const mgTLData& tld){
	if(m_rects) delete m_rects;
	m_rects=tld.m_rects;

	if(m_triangles) delete m_triangles;
	m_triangles=tld.m_triangles;

	if(m_tlpoints) delete m_tlpoints;
	m_tlpoints=tld.m_tlpoints;

	if(m_tlparam) delete m_tlparam;
	m_tlparam=tld.m_tlparam;

	if(m_tex_coords) delete m_tex_coords;
	m_tex_coords=tld.m_tex_coords;
	m_textured=tld.m_textured;

	if(m_tex_distorted) delete m_tex_distorted;
	m_tex_distorted=tld.m_tex_distorted;
	m_tex_distortion_ratio=tld.m_tex_distortion_ratio;

	mgTLData* ptld=const_cast<mgTLData*>(&tld);
	ptld->m_rects=0;
	ptld->m_triangles=0;
	ptld->m_tlpoints=0;
	ptld->m_tlparam=0;
	ptld->m_tex_coords=0;
	ptld->m_tex_distorted=0;
	return *this;
}

///////// member function ////////////

//Find the rects of level number (level-1)(which are supposed to be textured
//already) and texture their not-textured neighbor rects.
//Function's return value is
//false if no rects of (level-1) which were not textured were found.
//true if some rects of (level-1) which had non_textured neighbor rects
//were found and textured.
bool mgTLData::compute_level_texture(
	int level	//Level number to texture.
){
	assert(level>=2);
	bool textured=false;
	int levelm1=level-1;
	std::deque<mgTLRect*> rects;//modified rects will be prepended.

	mgTLRects::iterator i=m_rects->begin(), ie=m_rects->end();
	for(;i!=ie;i++){
		mgTLRect& recti=**i;
		if(recti.texture_level()==levelm1){
			if(recti.set_neighbor_tex_coord_data(*this,rects))
				textured=true;
		}
	}
	if(textured)
		m_textured=PARTIALLY_TEXTURED;
	return textured;
}

//compute all the non corner texture coordinates of this m_tlpoints.
//compute_texture_coordinates can be invoked only after all the rects are
//textured.
void mgTLData::compute_non_corner_texture_coordinates(){
	//set the texture coordinates of non-corner points.
	mgTLTriangles& tris=triangles();
	mgTLPoints& tluvs=tlpoints();
	mgTLTriangles::iterator j=tris.begin(), jend=tris.end();
	for(; j != jend; j++){
		mgTLTriangle& tri=**j;
		if(tri.getGeometryType()==mgTESTRIANG_STRIP)
			continue;//Since triangle has only rect corner points.
			
		mgTLRect* triRect=tri.rect();
		assert(triRect->is_textured());////////
		mgTLTriangle::iterator kid=tri.begin(), kidend=tri.end();
		for(; kid!=kidend; kid++){
			size_t id=*kid;
			MGPosition& st=get_texcoord(id);
			if(!st.is_null())
				continue;
			st=triRect->get_tex_coord(*this,tluvs[id]);
		}
		//tri.print(std::cout,*this);
	}
}

//compute all the texture coordinates of this m_tlpoints.
//compute_texture_coordinates can be invoked only after some of the rects
//are already textured.
void mgTLData::compute_texture_coordinates_from_neighbor(){
	//find a rect that is not textured and has some neighbors textured.
	mgTLRects::iterator i=m_rects->begin(), ie=m_rects->end();
	mgTLRect* rect=0;
	for(;i!=ie;i++){
		mgTLRect& recti=**i;
		if(recti.is_textured())
			continue;

		for(size_t iperi=0; iperi<4; iperi++){
			std::list<mgTLRect*> neibors=recti.neighbours(iperi);
			if(neibors.size()==0)
				continue;//continue to next perimeter.

			std::list<mgTLRect*>::iterator j=neibors.begin(), je=neibors.end();
			for(; j!=je; j++){
				mgTLRect* rectj=*j;
				if(!rectj->is_textured())
					continue;
				recti.compute_texture_by_triangles(*this);
				rect=&recti;
				break;
			}//end of j(neighbor) loop.
			if(rect)
				break;//break out of iperi loop.
		}//end of iperi loop.

		if(!rect)
			continue;
		std::deque<mgTLRect*> rects;	//modified rects will be prepended.
		rects.push_back(rect);
		while(rects.size()){
			mgTLRect* rect_modified=rects.back(); rects.pop_back();
			rect_modified->set_neighbor_tex_coord_data(*this,rects);	//modified rects will be prepended.
		}
	}
	m_textured=TOTALLY_TEXTURED;
}

//compute all the texture coordinates of this m_tlpoints.
//compute_texture_coordinates can be invoked only after triangulate().
void mgTLData::compute_texture_coordinates(
	 const MGPosition& uv,	//surface parameter of the point st.
	 const MGPosition& st,	//texture coordinate of uv.
	 const MGUnit_vector& texXaxis,//world coordinate x axis vector.
	 double distortion_ratio//When >1. distorted points(points whose ratio are more
							//than distortion_ratio or less than 1/distortion_ratio will be
							//stored in m_tex_distorted.
){
	MGHHisect isects;
	int next=0;
	compute_texture_coordinates_partially(next,isects,uv,st,texXaxis,distortion_ratio);
}

//compute all the texture coordinates of this m_tlpoints.
//compute_texture_coordinates can be invoked only after triangulate().
//Function's return value is true if textured.
bool mgTLData::compute_texture_coordinates_partially(
	int& id_hhis,	//id of hhis that define the current priority line of uv is input.
					//When id_hhis>=hhis.size(), no priority lines are input, and
					//all the rects will be textured from the point uv.
					//id of hhis to process next face(mgTLData) will be output.
	const MGHHisect& hhis,
	const MGPosition& uv,
	const MGPosition& st,	//texture coordinate of uv.
	const MGUnit_vector& texXaxis,//world coordinate x axis vector.
 	double distortion_ratio//When >1. distorted points(points whose ratio are more
							//than distortion_ratio or less than 1/distortion_ratio will be
							//stored in m_tex_distorted.
){
	int num_uvline=hhis.num_of_uvline();
	if(m_textured==TOTALLY_TEXTURED || !m_triangles || m_triangles->size()==0){
		id_hhis=num_uvline;
		return false;
	}

	//search the rect uv belongs to.
	mgTLRect* rectuvp=find_rect(uv);
	if(!rectuvp){
		id_hhis=num_uvline;
		return false;//when not found, return.
	}

	if(!m_tex_coords)
		initialize_tex_coords(distortion_ratio);
	const MGSurface& srf=surface();

	//1. set the texture coordinates of the found rect.
	if(!rectuvp->is_textured()){
		rectuvp->set_initial_tex_from_1point(srf.eval(uv),st,texXaxis,*this);
	}

	double errorSurf=srf.param_error();
	std::deque<mgTLRect*> rects;
	if(id_hhis<int(hhis.num_of_uvline())){
		int next=id_hhis;
		mgTLRect* rect_next=rectuvp;
		mgTLRect* rect_current;
		int perim_into, perim_going_out=-1;
		do{
			rect_current=rect_next;
			perim_into=perim_going_out;
			rect_current->set_texture_level(1);
			rect_current->set_neighbor_tex_coord_data(*this,rects);
			//Search the neiboring rect from fpoints.
			next=rect_current->find_neighbor_rect(
				errorSurf,next,hhis,perim_into,rect_next,perim_going_out);
		}while(rect_next);
		m_textured=PARTIALLY_TEXTURED;
		//std::cout<<(*this)<<std::endl;///////////
		id_hhis=next;
	}else{
		rects.push_back(rectuvp);
		while(rects.size()){
			mgTLRect* rect=rects.back(); rects.pop_back();
			rect->set_neighbor_tex_coord_data(*this,rects);	//modified rects will be prepended.
		}
		m_textured=TOTALLY_TEXTURED;
		id_hhis=num_uvline;
	}
	return true;
}

//compute the texture coordinates of the point uv from the data of 2 points,
//(uv1, st1) and (uv2, st2). Here uvx means the surface parameter and stx
//does texture coordinates known.
void mgTLData::compute_uv_texture_coordinate(
	const MGVector& normal,	//Normal closest plane of the objective mgTLRect.
	const MGPosition& uv,	//surface parameter of the point st to compute.
	const MGPosition& uv1,	//surface parameter of the point st1.
	const MGPosition& st1,	//texture coordinate of uv1.
	const MGPosition& uv2,	//surface parameter of the point st2.
	const MGPosition& st2,	//texture coordinate of uv2.
	MGPosition& st	//texture coordinate of uv will be output.
)const{
	const MGSurface& surf=surface();
	assert(!st1.is_null() && !st2.is_null());

	const MGPosition& xyz1=surf.eval(uv1);
	MGVector ST12=st2-st1;

	MGVector P1toP2=surf.eval(uv2)-xyz1;
	MGVector P1toP=surf.eval(uv)-xyz1;

	double L12=P1toP2.len();
	double stlen=ST12.len();
	double L=P1toP.len();
	
	double theta=P1toP2.angle2pai(P1toP,normal);
	MGMatrix rot; rot.set_rotate_2D(theta);
	MGVector stV=ST12*rot;
	stV*=L/L12;

	st=st1+stV;
}

//与えられたrectから開始してv=maxの辺よりストリップを作成する
//ストリップが求まったらm_trianglesに追加される
//Function's return value is true if pRect is triangulated,
//false if not.
bool mgTLData::createStrip(
	mgTLRect* pRect
){
	if(!pRect->hasStripCondition())
		return false;//はじめのrectがストリップにできない

	mgTLTriangle *pTriangle = new mgTLTriangle(pRect,mgTESTRIANG_STRIP);
	pTriangle->push_back(pRect->Pid(0)); pTriangle->push_back(pRect->Pid(1));
	int triad=m_triangles->size();
	pRect->set_triangled(triad);
	m_triangles->push_back(pTriangle);

	//はじめのrectのspan_uを取得する
	double u0=pRect->umin(), u1=pRect->umax();
	mgTLRect* pCurRect = pRect, *pNextRect = NULL;

	//v=maxの辺に隣接するrectを取得する。hasStripCondition=tureより隣接するrectは１つ
	//0になったらperimeterなので終了
	do{
		const std::list<mgTLRect*>& rectNeighbour = pCurRect->neighbours(2);//v=maxの辺
		if(rectNeighbour.size()==0){
			break;//隣がなくなったので終了
		}
		mgTLRect::CRecNItr iter = rectNeighbour.begin();
		pNextRect = *iter;	//隣のrectのポインタ
		if(!pNextRect->hasStripCondition())
			break;	//ストリップに出来ないので終了
		if((u0 != pNextRect->umin()) || (u1 != pNextRect->umax()))
			break;

		//ストリップになるときの処理
		pTriangle->push_back(pNextRect->Pid(0));
		pTriangle->push_back(pNextRect->Pid(1));
		pNextRect->set_triangled(triad);
		pCurRect = pNextRect;	//隣のをカレントにする
	}while(true);

	pTriangle->push_back(pCurRect->Pid(3));
	pTriangle->push_back(pCurRect->Pid(2));
	return true;
}

////////// private class //////////

//Private class to sort the fans in mgTLFans according to the vertex number.
class mgTLFanSize{
private:
	size_t m_center;	//id of mgTLFans(center id).
	size_t m_vnum;	//number of vertices of the mgTLFan.
public:
	size_t center()const{return m_center;};
	void set(size_t center, size_t vnum){m_center=center; m_vnum=vnum;};
	size_t vnum()const{return m_vnum;};
};

//mgTlTriangleをpolyの長さの順にソートするためのクラス
class mgTlfansizeSort{
public:
	bool operator()(const mgTLFanSize& tf1, const mgTLFanSize& tf2)
		const{return tf1.vnum() > tf2.vnum();}
};

//多角形からファンを作成して、m_trianglesに追加する
void mgTLData::createMultiFanSet(
	mgTLTriangles& polygons
){
	if(polygons.size()<=0)
		return;

	mgTLparameter& param=*m_tlparam;	//tessellation parameter.
	double uerror=param.isect_uerror(), verror=param.isect_verror();
	//各多角形のファンを生成する
	mgTLRect* rect=polygons.front()->rect();
	mgTLTriangles::triIterator iter = polygons.begin(), itend=polygons.end();
	mgTLPoints& tlpoints=*(param.tlpoints());
	for(; iter != itend; iter++){
		mgTLTriangle& polygon=**iter;
		if(polygon.size()<3)
			continue;
		mgTLFans fans(uerror,verror,tlpoints,polygon);//全体の三角形FANのベクトル

		//頂点と周辺の頂点リストのベクトルから3角形FANのベクトルを作成する
		size_t nfan=fans.size();
		std::vector<mgTLFanSize> fansizes(nfan);
		for(size_t i=0; i<nfan; i++){
			fansizes[i].set(i,fans[i]->size());
		}
		std::sort(fansizes.begin(), fansizes.end(), mgTlfansizeSort());
		std::vector<bool> vused(nfan,false);
			//Flag if the corresponding vertex is already processed to make fan.
			//vused[i] corrsponds to the vertex fans[i] is used or not.
		for(size_t ifan=0; ifan<nfan; ifan++){
			//Loop over fansizes vector. ifan is an id of fansizes.
			size_t center=fansizes[ifan].center();
			mgTLFan& fan=*(fans[center]);
			size_t nvert=fan.size();
			if(nvert<=1) continue;//To process the next fan.

			mgTLTriangle* cfan=new mgTLTriangle(rect,mgTESTRIANG_FAN);
			cfan->push_back(fans.pointID(center));

			//周辺の頂点のうち未使用点を抽出する
			for(size_t m=0; m<nvert; m++){
				size_t v1=fan[m];
				//現在頂点v1がすでに作成されたfanの一部であれば元のfanを分割する
				if(vused[v1]){
					if(cfan->size()>2){
						int triad=m_triangles->size();
						m_triangles->push_back(cfan);
						rect->set_triangled(triad);
					}else delete cfan;
					cfan=new mgTLTriangle(rect,mgTESTRIANG_FAN);
					cfan->push_back(fans.pointID(center));
				}else cfan->push_back(fans.pointID(v1));
			}
			if(cfan->size() > 2){
				int triad=m_triangles->size();
				m_triangles->push_back(cfan);
				rect->set_triangled(triad);
			}else delete cfan;
			vused[center]=true;
		}	

	}
}

//多角形からファンを作成して、m_trianglesに追加する
void mgTLData::createSingleFanSet(
	mgTLTriangles& polygons
){
	if(polygons.size()<=0)
		return;

	mgTLparameter&  param=*(m_tlparam);
	mgTLPoints& tlpoints=*(param.tlpoints());	//cout<<tlpoints<<endl;
	double uerror=param.isect_uerror(), verror=param.isect_verror();

	//各多角形のファンを生成する
	mgTLTriangles::triIterator iter = polygons.begin(), itend=polygons.end();
	int i=0; //cout<<"ploygon size="<<polygons.size()<<endl;
	for(; iter != polygons.end(); iter++){
		mgTLTriangle& polygon=**iter;
		mgTLRect* rect=polygon.rect();
		//cout<<endl<<i<<"::***Polygon=";polygon.print(cout,tlpoints);i++;////////////////
		mgTLFans fans(uerror,verror,tlpoints,polygon);//全体の三角形FANのベクトル
		//cout<<fans;/////////
		//頂点と周辺の頂点リストのベクトルから3角形FANのベクトルを作成する
		size_t nfan=fans.size();
//		size_t j=0;
		for(size_t center=0; center<nfan; center++){
			mgTLFan& fan=*(fans[center]);
			size_t nvert=fan.size();
			if(nvert<=1) continue;//To process the next fan.

			//周辺の頂点のうち未使用点を抽出する
			size_t nvm1=nvert-1;
			for(size_t m=0; m<nvm1; m++){
				size_t v1=fan[m];
				if(v1<center) continue;
				size_t v2=fan[m+1];
				if(v2<center){m++;continue;}
				size_t idc=fans.pointID(center), id1=fans.pointID(v1), id2=fans.pointID(v2);
				mgTLTriangle* cfan=	new mgTLTriangle(rect,idc,id1,id2);
				int triad=m_triangles->size();
				m_triangles->push_back(cfan);
				rect->set_triangled(triad);
//				cout<<"("<<i++<<","<<j++<<")::";
//				cfan->print(cout,tlpoints);
			}
		}	
		//cout<<")";
	}
}

//Compute minimum and maximum curvature of MGCL::SURFACE_CURVATURE_KIND kind.
void mgTLData::get_min_max_curvatures(
	MGCL::SURFACE_CURVATURE_KIND kind,
	double& minimum,
	double& maximum
)const{
	const MGSurface& surf=surface();
	const mgTLPoints& tlpts=tlpoints();
	size_t npoint=tlpts.size();

	double kappa;
	minimum=maximum=0.;
	mgTLPoints::const_iterator l=tlpts.begin(), lend=tlpts.end();
	if(l==lend)
		return;

	double curvature[4];
	MGUnit_vector N;
	surf.curvatures(*l++,curvature,N);
	minimum=maximum==curvature[kind];
	for(; l!=lend; l++){
		const MGPosition& uv = *l;
		surf.curvatures(uv,curvature,N);
		kappa=curvature[kind];
		if(maximum<kappa)
			maximum=kappa;
		else if(minimum>kappa)
			minimum=kappa;
	}
}

//Get texture coordinates, given m_tex_coords id i.
const MGPosition& mgTLData::get_texcoord(size_t i)const{
	return (*m_tex_coords)[i];
}
MGPosition& mgTLData::get_texcoord(size_t i){
	return (*m_tex_coords)[i];
}

//Assumed this is already textured, compute texture coordinate at the point uv,
//and the world coordinate direction of the s-axis.
//Function's return value is true if data is obtained successfully,
//false if the rect at uv is not textured and data was not obtained.
bool mgTLData::get_texXaxis(
	const MGPosition& uv,	//parameter value of the face to compute at.
	MGPosition& st,			//texture coordinate will be output.
	MGUnit_vector& texXaxis2//the s-axis direction at uv in world coordinate.
)const{
	if(textured()==NOT_TEXTURED)
		return false;
	mgTLRect* rect=find_rect(uv);
	if(!rect)
		return false;
	if(!rect->is_textured())
		return false;

	st=rect->get_tex_coord(*this,uv);
	const MGSurface& srf=surface();
	MGPlane plane=rect->compute_plane(srf);
	const MGUnit_vector& N=plane.normal();
	MGPosition uv0=rect->uv(0), uv2=rect->uv(2);
	texXaxis2=srf.eval(uv2)-srf.eval(uv0);

	const MGPosition& st0=get_texcoord(rect->Pid(0));
	const MGPosition& st2=get_texcoord(rect->Pid(2));
	MGVector STV02=st2-st0;

	double angl=STV02.angle2pai(MGVector(1.,0.),MGVector(0.,0.,1.));
	MGMatrix rot; rot.set_rotate_3D(N,angl);
	texXaxis2*=rot;
	return true;
}

//Get surface parameters, given m_tlpoints id i.
const MGPosition& mgTLData::get_uv(size_t i)const{
	return (*m_tlpoints)[i];
}

//Get surface parameters, given m_tlpoints id i.
MGPosition& mgTLData::get_uv(size_t i){
	return (*m_tlpoints)[i];
}
//Get texture coordinates and world coordinates, given m_tlpoints id i.
void mgTLData::get_texcoord_world(
	size_t i,
	MGPosition& st,	//(s,t) will be output, which is biased by base.
	MGPosition& world,//world (x,y,x) value of the point (s,t) will be output.
	const MGVector& base
)const{
	st=(*m_tex_coords)[i];
	st-=base;
	world= m_tlparam->get_surface().eval((*m_tlpoints)[i]);
}

//Obtain the number of points of the i-th triangle.
size_t mgTLData::number_of_points(size_t i)const{
	return triangle(i).size();
}

//Search neighboring faces that is alreday texture computed,
//then after neighboring rectagles texture computation, compute the whole recttangle
//in this tlData.
//This is assumed not to be texture-computed yet.
//Function's return value is true if neighboring texture-computed face is
//found and texture-computed.
//false if not texture-computed.
bool mgTLData::find_textured_and_compute(
	mgTLDataVector& tldvec,
 	double distortion_ratio//When >1. distorted points(points whose ratio are more
							//than distortion_ratio or less than 1/distortion_ratio will be
							//stored in m_tex_distorted.
){
	const MGFace& f=*(face_pointer());
	if(!(&f)) return false;
	const MGShell* shell=f.parent_shell();
	if(!shell) return false;

	if(!m_tex_coords)
		initialize_tex_coords(distortion_ratio);

	//1. compute texture coords of the neighboring edges.
	bool edge_computed=false;
	size_t nloop=f.number_of_loops();
	for(size_t i=0; i<nloop; i++){
		const MGLoop* loopi=f.loop(i);
		size_t nedge=loopi->number_of_edges();
		for(size_t j=0; j<nedge; j++){
			const MGEdge* edgej=loopi->edge(j);
			std::vector<const MGEdge*> patrner_edges=edgej->partner_edges();
			size_t npatrner_edges=patrner_edges.size();
			if(!npatrner_edges)
				continue;

			for(size_t k=0; k<npatrner_edges; k++){
				const MGEdge* pedgek=patrner_edges[k];
				const MGFace* fk=pedgek->face();
				mgTLData* tld=tldvec.find_tldata(fk);
				if(!tld)
					continue;
				if(!tld->textured())
					continue;
				compute_edge_texture_coordinates(*edgej,*pedgek,*tld);
				edge_computed=true;
			}
		}
	}
	if(!edge_computed)
		return false;

	//2. compute texture coords of all the rects of 1.
	// search the rect whose two vertices are texture computed in 1. process.
	//and set the texture data.
	mgTLRect *rect, *rect2;
	set_up_texture_data(true,rect);

	//3. compute texture coords of all the rects that shares vertices of 2.
	if(!rect) return false;
	set_up_texture_data(false,rect2);

	//4. compute texture coords of all the rects neiboring to 2 and 3's 
	//rects.
	compute_texture_coordinates_from_neighbor();

	m_textured=TOTALLY_TEXTURED;
	return true;
	
	//std::cout<<(*this)<<std::endl;////////////////*************
}

//Initialize m_tex_coords and all the rects' m_center and m_normal.
void mgTLData::initialize_tex_coords(
	double distortion_ratio
){
	//initialize the texture coordinates.
	m_tex_coords=new mgTLTexPoints(m_tlpoints->size());//,image_width,image_height);
	if(distortion_ratio>1.){
		m_tex_distorted=new std::vector<MGPosition>;
		m_tex_distortion_ratio=distortion_ratio;
	}
}

//set up texture data of this data from already set tlTexPoints data.
void mgTLData::set_up_texture_data(
	bool from_not_modify,//true if only tlTexPoints data set_as_not_modify is to use.
	mgTLRect*& rect	//1st modified rect will be output.
){
	rect=0;
	mgTLRects& rects=tlrects();
	mgTLTexPoints& sts=tlTexPoints();
	mgTLPoints& uvs=tlpoints();
	const MGSurface& surf=surface();

	mgTLRect* rect2;
	do{
	rect2=0;
	// search the rect whose two vertices are texture computed.
	mgTLRects::iterator irect=rects.begin(), irecte=rects.end();
	for(; irect!=irecte; irect++){
		mgTLRect& recti=**irect;
		if(recti.is_textured())
			continue;

		size_t iperi1, iperi2;
		size_t id1, id2;

		if(from_not_modify){
//When vertices of shared edges of a shell computation
		for(iperi1=0; iperi1<4; iperi1++){
			id1=recti.Pid(iperi1);
			if(sts.was_set_as_not_modify(id1))
				break;
		}
		if(iperi1==4) continue;

		for(iperi2=iperi1+1;iperi2<4; iperi2++){
			id2=recti.Pid(iperi2);
			if(sts.was_set_as_not_modify(id2))
				break;
		}
		if(iperi2==4) continue;

		}else{
//When general vertices computation
		for(iperi1=0; iperi1<4; iperi1++){
			id1=recti.Pid(iperi1);
			if(!(sts[id1].is_null()))
				break;
		}
		if(iperi1==4) continue;
		

		iperi2=iperi1+1;
		for(;iperi2<4; iperi2++){
			id2=recti.Pid(iperi2);
			if(!(sts[id2].is_null()))
				break;
		}
		if(iperi2==4) continue;
		
		}

		recti.compute_texture_by_triangles(*this);
		if(!rect)
			rect=&recti;
		else{
			if(rect->uspan()<recti.uspan() || rect->vspan()<recti.vspan())
				rect=&recti;
		}
		rect2=rect;
	}
	if(from_not_modify)break;
	}while(rect2);
}

//compute texture coordinates of this_edge, given parameter edge partnerE connected to already
//computed face.
void mgTLData::compute_edge_texture_coordinates(
	const MGEdge& this_edge,	//edge of this TLData connected to already computed face.
	const MGEdge& partnerE,		//a partner edge of this_edge. The partner edge is already
		//computed face's edge.
	mgTLData& partner_tld		//partnerE's TLData which is already texture-computed.
){
	if(!m_triangles)
		return;
	if(m_triangles->size()==0)
		return;

	std::auto_ptr<const MGCurve> common_curve((this_edge.make_binder_with_curve())->curve_limitted());
	const MGBox& edge_box=this_edge.box();
	mgTLPoints& tluvs=tlpoints();

	double crvError=parameter().get_tess_crvError();
	double errorSave=MGTolerance::set_wc_zero(crvError*2.);
	//Set texture coordinates that share the common edge(common_curve).
	//Modified points' rect will be sotred in rects.
	const MGSurface& surf=surface();//this_edge's surface.
	size_t nuvs=tluvs.size();
	for(size_t i=0; i<nuvs; i++){
		const MGPosition& uv=tluvs[i];
		if(edge_box.includes(uv)){
			MGPosition P=surf.eval(uv);
			double t;
			if(common_curve->on(P,t)){
				double t_partner=partnerE.param_pcell(t);
				MGPosition uv_partner=partnerE.eval(t_partner);
				MGPosition st_partner=partner_tld.tex_coord(uv_partner);
				set_texcoord_as_not_modify(i,st_partner);
			}
		}
	}
	MGTolerance::set_wc_zero(errorSave);
}

//Record the distortion point.
//uvr is (uvr[0], uvr[1]) is the surface parameter,
//and uvr[2] is the distortion ratio.
void mgTLData::record_tex_distortion(const MGPosition& uvr){
	m_tex_distorted->push_back(uvr);
}

void mgTLData::set_texcoord(size_t i, double s, double t){
	m_tex_coords->set_texcoord(i,s,t);
}
void mgTLData::set_texcoord(size_t i, const MGPosition &texc){
	m_tex_coords->set_texcoord(i,texc);
}
void mgTLData::set_texcoord_as_not_modify(size_t i, const MGPosition &texc){
	m_tex_coords->set_texcoord_as_not_modify(i,texc);
}

//Obtain the number of triangles obtained.
size_t mgTLData::number_of_triangles()const{
	return m_triangles->size();
}

//Obtain the MGFace pointer if this is the data of MGFace,
//else null will be returned.
const MGFace* mgTLData::face_pointer()const{
	return &(m_tlparam->get_face());
}

const MGSurface& mgTLData::surface()const{
	return m_tlparam->get_surface();
}

void mgTLData::triangulate(
	MGCL::fan_kind fkind			//実行パラメータデフォルトは1Fan,1Triangle
){
	//ストリップセットを作成し、残りからファンを作成するための多角形を作成する
	mgTLparameter& param=parameter();
	mgTLRects::RecIterator i=m_rects->begin(), ie=m_rects->end();
	for(; i != ie; i++){
		mgTLRect& rect=**i;
		//if(rect.umin()>53. && rect.umax()<107. &&
		//	rect.vmin()>104. && rect.vmax()<177.)
		//	cout<<rect;//////////////

		//OUT or ストリップ作成済み のときは飛ばす
		if(rect.status()==MGRECT_OUT || rect.is_triangled())
			continue;

		//ストリップの作成を行う
		if(fkind==MGCL::SINGLE_TRIANGLE || fkind==MGCL::MULTIPLE_TRIANGLES){
			if(createStrip(&rect))
				continue;
		}

		//ファンの作成を行う
		//rectからポリゴンを作成する
		mgTLTriangles polygons;
		rect.createPolygon(param,polygons);

		//Generate trianlges from polygons.
		if(fkind==MGCL::SINGLE_TRIANGLE || fkind==MGCL::SINGLE_TRIANGLE_NO_STRIP)
			//1Fan,1Triangle
			createSingleFanSet(polygons);
		else
			createMultiFanSet(polygons);
			//polygons.print(cout,*m_tlparam.tlpoints());
	}
}

//Compute uv's texture coordinates.
MGPosition mgTLData::tex_coord(const MGPosition& uv){
	mgTLRect* rect=find_rect(uv);
	return rect->get_tex_coord(*this,uv);
}

//Compute maximum  width and height out of all the tex_planes.
void mgTLData::tex_plane_maximum(double& width, double& height)const{
	width=height=0.;
	if(!textured())
		return;

	const std::vector<mgTLTexPlane*>& planes=tex_planes();
	std::vector<mgTLTexPlane*>::const_iterator i=planes.begin(), ie=planes.end();
	for(; i!=ie; i++){
		const mgTLTexPlane& pli=**i;
		double wi=pli.width(), hi=pli.height();
		if(wi>width)
			width=wi;
		if(hi>height)
			height=hi;
	}
}

//Obtain i-th triangle in this mgTLData.
const mgTLTriangle& mgTLData::triangle(size_t i)const{
	return (*m_triangles)[i];
}
mgTLTriangle& mgTLData::triangle(size_t i){
	return (*m_triangles)[i];
}

//Obtain the i-th triangle's type.
//=mgTESTRIANG_FAN:fan, =mgTESTRIANG_STRIP:strip.
mgTESTRIANG mgTLData::triangle_type(size_t i)const{
	return triangle(i).getGeometryType();
}

//Obtain the i-th traiangle's j-th point's surface parameter value.
const MGPosition& mgTLData::uv(size_t i, size_t j)const{
	return m_triangles->m_triangles[i]->uv(j,*m_tlpoints);
}
MGPosition& mgTLData::uv(size_t i, size_t j){
	return m_triangles->m_triangles[i]->uv(j,*m_tlpoints);
}

std::ostream& operator<<(std::ostream& out, const mgTLData& data){
	out<<"mgTLData="<<(&data)<<", m_textured="<<data.m_textured<<std::endl;
	if(data.m_rects)
		out<<std::endl<<*(data.m_rects);
	if(data.m_triangles)
		out<<std::endl<<*(data.m_triangles);
	if(data.m_tlpoints)
		out<<*(data.m_tlpoints);
	if(data.m_tex_coords)
		out<<*(data.m_tex_coords);
	out<<"  END OF mgTLData."<<std::endl;
	return out;
}
