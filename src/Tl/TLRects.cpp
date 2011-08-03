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
#include "Tl/TLTexPlane.h"
#include "Tl/TLPoints.h"
#include "Tl/TLparameter.h"
#include "Tl/TLisectsList.h"
#include "Tl/TLRect.h"
#include "Tl/TLRects.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//////////// mgTLRects ///////////
//mgTLRects is a proprietry class for Face tessellation.
//mgTLRects holds all the subdivided rectangles for the triangulation.

bool u_is_long(const MGSurface& srf, mgTLRect* rect){
	double u0=rect->umin(), u1=rect->umax(), v0=rect->vmin(), v1=rect->vmax();
	MGVector P00=srf.eval(u0,v0), P01=srf.eval(u0, v1), P10=srf.eval(u1,v0);
	MGVector alongU0=P10-P00;
	double lenu=alongU0%alongU0;
	MGVector alongV0=P01-P00;
	double lenv=alongV0%alongV0;
	return lenu>=lenv;
}

////////// constructor /////////////

mgTLRects::mgTLRects(
	mgTLparameter& param,	//tessellation parameter.
		//Point((u,v) of the face) array will be returned into param.m_points.
		//mgTLRect's m_Pid[.] is the id of this points.
		//points[m_Pid[.]] is the MGPosition data (u,v) of the face.
	size_t divnum	//Minimum number of the devided rectangles.
):m_param(&param),m_points(param.tlpoints()){
	const MGSurface& srf=param.get_surface();
	double tesSrfTol=param.get_tess_srfError();

	mgTLRect* rect2;
	mgTLRect* rect=init();//rect is the initial rectangle.
	int recStatus=rect->status();//Save the status.
	//cout<<(*rect)<<endl; cout<<(*m_param)<<endl;

	std::stack<mgTLRect*> rectStack, prerects;
	bool divideU;
	int st;

	double log2=0.6931472;
	int max_divnum=int(log(double(divnum))/log2);
	if(divnum>0) max_divnum++;
	while(true){
		st=rect->status();
			//	2:MGRECT_OUT,	//Whole rectangle is outside the Face.
		if(st==MGRECT_OUT) m_rects.push_back(rect);
		else{
			if(rect->m_divnum<max_divnum){
				prerects.push(rect2=rect->subdivide(*this,u_is_long(srf,rect)));
					//cout<<"2="<<*(rect2)<<endl<<" 1="<<*rect<<endl;///////
				continue;
			}
			rectStack.push(rect);//cout<<(*rect)<<endl;
		}
		if(prerects.empty()) break;
		rect=prerects.top(); prerects.pop();
	}

	if(rectStack.empty()) goto end_process;
	rect=rectStack.top(); rectStack.pop();
	while(true){
		st=rect->status();
		//	1:MGRECT_IN,	//Whole rectangle is inside the Face.
		//	2:MGRECT_OUT,	//Whole rectangle is outside the Face.
		//	3:MGRECT_ON,	//Some trimming curve is crossing the rectangle.
		//					//So part is inside, and part is outside the face.
		//	4:MGRECT_OVER,	//In the case that whole inner boundary loop is included
							//in the rectangle, perimeter loops exists, or at least
							//one edge of the outer loop is not on the surface perimeter.
		//	5:MGRECT_ONANDOVER
		if(st==MGRECT_OUT) goto PushandNext;
		if(st==MGRECT_IN || (st==MGRECT_ON && rect->has_trim())){
			if(is_flat(*rect,divideU)){
				rect_perimeter_subdivide(rect);
				goto nextRect;
			}else{
				if(divideU){if(rect->uspan()<=param.get_UError()) goto PushandNext;}
				else{	    if(rect->vspan()<=param.get_VError()) goto PushandNext;}
			}
		}else{
			divideU=rect->ratio_sqr(param)>=1.;
		}
		//If rect is not planar, the status is OVER, or rect has no trims,
		//we subdivide rect and check again.
		//if(u0>.1.5 && v0>75.)
		//	rect=rect;///////cout
		rectStack.push(rect2=rect->subdivide(*this,divideU));
			//cout<<"2="<<*(rect2)<<endl<<" 1="<<*rect<<endl;///////
		continue;

PushandNext:
		m_rects.push_back(rect);
nextRect:
		if(rectStack.empty()) break;
		rect=rectStack.top(); rectStack.pop();
	}

end_process:
	//Normalize each rectangle's trim points.
	RecIterator i=m_rects.begin(), ie=m_rects.end();
	for(; i!=ie; i++){
		mgTLRect& recti=**i;
		int sta=recti.status();
		if(sta!=MGRECT_OUT && sta!=MGRECT_IN) recti.normalize_trim_points();
	}
	//cout<<(*this)<<endl;
}

mgTLRects::mgTLRects(const mgTLRects& rects2)	//Copy constructor.
:m_param(rects2.m_param), m_points(rects2.m_param->tlpoints()){
	RecIterator j=begin(), je=end();
	for(; j!=je; j++) delete *j;

	size_t n=rects2.size();
	m_rects.resize(n);

	const_RecIterator i=rects2.begin(), ie=rects2.end();
	j=begin(), je=end();
	for(; i!=ie; i++, j++) (*j)=new mgTLRect(**i);
}

///////// destructor /////////////
mgTLRects::~mgTLRects(){
	RecIterator i=begin(), ie=end();
	for(; i!=ie; i++) delete *i;

	size_t n=m_texture_planes.size();
	for(size_t j=0; j<n; j++){
		delete m_texture_planes[j];
	}
}

////////// operator overload //////////
mgTLRects& mgTLRects::operator= (const mgTLRects& rects2){
	m_param=rects2.m_param;
	m_points=rects2.m_param->tlpoints();

	RecIterator j=begin(), je=end();
	for(; j!=je; j++) delete *j;

	size_t n=rects2.size();
	m_rects.resize(n);

	const_RecIterator i=rects2.begin(), ie=rects2.end();
	j=begin(), je=end();
	for(; i!=ie; i++, j++) (*j)=new mgTLRect(**i);
	return *this;
}

class mgTLperiPoint{//Local class for mgTLRects::init().
public:
	size_t m_pnum;	//Perimeter number of the surface.
	mgTLisect m_ip;	//End point of non-perimeter edge.
	mgTLperiPoint(size_t pnum, const mgTLisect& ip):m_pnum(pnum),m_ip(ip){;};
};

//Get a loop's increase or decrease flag for mgTLisect from the start and punm.
int mgTLprInc(
	bool start,		//start point or not. =true: start.
	size_t pnum){	//perimeter number.
	int increase;
	if(pnum==1 || pnum==2) increase=-1;
	else increase=1;
	if(!start) increase*=-1;
	return increase;
}

//Store mgTLperiPoint of ei' start(start=true) or end to isvec.
//This is for an exclusive use of mgTLRects::init().
void mgTLprPpush(
	std::vector<mgTLperiPoint>& isvec,//Local variable of init().
	const MGLoop& lp,	//Loop of the edge ei.
	size_t ei,			//Edge id that indicates the target edge,
						//which must be a perimeter edge.
	bool start			//true if starting point of the ei, false if ending.
){
	const MGEdge* edg=lp.edge(ei);
	MGPosition uv;
	if(start) uv=edg->start_point(); else uv=edg->end_point();

	const MGSurface& srf=*(lp.surface());
	size_t pnum; srf.on_a_perimeter(uv(0), uv(1), pnum);
	double t=uv[pnum%2];
	int increase=mgTLprInc(!start,pnum);
	double s; if(start) s=edg->param_s(); else s=edg->param_e();
	isvec.push_back(
		mgTLperiPoint(pnum,	mgTLisect(t,MGLEPoint(lp,ei,s), increase)) );
}

////////Member function//////////

//Add position data P into m_position.
//Function's return value i is the id of the stored position.
//m_position[i] will be P.
size_t mgTLRects::add_point(double u, double v){
	return m_points->add(u,v);
}

//Find mgTLRect that includes the surface parameter uv.
mgTLRect* mgTLRects::find_rect(const MGPosition& uv){
	// search the rect uv belongs to.
	mgTLRects::iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		if((**i).includes(uv))
			break;
	}
	if(i==ie)
		return 0;//when not found.
	else
		return *i;
}

//init() builds the 1st rectangle information of the whole face.
//Function's return value is mgTLRect* of the 1st rectangle.
//init() also builds u=min and max, v=min and max mgTLisects of
//mgTLisectsList uList and vList in the m_param.
mgTLRect* mgTLRects::init(){

	mgTLparameter& param=*m_param;
	const MGSurface& srf=param.get_surface();
	const MGFace& f=param.get_face();
	double u0,u1, v0,v1;
	if(param.is_face()){
/*		MGBox uvbox=f.box_param();
		uvbox*=1.005;
		const MGInterval& urng=uvbox[0];
		const MGInterval& vrng=uvbox[1];

		double u0s=srf.param_s_u();
		u0=urng.low_point();
		if(u0<u0s) u0=u0s;

		double u1s=srf.param_e_u();
		u1=urng.high_point();
		if(u1>u1s) u1=u1s;

		double v0s=srf.param_s_v();
		v0=vrng.low_point();
		if(v0<v0s) v0=v0s;

		double v1s=srf.param_e_v();
		v1=vrng.high_point();
		if(v1>v1s) v1=v1s;
*/
		const MGBox& uvbox=f.box_param();

		const MGInterval& urng=uvbox[0];
		const MGInterval& vrng=uvbox[1];
		u0=urng.low_point(); u1=urng.high_point();
		v0=vrng.low_point(); v1=vrng.high_point();
	}else{
		const double* prang=param.suf_param_range();
		u0=prang[0]; u1=prang[1];
		v0=prang[2]; v1=prang[3];
	}
	mgTLRect* rect=new mgTLRect(u0,u1,v0,v1);

	if(!param.is_face()){
		rect->m_Pid[0]=add_point(u0,v0);
		rect->m_Pid[1]=add_point(u1,v0);
		rect->m_Pid[2]=add_point(u1,v1);
		rect->m_Pid[3]=add_point(u0,v1);
		for(size_t i=0; i<4; i++) rect->m_Pin[i]=true;
		rect->m_status=MGRECT_IN;
		return rect;
	}

	if(f.no_outer_boundaries()){
		for(size_t i=0; i<4; i++) rect->m_Pin[i]=true;
		if(f.number_of_inner_boundaries()) rect->m_status=MGRECT_OVER;
		else rect->m_status=MGRECT_IN;
		rect->m_Pid[0]=add_point(u0,v0);
		rect->m_Pid[1]=add_point(u1,v0);
		rect->m_Pid[2]=add_point(u1,v1);
		rect->m_Pid[3]=add_point(u0,v1);
		return rect;
	}

	rect->m_Pin[0]=f.in_range(u0,v0);
	rect->m_Pin[1]=f.in_range(u1,v0);
	rect->m_Pin[2]=f.in_range(u1,v1);
	rect->m_Pin[3]=f.in_range(u0,v1);
/*
	//When corner points are on a non perimeter-edge, it should be out of the
	//face parameter range.
	size_t nloop=f.number_of_perimeter_boundaries();
	if(!nloop) nloop=1;//In this cadse, loop is only the one outer boundary.
	const double uspan=u1-u0; const double vspan=v1-v0;
	for(size_t j=0; j<nloop; j++){
	std::vector<bool>& edgPerimj=m_param->get_edgPerim()[j];
	const MGLoop& lp=*(f.loop(j));
	size_t nedge=lp.number_of_pcells();
	for(size_t i=0; i<nedge; i++){
	if(!edgPerimj[i]){//If this edge is not on a surface perimeter.
		const MGEdge& edg=*(lp.edge(i));
		MGPosition ES[2]={edg.start_point(), edg.end_point()};
		for(size_t m=0; m<2; m++){
		if(MGREqual_base(ES[m][0],u0,uspan)){
			if(MGREqual_base(ES[m][1],v0,vspan)){
				rect->m_Pin[0]=false; continue;
			}
			if(MGREqual_base(ES[m][1],v1,vspan)){
				rect->m_Pin[3]=false; continue;
			}
		}
		if(MGREqual_base(ES[m][0],u1,uspan)){
			if(MGREqual_base(ES[m][1],v0,vspan)){
				rect->m_Pin[1]=false; continue;
			}
			if(MGREqual_base(ES[m][1],v1,vspan)){
				rect->m_Pin[2]=false; continue;
			}
		}
		}
	}
	}
	}
*/

	rect->m_Pid[0]=add_point(u0,v0);
	rect->m_Pid[1]=add_point(u1,v0);
	rect->m_Pid[2]=add_point(u1,v1);
	rect->m_Pid[3]=add_point(u0,v1);

	const size_t nploop=f.number_of_perimeter_boundaries();
	std::vector<mgTLperiPoint> isvec;
		//All the end point of non-perimeter edge will be stored in this isvec.
	if(nploop){

//When perimeter boudaries exist.
	for(size_t j=0; j<nploop; j++){
		const MGLoop& lp=*(f.loop(j));
		size_t nedge=lp.number_of_pcells();
		size_t pnum1, pnum2;
		lp.both_end_on_perimeter(pnum1, pnum2);

		int increase=mgTLprInc(true,pnum1);
		double t=lp.start_point()[pnum1%2];
		isvec.push_back(
			mgTLperiPoint(pnum1,mgTLisect(t,lp.start_LPoint(),increase)));

		increase=mgTLprInc(false,pnum2);
		t=lp.end_point()[pnum2%2];
		isvec.push_back(
			mgTLperiPoint(pnum2,mgTLisect(t,lp.end_LPoint(),increase)));
	}

	}else{

//When one outer boudary exists.
	const MGLoop& olp=*(f.loop(size_t(0)));//cout<<olp<<endl;////////////////
	size_t nedge=olp.number_of_edges();
	std::vector<bool>& edgPerim0=m_param->get_edgPerim()[0];

	//Find 1st id of edgPerim0 which is not a perimeter edge.
	size_t i,i1, i1Pre;
	for(i1=0; i1<nedge; i1++){ if(!edgPerim0[i1]) break;};

	if(i1<nedge){//If a non-perimeter edge found.
		bool wasPeri=false;
		i1Pre=i1++;
		for(i=1; i<nedge; i++){
			i1=i1%nedge;
			if(edgPerim0[i1]){//If this edge is on a surface perimeter.
				if(!wasPeri) mgTLprPpush(isvec,olp,i1,true);
				wasPeri=true;
			}else{//If this edge is not on a surface perimeter.
				if(wasPeri) mgTLprPpush(isvec,olp,i1Pre,false);
				wasPeri=false;
			}
			i1Pre=i1++;
		}
		if(wasPeri) mgTLprPpush(isvec,olp,i1Pre,false);
			//When last edge was on a perimeter.
	}

	}

	const int nis=isvec.size();
	if(nis){

//Change isvec to mgTLRecisects format. Also, build u=min and max, v=min and max
//mgTLisects of mgTLisectsList uList and vList.
	mgTLisectsList& uList=param.get_uList();
	mgTLisectsList& vList=param.get_vList();
	mgTLisects& u0isects=*(uList.m_isectsList.begin());
	mgTLisects& u1isects=*(uList.m_isectsList.rbegin());
	mgTLisects& v0isects=*(vList.m_isectsList.begin());
	mgTLisects& v1isects=*(vList.m_isectsList.rbegin());

//Find the location of isvec whose pnum is first minimum one.
	int m, sid=0, nism1=nis-1;
	for(m=0; m<nism1; m++){
		size_t mp1=m+1;
		if(isvec[m].m_pnum>isvec[mp1].m_pnum){sid=mp1; break;}
	}

	int idu0,idv1, nu0, nu1, nv0, nv1;
	nu0=nu1=nv0=nv1=0;
	for(m=0; m<nis; m++){
		size_t pnum=isvec[sid].m_pnum;
		if(pnum==0){
			v0isects.push_back(isvec[sid].m_ip); nv0++;
		}else if(pnum==1){
			u1isects.push_back(isvec[sid].m_ip); nu1++;
		}else if(pnum==2){
			if(!nv1) idv1=sid;
			nv1++;
		}else{
			if(!nu0) idu0=sid;
			nu0++;
		}
		sid++; if(sid==nis) sid=0;
	}

	mgTLRecisects& recisects=rect->m_Recisects;
	recisects.resize(nis);
	mgTLRecisects::ReciItr recisItr=recisects.begin();
	mgTLisects::iterator isi;
	isi=v0isects.begin();
	for(m=0; m<nv0; m++) *(recisItr++)=mgTLRecisect(0,isi++);
	isi=u1isects.begin();
	for(m=0; m<nu1; m++) *(recisItr++)=mgTLRecisect(1,isi++);
	if(nv1){
		size_t idv11=idv1+nv1;
		for(m=0; m<nv1; m++) v1isects.push_back(isvec[--idv11%nis].m_ip);
		isi=v1isects.begin()+nv1-1;
		for(m=0; m<nv1; m++) *(recisItr++)=mgTLRecisect(2,isi--);
	}
	if(nu0){
		size_t idu00=idu0+nu0;
		for(m=0; m<nu0; m++) u0isects.push_back(isvec[--idu00%nis].m_ip);
		isi=u0isects.begin()+nu0-1;
		for(m=0; m<nu0; m++) *(recisItr++)=mgTLRecisect(3,isi--);
	}
	//cout<<(*rect)<<endl;
	//Check rap-around trim points(last and 1st status).
	int firstS=recisects[0].status(),
		lastS=recisects[nism1].status();
	if((lastS==1 && (firstS==2 || firstS==0)) || (lastS==0 && firstS==2)){
		mgTLRecisect is=recisects[nism1];
		for(int ii=nism1-1; ii>=0; ii--)
			recisects[ii+1]=recisects[ii];
		recisects[0]=is;
	}

	//Arrange rect->m_Pin[.] according to the result of mgTLRecisect.
	for(m=0; m<nism1; m+=2){
		int pm=recisects[m].perim(), pmp1=recisects[m+1].perim();
		if(pm != pmp1){
			int pnext=pm;
			do{
				pnext=(pnext+1)%4;
				rect->m_Pin[pnext]=false;
			}while(pnext!=pmp1);
		}
	}

	}else{
		if(!rect->m_Pin[0] || !rect->m_Pin[1]
			|| !rect->m_Pin[2] || !rect->m_Pin[3])
		//When no non-perimeter edges exist and one of the corner points are
		//out of face, all of the corner points should be out of the face.
		rect->m_Pin[0]=rect->m_Pin[1]=rect->m_Pin[2]=rect->m_Pin[3]=false;
	}

	//Set m_status of initial rect.
	if(f.number_of_inner_boundaries()){
		if(nis) rect->m_status=MGRECT_ONANDOVER;
		else rect->m_status=MGRECT_OVER;
	}else{
		if(nis) rect->m_status=MGRECT_ON;
		else if(!rect->m_Pin[0] || !rect->m_Pin[1] || 
			!rect->m_Pin[2] || !rect->m_Pin[3]) rect->m_status=MGRECT_OVER;
		else{
			if(nploop) rect->m_status=MGRECT_OVER;//Case that perimeter loops exists
			else{//Case of one outer boundary loop.
				std::vector<bool>& pedge=(param.get_edgPerim())[0];
				size_t i, npe=pedge.size();
				for(i=0; i<npe; i++){if(!pedge[i]) break;};
				if(i<npe) rect->m_status=MGRECT_OVER;
				else rect->m_status=MGRECT_IN;
			}
		}
	}
	return rect;
}

double get_length(
	MGPosition Pn[9],
	bool& direction
){
	MGVector alongU0=Pn[2]-Pn[0];
	MGVector alongU1=Pn[8]-Pn[6];
	double lenu=alongU0%alongU0+alongU1%alongU1;
	MGVector alongV0=Pn[6]-Pn[0];
	MGVector alongV1=Pn[8]-Pn[2];
	double lenv=alongV0%alongV0+alongV1%alongV1;
	direction=lenu>=lenv;

	double len;
	if(direction){
		len=lenu;
	}else{
		len=lenv;
	}
	len*=0.5;
	return len;
}

//Test if surface of the rect is flat, that is, if within the tolerance.
//When m_param.texture() is true, test if rect is within the texture tolerance
//and if so, the texture plane will be pushed back to m_texture_planes.
bool mgTLRects::is_flat(
	mgTLRect& rect,
	bool& direction//  true: u-direction is more non flat.
					// false: v-direction is more non flat.
){
	size_t i;	//id of Pn[].
	MGPosition Pn[9];
	MGVector Nn[9];
	MGPosition P;
	MGUnit_vector N;
	const MGSurface& surf=m_param->get_surface();
	double surftol=m_param->get_tess_srfError();
	double u0=rect.umin(), u1=rect.umax(), v0=rect.vmin(), v1=rect.vmax();
	surf.compute_sample_point(u0,u1,v0,v1,Pn,P,N,Nn);

	double x,xmax=-1., d=P%N;
	double dist[9];
	bool flat=true;
	for(i=0; i<9; i++){
		x=dist[i]=d-Pn[i]%N;
		if(x<0.) x=-x;
		if(x>xmax) xmax=x;
		if(x>surftol) flat=false;
	}

	double max_edge_len=m_param->get_max_edge_len();
	bool texture=m_param->texture();
	double tex_stol,tex_atol,tex_maxedge_len;//Texture parameter.
	MGPlane tex_plane;
	double tex_width, tex_height;
	if(texture)
		m_param->get_texture_parameter(tex_stol,tex_atol,tex_maxedge_len);

	if(flat){
		double len=get_length(Pn,direction);//length of the max perimeter.
		if(max_edge_len>0.){
			double melen2=max_edge_len;
			melen2*=melen2;
			if(len>melen2) flat=false;
		}
		if(!texture) return flat;

		//test and get texture mapping plane.
		if(!rect.tex_within_tol()){
			if(tex_maxedge_len>0.){
				if(len>tex_maxedge_len*tex_maxedge_len)
					return false;
			}

			if(test_to_get_approximate_plane(Pn,Nn,P,N,
				tex_stol,tex_atol,tex_plane,tex_width,tex_height)){
				mgTLTexPlane* tpl=new mgTLTexPlane(tex_plane,tex_width,tex_height,m_rects.size());
				m_texture_planes.push_back(tpl);
				rect.set_tex_within_tol();
			}else
				return false;
		}
		return flat;
	}

	double um=(u0+u1)*0.5, vm=(v0+v1)*0.5;
	MGVector dfdu=surf.eval(um,vm,1,0), dfdv=surf.eval(um,vm,0,1);
	double ulen=dfdu.len()*(u1-u0), vlen=dfdv.len()*(v1-v0);
	if(ulen*5.<vlen) direction=false;
	else if(ulen>vlen*5.) direction=true;
	else{
		double udif=MGCL::Max3(dist[0],dist[1],dist[2])
					+MGCL::Max3(dist[3],dist[4],dist[5])
					+MGCL::Max3(dist[6],dist[7],dist[8]);
		double vdif=MGCL::Max3(dist[0],dist[3],dist[6])
					+MGCL::Max3(dist[1],dist[4],dist[7])
					+MGCL::Max3(dist[2],dist[5],dist[8]);
		double error=MGTolerance::wc_zero_sqr()*15.;

		double difmax;
		if(udif>=vdif){
			difmax=udif;
			direction=true;
		}else{
			difmax=vdif;
			direction=false;
		}
		if(difmax<=error){
			direction=ulen>=vlen;
		}
	}

	if(!texture) return false;
	//test and get texture mapping plane.
	bool dir2;
	if(!rect.tex_within_tol()){
		double len=get_length(Pn,dir2);//length of the max perimeter.
		if(tex_maxedge_len>0.){
			if(len>tex_maxedge_len*tex_maxedge_len)
				return false;
		}
		if(test_to_get_approximate_plane(Pn,Nn,P,N,
			tex_stol,tex_atol,tex_plane,tex_width,tex_height)){
			mgTLTexPlane* tpl=new mgTLTexPlane(tex_plane,tex_width,tex_height,m_rects.size());
			m_texture_planes.push_back(tpl);
			rect.set_tex_within_tol();
		}
	}
	return false;
}

ostream& operator<< (ostream& out, const mgTLRects& rects){
	size_t n, j;
	out<<*(rects.m_param)<<endl;
	n=rects.m_rects.size();
	out<<"TLRects::num of Rects="<<n<<endl;
	mgTLRects::const_RecIterator i=rects.begin(), ie=rects.end();
	j=0;
	for(; i!=ie; i++) out<<j++<<"="<<(**i)<<endl;
//	return out;///////

	n=rects.points().size();
	out<<"TLRects::num of points="<<n;
	j=0;
	for(j=0; j<n; j++){
		if(j%4) out<<","; else out<<endl;
		out<<j<<"="<<(*rects.m_points)[j];
	}
	n=rects.m_texture_planes.size();
	out<<std::endl<<"TLRects::num of texture_planes="<<n;
	for(j=0; j<n; j++){
		if(j%4) out<<","; else out<<endl;
		out<<j<<"="<<rects.m_texture_planes[j];
	}
	return out;
}

//Check if rect perimeter(that is rect's u=min,max or v=min,max parameter line)
//is within surface tolerance.
//If a rect  perimeter is not in surface tolerance,
//subdivide uvrect and push back to m_rects.
void mgTLRects::rect_perimeter_subdivide(
	mgTLRect* rect
){
	mgTLparameter& param=parameter();
	const MGSurface& srf=param.get_surface();
	double error2=param.get_tess_srfError(); error2*=error2;	

	std::stack<mgTLRect*> rectStack;
	while(true){
		if(rect->status()==MGRECT_OUT){
			m_rects.push_back(rect);
		}else{
			if(rect->perim_within_tol(true, srf, error2)){//Check along u.
				if(rect->perim_within_tol(false, srf, error2)){//Check along v.
					surface_perimeter_subdivide(rect);
				}else{
					rectStack.push(rect->subdivide(*this,false));
					continue;
				}
			}else{
				rectStack.push(rect->subdivide(*this,true));
				continue;
			}
		}
		if(rectStack.empty()) break;
		rect=rectStack.top(); rectStack.pop();
	}
}

//Check if surface perimeter is within curve tolerance.
//If the perimeter is not in curve tolerance, subdivide uvrect
// and push back to m_rects.
void mgTLRects::surface_perimeter_subdivide(mgTLRect* rect){
	mgTLparameter& param=parameter();
	const MGSurface& srf=param.get_surface();
	double tesCrvTol=param.get_tess_crvError();

	bool bSubDivU;
	//u方向に分割するかどうかのフラグ、checkPerimeterBoundary=trueのときは終わり

	std::stack<mgTLRect*> rectStack;
	while(true){
		if(rect->status()==MGRECT_OUT){
			m_rects.push_back(rect);
		}else if(rect->checkPerimeterBoundary(m_points,tesCrvTol,srf,bSubDivU)){
		//If perimeter is ok, check maximum ratio.
			ratio_subdivide(rect);
		}else{
			bool toosmall;
			if(bSubDivU) toosmall=rect->uspan()<=param.get_UError();
			else		 toosmall=rect->vspan()<=param.get_VError();
			if(toosmall){
				m_rects.push_back(rect);
			}else{
				rectStack.push(rect->subdivide(*this,bSubDivU));
				continue;
			}
		}
		if(rectStack.empty()) break;
		rect=rectStack.top(); rectStack.pop();
	}
}

//Check if ratio of u span and v span of uvrect is within maximum.
//If the ratio exceeds maximum, subdivide uvrect and push back to m_rects.
void mgTLRects::ratio_subdivide(mgTLRect* rect){
	mgTLparameter& param=parameter();
	double max_ratio_sqr=param.get_max_ratio_sqr();
	double one_max_ratio_sqr=1./max_ratio_sqr;

	double r=rect->ratio_sqr(param);// = du/dv
	const bool udiv=(r>max_ratio_sqr);
	if(one_max_ratio_sqr<=r && !udiv) m_rects.push_back(rect);
	else{
		std::stack<mgTLRect*> rectStack;
		while(true){
			if(rect->status()==MGRECT_OUT){
				m_rects.push_back(rect);
			}else{
				r=rect->ratio_sqr(param);// = du/dv
				if(udiv && r<=max_ratio_sqr) m_rects.push_back(rect);
				else if(!udiv && r>=one_max_ratio_sqr) m_rects.push_back(rect);
				else{
					bool toosmall;
					if(udiv) toosmall=rect->uspan()<=param.get_UError();
					else	 toosmall=rect->vspan()<=param.get_VError();
					if(toosmall){
						m_rects.push_back(rect);
					}else{
						rectStack.push(rect->subdivide(*this,udiv));
						continue;
					}
				}
			}
			if(rectStack.empty()) break;
			rect=rectStack.top(); rectStack.pop();
		}
	}

}

std::ostream& operator<< (std::ostream& out, const mgTLTexPlane& texplane){
	out<<"mgTLTexPlane::m_plane="<<texplane.m_plane<<std::endl;
	out<<"mgTLTexPlane::m_width="<<texplane.m_width<<", m_height="<<texplane.m_height;
	out<<", m_index="<<texplane.m_index;
	return out;
}
