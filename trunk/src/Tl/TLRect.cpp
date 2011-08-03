/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Interval.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/LBRep.h"
#include "mg/CCisect_list.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/HHisect.h"
#include "Tl/TLparameter.h"
#include "Tl/TLisectsList.h"
#include "Tl/TLRect.h"
#include "Tl/TLRects.h"
#include "Tl/TLTriangle.h"
#include "Tl/TLTriangles.h"
#include "Tl/TLData.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//mgTLRect is a proprietry class for Face tessellation.
//mgTLRect holds all the necessary information for the triangulation
//of a retangle face parameter space.

////////// Constructor ///////////////

////////// constructor /////////////
mgTLRect::mgTLRect():m_status(0),m_divnum(0),
m_textured(0),m_texcoord_ratio(1.),m_triangle(-1),m_tex_within_tol(false){
	m_Pid[0]=m_Pid[1]=m_Pid[2]=m_Pid[3]=0;
	m_Pin[0]=m_Pin[1]=m_Pin[2]=m_Pin[3]=0;
}

//Construct from the parameter box (u0, v0) - (u1, v1)
mgTLRect::mgTLRect(double u0, double u1, double v0, double v1)
: m_u0(u0), m_u1(u1), m_v0(v0),m_v1(v1)
,m_status(1),m_divnum(0),
m_textured(0),m_texcoord_ratio(1.),m_triangle(-1),m_tex_within_tol(false){
	m_Pid[0]=m_Pid[1]=m_Pid[2]=m_Pid[3]=0;
	m_Pin[0]=m_Pin[1]=m_Pin[2]=m_Pin[3]=0;
}

//Construct a higher parameter value mgTLRect by subdividing the input rect
// along u(uDirection=ture) or v direction. This constructor builds
//the following data:
//m_urange, m_vrange, m_Pid, m_Pin, m_status,m_divnum.
//m_Pin[0], [3] (when uDirection=true) or m_Pin[0], [1] are not set.
mgTLRect::mgTLRect(bool uDirection, mgTLRect& rect)
:m_status(rect.m_status),
m_textured(0),m_texcoord_ratio(1.),m_triangle(-1),
m_tex_within_tol(rect.m_tex_within_tol){
	if(uDirection){
		m_u0=(rect.m_u0+rect.m_u1)*0.5; m_u1=rect.m_u1;
		m_v0=rect.m_v0; m_v1=rect.m_v1;

		m_Pid[0]=0; m_Pin[0]=0;
		m_Pid[1]=rect.m_Pid[1]; m_Pin[1]=rect.m_Pin[1];
		m_Pid[2]=rect.m_Pid[2]; m_Pin[2]=rect.m_Pin[2];
		m_Pid[3]=0; m_Pin[3]=0;
	}else{
		m_u0=rect.m_u0; m_u1=rect.m_u1;
		m_v0=(rect.m_v0+rect.m_v1)*0.5; m_v1=rect.m_v1; 

		m_Pid[0]=0; m_Pin[0]=0;
		m_Pid[1]=0; m_Pin[1]=0;
		m_Pid[2]=rect.m_Pid[2]; m_Pin[2]=rect.m_Pin[2];
		m_Pid[3]=rect.m_Pid[3]; m_Pin[3]=rect.m_Pin[3];
	}
	rect.m_divnum++;
	m_divnum=rect.m_divnum;
}

////////// Member Function /////////////

//Add vertices of the trimming line(boundary line) from isects()[iss] to
// [iss+1], converting the boundary line to a polyline.
void mgTLRect::add_boundary_vertices(
	double crvTol,
	size_t iss,
	mgTLTriangle& polygon,
	mgTLPoints& Tlpoints
)const{
//	cout<<"******"<<(*this);//////////////
	double u0=umin(), u1=umax(), v0=vmin(), v1=vmax();
	//if(u0>0.21 && u1<.32 && v0>.42 && v1<.65)
	//	u0=u0;
	const mgTLRecisects& trims = isects();
	const MGLEPoint& trim1 = trims[iss].is()->isect();
	const MGLEPoint& trim2 = trims[iss+1].is()->isect();
	const MGLoop* ploop = trim1.loop();
	const MGSurface* surf=ploop->surface();
	const int numEdge=ploop->number_of_pcells();
	const int edge1=trim1.edge_num(), edge2=trim2.edge_num();
	double ctol2=crvTol*crvTol;

	polygon.push_back(get_point_id(trim1.eval(),Tlpoints));//Starting point.

	int num=edge2-edge1, curEdge=edge1;
	if(trim2<trim1) num+=numEdge;

	for(int i=0; i<=num; i++, curEdge++){
		if(curEdge == numEdge){
			if(!ploop->closed()) break;	//開いたループのとき終わり
			curEdge = 0;	//エッジが終わったら最初から始める
		}
		//パラメータエッジカーブpcrv,バインダエッジカーブbcrvを作成する
		const MGEdge* pedge = ploop->edge(curEdge);
		MGTrimmedCurve pcrv(pedge->trimmed_curve());//cout<<"pcrv="<<pcrv<<endl;/////////
		const MGEdge* bedge = pedge->make_binder_with_curve();//cout<<(*pedge)<<(*bedge)<<endl;
		const MGSurface* surf=pedge->star_surface();

		if(bedge->number_of_partner_members()<=1){
		//When binder edge has only one partner member.********

		MGSurfCurve pcrvS(*surf,pcrv);
		//pcrvのそれぞれのパラメータ範囲を求める
		double ss, se;	//parameter of pcell
		if(i==0) ss=trim1.param(); else ss=pcrv.param_s();
		if(i==num) se=trim2.param(); else se=pcrv.param_e();

		MGPosition P0=pcrvS.eval(ss);
		std::stack<double> sstk;
		sstk.push(se);
		while(sstk.size()){
			double s=sstk.top();
			MGPosition P1=pcrvS.eval(s);
			double sm=(ss+s)*.5;
			MGPosition PmS=pcrvS.eval(sm);
			MGVector dif1=PmS-P0, dif2=P1-P0;
			double a=dif1%dif2, l1=dif1%dif1, l2=dif2%dif2;
			if((l1*l2-a*a)<=ctol2*l2){
				polygon.push_back(get_point_id(pcrv.eval(s),Tlpoints));
				sstk.pop();
				P0=P1;ss=s;
			}else{
				sstk.push(sm);
			}
		}

		}else{
		//When binder edge has more than one partner members.*******
		//cout<<"pcrv start="<<pcrv.start_point()<<", end="<<pcrv.end_point()<<endl;

		MGTrimmedCurve bcrv(bedge->trimmed_curve());
		//cout<<"bcrv="<<bcrv<<endl;/////////
		bool SameDirection=pedge->equal_direction_to_binder();

		//pcrv,bcrvのそれぞれのパラメータ範囲を求める
		double ss, se;	//parameter of pcell
		double ts, te;	//parameter of bcell
		if(i==0){
			ss=trim1.param();
			ts=pedge->param_bcell(ss);
		}else{
			ss=pcrv.param_s();
			if(SameDirection) ts=bcrv.param_s();
			else ts=bcrv.param_e();
		}
		if(i==num){
			se=trim2.param();
			te=pedge->param_bcell(se);
		}else{
			se=pcrv.param_e();
			if(SameDirection) te=bcrv.param_e();
			else te=bcrv.param_s();
		}
		//cout<<"ts"<<ts<<",te="<<te<<endl;///////////			

		//BinderEdgeをcrvTolでポリゴン化する
//		cout<<"**pre "<<(bcrv)<<ts<<","<<te<<endl;
//		if(ts>=54.3477 && te<=54.3479)
//				ts=ts;
		MGInterval trimI;
		if(SameDirection){
//			if(ts>=te) continue;
			trimI=MGInterval(ts, te);
		}else{
//			if(te>=ts) continue;
			trimI=MGInterval(te, ts);
		}
		if(trimI.length()<=bcrv.param_error()) continue;
		bcrv.limit(trimI); //cout<<"aft "<<(bcrv)<<ts<<te<<endl;
		double dPreLineTol = MGTolerance::set_line_zero(crvTol);
		MGLBRep trimLbrep(bcrv, 2);//cout<<trimLbrep;////////////
		MGTolerance::set_line_zero(dPreLineTol);

		//BinderEdgeをポリゴン化した結果をpolygonに追加する
		int npoint=trimLbrep.bdim();
		const MGKnotVector& t=trimLbrep.knot_vector();
		MGPosition pnt;
		double s=ss, sguess;
		double* sguessp=&sguess;
		double ratio=(se-ss)/(te-ts);
		if(!SameDirection){	//pedge,bedgeの向きが逆の時の処理
			int j=npoint-1,jend=1;
			//if(i==num) jend=2;
			for(; j>=jend; j--){
				sguess=s+(t[j]-t[j+1])*ratio;
				s=pedge->param_pcell(t[j], sguessp);
				if(s<ss) s=ss; if(s>se) s=se;
				pnt=pcrv.eval(s);
				if(pnt[0]<u0) pnt(0)=u0; if(pnt[0]>u1) pnt(0)=u1;
				if(pnt[1]<v0) pnt(1)=v0; if(pnt[1]>v1) pnt(1)=v1;
				//cout<<"t="<<t[j]<<"BP="<<bcrv.eval(t[j])<<"s="<<s<<","<<pnt<<" SP="<<surf->eval(pnt)<<endl;//////////
				polygon.push_back(get_point_id(pnt,Tlpoints));
			}
		}else{
			int j=2, jend=npoint;
			//if(i==num) jend=npoint-1;
			for(; j<=jend; j++){
				sguess=s+(t[j]-t[j-1])*ratio;
				s=pedge->param_pcell(t[j], sguessp);
				if(s<ss) s=ss; if(s>se) s=se;
				pnt=pcrv.eval(s);
				if(pnt[0]<u0) pnt(0)=u0; if(pnt[0]>u1) pnt(0)=u1;
				if(pnt[1]<v0) pnt(1)=v0; if(pnt[1]>v1) pnt(1)=v1;
				//cout<<"t="<<t[j]<<"BP="<<bcrv.eval(t[j])<<"s="<<s<<","<<pnt<<" SP="<<surf->eval(pnt)<<endl;//////////
				polygon.push_back(get_point_id(pnt,Tlpoints));
			}
		}
	}

	}

	//polygon.push_back(Tlpoints.add(trim2.eval()));//Ending point.
}

//Add rect2 after rect1 in the neibour sequence of this perimeter peri.
void mgTLRect::add_neighbour(
	size_t peri,		//Perimeter num of this.
	mgTLRect* rect1,
	mgTLRect* rect2
){
	container_type& neibrs=m_neighbours[peri];
	RecNItr i=neibrs.begin(), ie=neibrs.end();
	for(;i!=ie;i++)
		if((*i)==rect1){
			i++; neibrs.insert(i,rect2);
			return;
		}
}

//Add vertices of the neighboring rectangles' vertices at the perimeter peri
//from parameter s1 to s2. s1 and s2 are in the order of the rectangle. That is
//if peri==0 or 1, in the order of the u or v parameter value, and
//if peri==2 or 3, in the reverse order or the u or v parameter value.
void mgTLRect::add_neighbor_vertices(
	size_t peri,
	double s1,
	double s2,
	mgTLTriangle& polygon
)const{//隣のrectの頂点を求める。頂点はパラメータの小さい順に並んでいる。
	const std::list<mgTLRect*>& nrects = neighbours(peri);
	if(peri==0){
		CRecNItr i=nrects.begin(), ie=nrects.end();
		for(; i!=ie; i++){
			double u=(*i)->umax();
			if(s1<u && u<s2) polygon.push_back((*i)->Pid(2));
		}
	}else if(peri==1){
		CRecNItr i=nrects.begin(), ie=nrects.end();
		for(; i!=ie; i++){
			double v = (*i)->vmax();
			if(s1<v && v<s2) polygon.push_back((*i)->Pid(3));
		}
	}else if(peri==2){
		std::list<mgTLRect*>::const_reverse_iterator
			i=nrects.rbegin(), ie=nrects.rend();
		for(; i!=ie; i++){
			double u = (*i)->umin();
			if(s2<u && u<s1) polygon.push_back((*i)->Pid(0));
		}
	}else{
		std::list<mgTLRect*>::const_reverse_iterator
			i=nrects.rbegin(), ie=nrects.rend();
		for(; i!=ie; i++){
			double v = (*i)->vmin();
			if(s2<v && v<s1) polygon.push_back((*i)->Pid(1));
		}
	}
}

//Test if all corner points are in the face or not.
//Return true if all corner points are in.
bool mgTLRect::all_corners_in()const{
	if(!Pin(0)) return false;
	if(!Pin(1)) return false;
	if(!Pin(2)) return false;
	if(!Pin(3)) return false;
	return true;
}

//get center(surface parameter) of this rect.
MGPosition mgTLRect::center_uv()const{
	return MGPosition((umin()+umax())*.5,(vmin()+vmax())*.5);
}

//Change neighbours of this rect's neibrs1 from N1Itr to neibrs1.end() to neibrs2's
//neighbours. Here, neibrs2 is rect2->m_neighbours[peri], and
//neibrs1 is this->m_neighbours[peri].
void mgTLRect::change_neighbours(
	mgTLRect* rect2,
	size_t peri,
	mgTLRect::RecNItr N1Itr
){
	container_type& neibrs1=m_neighbours[peri];
	mgTLRect::RecNItr N1ItrE=neibrs1.end();
	if(N1ItrE==N1Itr) return;

	container_type& neibrs2=rect2->m_neighbours[peri];
	mgTLRect::RecNItr N2Itr=neibrs2.end(), N2ItrE;
	neibrs2.splice(N2Itr,neibrs1,N1Itr,N1ItrE);
	N2Itr=neibrs2.begin();
	N2ItrE=neibrs2.end();
	size_t periOpo=(peri+2)%4;
	N2ItrE--;
	if(N2Itr==N2ItrE){//When opposite perim is only one.
		container_type& neibrsOpo=(**N2Itr).m_neighbours[periOpo];
		//Search this pointer in rect2->m_neighbours[peri2opo].
		mgTLRect::RecNItr opoItr=neibrsOpo.begin(), opoItrE=neibrsOpo.end();
		for(;opoItr!=opoItrE; opoItr++)
			if(*opoItr==this){*opoItr=rect2; break;}
	}else{//When opposite perims are more than one.
		mgTLRect* nrect;
		for(;N2Itr!=N2ItrE; N2Itr++){
			nrect=*N2Itr;
			*(nrect->m_neighbours[periOpo].rbegin())=rect2;
		}
		nrect=*N2Itr;
		*(nrect->m_neighbours[periOpo].begin())=rect2;
	}
}

//Test if this rect has perimeter boundary, check the perimeter boundary
//is within the curve tolerance error.
bool mgTLRect::checkPerimeterBoundary(
	mgTLPoints* points, double crvTol, const MGSurface &srf, bool &bSubDivU){
	double umid=(m_u0+m_u1)*0.5, vmid=(m_v0+m_v1)*0.5;
	bool bPeri[4] = {
		m_v0 == srf.param_s_v(), m_u1 == srf.param_e_u(),
		m_v1 == srf.param_e_v(), m_u0 == srf.param_s_u()};
	for(int i = 0; i < 4; i++){
		//check if is perimeter boundary
		if(!bPeri[i]) continue;	//面のペリメータかどうかのチェック
		int next = i+1;
		if(next == 4) next = 0;
		MGPosition uv1, uv2;
		if(!Pin(i) && !Pin(next)) continue;	//ペリメータが有効かどうか
		MGPosition posMidRect = (srf.eval(uv(i)) + srf.eval(uv(next))) / 2.0;
		MGPosition posMidParam(2);
		if(i == 0 || i == 2){
			bSubDivU = true;
			posMidParam(0) = umid;
			posMidParam(1) = (i == 0) ? m_v0 : m_v1;
		}else{
			bSubDivU = false;
			posMidParam(1) = vmid;
			posMidParam(0) = (i == 3) ? m_u0 : m_u1;
		}
		MGPosition posMid = srf.eval(posMidParam);
		if((posMidRect - posMid).len() > crvTol){
			return false;
		}
	}
	return true;
}

//compare which is greater, t1 or t2 at the perimeter perim.
//latter position around the anti-clockwise round of the rectangle is greater.
//If t1 is greater or equal to t2 return true.
bool mgTLRect::compare(size_t perim,double t1,double t2)const{
	if(perim<=1) return t1>=t2;
	return t1<=t2;
}

//Compute the closeset plnae of this rect.
MGPlane mgTLRect::compute_plane(const MGSurface& surf)const{
	MGPlane plane;
	surf.get_approximate_plane(m_u0,m_u1,m_v0,m_v1,plane);
	return plane;
}

//矩形交点からトリムポリゴンを作成し polygons に追加する
//トリムされている場合の処理を行う。頂点は mgTLPoints に追加する
void mgTLRect::createPolygon(
	mgTLparameter& param,
	mgTLTriangles& polygons
){
	//ファンの作成を行う
	//rectからポリゴンを作成する
	if(status()==MGRECT_ON){
		//トリムされている場合の処理
		createTrimPolygon(param,polygons);
//			if(rect.vmax()<.032)
//				(polygons.m_triangles.back())->print(cout,tlpoints);//////////////
	}else{
		//トリムされていない場合の処理
		createNonTrimPolygon(polygons);
	}
}

//矩形からポリゴンを作成し polygon を作成し、polygons に追加する
//トリムされていない場合の処理を行う
void mgTLRect::createNonTrimPolygon(mgTLTriangles& polygons){
	mgTLTriangle& polygon = *(new mgTLTriangle(this));
	//各辺毎に頂点を調べていく
	for(int peri = 0; peri < 4; peri++){
		//頂点をポリゴンに追加する
		polygon.push_back(Pid(peri));

		//隣のrectの頂点を追加する
		std::vector<size_t> vecPolygon;
		getNeighborPoint(peri, vecPolygon);
		if(peri == 0 || peri == 1){	//辺番号0,1は順番そのままでよい
			std::vector<size_t>::iterator iter = vecPolygon.begin();
			for(; iter != vecPolygon.end(); iter++){
				polygon.push_back(*iter);
			}
		}else{	//辺番号2,3は頂点インデックスの順番を逆にする
			std::vector<size_t>::reverse_iterator iter = vecPolygon.rbegin();
			for(; iter != vecPolygon.rend(); iter++){
				polygon.push_back(*iter);
			}
		}
	}
	polygons.push_back(&polygon);
}

//矩形交点からトリムポリゴンを作成し polygons に追加する
//トリムされている場合の処理を行う。頂点は mgTLPoints に追加する
void mgTLRect::createTrimPolygon(
	mgTLparameter& param,
	mgTLTriangles& polygons
){
	double crvTol=param.get_tess_crvError();
	mgTLPoints& tlpoints=*(param.tlpoints());

//	double u0=umin(), u1=umax(), v0=vmin(), v1=vmax();////////////*****
//	if(v1<.032)///////////////************
		//cout<<endl<<(*this);///////////////************
	mgTLRecisects& iscts=isects();
	size_t trimNum=iscts.size();//トリム数iscts.size()は偶数になっている
	if(trimNum<=1){
		if(!all_corners_in()) return;
		polygons.push_back(new mgTLTriangle(this,4,Pid()));
		return;
	}

	//トリム開始点よりポリゴンを求める
	for(size_t i = 0; i < trimNum; i+=2){
		//トリム多角形のポリゴンを求めてpolygonに追加する
		mgTLRecisect& rectiSS= iscts[i];
		int periS = rectiSS.perim();//Save the 1st trim point's perimeter number.
		if(periS<0) continue;		//If this rectiSS was processed already.

		double sS = rectiSS.is()->t();//Save the 1st trim point's parameter value.
		if(periS%2){
			if(sS<m_v0) sS=m_v0;
			if(sS>m_v1) sS=m_v1;
		}else{
			if(sS<m_u0) sS=m_u0;
			if(sS>m_u1) sS=m_u1;
		}
		size_t iNow=i;
		mgTLTriangle *polygon = new mgTLTriangle(this);
boundary_again:
		add_boundary_vertices(crvTol, iNow, *polygon, tlpoints);
		iscts[iNow].set_processed();

		//トリム多角形の終点から隣接rectの頂点があるかどうか調べる
		const mgTLRecisect& rectiE = iscts[iNow+1];
		int peri = rectiE.perim();	//終点のペリメータ番号を取得する
		double s = rectiE.is()->t();	//終点のパラメータ値を取得する
		size_t iNxt;
		if(NextStartIsect(peri,s,periS,sS,iNxt)){
			double e=iscts[iNxt].is()->t();
			add_neighbor_vertices(peri,s,e,*polygon);
			iNow=iNxt;
			goto boundary_again;
				//to repeat the process trim lines at the same perimeter.
		}

no_more_isect_at_peri:
		if(peri==periS && compare(peri,sS,s)){
		//When sS is greater than s
			double e=sS;
			add_neighbor_vertices(peri,s,e,*polygon);
			polygons.push_back(polygon);
			//polygon->print(std::cout,tlpoints);//cout
			continue;
		}else{
			double e=end_point(peri);
			add_neighbor_vertices(peri,s,e,*polygon);
			peri=(peri+1)%4;
			polygon->push_back(Pid(peri));
		}

		s=start_point(peri);
		if(NextStartIsect(peri,s,periS,sS,iNxt)){
			double e=iscts[iNxt].is()->t();
			add_neighbor_vertices(peri,s,e,*polygon);
			iNow=iNxt;
			goto boundary_again;
				//to repeat the process trim lines at the same perimeter.
		}else goto no_more_isect_at_peri;
	}
	
}

//divide neighbours at the parameter value t of the perimeter num peri,
//and put the latter part to rect2.
void mgTLRect::divide_neighbours(size_t peri, double t, mgTLRect& rect2){
	container_type& neibrs1=m_neighbours[peri];
	RecNItr N1Itr;
	bool t_is_on=find_neighbour(peri,t,N1Itr);
		//Find the 1st neighbour whose minimum value of the perimeter
		//is greater or equal to t.

	mgTLRect::RecNItr N1ItrE=neibrs1.end();
	if(N1Itr != N1ItrE){
		if(t_is_on) change_neighbours(&rect2,peri,N1Itr);
		else{
			mgTLRect* neibr=*N1Itr; N1Itr++;//cout<<*neibr;
			change_neighbours(&rect2,peri,N1Itr);
			rect2.m_neighbours[peri].push_front(neibr);
			size_t peri2=(peri+2)%4;
			neibr->add_neighbour(peri2,this,&rect2);//cout<<*neibr;
				//Add rect2 after neibour sequence of neibr's perimeter peri2.
		}
	}
}

//Get the end point parameter value of the perimeter peri.
//End point means in the order of the rectangle' anti-cloclwise order.
double mgTLRect::end_point(size_t peri)const{
	switch(peri){
	case 0: return umax();
	case 1: return vmax();
	case 2: return umin();
	default: return vmin();
	}
}

//Find the 1st neighbour iterator i whose minimum value at the perimeter peri
//is less or equal to t and the maximum value is greater than t.
//Function's return value is true if t is equal to minimum value of the perimeter.
bool mgTLRect::find_neighbour(size_t peri, double t, RecNItr& i){
	i=m_neighbours[peri].begin();
	RecNItr iE=m_neighbours[peri].end();
	bool t_isU=(peri%2)==0, t_is_on=false;
	double min, max;
	for(;i!=iE; i++){
		if(t_isU){
			min=(**i).umin(); max=(**i).umax();
		}else{
			min=(**i).vmin(); max=(**i).vmax();
		}
		if(min<=t && t<max){ t_is_on=t==min; break;}
	}
	return t_is_on;
}

//Find the neighbor rect or face of rect_next following the line of fpoints.
//Function's return value is id of hhis to process the next face(mgTLData).
//When hhis moves to a different face, rect_next=0 and the proper next hhis id will
//be output.
int mgTLRect::find_neighbor_rect(
	double errorSurf,//surface parameter space error.
	int id_hhis,	//id of hhis that defines the current priority line of uv is input.
					//When id_hhis>=hhis.size(), no priority lines are input.
	const MGHHisect& hhis,
	int perim_into,	//The perimeter number(of this) that came into this rect.
					//perim_into=-1 when no previous rect.
	mgTLRect*& rect_next,//next rect will be output.
						//When rect_next=0 is output, no next in this face.
	int& perim_going_out//rect_next's perimeter number.
){
	perim_going_out=-1;//number of perimeter of no next perimeter.
	rect_next=0;

	int next_id_hhis=id_hhis;
	int nhhi=hhis.num_of_uvline();
	if(id_hhis>=nhhi)
		return nhhi;

	while(includes2(hhis.uvline1(id_hhis).uvline().end_point())){
		//If the end point of uvlinei is inside this rect.
		next_id_hhis=id_hhis+1;
		if(next_id_hhis>=nhhi)
			return next_id_hhis;

		if(hhis.uvline1(id_hhis).face()!=hhis.uvline1(next_id_hhis).face())
			return next_id_hhis;
		id_hhis=next_id_hhis;
	}

	const MGCurve& uvlinei=hhis.uvline1(next_id_hhis).uvline();
	//std::cout<<uvlinei<<std::endl;///////////********
	MGPosition uvE=uvlinei.end_point(), uvS=uvlinei.start_point();
	double errorSave=MGTolerance::set_wc_zero(errorSurf);
	for(int i=0; i<4; i++){
		if(i==perim_into)
			continue;

		MGPosition P0=uv(i), P1=uv((i+1)%4);
		MGStraight perimSL(P1,P0);
		MGCCisect_list islist=uvlinei.isect(perimSL);
		size_t nis=islist.size();
		if(nis==0 || nis>1)
			continue;//When intersection is more than 1, line goes in and out.

		const MGPosition& isuv=islist.front().point();
		if(isuv==uvS)
			continue;
		if(isuv==uvE){
			int next_id_hhis2=next_id_hhis+1;
			if(next_id_hhis2>=nhhi){
				MGTolerance::set_wc_zero(errorSave);//Restore.
				return nhhi;
			}
			if(hhis.uvline1(next_id_hhis).face()!=hhis.uvline1(next_id_hhis2).face()){
				MGTolerance::set_wc_zero(errorSave);//Restore.
				return next_id_hhis2;
			}
			next_id_hhis=next_id_hhis2;
		}

		mgTLRect::RecNItr ineibor_rect;
		double t=isuv[i%2];
		find_neighbour(i,t,ineibor_rect);
		if(ineibor_rect!=m_neighbours[i].end()){
			rect_next=*ineibor_rect;
			perim_going_out=opposite_perim(i);
			break;
		}
	}

	MGTolerance::set_wc_zero(errorSave);//Restore.
	return next_id_hhis;
}

//Get the point id of uv.
//get_point_id will() tests if uv is equal to a corner data of this rect, and
//if the same, the corner point id will be returned, else uv will be
//added to Tlpoints.
size_t mgTLRect::get_point_id(
	const MGPosition& uv,	//parameter value to get the texture, which
			//may not be on the perimeter i.
	mgTLPoints& Tlpoints
)const{
	int cid=get_corner_id(uv);
	if(cid<0)
		return Tlpoints.add(uv);
	return m_Pid[cid];
}

//Compute texture coordinates of the point world coord point xyz.
//This is valid only after four corner points are textured.
MGPosition mgTLRect::get_tex_coord(
	const mgTLData& tldata,
	const MGPosition& uv	//parameter value to get the texture, which
			//may not be on the perimeter i.
)const{
/*
	int pnum=is_on_perimeter(uv);
	if(pnum>=0){//If on a perimeter
		return compute_perimeter_tex_coord(tldata,size_t(pnum),uv);
	}

	const MGPosition& st0=tldata.get_texcoord(m_Pid[0]);
	const MGPosition& st1=tldata.get_texcoord(m_Pid[1]);
	const MGPosition& st2=tldata.get_texcoord(m_Pid[2]);
	const MGPosition& st3=tldata.get_texcoord(m_Pid[3]);

	double a=(uv[0]-m_u0)/uspan(), b=(uv[1]-m_v0)/vspan();
	MGPosition st01=st0+(st1-st0)*a;
	MGPosition st32=st3+(st2-st3)*a;
	return st01+(st32-st01)*b;
*/

	MGPlane plane=compute_plane(tldata.surface());
	const MGVector& nrml=plane.normal();

	const MGPosition& st0=tldata.get_texcoord(m_Pid[0]);
	assert(!st0.is_null());
	const MGPosition& uv0=tldata.get_uv(m_Pid[0]);

	const MGPosition& st2=tldata.get_texcoord(m_Pid[2]);
	assert(!st2.is_null());
	const MGPosition& uv2=tldata.get_uv(m_Pid[2]);

	MGPosition st;
	tldata.compute_uv_texture_coordinate(nrml,uv,uv0,st0,uv2,st2,st);
	return st;

}

//Compute trimming points(intersections with trimming loops with the rectangle)
//at the newly generated adjacent perimeters when old rectangle is subdivided
//into half at the parameter t(u or v according to the coordinate kind kcod).
//Trimming points will be recomputed and stored in this and rect2's m_Recisects.
void mgTLRect::get_trim_point(
	mgTLparameter& param,	//Tessellation parameter
	double t,				//Subdividing parameter value u or v.
	size_t kcod,			//t's coordinate kind. 0:u, 1:v.
	double s0, double s1,	//parameter range of the other coordinate than t.
	mgTLRect& rect2			//2nd rectangle after subdivided.
){
	//Devide trim points into the upper and lower parts
	//and put them into istemp1,2.
	mgTLRecisects istemp1, &istemp2=rect2.m_Recisects;
	mgTLRecisects::ReciItr j=m_Recisects.begin(), je=m_Recisects.end();
	bool uDir= kcod==0;
	for(; j!=je; j++){
		if(j->less_than(uDir,t)){
			istemp1.push_back(*j);
		}else{
			istemp2.push_back(*j);
		}
	}

	//Get intersections with a perimeter and trimming loops.
	mgTLisectsList::issListItr i;
	if(kcod) i=param.get_vList().get(param,t,kcod);
	else     i=param.get_uList().get(param,t,kcod);
	mgTLisects::iterator i0,i1, i0save,i1save;//cout<<(*i)<<endl;//////
	i->locate(s0,s1,i0,i1);
		//i0 and i1 are new trim points range in i(mgTLisects of t).

	//Put the new trim points into istemp1,2.
	i0save=i0; i1save=i1;
	if(!kcod){//When dividing u direction
		while(i1!=i0){ i1--; istemp2.push_back(mgTLRecisect(3,i1));}
		//Build this m_Recisects.
		while(i0!=i1save){ istemp1.push_back(mgTLRecisect(1,i0++));}
	}else{//When dividing v direction
		//Build rect2's m_Recisects.
		while(i0!=i1){ istemp2.push_back(mgTLRecisect(0,i0++));}
		//Build this m_Recisects.
		while(i1!=i0save){ i1--; istemp1.push_back(mgTLRecisect(2,i1));}
	}
	m_Recisects=istemp1;
}

//Make the uv included in this rectangle's parameter space.
void mgTLRect::make_included(MGPosition& uv)const{
	if(uv[0]<m_u0)
		uv(0)=m_u0;
	if(uv[0]>m_u1)
		uv(0)=m_u1;
	if(uv[1]<m_v0)
		uv(1)=m_v0;
	if(uv[1]>m_v1)
		uv(1)=m_v1;
}

//ペリメタ番号periとパラメータ値pを与えて次のスタート交点を求める
//求まったときtrueが返却されindexに交点インデックスが入る
bool mgTLRect::NextStartIsect(
	size_t peri,//perimeter number to get.
	double p,	//last isect's parameter value.
	size_t periS,//polygon's starting perimeter number and 
	double sS,	//the parameter value.
	size_t& index
)const{
	bool bGet = false;	//求まったかどうかのフラグ
	size_t n=m_Recisects.size();
	double spre = 0.0;
	for(size_t i=0; i<n; i+=2){//交点の開始点だけを調べるので偶数のみを調べる
	if(m_Recisects[i].perim() == peri){

	double s=m_Recisects[i].is()->t();
	if(peri<=1){
		if(s>=p){
			//初めて求まったときの処理
			if(!bGet){
				index=i; spre=s; bGet=true;
			}else{ //２回目以降のときの処理
				if(s<spre){	//今までよりも小さかったら新しく求まったのを使用する
					index=i; spre=s;
				}
			}
		}
	}else{
		if(s<=p){
			//初めて求まったときの処理
			if(!bGet){
				index=i; spre=s; bGet=true;
			}else{ //２回目以降のときの処理
				if(s>spre){ //今までよりも大きかったら新しく求まったのを使用する
					index=i; spre=s;
				}
			}
		}
	}

	}
	}
	if(bGet){
		if(peri!=periS) return true;
		else if(compare(peri,p,sS)) return true;
		return compare(peri,sS,spre);
	}
	return false;
}

//rectがストリップになる条件を持っているかどうか調べる
//条件 ・隣り合うrectが1以下 ・トリムされていない ・triangulationされていない
bool mgTLRect::hasStripCondition()const{
	if(status() != MGRECT_IN || is_triangled())
		return false;
	//隣接するRectが各辺１以下かどうか調べる
	for(size_t i=0; i<4; i++){
		const std::list<mgTLRect*>& rectNeighbour = neighbours(i);
		if(rectNeighbour.size()>1)
			return false;
	}
	return true;
}

//Test if this rect has trimming intersection or not on the perimeter peri.
bool mgTLRect::has_trim(size_t peri)const{
	mgTLRecisects::CReciItr i=m_Recisects.begin(), ie=m_Recisects.end();
	for(; i!=ie; i++){
		if(peri==i->perim()) return true;
	}
	return false;
}

//Test if this rectangle includes the parameter value uv.
//Returns true if includes uv, false if not.
bool mgTLRect::includes(const MGPosition& uv)const{
	if(uv[0]<m_u0) return false;
	if(uv[0]>m_u1) return false;
	if(uv[1]<m_v0) return false;
	if(uv[1]>m_v1) return false;
	return true;
}

//Test if this rectangle includes the parameter value uv.
//Returns true if includes uv, false if not.
bool mgTLRect::includes2(const MGPosition& uv)const{
	double tol=MGTolerance::rc_zero()*2.;
	double uerror=uspan()*tol;
	if(uv[0]<m_u0-uerror) return false;
	if(uv[0]>m_u1+uerror) return false;
	double verror=vspan()*tol;
	if(uv[1]<m_v0-verror) return false;
	if(uv[1]>m_v1+verror) return false;
	return true;
}

//test if the parameter uv is on a corner of this rect.
//returned is:
//=-1: when not on any corner,
//>=0: when on a corner, the corner id will be returned.
int mgTLRect::get_corner_id(const MGPosition& uv)const{
	double error=MGTolerance::rc_zero();
	double utol=uspan()*error;
	double vtol=vspan()*error;

	if(fabs(uv[0]-m_u0)<=utol){
		if(fabs(uv[1]-m_v0)<=vtol)
			return 0;
		else if(fabs(uv[1]-m_v1)<=vtol)
			return 3;
	}

	if(fabs(uv[0]-m_u1)<=utol){
		if(fabs(uv[1]-m_v0)<=vtol)
			return 1;
		else if(fabs(uv[1]-m_v1)<=vtol)
			return 2;
	}
	return -1;
}

//test if the parameter uv is on a perimeter of this rect.
//returned is:
//=-1: when not on any perimeter,
//>=0: when on a perimeter, perimeter id will be returned.
int mgTLRect::is_on_perimeter(const MGPosition& uv)const{
	double error=MGTolerance::rc_zero();
	double utol=uspan()*error;
	if(fabs(uv[0]-m_u0)<=utol)
		return 3;
	if(fabs(uv[0]-m_u1)<=utol)
		return 1;

	double vtol=vspan()*error;
	if(fabs(uv[1]-m_v0)<=vtol)
		return 0;
	if(fabs(uv[1]-m_v1)<=vtol)
		return 2;

	return -1;
}

//Normalize trimming points.
//Re-order the recisect seq 
void mgTLRect::normalize_trim_points(){
	size_t n=m_Recisects.size(); if(!n) return;
	n/=2; if(!n){m_Recisects.clear(); return;}
	mgTLRecisects recis(2*n);

	size_t id=0;
	for(size_t i=0; i<n; i++){
		//Find minimum(in the same loop order) going_in reci inRi.
		mgTLRecisects::ReciItr
			inRip,outRi,j=m_Recisects.begin(),je=m_Recisects.end();
		inRip=j;
		for(; j!=je; j++){if((*j)<(*inRip)) inRip=j;}
		if(inRip->going_out()) break;
		mgTLRecisect inRi=*inRip;
		m_Recisects.erase(inRip);

		//Find going out reci next to inRi.
		mgTLRecisects::ReciItr outRiSv;
		const MGLEPoint& inRiLE=inRi.is()->isect();
		const MGLEPoint* outRiLE, *outRiLEsv;
		outRiSv=outRi=je=m_Recisects.end();
		for(j=m_Recisects.begin(); j!=je; j++){
			const MGLEPoint& jLE=j->is()->isect();
			if(j->same_loop(inRi)){
				if(j->going_out()){
					if(jLE>inRiLE){
						if(outRi==je) {outRi=j; outRiLE=&jLE;}
						else if(jLE<*outRiLE) {outRi=j; outRiLE=&jLE;}
					}else{
						if(outRiSv==je) {outRiSv=j; outRiLEsv=&jLE;}
						else if(jLE<*outRiLEsv) {outRiSv=j; outRiLEsv=&jLE;}
					}
				}
			}
		}
		if(outRi==je){//When outRiLE that is greater than inRiLE was not found.
			outRi=outRiSv;
		}
		if(outRi!=je){//When outRi found.
			recis[id*2]=inRi;
			recis[id*2+1]=*outRi; m_Recisects.erase(outRi);
			id++;
		}
	}
	if(id!=n) recis.resize(id*2);
	m_Recisects=recis;
}

//Check if this rect's perimeter(that is rect's v=min,max, when along_u is true
//or u=min,max parameter line, when aling_u is false) is within error .
bool mgTLRect::perim_within_tol(
	bool along_u,
	const MGSurface& srf,
	double error2		//square of error allowed
)const{
	double t0, tm ,t1, s1;
	MGPosition uv0(2), uv1(2), uv2(2);
	size_t id;
	if(along_u){
		uv0(1)=uv1(1)=uv2(1)=m_v0; s1=m_v1;
		t0=m_u0; t1=m_u1;
		id=0;
	}else{
		uv0(0)=uv1(0)=uv2(0)=m_u0; s1=m_u1;
		t0=m_v0; t1=m_v1;
		id=1;
	}
	tm=(t0+t1)*0.5;
	uv0(id)=t0; uv1(id)=tm; uv2(id)=t1;
	MGVector P0=srf.eval(uv0);
	MGVector P1=srf.eval(uv1);
	MGVector P2=srf.eval(uv2);

	MGVector P01(P1-P0), P02(P2-P0);
	double a=P01%P02, len1= P01%P01, len2=P02%P02;
	if((len1*len2-a*a)>len2*error2) return false;

	id=(id+1)%2;
	uv0(id)=uv1(id)=uv2(id)=s1;
	P0=srf.eval(uv0);
	P1=srf.eval(uv1);
	P2=srf.eval(uv2);

	P01=MGVector(P1-P0); P02=MGVector(P2-P0);
	a=P01%P02; len1=P01%P01; len2=P02%P02;
	if((len1*len2-a*a)>len2*error2) return false;
	return true;
}

void mgTLRect::print_neighbours(ostream& out) const{
	out<<this<<" perim";
	for(size_t i=0; i<4; i++){
		const container_type& neibrs=m_neighbours[i];
		mgTLRect::CRecNItr jj=neibrs.begin(), jE=neibrs.end();
		out<<i<<"=";
		for(; jj!=jE; jj++) out<<(*jj)<<" ";
		if(i!=3) out<<",";
	}
	out<<endl;
}

void mgTLRect::print_triangles(std::ostream& out, mgTLData& tldata) const{
	out<<"mgTLRect="<<this<<" Triangles::"<<endl;
	int idtri=m_triangle;
	if(idtri==-1)
		return;
	mgTLTriangles& tris=tldata.triangles();
	int ntri=tris.size();
	do{
		mgTLTriangle& tri=tris[idtri];
		tri.print(out,tldata);
		idtri++;
	}while(idtri<ntri && tris[idtri].rect()==this);
	out<<endl;
}

//Compute ratio square of du/dv of the rectangle.
double mgTLRect::ratio_sqr(const mgTLparameter& param)const{
	double u0=umin(), u1=umax(), v0=vmin(), v1=vmax();
	const MGSurface& f=param.get_surface();
	MGVector f00=f.eval(u0,v0);
	MGVector f1000=f.eval(u1,v0)-f00;
	MGVector f0100=f.eval(u0,v1)-f00;
	double ulen=f1000%f1000, vlen=f0100%f0100;
	double error=param.get_ratio_error();
	if(ulen<=error || vlen<=error)
		return 1.;
	return ulen/vlen;
}

//Obtain previous or aft ReciItr of i in the ReciItr loop of this
//rectangle. Input i must not be end().
mgTLRecisects::ReciItr mgTLRect::reci_aft(mgTLRecisects::ReciItr i){
	i++;
	if(i==m_Recisects.end()) i=m_Recisects.begin();
	return i;
}
mgTLRecisects::ReciItr mgTLRect::reci_pre(mgTLRecisects::ReciItr i){
	if(i==m_Recisects.begin()){
		i=m_Recisects.begin()+m_Recisects.size()-1;
	}else{
		i--;
	}
	return i;
}

//Set the new m_status of mgTLRect after subdivided.
//m_box, m_Rectisects, and m_Pin[0] must be updated before use.
//These are referenced in set_status(). old m_status is referenced also.
void mgTLRect::set_status(
	const MGFace& f	//Original face.
	){
	if(m_status==MGRECT_OVER){
		if(has_trim()){
			if(f.hasLoop(MGBox(MGInterval(m_u0,m_u1), MGInterval(m_v0,m_v1))))
				m_status=MGRECT_ONANDOVER;
			else m_status=MGRECT_ON;
		}else{
			if(!f.hasLoop(MGBox(MGInterval(m_u0,m_u1), MGInterval(m_v0,m_v1)))){
				if(m_Pin[0]) m_status=MGRECT_IN;
				else m_status=MGRECT_OUT;
			}
		}
	}else if(m_status==MGRECT_ON){
		if(has_trim()){
			m_status=MGRECT_ON;
		}else{
			if(!m_Pin[0] || !m_Pin[1] || !m_Pin[2] || !m_Pin[3]) m_status=MGRECT_OUT;
			else m_status=MGRECT_IN;
		}
	}else if(m_status==MGRECT_ONANDOVER){
		if(has_trim()){
			if(f.hasLoop(MGBox(MGInterval(m_u0,m_u1), MGInterval(m_v0,m_v1))))
				m_status=MGRECT_ONANDOVER;
			else m_status=MGRECT_ON;
		}else if(m_Pin[0]){
			if(f.hasLoop(MGBox(MGInterval(m_u0,m_u1), MGInterval(m_v0,m_v1))))
				m_status=MGRECT_OVER;
			else m_status=MGRECT_IN;
		}else m_status=MGRECT_OUT;
	}
	//When status() was MGRECT_IN, new status is also MGRECT_IN.

	return;
}

//Assuming this is already textured, 
//set all the neighboring rect's texture coordinates data
//invoking compute_texture_by_triangles().
bool mgTLRect::set_neighbor_tex_coord_data(
	mgTLData& tldata,
	std::deque<mgTLRect*>& rects	//modified rects will be prepended.
){
	//if(m_u1<.251&& m_v1<.5626 &&m_v0>.499)
	//	m_u1=m_u1;//****************************
	bool textured=false;
	int level=m_textured+1;
	for(size_t i=0; i<4; i++){
		std::list<mgTLRect*> neibors=neighbours(i);
		std::list<mgTLRect*>::iterator j=neibors.begin(), je=neibors.end();
		for(; j!=je; j++){
			mgTLRect* rectj=*j;
			if(rectj->status()==MGRECT_OUT)
				continue;
			if(rectj->is_textured())
				continue;
			rectj->compute_texture_from_neighbor(tldata,*this,i);
			//rectj->compute_texture_by_triangles(tldata);
			rectj->set_texture_level(level);
			rects.push_front(rectj);
			textured=true;
		}
	}
	return textured;
}

//Set all the texture coordinates data of this rect's vertex.
void mgTLRect::set_initial_tex_from_1point(
	const MGPosition& xyz1,	//parameter value of the surfaceof texture st1.
	const MGPosition& st1,	//texture coordinates of uv1
	const MGUnit_vector& xaxis,	//texture x-axis.
	mgTLData& tldata
){
	//define
	MGPlane plane=compute_plane(tldata.surface());
	const MGUnit_vector& normal=plane.normal();
	MGUnit_vector xaxis2=normal*xaxis*normal;//adjust the xaxis.
	MGPosition P1=plane.eval(plane.uv(xyz1));//normalize xyz1 so as to lie on m_plane.
	
	mgTLPoints& pointsuv=tldata.tlpoints();
	const MGSurface& surf=tldata.surface();
	for(size_t i=0; i<4; i++){
		size_t idi=m_Pid[i];
		MGPosition P2=surf.eval(pointsuv[idi]);
		P2=plane.eval(plane.uv(P2));//normalize P2 so as to lie on plane.
		MGVector P1toP2=P2-P1;
		double L=P1toP2.len();
		//L/=m_texcoord_ratio;
	
		double theta=xaxis2.angle2pai(P1toP2,normal);
		MGPosition st2=st1+MGPosition(L*cos(theta), L*sin(theta));
		//tldata.get_texcoord(idi)=st2;
		tldata.set_texcoord(idi,st2);
	}
	compute_texture_of_non_corner(tldata);

	m_textured=1;//Level 1
	//print_triangles(std::cout,tldata);///////////***************
}

//neighbor is already textured, and neighbor  and this is connected
//along neghbor's perimeter number nperim. Then
//compute_texture_from_neighbor will texture-compute this rect.
void mgTLRect::compute_texture_from_neighbor(
	mgTLData& tldata,
	mgTLRect& neighbor,
	size_t nperim
){
	MGPlane plane=compute_plane(tldata.surface());
	const MGVector& nrml=plane.normal();

	//1. Compute corner texture.
	size_t v1=nperim, v2=(nperim+1)%4;
	const MGPosition& uv1=tldata.get_uv(neighbor.m_Pid[v1]);
	const MGPosition& st1=tldata.get_texcoord(neighbor.m_Pid[v1]);
	const MGPosition& uv2=tldata.get_uv(neighbor.m_Pid[v2]);
	const MGPosition& st2=tldata.get_texcoord(neighbor.m_Pid[v2]);
	for(size_t i=0; i<4; i++){
		MGPosition& st=tldata.get_texcoord(m_Pid[i]);
		if(st.is_null()){
			MGPosition& uv=tldata.get_uv(m_Pid[i]);
			tldata.compute_uv_texture_coordinate(nrml,uv,uv1,st1,uv2,st2,st);
		}
	}

	//1. Compute non-corner texture.
	compute_texture_of_non_corner(tldata);
}

//Assume four corner vertices are textured, texture other vertices.
void mgTLRect::compute_texture_of_non_corner(
	mgTLData& tldata
){
	int idtri=m_triangle;
	if(idtri==-1)
		return;

	mgTLTriangles& tris=tldata.triangles();
	int ntri=tris.size();
	do{
		mgTLTriangle& tri=tris[idtri];
		size_t nvrtx=tri.size();
		for(size_t j=0; j<nvrtx; j++){
			size_t idj=tri[j];
			MGPosition& stj=tldata.get_texcoord(idj);
			if(stj.is_null()){
				//stj=get_tex_coord(tldata,tldata.get_uv(idj));
				tldata.set_texcoord(idj,get_tex_coord(tldata,tldata.get_uv(idj)));
			}
		}
		idtri++;
	}while(idtri<ntri && tris[idtri].rect()==this);
}

//Compute texture coordinates of the point world coord point xyz.
//This is valid only after texutured.
MGPosition mgTLRect::compute_perimeter_tex_coord(
	const mgTLData& tldata,
	size_t i,//perimeter number where uv lies on.
	const MGPosition uv	//parameter value to get the texture, which
			//must be on the perimeter i.
)const{
	size_t id1,id2;
	double span;
	double delta;
	if(i%2){
		//perimeter 1,3=v variable periemter.
		span=vspan();
		delta=uv[1]-m_v0;
		if(i==3){
			id1=0; id2=3;
		}else{
			id1=1; id2=2;
		}
	}else{
		//perimeter 0,2=u variable periemter.
		span=uspan();
		delta=uv[0]-m_u0;
		if(i==0){
			id1=0; id2=1;
		}else{
			id1=3; id2=2;
		}
	}
	const MGPosition& st1=tldata.get_texcoord(m_Pid[id1]);
	const MGPosition& st2=tldata.get_texcoord(m_Pid[id2]);
	MGPosition st=st1+(st2-st1)*(delta/span);
	return st;
}

//Assurme this is not textured yet, and at least one of the
//neighboring perimeter is textured. Then compute_texture_by_triangles
//will perform texture computation by invoking mgTriangle::texture.
void mgTLRect::compute_texture_by_triangles(
	mgTLData& tldata
){
	mgTLTriangles& tris=tldata.triangles();
	int ntris=tris.size();

	//if(m_u1<.251&& m_v0>.49&& m_v1<.5626)
	//	m_u1=m_u1;/////////****************
	int ntextured;
	do{
		ntextured=0;
		int i=m_triangle;
		do{
			mgTLTriangle& tri=tris[i];
			if(tri.texture_triangle(tldata)){
				ntextured++;
			}
			i++;
		}while(i<ntris && tris[i].rect()==this);

	}while(ntextured);
	compute_texture_of_non_corner(tldata);
	set_texture_level(10000);//10000 is unknown texture level.

	//print_triangles(std::cout,tldata);///////////***************
}

//set triangled flag.
void mgTLRect::set_triangled(int tri){
	if(m_triangle==-1)
		m_triangle=tri;
}

//Get the start point parameter value of the perimeter peri.
//Start point means in the order of the rectangle' anti-cloclwise order.
double mgTLRect::start_point(size_t peri)const{
	switch(peri){
	case 0: return umin();
	case 1: return vmin();
	case 2: return umax();
	default: return vmax();
	}
}

//add (u,v) to by add_point of rects.
size_t mgTLRect::add_rect_point(
	size_t pnum,		//perimeter number of this rect.
	double u,
	double v,
	mgTLRects& rects	//mgTLRects into which all the rects are stored.
){
	const size_t vnum[4]={2,3,1,2};
	container_type& neibrs=m_neighbours[pnum];
	if(neibrs.size()==0) 
		return rects.add_point(u,v);
	else{
		mgTLRect& neibor_rect=**(neibrs.rbegin());
		bool already_registered;
		if((pnum+1)%2)
			already_registered=neibor_rect.umax()==u;
		else
			already_registered=neibor_rect.vmax()==v;

		if(already_registered)
			return neibor_rect.Pid(vnum[pnum]);
		else 
			return rects.add_point(u,v);
	}
}

//////////////////////
//Divide into two parts at the middle of u(uDirection=true) or
//v(uDirection=false) parameter.
//The divided two are put in this and function's return value.
//The new this will be located first in the sequence of
//the uv rectangle list mgTLRects.
mgTLRect* mgTLRect::subdivide(
	mgTLRects& rects,	//mgTLRects into which all the rects are stored.
	bool uDirection		//direction of sibdivision. =true:udirection.
){
	assert(valid());
	//std::cout<<(*this);

	mgTLparameter& param=rects.parameter();
	const MGFace& f=param.get_face();
	double u0=umin(), u1=umax(), v0=vmin(), v1=vmax();
/*	if(u0>240. && u1<255.)
		u0=u0;///////cout*/

	mgTLRect* rect2=new mgTLRect(uDirection, *this);
	if(uDirection){
		double u=(u0+u1)*0.5;
//		if(u<734. && u>720.)
//			u=u;

		//Set neighbours.
//			print_neighbours(cout);///////////
		container_type& neibrs11=m_neighbours[1];
		change_neighbours(rect2,1,neibrs11.begin());
		neibrs11.push_back(rect2);
		rect2->m_neighbours[3].push_back(this);
		divide_neighbours(0,u,*rect2);
		divide_neighbours(2,u,*rect2);
//			cout<<"u::"; print_neighbours(cout);///////
//			rect2->print_neighbours(cout);cout<<endl;/////

		container_type& neibrs10=m_neighbours[0];
		size_t P0=0, P1=0;
		m_u1=u;
		if(m_status==MGRECT_IN){
			P0=add_rect_point(0,u,v0,rects);
			P1=add_rect_point(2,u,v1,rects);
			rect2->m_Pid[0]=m_Pid[1]=P0; rect2->m_Pid[3]=m_Pid[2]=P1;
			m_Pin[1]=m_Pin[2]=rect2->m_Pin[0]=rect2->m_Pin[3]=true;
			assert(valid());
			return rect2;
		}
		P0=add_rect_point(0,u,v0,rects);
		P1=add_rect_point(2,u,v1,rects);
		rect2->m_Pid[0]=m_Pid[1]=P0; rect2->m_Pid[3]=m_Pid[2]=P1;
		//Build m_Recisects of this and rect2.
		get_trim_point(param,u,0,v0,v1,*rect2);
		if(m_status==MGRECT_OVER){
			m_Pin[1]=m_Pin[2]=rect2->m_Pin[0]=rect2->m_Pin[3]=m_Pin[0];
		}else if(m_status==MGRECT_ON || m_status==MGRECT_ONANDOVER){
			if(m_Pin[0]!=m_Pin[1] || has_trim(0)){
				m_Pin[1]=rect2->m_Pin[0]=f.in_range(u,v0);
			}else{
				m_Pin[1]=rect2->m_Pin[0]=m_Pin[0];
			}
			if(m_Pin[3]!=m_Pin[2] || has_trim(2)){
				m_Pin[2]=rect2->m_Pin[3]=f.in_range(u,v1);
			}else{
				m_Pin[2]=rect2->m_Pin[3]=m_Pin[3];
			}
		}
	}else{
		double v=(v0+v1)*0.5;
		//if(140.<v && v<142. && u1<55.)
		//	v=v;//cout
		//Set neighbours.
//			print_neighbours(cout);//////////
		container_type& neibrs12=m_neighbours[2];
		change_neighbours(rect2,2,neibrs12.begin());
		neibrs12.push_back(rect2);
		rect2->m_neighbours[0].push_back(this);
		divide_neighbours(1,v,*rect2);
		divide_neighbours(3,v,*rect2);
			//cout<<"this::"; print_neighbours(cout);///////
			//cout<<"2::";rect2->print_neighbours(cout);cout<<endl;/////
		container_type& neibrs13=m_neighbours[3];
		size_t P0=0,P1=0;
		m_v1=v;
		if(m_status==MGRECT_IN){
			P0=add_rect_point(3,u0,v,rects);
			P1=add_rect_point(1,u1,v,rects);
			m_Pid[3]=rect2->m_Pid[0]=P0; m_Pid[2]=rect2->m_Pid[1]=P1;
			m_Pin[2]=m_Pin[3]=rect2->m_Pin[0]=rect2->m_Pin[1]=true;
			assert(valid()); return rect2;
		}
		P0=add_rect_point(3,u0,v,rects);
		P1=add_rect_point(1,u1,v,rects);
		m_Pid[3]=rect2->m_Pid[0]=P0; m_Pid[2]=rect2->m_Pid[1]=P1;
		//Build m_Recisects of this and rect2.
		get_trim_point(param,v,1,u0,u1,*rect2);
		if(m_status==MGRECT_OVER){
			m_Pin[2]=m_Pin[3]=rect2->m_Pin[0]=rect2->m_Pin[1]=m_Pin[0];
		}else if(m_status==MGRECT_ON || m_status==MGRECT_ONANDOVER){
			if(m_Pin[3]!=m_Pin[0] || has_trim(3)){
				m_Pin[3]=rect2->m_Pin[0]=f.in_range(u0,v);
			}else{
				m_Pin[3]=rect2->m_Pin[0]=m_Pin[0];
			}
			if(m_Pin[1]!=m_Pin[2] || has_trim(1)){
				m_Pin[2]=rect2->m_Pin[1]=f.in_range(u1,v);
			}else{
				m_Pin[2]=rect2->m_Pin[1]=m_Pin[1];
			}
		}
	}
	set_status(f); rect2->set_status(f);
	assert(valid());
	assert(rect2->valid());
	return rect2;
}

/*
//Test if the distance between m_Recisects[j1] and m_Recisects[j2] is too small
//(within the tolerance error) and we can neglect these two points.
bool mgTLRect::too_near(double error, size_t j1, size_t j2)const{
	MGPosition uv1=m_Recisects[j1].is()->isect().eval();
	MGPosition uv2=m_Recisects[j2].is()->isect().eval();
	MGVector diff=uv1-uv2; double len=diff%diff;
	return len<=error;
}
*/

//Return this rectangle's corner point i's (u,v).
MGPosition mgTLRect::uv(size_t i)const{
	switch(i){
	case 0: return MGPosition(m_u0, m_v0);
	case 1: return MGPosition(m_u1, m_v0);
	case 2: return MGPosition(m_u1, m_v1);
	case 3: return MGPosition(m_u0, m_v1);
	}
	return MGPosition(m_u0, m_v0);
}

//隣のrectの頂点を求める。頂点はパラメータの小さい順に並んでいる。
void mgTLRect::getNeighborPoint(size_t peri, std::vector<size_t>& vecPoint) const{
	const std::list<mgTLRect*>& listRectNeighbour = neighbours(peri);
	CRecNItr iter = listRectNeighbour.begin();
	double dUmax = umax(), dVmax = vmax();
	for(; iter != listRectNeighbour.end(); iter++){
		if(peri == 0 || peri == 2){
			double u = (*iter)->umax();
			if(dUmax > u){	//自身の矩形範囲内の時に、vecPointに追加する
				if(peri == 0){	//頂点がOUTの時は追加しない
					if((*iter)->Pin(2))vecPoint.push_back((*iter)->Pid(2));
				}else{
					if((*iter)->Pin(1))vecPoint.push_back((*iter)->Pid(1));
				}
			}else break;
		}else{
			double v = (*iter)->vmax();
			if(dVmax > v){	//自身の矩形範囲内の時に、vecPointに追加する
				if(peri == 1){	//頂点がOUTの時は追加しない
					if((*iter)->Pin(3))vecPoint.push_back((*iter)->Pid(3));
				}else{
					if((*iter)->Pin(2))vecPoint.push_back((*iter)->Pid(2));
				}
			}else break;
		}
	}
}

bool mgTLRect::valid()const{
	if(m_status==MGRECT_OUT)
		return true;
	if(m_status!=MGRECT_IN || 
		(m_Pin[0]==true  && m_Pin[1]==true && m_Pin[2]==true  && m_Pin[3]==true)){
		size_t n=m_Recisects.size();
		if(n<=1) return true;
		for(size_t i=0; i<n; i++){
			double s=m_Recisects[i].is()->t();
			if((m_Recisects[i].m_perim)%2){
				if(s<vmin() || s>vmax()){
					double vspn=vspan()*double(m_divnum+1);
					if(!MGRZero2(s-vmin(),vspn)&&!MGRZero2(s-vmax(),vspn)){
					cout<<" mgTLRect::valid, this="<<this<<", m_Recisects["<<i<<"].is()->t() is invalid"<<endl;
					return true;
					}
				}
			}else{
				if(s<umin() || s>umax()){
					double uspn=uspan()*double(m_divnum+1);
					if(!MGRZero2(s-umin(),uspn)&&!MGRZero2(s-umax(),uspn)){
					cout<<" mgTLRect::valid, this="<<this<<", m_Recisects["<<i<<"].is()->t() is invalid"<<endl;
					return true;
					}
				}
			}
		}
		return true;
		}else{
			cout<<" mgTLRect::valid, this="<<this<<", m_status or m_Pin[.] is invalid."<<endl;
			return true;
		}
}

ostream& operator<< (ostream& out, const mgTLRect& rect){
	double u0=rect.umin(), u1=rect.umax(), v0=rect.vmin(), v1=rect.vmax();
	out<<"TLRect="<<(&rect)<<" u("<<u0<<","<<u1<<")"
		<<"-v("<<v0<<","<<v1<<"),"
		<<"m_status=";
	char* type[6]={"UNKNOWN","IN","OUT","ON","OVER","ONANDOVER"};
	out<<type[rect.m_status];
	out<<",triangle="<<rect.m_triangle<<endl;
	for(size_t i=0; i<4; i++){
		out<<i<<"::Pin="<<rect.m_Pin[i];
		out<<",Pid="<<rect.m_Pid[i]<<",";
		out<<" neibours=";
		const mgTLRect::container_type& neibrs=rect.m_neighbours[i];
		mgTLRect::CRecNItr j=neibrs.begin(), jE=neibrs.end();
		for(; j!=jE;){out<<(*j++); if(j!=jE) out<<",";}
		if(i%2) out<<endl; else out<<", ";
	}
	mgTLRecisects::CReciItr itr=rect.m_Recisects.begin(),
		itre=rect.m_Recisects.end();
	for(; itr!=itre; itr++) out<<(*itr)<<endl;
	return out;
}
