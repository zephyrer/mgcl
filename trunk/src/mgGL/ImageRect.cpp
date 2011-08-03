/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include "mg/Tolerance.h"
#include "mg/Box.h"
#include "mg/Straight.h"
#include "mg/CCisect.h"
#include "mg/CParam_list.h"
#include "mgGL/ImageRect.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

mgImageRect::mgImageRect(
	const MGBox& box,//The box (left,bottom, right, top)
	const MGPosition st_tri[3],
	const MGPosition world_tri[3]
){
	MGVector v01=st_tri[1]-st_tri[0];
	MGVector v02=st_tri[2]-st_tri[0];
	if((v01*v02)[2]>=0){
		for(size_t i=0; i<3; i++){
			m_triangle[i]=st_tri+i;
			m_world_tri[i]=world_tri+i;
		}
	}else{
		m_triangle[0]=st_tri;
		m_triangle[2]=st_tri+1;
		m_triangle[1]=st_tri+2;
		m_world_tri[0]=world_tri;
		m_world_tri[2]=world_tri+1;
		m_world_tri[1]=world_tri+2;
	}
	const MGInterval& urng=box[0];
	const MGInterval& vrng=box[1];
	m_t[3]=urng.low_point(); m_t[1]=urng.high_point();
	m_t[0]=vrng.low_point(); m_t[2]=vrng.high_point();

	for(int tri_pn=0; tri_pn<3; tri_pn++){//tri_pn is the perimeter number of the triangle.
		int ntri_pn=(tri_pn+1)%3;
		const MGPosition& st=*(m_triangle[tri_pn]);
		const MGPosition& stn=*(m_triangle[ntri_pn]);
		MGStraight sl01(stn,st);
		for(int pn=0; pn<4; pn++){//pn is the perimeter number of the rectangle.
			int cid=(pn+1)%2;//coordinate id, 0:u-value, 1:v-value.
			MGCParam_list isects=sl01.isect_1D(m_t[pn],cid);
			if(isects.empty())
				continue;

			MGPosition uv(sl01.eval(isects.front()));
			double tmin, tmax;
			if(cid){
				tmin=m_t[3]; tmax=m_t[1];
			}else{
				tmin=m_t[0]; tmax=m_t[2];
			}
			double p=uv[(cid+1)%2];
			if(tmin<p && p<tmax)
				m_isects.push_back(mgIRisect(pn,p,tri_pn));			
		}
	}
	if(!m_isects.empty()){
		std::sort(m_isects.begin(), m_isects.end());
		if(in_triangle(m_t[3],m_t[0])){
			mgIRisect lastis=m_isects.back();
			m_isects.pop_back();
			m_isects.push_front(lastis);
		}
	}
	/*std::cout<<"m_isects=";
	for(size_t i=0; i<m_isects.size(); i++)	std::cout<<m_isects[i];
	std::cout<<std::endl;*/
}

//Add vertices of the triangle to stpoly from the isect i to ip1.
//When ip1==m_isects.end(), ip1 means isects.begin().
void mgImageRect::add_trim_point(
	iterator i,		//m_isects' iterator.
	iterator ip1,	//m_isects' iterator.
	std::vector<MGPosition>& stpoly
){
	if(ip1==m_isects.end())
		ip1=m_isects.begin();

	int tid0=i->tri_id(), tid1=ip1->tri_id();
	for(int tid=tid0; tid!=tid1;){
		int ntid=(tid+1)%3;
		stpoly.push_back(*(m_triangle[ntid]));
		tid=ntid;
	}
}

//Test if (u,v) is in the triangle.
bool mgImageRect::in_triangle(double u, double v)const{
	MGPosition uv(u,v);
	MGVector v01=*(m_triangle[1])-*(m_triangle[0]);
	MGVector v02=*(m_triangle[2])-*(m_triangle[0]);
	MGVector v0uv=uv-*(m_triangle[0]);
	if((v02*v0uv)%(v0uv*v01)<0.)
		return false;

	MGVector v12=*(m_triangle[2])-*(m_triangle[1]);
	MGVector v1uv=uv-*(m_triangle[1]);
	if((v12*v1uv)%(v1uv*v01)>0.)
		return false;

	return true;
}

//矩形交点からトリムポリゴンを作成し stpoly に追加する
void mgImageRect::createTrimPolygon(
	std::vector<MGPosition>& stpoly
		//The vertices of the trimmed polygon of (s,t) coordinates will be output.
){
	//std::cout<<(*this);

	size_t trimNum=m_isects.size();//トリム数m_isects.size()は偶数になっている
	if(trimNum<=1){
		if(in_triangle(m_t[3],m_t[0])){
			for(size_t iv=0; iv<4; iv++)
				stpoly.push_back(start_point(iv));//1st point.
		}
	}else{
		iterator i=m_isects.begin(), ip1, iend=m_isects.end();
		while(i!=iend){
			int periS=i->perimeter(), periE;
			stpoly.push_back(isect_uv(*i++));//1st trim point.

			periE=i->perimeter();
				//Include rect's corner point into the polygon.
			int nperi=periS;
			while(nperi!=periE){
				nperi=(nperi+1)%4;
				stpoly.push_back(start_point(nperi));//1st point.
			}
			stpoly.push_back(isect_uv(*i));//2nd trim point.
			ip1=i+1;// ip1++;
			add_trim_point(i,ip1,stpoly);
			i=ip1;
		}
	}
	//for(std::vector<MGPosition>::iterator ii=stpoly.begin(); ii!=stpoly.end(); ii++)
	//	std::cout<<(*ii);
	//std::cout<<std::endl;
}

MGPosition mgImageRect::isect_uv(const mgIRisect& is)const{
	int peri=is.perimeter();
	if(peri%2)
		return MGPosition(m_t[peri],is.t());
	else
		return MGPosition(is.t(),m_t[peri]);
}

//Get the start point parameter value (u,v) of the perimeter peri.
//Start point means in the order of the rectangle' anti-clockwise order.
MGPosition mgImageRect::start_point(int peri)const{
	int pre_peri=(peri+3)%4;
	if(peri%2)
		return MGPosition(m_t[peri],m_t[pre_peri]);
	else
		return MGPosition(m_t[pre_peri],m_t[peri]);
}

//Convert (s,t) coordinates to world coordinates using m_triangle and m_world_tri.
void mgImageRect::convert_to_world(
	const MGPosition& st,	//(s,t) coordinates
	MGPosition& world		//world coordinates
)const{
	MGVector diff;
	double error=MGTolerance::wc_zero();
	double lenmax=-1.;
	size_t i0=0, i1,i2;
	for(size_t i=0; i<3; i++){
		diff=st-*(m_triangle[i]);
		double len=fabs(diff[0])+fabs(diff[1]);
		if(len<=error){
			world=*(m_world_tri[i]);
			return;
		}
		if(lenmax<len){
			i0=i;
			lenmax=len;
		}
	}
	i1=(i0+1)%3; i2=(i0+2)%3;
	const MGPosition& st0=*(m_triangle[i0]);
	const MGPosition& st1=*(m_triangle[i1]);
	const MGPosition& st2=*(m_triangle[i2]);
	MGVector v0st=st-st0;
	MGStraight sl0st(MGSTRAIGHT_UNLIMIT,v0st,st0);
	MGVector v12=st2-st1;
	MGStraight sl12(MGSTRAIGHT_UNLIMIT,v12,st1);
	MGCCisect isect;
	sl0st.relation(sl12,isect);//Get the intersection.
	const MGPosition& st12m=isect.point();

	double l12=v12.len();
	double a=v0st.len(), b=(st12m-st).len(), c=(st12m-st1).len();
	const MGPosition& P0=*(m_world_tri[i0]);
	const MGPosition& P1=*(m_world_tri[i1]);
	const MGPosition& P2=*(m_world_tri[i2]);
	MGPosition P12m=(P2*c+P1*(l12-c))/l12;
	world=(P12m*a+P0*b)/(a+b);
}

void mgImageRect::convert_to_world(
	const std::vector<MGPosition>& stpoly,
	std::vector<MGPosition>& worldpoly
)const{
	worldpoly.clear();
	size_t n=stpoly.size();
	MGPosition P;
	for(size_t i=0; i<n; i++){
		const MGPosition& st=stpoly[i];
		convert_to_world(st,P);
		worldpoly.push_back(P);
	}
}

std::ostream& operator<<(std::ostream& out, const mgImageRect& rect){
	out<<std::endl<<"mgImageRect::rectangle(";
	out<<rect.m_t[0]<<",";out<<rect.m_t[1]<<",";out<<rect.m_t[2]<<",";
	out<<rect.m_t[3]<<")";
	out<<", m_isects::n="<<rect.m_isects.size()<<"::"<<std::endl;
	int j=0;
	mgImageRect::const_iterator i=rect.m_isects.begin(), ip1, iend=rect.m_isects.end();
	while(i!=iend){
		out<<++j<<": "<<(*i++);
	}
	return out;
}
