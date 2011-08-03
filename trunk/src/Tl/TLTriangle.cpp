/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/Surface.h"
#include "topo/Face.h"
#include "Tl/TLPoints.h"
#include "Tl/TLTriangle.h"
#include "Tl/TLTexPoints.h"
#include "Tl/TLData.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//construct a triangle whose type is mgTESTRIANG_FAN, and number of
//vertices are 3 including the center.
mgTLTriangle::mgTLTriangle(
	mgTLRect* rect,	//mgTLRect this belongs to.
	size_t center, size_t v1, size_t v2	//ÇRäpå`í∏ì_
):m_rect(rect), m_type(mgTESTRIANG_FAN), m_indices(3){
	m_indices[0]=center;
	m_indices[1]=v1;
	m_indices[2]=v2;
}

//construct a triangle whose type is mgTESTRIANG_FAN, and number of
//vertices is n.
mgTLTriangle::mgTLTriangle(
	mgTLRect* rect,//mgTLRect this belongs to.
	size_t n,			//number of vertices.
	const size_t* vertices	//vertices of array length n.
):m_rect(rect), m_type(mgTESTRIANG_FAN), m_indices(n){
	for(size_t i=0; i<n; i++) m_indices[i]=vertices[i];
}

//Erase i-th element of m_indices.
mgTLTriangle::IndexItr mgTLTriangle::erase(size_t i){
	IndexItr itr=begin()+i;
	return m_indices.erase(itr);
}

const char* type[3]={"mgTESTRIANG_UNKNOWN","mgTESTRIANG_FAN","mgTESTRIANG_STRIP"};

std::ostream& operator<< (std::ostream& out, const mgTLTriangle& triangle){
	out<<"TLTriangle="<<(&triangle)<<" num of points="<<triangle.m_indices.size()
		<<",type="<<type[triangle.getGeometryType()]<<",rect="<<triangle.m_rect
		<<":: "<<std::endl;
	mgTLTriangle::CIndexItr iter = triangle.begin();
	for(; iter != triangle.end(); iter++){out<<*iter<<",";}
	out<<std::endl;
	return out;
}

void mgTLTriangle::print(std::ostream& out, const mgTLData& tld)const{
	out<<"TLTriangle="<<this<<" num of points="<<this->m_indices.size()
		<<",type="<<type[this->getGeometryType()]<<",rect="<<m_rect
		<<":: "<<std::endl;
	mgTLTriangle::CIndexItr iter = this->begin(), iend=this->end();
	const mgTLPoints& points=tld.tlpoints();
	const mgTLTexPoints& stpoints=tld.tlTexPoints();
	const MGSurface& srf=tld.surface();
	for(; iter != iend; iter++){
		size_t id=*iter;
		const MGPosition& uv=points[id];
		out<<id<<":uv="<<uv;
		if(&stpoints) out<<",st="<<stpoints[id];
		out<<",xyz="<<srf.eval(uv)<<",";
		out<<std::endl;
	}
}

//propagate texture computation from the perimeter (0, textured_vertex)
//to (0,1) perimeter(when increase=false) or to (0,n-1) perimeter
//(when increase=true). The vertex 0 must be already textured.
void mgTLTriangle::propagate_texture(
	mgTLData& tldata,
	int textured_vertex,
	bool increase
){
	int i, n=size();
	mgTLTexPoints& sts=tldata.tlTexPoints();

	if(increase){
		
	//propagate texture computation to perimeter 0-(n-1).
	for(i=textured_vertex+1; i<n; i++){
		size_t idi=m_indices[i];
		if(!sts[idi].is_null())
			return;
		tldata.set_texcoord(idi,third_texture(tldata,increase,0,i-1,i));
		//sts[idi]=third_texture(tldata,increase,0,i-1,i);
	}

	}else{

	//propagate texture computation to perimeter 0-1.
	for(i=textured_vertex-1; i>=1; i--){
		size_t idi=m_indices[i];
		if(!sts[idi].is_null())
			return;
		tldata.set_texcoord(idi,third_texture(tldata,increase,0,i+1,i));
		//sts[idi]=third_texture(tldata,increase,0,i+1,i);
	}

	}

}

//Assumed this type is mgTESTRIANG_FAN and at least one of the
//perimeter is textued, compute all the other perimeter's texture.
//Function's return value is true: if som of the triangles are textured,
//false if no triangles are textured.
bool mgTLTriangle::texture_triangle(
	mgTLData& tldata
){
	assert(m_type==mgTESTRIANG_FAN);

	const MGSurface& surf=tldata.surface();
	mgTLPoints& uvs=tldata.tlpoints();
	mgTLTexPoints& sts=tldata.tlTexPoints();
	int n=size();
	assert(n>=3);
	int nm1=n-1;

	bool textured=false;
//1. Texture node 0 if not textured yet.
	int id0=m_indices[0], id1=m_indices[1];
	if(sts[id0].is_null()){
		//find the most adquate perimeter to propagate the texture.
		MGPosition Pi, Pim1;
		if(!sts[id1].is_null())
			Pim1=surf.eval(uvs[id1]);

		double len=-1.;
		int isave=-1;
		for(int i=2; i<n; i++){
			int idi=m_indices[i];
			if(sts[idi].is_null())
				Pi.set_null();
			else
				Pi=surf.eval(uvs[idi]);
			if(!Pi.is_null() && !Pim1.is_null()){
				MGVector V=Pi-Pim1;
				double len1=V%V;
				if(len1>len){
					isave=i;
					len=len1;
				}
			}
			Pim1=Pi;
		}
		if(isave==-1)
			return false;

		//Set the texture of id0.
		tldata.set_texcoord(id0,third_texture(tldata,true,isave-1,isave,0));
		//sts[id0]=third_texture(tldata,true,isave-1,isave,0);
		textured=true;
	}

	MGPosition P0=surf.eval(uvs[id0]);
//2. find the vertex range that are not textured.
	int i1, i2;//from i1 to i2-1 will be the vertex range to texture.
	for(i1=1; i1<n; i1++){
		int idi1=m_indices[i1];
		if(!sts[idi1].is_null())
			continue;

		for(i2=i1+1;i2<n; i2++){
			size_t idi2=m_indices[i2];
			if(!sts[idi2].is_null())
				break;
		}
		if(i2>=n){
		//the case that no textured vertices found until the end.
			if(i1==1)
				return false;//no textured vertices found from the begining to the end.
			propagate_texture(tldata,i1-1,true);
			return true;
		}

		//i2 is textured here.
		if(i1==1){
			propagate_texture(tldata,i2,false);//texture from i2-1 to 1
			return true;
		}

	//Here, from i1 to i2-1 is not textured and i1-1 and i2 are textured.
		MGPosition Pi1m1=surf.eval(uvs[m_indices[i1-1]]);
		MGPosition Pi2=surf.eval(uvs[m_indices[i2]]);
		MGVector V1=Pi1m1-P0;
		MGVector V2=Pi2-P0;
		if(V1%V1>=V2%V2){
			propagate_texture(tldata,i1-1,true);
		}else{
			propagate_texture(tldata,i2,false);
		}
		i1=i2;
	}
	return textured;
}

//Compute 3rd texture, given ids of the points.
MGPosition mgTLTriangle::third_texture(
	mgTLData& tldata,
	bool right_hand,
	size_t i1,	//texture-known point 1.
	size_t i2,	//texture-known point 2.
	size_t i	//point to compute texture.
)const{
	const MGSurface& surf=tldata.surface();
	size_t id1=m_indices[i1], id2=m_indices[i2], id=m_indices[i];
	const MGPosition& st1=tldata.get_texcoord(id1);
	const MGPosition& st2=tldata.get_texcoord(id2);
	assert(!st1.is_null() && !st2.is_null());

	MGPosition uv=tldata.get_uv(id), uv1=tldata.get_uv(id1), uv2=tldata.get_uv(id2);
	const MGPosition& xyz1=surf.eval(uv1);
	const MGPosition& xyz2=surf.eval(uv2);
	const MGPosition& xyz=surf.eval(uv);

	MGVector ST12=st2-st1;
	MGUnit_vector Xst(ST12);
	MGUnit_vector Yst=MGVector(-Xst[1],Xst[0]);

	MGVector P1toP2=xyz2-xyz1;
	MGVector P1toP=xyz-xyz1;
	double L12=P1toP2.len();
	double stlen=ST12.len();
	double L=P1toP.len();
	if(tldata.distortion_record_required()){
		double ratio;
		if(stlen>L12) ratio=stlen/L12;
		else ratio=L12/stlen;
		if(tldata.tex_distorted(ratio)){
			MGPosition uvmid=(uv+uv1+uv2)/3.;
			tldata.record_tex_distortion(MGPosition(uvmid[0],uvmid[1],ratio));
		}
	}
	
	double theta=P1toP2.angle(P1toP);
	if(!right_hand)
		theta*=-1.;
	return st1+L*cos(theta)*Xst+L*sin(theta)*Yst;

}

//Get parameter (u,v) of the surface from the id of this triangle.
const MGPosition& mgTLTriangle::uv(size_t id, const mgTLPoints& tlpoints)const{
	return tlpoints[operator[](id)];
}
MGPosition& mgTLTriangle::uv(size_t id, mgTLPoints& tlpoints){
	return tlpoints[operator[](id)];
}

//Get world coordinates of the surface from the id of this triangl.
//The output is surf.eval(uv(id,tlpoints));
MGPosition mgTLTriangle::world(
	size_t id,					//id of this triangle.
		//0<= id <= size().
	const mgTLPoints& tlpoints,//tlpoints[(*this)[id]] is (u,v),
		// or uv(id, tlpoints) is (u,v).
	const MGSurface& surf		//Tessellated surface. When the object was
		//a face f, surf=*(f.surface());
)const{
	return surf.eval(uv(id,tlpoints));
}

//Get world coordinates of the surface from the id of this triangl.
//The output is surf.eval(uv(id,tlpoints));
MGPosition mgTLTriangle::world(
	size_t id,	//id of this triangle. 0<= id <= size().
	const mgTLPoints& tlpoints,//tlpoints[(*this)[id]] is (u,v),
		// or uv(id, tlpoints) is (u,v).
	const MGFace& f		//Tessellated surface. When the object was
)const{
	return f.surface()->eval(uv(id,tlpoints));
}
