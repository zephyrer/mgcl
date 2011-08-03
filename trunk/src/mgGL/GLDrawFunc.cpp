/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"

#include "mg/Tolerance.h"
#include "mg/DrawFunc.h"
#include "mg/Object.h"
#include "mg/Box.h"
#include "mg/LBRep.h"
#include "mg/SPointSeq.h"
#include "mg/FSurface.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/MGStl.h"
#include "mg/PickObjectCB.h"
#include "mg/PickObjectFB.h"
#include "mg/PickObjectSB.h"
#include "topo/Edge.h"
#include "topo/Face.h"
#include "Tl/TLTriangles.h"
#include "Tl/TLTexPoints.h"
#include "Tl/TLRect.h"
#include "Tl/TLData.h"
#include "Tl/TLDataVector.h"
#include "Tl/TLTexPlane.h"
#include "mgGL/GLDrawFunc.h"
#include "mgGL/OpenGLView.h"

//#include <algorithm>

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

namespace mgGDL{

struct vertex3dv{
	void operator()(const MGPosition& pos) const{::glVertex3dv(pos.data());}
};

//OpenGL declaration for shading. setup_easy_shade calls the follows:
//1. glEnable for GL_LIGHT0, GL_LIGHTING, GL_COLOR_MATERIAL.
//2. glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
//3. glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
//4. glDisable(GL_CULL_FACE);// 両面を描画
//5. glColorMaterial(GL_FRONT_AND_BACK) for GL_AMBIENT & GL_DIFFUSE.
void setup_easy_shade(
){
	// シェーディング描画を行う
	glEnable(GL_LIGHT0);// 光源を有効化
	glEnable(GL_LIGHTING);// ライティング開始
	glDisable(GL_CULL_FACE);// 両面を描画
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);// 図形の両面を照光
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT);
	glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);// 図形の両面を塗りつぶす
}

//OpenGL declaration for shading. setup_wire_drawing calls the follows:
void setup_wire_drawing(){
	glDisable(GL_COLOR_MATERIAL);// ライトをオフ
	glDisable(GL_LIGHTING);// ライティングオフ	

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);// ワイヤーフレームを描画		
}

//Draw an arrow symbol with implementation of OpenGL.
//data[0] is the origin of the arrow, data[1] is the top of the arrow,
//data[2], [3] are two bottoms of arrowhead.
void MGDrawArrow(MGPosition pos[4]){
	glBegin(GL_LINE_STRIP);
	glVertex3dv(pos[2].data());
	glVertex3dv(pos[1].data());
	glVertex3dv(pos[3].data());
	glEnd();
	glBegin(GL_LINES);
	glVertex3dv(pos[1].data());
	glVertex3dv(pos[0].data());
	glEnd();
}

void MGDrawArrow(const MGCurve& curve, int ndiv){
	double ts=curve.param_s(), te=curve.param_e();
	MGPosition pos[4];
	for(int i = 0; i <= ndiv; i++){
		double param = (ts*(ndiv-i)+te*i)/ndiv;
		curve.arrow(param, pos);
		mgGDL::MGDrawArrow(pos);
	}
}

//const float ucolor[3]={0.895f, 0.6f, 0.6f};
//const float vcolor[3]={0.6f, 0.895f, 0.6f};
//const float white[3] = {1.0f, 1.0f, 1.0f};

void MGDrawArrow(const MGFSurface& face, int udiv, int vdiv){
//	float saved[3];
//	::glGetFloatv(GL_CURRENT_COLOR, &saved[0]);

	MGPosition uv(2), pos[10];
	const MGBox& box = face.box_param2();
	double us = box[0].low_point(), ue = box[0].high_point(),
		vs = box[1].low_point(), ve = box[1].high_point();
	const MGColor& ucolor=MGColor::get_instance(MGColor::Red);
	const MGColor& vcolor=MGColor::get_instance(MGColor::Green);
	const MGColor& white=MGColor::get_instance(MGColor::White);
	for(int i = 0; i <= udiv; i++){
		uv(0) = (us*(udiv-i)+ue*i)/udiv;
		for(int j = 0; j <= vdiv; j++){
			uv(1) = (vs*(vdiv-j)+ve*j)/vdiv;
			face.arrow(uv, pos);

			ucolor.exec();
			mgGDL::MGDrawArrow(&pos[0]);
			pos[3] = pos[0];
			
			vcolor.exec();
			mgGDL::MGDrawArrow(&pos[3]);
			pos[6] = pos[0];
			
			white.exec();
			mgGDL::MGDrawArrow(&pos[6]);
		}
	}
//	::glColor3fv(saved);
	::glColor3f(0.0f, 0.0f, 0.0f);
}

// Draw an object of class MGBox, by wireframe.
void MGDrawBox(const MGBox& box){
	size_t sd = box.sdim();
	if(sd == 2){
		double u0=box[0].low_point(), u1=box[0].high_point();
		double v0=box[1].low_point(), v1=box[1].high_point();
		glBegin(GL_LINE_LOOP);
			glVertex3d(u0, v0, 0.); glVertex3d(u1, v0, 0.);
			glVertex3d(u1, v1, 0.); glVertex3d(u0, v1, 0.);
		glEnd();
		return;
	}else if(sd == 3){
		const MGInterval& x = box[0];
		const MGInterval& y = box[1];
		const MGInterval& z = box[2];
	
		glBegin(GL_LINE_LOOP);
		{
			glVertex3d(x.low_point(), y.low_point(), z.low_point());
			glVertex3d(x.high_point(), y.low_point(), z.low_point());
			glVertex3d(x.high_point(), y.high_point(), z.low_point());
			glVertex3d(x.high_point(), y.low_point(), z.low_point());
		}
		glEnd();
		glBegin(GL_LINE_LOOP);
		{
			glVertex3d(x.low_point(), y.low_point(), z.high_point());
			glVertex3d(x.high_point(), y.low_point(), z.high_point());
			glVertex3d(x.high_point(), y.high_point(), z.high_point());
			glVertex3d(x.high_point(), y.low_point(), z.high_point());
		}
		glEnd();
		glBegin(GL_LINES);
		{
			glVertex3d(x.low_point(), y.low_point(), z.low_point());
			glVertex3d(x.low_point(), y.low_point(), z.high_point());

			glVertex3d(x.high_point(), y.low_point(), z.low_point());
			glVertex3d(x.high_point(), y.low_point(), z.high_point());

			glVertex3d(x.high_point(), y.high_point(), z.low_point());
			glVertex3d(x.high_point(), y.high_point(), z.high_point());

			glVertex3d(x.high_point(), y.low_point(), z.low_point());
			glVertex3d(x.high_point(), y.low_point(), z.high_point());
		}
		glEnd();
	}
}

// Draw a control points, dotted lines shall be drawn
// between point[i-1] and point[i], for i = 1, .., length()-1.
void MGDrawPointSeq(
	const MGBPointSeq& bp,
	bool draw_points,		//True if points be drawn.
	bool dotted				//true if dotted lines be drawn.
){
//	if(bp.sdim() < 3) return;

	if(dotted){
		glEnable(GL_LINE_STIPPLE);
		short pattern=0x5555;
		glLineStipple(2,pattern);
	}
	MGDrawPolyline(bp);
	if(dotted){
		glDisable(GL_LINE_STIPPLE);
	}
	if(draw_points){
		size_t n=bp.length();
		for(size_t i = 0; i < n; i++){
			MGDrawPoint(bp(i));
		}
	}
}

void MGDrawPointSeq(
	const MGSPointSeq& sp,
	bool draw_points,		//True if points be drawn.
	bool dotted				//true if dotted lines be drawn.
){
//	if(sp.sdim() < 3) return;

	if(dotted){
		glEnable(GL_LINE_STIPPLE);
		short pattern=0x5555;
		glLineStipple(2,pattern);
	}
	
	size_t i, nu = sp.length_u(), nv = sp.length_v();
	for(i = 0; i < nu; i++){
		glBegin(GL_LINE_STRIP);
		for(size_t j = 0; j < nv; j++){
			glVertex3d(sp(i,j,0),sp(i,j,1),sp(i,j,2));
		}
		glEnd();
	}

	for(i = 0; i < nv; i++){
		glBegin(GL_LINE_STRIP);
		for(size_t j = 0; j < nu; j++){
			glVertex3d(sp(j,i,0),sp(j,i,1),sp(j,i,2));
		}
		glEnd();
	}

	if(dotted){
		glDisable(GL_LINE_STIPPLE);
	}
	if(draw_points){
		for(size_t i = 0; i < nu; i++){
			for(size_t j = 0; j < nv; j++){
				MGDrawPoint(sp(i, j));
			}
		}
	}
}

//Draw curvature variation graph so-called Hige.
void MGDrawCurvaGraph(const MGCurve& crv, double scale, int density, bool use_radius){
	assert(density > 0);

	MGUnit_vector T, N, B;	
	double curva, torsion;
	double revdens = 1. / density;
	MGPosition pos(3);
	MGVector v(crv.start_point());  // v : vertex of the hige's polyline
	const MGKnotVector& kv = crv.knot_vector();

	glBegin(GL_LINES);
	for(size_t i = kv.order()-1; i < kv.bdim(); i++){
		MGInterval iv(kv[i], kv[i+1]);
		double dt = (kv[i+1]-kv[i])*revdens;
		double t = kv[i];
		for(int j = 0; j < density; j++, t += dt){
			glVertex3dv(v.data());
			pos = crv.eval(t);     // pos : on the curve
			crv.Frenet_frame(t, T, N, B, curva, torsion);
			// N is from pos to the center of the osculating circle at pos.

			if(!use_radius){
				v = pos - scale * (N * curva); // curvature
			}else{
				v = pos + scale * (N / curva); // radius of curvature
			}

			glVertex3dv(v.data());
			glVertex3dv(pos.data());
			glVertex3dv(v.data());
		}
	}
	glVertex3dv(v.data());
	pos = crv.end_point();
	crv.Frenet_frame(crv.param_e(), T, N, B, curva, torsion);
	if(!use_radius){
		v = pos - scale * (N * curva); // curvature
	}else{
		v = pos + scale * (N / curva); // radius of curvature
	}
	glVertex3dv(v.data());
	glVertex3dv(pos.data());
	glVertex3dv(v.data());
	glEnd();
}

//Draw a point using openGL functions.
void MGDrawPoint(double x,double y,double z){
//	glColor3fv(PColor1);//Boundary color
	float psize=MGDrawFunc::point_size();
	glPointSize(psize);
	glBegin(GL_POINTS);
	glVertex3d(x,y,z);
	glEnd();

	glPushAttrib(GL_CURRENT_BIT);
	MGDrawFunc::PColor2().exec();//Inner color.
	glPointSize(float(psize-2.));
	glBegin(GL_POINTS);
	glVertex3d(x,y,z);
	glEnd();
	glPopAttrib();
}

void MGDrawPoint(const MGPosition& pos){
	MGDrawPoint(pos[0], pos[1], pos[2]);
}

void MGGLBeginLINE_STRIP(){glBegin(GL_LINE_STRIP);}

void MGGLVertex(double x, double y, double z){glVertex3d(x,y,z);}

void MGGLEnd(){	glEnd();}

//Generate the OpenGL display list for GL3ParamView's parameter rectangle.
void param_rectangle(
	const MGBox& param_box	//Parameter range.
){
	double u0=param_box[0].low_point(), u1=param_box[0].high_point();
	double v0=param_box[1].low_point(), v1=param_box[1].high_point();
	glEnable(GL_LINE_STIPPLE);
	short pattern=0x5555;
	glLineStipple(2,pattern);
	glBegin(GL_LINE_LOOP);
		glVertex3d(u0, v0, 0.); glVertex3d(u1, v0, 0.);
		glVertex3d(u1, v1, 0.); glVertex3d(u0, v1, 0.);
	glEnd();
	glDisable(GL_LINE_STIPPLE);
}

//OpenGL display for the tessellation lines drawn in world view.
void mgGDL_WTess(
	const mgTLData& tld	//tessellation data.
){
	const mgTLTriangles& tris=tld.triangles();
	const MGSurface& surf=tld.surface();
	const mgTLPoints& tlpoints=tld.tlpoints();
	mgTLTriangles::const_triIterator i=tris.begin(), ie=tris.end();
	for(; i!=ie; i++){
		const mgTLTriangle& tri=**i;
		if(tri.size()<3) continue;//If number of vertices are less than 3.
		if(tri.getGeometryType()==mgTESTRIANG_FAN){
			mgTLTriangle::CIndexItr j=tri.begin(), je=tri.end();
			MGPosition pivot=surf.eval(tlpoints[*j++]);
			MGPosition pre=surf.eval(tlpoints[*j++]);
			MGPosition aft=surf.eval(tlpoints[*j++]);
			glBegin(GL_LINE_LOOP);
				glVertex3dv(pivot.data());
				glVertex3dv(pre.data());
				glVertex3dv(aft.data());
			glEnd();
			for(;j!=je; j++){
				pre=aft;
				aft=surf.eval(tlpoints[*j]);
				glBegin(GL_LINE_STRIP);
					glVertex3dv(pre.data());
					glVertex3dv(aft.data());
					glVertex3dv(pivot.data());
				glEnd();
			}
		}else{
			mgTLTriangle::CIndexItr j=tri.begin(), je=tri.end();
			MGPosition pre=surf.eval(tlpoints[*j++]);
			MGPosition current=surf.eval(tlpoints[*j++]);
			MGPosition aft=surf.eval(tlpoints[*j++]);
			glBegin(GL_LINE_LOOP);
				glVertex3dv(current.data());
				glVertex3dv(pre.data());
				glVertex3dv(aft.data());
			glEnd();
			for(;j!=je; j++){
				pre=current;
				current=aft;
				aft=surf.eval(tlpoints[*j]);
				glBegin(GL_LINE_STRIP);
					glVertex3dv(current.data());
					glVertex3dv(aft.data());
					glVertex3dv(pre.data());
				glEnd();
			}
		}
	}
}

//OpenGL display for the tessellation lines drawn in parameter view.
void mgGDL_PTess(
	const mgTLData& tld	//tessellation data.
){
	const mgTLTriangles& tris=tld.triangles();
	const mgTLPoints& tlpoints=tld.tlpoints();
	mgTLTriangles::const_triIterator i=tris.begin(), ie=tris.end();
	for(; i!=ie; i++){
		const mgTLTriangle& tri=**i;
		if(tri.size()<3) continue;//If number of vertices are less than 3.
		if(tri.getGeometryType()==mgTESTRIANG_FAN){
			mgTLTriangle::CIndexItr j=tri.begin(), je=tri.end();
			const double* pivot=tlpoints[*j++].data();
			const double* pre=tlpoints[*j++].data();
			const double* aft=tlpoints[*j++].data();
			glBegin(GL_LINE_LOOP);
				glVertex2dv(pivot);
				glVertex2dv(pre);
				glVertex2dv(aft);
			glEnd();
			for(;j!=je; j++){
				pre=aft;
				aft=tlpoints[*j].data();
				glBegin(GL_LINE_STRIP);
					glVertex2dv(pre);
					glVertex2dv(aft);
					glVertex2dv(pivot);
				glEnd();
			}
		}else{
			mgTLTriangle::CIndexItr j=tri.begin(), je=tri.end();
			const double* pre=tlpoints[*j++].data();
			const double* current=tlpoints[*j++].data();
			const double* aft=tlpoints[*j++].data();
			glBegin(GL_LINE_LOOP);
				glVertex2dv(current);
				glVertex2dv(pre);
				glVertex2dv(aft);
			glEnd();
			for(;j!=je; j++){
				pre=current;
				current=aft;
				aft=tlpoints[*j].data();
				glBegin(GL_LINE_STRIP);
					glVertex2dv(current);
					glVertex2dv(aft);
					glVertex2dv(pre);
				glEnd();
			}
		}
	}
}

//Draw a polyline using openGL functions.
//When clodes=true, 1st and last points will be connected.
void MGDrawPolyline(const MGBPointSeq& line, bool closed){
	size_t kind=GL_LINE_STRIP;
	if(closed) kind=GL_LINE_LOOP;
	glBegin(kind);
	size_t n = line.length();
	for(size_t i = 0; i < n; i++){
		glVertex3d(line(i,0),line(i,1),line(i,2));
	}
	glEnd();
}

//Draw a polyline using openGL functions.
//Draw a polyline using openGL functions.
//When cloded=true, 1st and last points will be connected.
void MGDrawPolyline(const std::vector<MGPosition>& line, bool closed){
	size_t kind=GL_LINE_STRIP;
	if(closed) kind=GL_LINE_LOOP;
	glBegin(kind);
	size_t n = line.size();
	for(size_t i = 0; i < n; i++){
		const MGPosition& Pi=line[i];
		glVertex3d(Pi[0],Pi[1],Pi[2]);
	}
	glEnd();
}
void MGDrawStraight(const MGPosition& end, const MGPosition& start){
	glBegin(GL_LINE_STRIP);
	glVertex3d(start[0],start[1],start[2]);
	glVertex3d(end[0],end[1],end[2]);
	glEnd();
}

// Renders curvatures mapping that comes into colorful image.
// A point whose curvature is within [lower, upper], the color varies.
void draw_surface_curvature_data(
	const mgTLData& tld,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper, //minimum and maximum value of the curvatures of the kind.
		//Whne lower>=upper, lower is set as the minimum value and upper is set
		//as the maximum value out of all the curvatures.
	double* real_lower,	//When not null, minimum curvature of the kind will be output.
	double* real_upper	//When not null, maximum curvature of the kind will be output.
){	
	const MGSurface& surf=tld.surface();
	const mgTLPoints& tlpoints=tld.tlpoints();
	size_t npoint=tlpoints.size();

	double kappa, max_curvature=0., min_curvature=0.;
	std::vector<MGPosition> curvaData(npoint);
		//In MGPosition, (0)=curvature, (1-3)=normal, (4-6)=point data.
	mgTLPoints::const_iterator l=tlpoints.begin();
	for(size_t k=0; k<npoint; k++,l++){
		const MGPosition& uv = *l;
		MGPosition& curvaDatak=curvaData[k];
		double curvature[4];
		MGUnit_vector N;
		surf.curvatures(uv,curvature,N);
		curvaDatak.resize(7);
		curvaDatak(0)=kappa=curvature[kind];
		if(k){
			if(max_curvature<kappa)
				max_curvature=kappa;
			else if(min_curvature>kappa)
				min_curvature=kappa;
		}else{
			max_curvature=min_curvature=kappa;
		}
		curvaDatak.store_at(1,N,0,3);
		curvaDatak.store_at(4,surf.eval(uv),0,3);
	}
	if(real_lower)
		*real_lower=min_curvature;
	if(real_upper)
		*real_upper=max_curvature;

	if(lower >= upper){
		lower = min_curvature;
		upper = max_curvature;
	}
	double mzero=MGTolerance::mach_zero();
	if(upper-lower<=2.*mzero){
		upper+=mzero;
		lower-=mzero;
	}

	const double R[3] = {1.0, 0.0, 0.0};
	const double Y[3] = {1.0, 1.0, 0.0};
	const double G[3] = {0.0, 1.0, 0.0};
	const double C[3] = {0.0, 1.0, 1.0};
	const double B[3] = {0.0, 0.0, 1.0};

	MGBPointSeq bp(5, 3);
	bp.store_at(4, R);
	bp.store_at(3, Y);
	bp.store_at(2, G);
	bp.store_at(1, C);
	bp.store_at(0, B);
	int err = 0;
	MGLBRep color(bp, err, 2);
	color.change_range(lower, upper);

	const mgTLTriangles& tris=tld.triangles();
	mgTLTriangles::const_triIterator i=tris.begin(), ie=tris.end();
	for(; i!=ie; ++i){
		const mgTLTriangle& tri=**i;
		if(tri.size()<3) continue;

		mgTLTriangle::CIndexItr j=tri.begin(), je=tri.end();
		(tri.getGeometryType()==mgTESTRIANG_FAN)
			? glBegin(GL_TRIANGLE_FAN)
			: glBegin(GL_TRIANGLE_STRIP);
		for(;j!=je; ++j){
			MGPosition& curvaDatak=curvaData[*j];

			//Nomal.
			glNormal3d(curvaDatak[1],curvaDatak[2],curvaDatak[3]);

			//Curvature color.
			double curvature = curvaDatak[0];
			if(curvature > upper) curvature = upper;
			else if(curvature < lower) curvature = lower;
			glColor3dv(color.eval(curvature).data());

			//Position data.
			glVertex3d(curvaDatak[4],curvaDatak[5],curvaDatak[6]);
		}
		glEnd();
	}
}

// Renders curvatures mapping that comes into colorful image.
// A point whose curvature is within [lower, upper], the color varies.
void draw_surface_curvature(
	const mgTLDataVector& tldvec,
	MGCL::SURFACE_CURVATURE_KIND kind,
	double lower, double upper, //minimum and maximum value of the curvatures of the kind.
		//Whne lower>=upper, lower is set as the minimum value and upper is set
		//as the maximum value out of all the curvatures.
	double* real_lower,	//When not null, minimum curvature of the kind will be output.
	double* real_upper	//When not null, maximum curvature of the kind will be output.
){
	mgTLDataVector::const_iterator datai=tldvec.begin(), dataie=tldvec.end();
	for(; datai!=dataie; datai++){
		draw_surface_curvature_data(*datai,kind,lower,upper,real_lower,real_upper);
	}
}

//OpenGL shading display of a tesselated data tld.
void shade(
	const mgTLData& tld
){
	const mgTLTriangles& tris=tld.triangles();
	const MGSurface& surf=tld.surface();
	const mgTLPoints& tlpoints=tld.tlpoints();

	for(mgTLTriangles::const_triIterator i=tris.begin(); i!=tris.end(); ++i){
		const mgTLTriangle& tri=**i;
		if(tri.size()<3) continue;

		if(tri.getGeometryType()==mgTESTRIANG_FAN){
			glBegin(GL_TRIANGLE_FAN);
			//std::cout<<"start of tri:GL_TRIANGLE_FAN"<<std::endl;
		}else{
			glBegin(GL_TRIANGLE_STRIP);
			//std::cout<<"start of tri:GL_TRIANGLE_STRIP"<<std::endl;
		}
		mgTLTriangle::CIndexItr j=tri.begin(), je=tri.end(); // size_t
		for(;j!=je; ++j){
			const MGPosition& uv = tlpoints[*j];
			glNormal3dv(surf.unit_normal(uv).data());
			glVertex3dv(surf.eval(uv).data());
			//std::cout<<"uv="<<uv<<",pos="<<surf.eval(uv);///////////
			//std::cout<<",normal="<<surf.unit_normal(uv)<<std::endl;///////////
		}
		glEnd();
	}
}

// MGStlオブジェクトを描画する
// glBiginからglEndまでを実行する
void draw_STL(
	 const MGStl& stl, // 描画するMGStlオブジェクト
	 bool invokeNormal // glNormalを呼ぶかどうか
){
	const std::vector<MGPosition>& vertices = stl.positions();
	const std::vector<MGUnit_vector>& normals = stl.normals();

	glBegin(GL_TRIANGLES);// 描画を開始
	int nTriang = normals.size(); // 三角形の個数だけループ
	size_t indices[3];
	for(int j = 0; j < nTriang; j++){
		stl.GetVertIndices(j, indices);
		// 頂点のインデックスを元に頂点の座標を取得する
		const MGPosition& position1 = vertices[indices[0]];
		const MGPosition& position2 = vertices[indices[1]];
		const MGPosition& position3 = vertices[indices[2]];

		if(invokeNormal)
			// シェーディングを行う場合
			glNormal3dv(normals[j].data());
		glVertex3dv(position1.data());
		glVertex3dv(position2.data());
		glVertex3dv(position3.data());
	}
	glEnd();
}

//draw poits sequence ipos with 2 colors, inner and outer.
void draw_points(
	const MGColor& boundary_color,
	const MGColor& inner_color,
	const std::vector<MGPosition>& ipos
){
	boundary_color.exec();//Boundary color
	size_t n = ipos.size();
	float psize=MGDrawFunc::point_size();
	glPointSize(psize);
	glBegin(GL_POINTS);
		std::for_each(ipos.begin(), ipos.end(), vertex3dv());
	glEnd();

	inner_color.exec();//Inner color.
	glPointSize(float(psize-2.));
	glBegin(GL_POINTS);
		std::for_each(ipos.begin(), ipos.end(), vertex3dv());
	glEnd();
}

}//End of namespace mgGDL.

//Make 2 types of display list of this gel(wire and shading).
//Return is the display list name.
size_t MGObject::make_display_list(
	double span_length,//span length to approximate by polyline.
	int line_density//line density to draw surface in wire mode.
)const{
	size_t name=dlist_name();
	size_t glname=name;

	//Draw wire mode display list.
	glNewList(glname, GL_COMPILE);
		glPushName(name);
		size_t mask=get_draw_attrib_mask();
		if(mask){
			glPushAttrib(mask);
			draw_attribute(false);
		}
		drawWire(span_length,line_density);
		if(mask) glPopAttrib();
		glPopName();
	glEndList();

	//Draw shading mode display list.
	glname=name+2;
	glNewList(glname, GL_COMPILE);
		glPushName(name);
		if(manifold_dimension()<2)
			glCallList(name);//Only wire is the only possible mode to draw.
		else{
			size_t mask=get_draw_attrib_mask();
			if(mask){
				glPushAttrib(mask);
				draw_attribute(false);
			}
			shade(span_length);
			if(mask)
				glPopAttrib();
		}
		glPopName();
	glEndList();
	return name;
}

///Make a display list without color of this gel.
///Return is the display list name.
size_t MGObject::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{
	size_t glname=dlist_name()+1;

	//Draw wire mode display list.
	glNewList(glname, GL_COMPILE);
		size_t mask=get_draw_attrib_mask();
		if(mask){
			glPushAttrib(mask);
			draw_attribute(true);
		}
		drawWire(span_length,line_density);
		if(mask) glPopAttrib();
	glEndList();
	return glname;
}

///Delete a display list of this gel.
void MGObject::delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	//mgSysGL pointer will be apppended.
)const{
	size_t nm=dlist_name();
	glDeleteLists(nm,3);
	glv.m_sysgllist.delete_lists_by_object_id(this,functions);
}

//Make a display list of this gel.
size_t MGGroup::make_display_list(
	double span_length,//span length to approximate by polyline.
	int line_density//line density to draw surface in wire mode.
)const{
	size_t name=make_only_call_list(false);
	MGGroup::const_iterator i,is=begin(), ie=end();	
	//Make display list of the gel that has an object.
	for(i=is; i!=ie; i++){
		if((*i)->no_display())
			continue;
		size_t name2=(*i)->dlist_name();
		if(name2)
			(*i)->make_display_list(span_length,line_density);
	}
	return name;
}

///Make a display list without color of this gel.
///Return is the display list name.
size_t MGFSurface::make_display_list_to_hilightFS(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{
	const MGObject* obj=object_pointer();
	size_t glname=obj->dlist_name()+1;

	//Draw wire mode display list.
	glNewList(glname, GL_COMPILE);
		size_t mask=obj->get_draw_attrib_mask();
		if(mask){
			glPushAttrib(mask);
			obj->draw_attribute(true);
		}
		drawWireFS_to_highlight(span_length,line_density);
		if(mask) glPopAttrib();
	glEndList();
	return glname;
}

///Make a display list without color of this gel.
///Return is the display list name.
size_t MGGroup::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{
	size_t name=make_only_call_list(true);
	MGGroup::const_iterator i,is=begin(), ie=end();	
	//Make display list of the gel that has an object.
	for(i=is; i!=ie; i++){
		const MGGel* gel=*i;
		if(gel->no_display())
			continue;
		size_t name2=gel->dlist_name();
		if(name2)
			gel->make_display_list_to_hilight(span_length,line_density);
	}
	return name;
}

//Delete a display list of this gel.
void MGGroup::delete_display_list(
	MGOpenGLView& glv,
	MGPvector<mgSysGL>& functions	//mgSysGL pointer will be apppended.
)const{
	MGGroup::const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		(*i)->delete_display_list(glv,functions);
	size_t nm=dlist_name();
	glDeleteLists(nm,1);
	glv.m_sysgllist.delete_lists_by_object_id(this,functions);
}

void MGPlane::display_arrows()const{
	MGVector U,V;
	get_uv_display_vector(U,V);
	const MGPosition& cen = center();
	MGBox box(cen-U,cen+U);
	box.expand(cen-V);
	box.expand(cen+V);
	MGPosition pos[10], uv=center_param();
	arrow(box,uv[0],uv[1],pos);

	const MGColor& ucolor=MGColor::get_instance(MGColor::Red);
	const MGColor& vcolor=MGColor::get_instance(MGColor::Green);
	const MGColor& white=MGColor::get_instance(MGColor::White);
	ucolor.exec();
	mgGDL::MGDrawArrow(&pos[0]);

	pos[3] = pos[0];			
	vcolor.exec();
	mgGDL::MGDrawArrow(&pos[3]);

	pos[6] = pos[0];			
	white.exec();
	mgGDL::MGDrawArrow(&pos[6]);
}

//Shade the object in world coordinates.
void MGStl::shade(
	double span_length	//Line segment span length.
)const{
	glShadeModel(GL_FLAT);// shading model
	mgGDL::draw_STL(*this,true);
}

//Draw 3D curve in world coordinates.
//The object is converted to curve(s) and is drawn.
void MGStl::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	mgGDL::draw_STL(*this);
}

// 三角形ごとの法線ベクトルを表示する
void MGStl::display_arrows()const{
	// 三角形の数を取得
	size_t nTriangle = m_vecNormlTriang.size();
	// 矢印描画のための座標値の配列
	MGPosition pos[4];
	// 色を示すインスタンスを生成
	const MGColor& white=MGColor::get_instance(MGColor::White);

	for(size_t i = 0; i < nTriangle; i++){
	// それぞれの三角形の法線方向に矢印を描画
		// 三角形の中点を取得
		size_t i3 = i*3;
		pos[0] = (m_vecPos[m_indices[i3]]+m_vecPos[m_indices[i3+1]]+m_vecPos[m_indices[i3+2]])/3;

		// 各辺の長さを取得
		double dist[3];
		dist[0] = m_vecPos[m_indices[i3]].distance(m_vecPos[m_indices[i3+1]]);
		dist[1] = m_vecPos[m_indices[i3+1]].distance(m_vecPos[m_indices[i3+2]]);
		dist[2] = m_vecPos[m_indices[i3]].distance(m_vecPos[m_indices[i3+2]]);
		
		// 三角形の中から最長の辺の長さを取得し
		// その1/2の値を矢印の軸の長さに用いる
		double max = dist[0];
		for(size_t j = 0; j < 2; j++){
			if(dist[j] < dist[j+1]){
				max = dist[j+1];
			}
		}
		double len = max/2;

		// 矢印の先端の座標を計算
		const MGVector& vecX = m_vecNormlTriang[i] * len;
		pos[1] = pos[0] + vecX;

		// 矢印の両端の座標を求める処理
		// 矢印の先端座標に加えるベクトルを計算
		const MGVector& head_rootx = vecX * .3;
		// 三角形の任意の辺のベクトルを取得し、面の法線べクトルとの積算を行う
		MGUnit_vector& arrowVec = (m_vecPos[m_indices[i3+1]]- m_vecPos[m_indices[i3]]).normalize();
		arrowVec *= m_vecNormlTriang[i];
		// 矢印の先端座標に加えるもう１つのベクトルを計算
		const MGVector& head_rooty = arrowVec*(.5*.3*vecX.len());
		// 矢印の両端の座標を計算
		pos[2]=pos[1]-head_rootx+head_rooty;
		pos[3]=pos[1]-head_rootx-head_rooty;

		// 矢印の描画を行う
		white.exec();
		mgGDL::MGDrawArrow(pos);
	}
}

//Shade the object in world coordinates.
void MGSurface::shade(
	double span_length	//Line segment span length.
)const{
	mgTLData tess(*this,mgTLInputParam(*this,span_length));
	glShadeModel(GL_SMOOTH);// shading model
	mgGDL::shade(tess);
}

//Draw 3D curve in world coordinates.
//The object is converted to curve(s) and is drawn.
void MGFSurface::drawWireFS(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	glPushAttrib(GL_LINE_BIT);
	glLineWidth(1.);
	MGPvector<MGCurve> ilines=inner_skeleton(line_density);
	MGPvector<MGCurve>::const_iterator i=ilines.begin(), ie=ilines.end();
	for(; i!=ie; i++)
		(*i)->drawWire(span_length,line_density);

	glLineWidth(2.);
	MGPvector<MGCurve> bndries=get_all_boundaries();
	MGPvector<MGCurve>::const_iterator j=bndries.begin(), je=bndries.end();
	for(; j!=je; j++)
		(*j)->drawWire(span_length,line_density);
	glPopAttrib();
}

//Draw 3D curve in world coordinates.
//The object is converted to curve(s) and is drawn.
void MGFSurface::drawWireFS_to_highlight(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	MGPvector<MGCurve> ilines=skeleton(line_density);
	MGPvector<MGCurve>::const_iterator i=ilines.begin(), ie=ilines.end();
	for(; i!=ie; i++)
		(*i)->drawWire(span_length,line_density);
}

//Make a display list of only call lists of this groupby the name to push for pick
//and display list name glname.
//When gels in this group are included in gels_to_delete, the display list
//will not be made.
void MGGroup::make_only_call_list_sub(
	size_t name,	//the name to glPushName for pick
	int mode,	//=0 when wire, =1 when wire and no color, =2 when SHADING.
	const std::vector<const MGGel*>* gels_to_delete
)const{
	std::vector<const MGGel*>::const_iterator gelsE,gelsS;
	if(gels_to_delete){
		gelsS=gels_to_delete->begin();
		gelsE=gels_to_delete->end();
	}

	bool no_color=false;//if true, color attribute will be neglected.
	if(mode==1)
		no_color=true;
	size_t glname=name+mode;
	glNewList(glname, GL_COMPILE);
		glPushName(name);
		size_t mask=get_draw_attrib_mask();
		if(mask){
			glPushAttrib(mask);
			draw_attribute(no_color);
		}
		MGGroup::const_iterator i,is=begin(), ie=end();
		for(i=is; i!=ie; i++){
			if((*i)->no_display())
				continue;

			if(gels_to_delete)
				if(std::find(gelsS,gelsE,*i)!=gelsE)
					continue;//If *i is found in gels_to_delete.
			size_t name2=(*i)->dlist_name();
			if(name2){
				name2+=mode;
				glCallList(name2);//call object.
			}else{
				(*i)->drawAttrib(no_color);//This will be attributes.
			}
		}
		if(mask)
			glPopAttrib();
		glPopName();
	glEndList();
}

///Make the display list of this object as a highlighted one.
void MGPickObject::make_display_list_to_hilight(
	double span_length,///<span length to approximate by polyline.
	int line_density///<line density to draw surface in wire mode.
)const{
	const MGGel* gl=gel();
	gl->make_display_list_to_hilight(span_length,line_density);
		//generate without color. When generating color neglecting display list
		//the display list name will be gel->dlist_name()+1.
}

static const float edgeColor[4]={1.,.5,.5,0.};//Edge color to hilight.
static const float white[4]={1.,1.,1.,0.};//Highlight back color.
static const float endPointColor[4]={.5,1.,.5,0.};//Start/End point color to hilight.

///Highlightthe object using the display list of this object.
void MGPickObject::hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const{
	glColor4fv(white);
	glLineWidth(2.);
	size_t nm=gel()->dlist_name()+1;
	glCallList(nm);

	glLineWidth(1.);
	Hcolor.exec();
	glCallList(nm);
}

///Highlight the object using the display list of this object.
void MGPickObjectFB::hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const{
	MGPickObject::hilight_using_display_list(Hcolor,span_length,line_density);
	glColor4fv(edgeColor);
	const MGEdge& e=*(edge());
	double slen=1.;
	e.drawWire(slen);
}

///Highlight the object using the display list of this object.
void MGPickObjectSB::hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const{
	MGPickObject::hilight_using_display_list(Hcolor,span_length,line_density);
	glColor4fv(edgeColor);
	std::auto_ptr<MGCurve> e(surface()->perimeter_curve(perimeter()));
	double slen=1.;
	e->drawWire(slen);
}

///Highlightthe object using the display list of this object.
void MGPickObjectCB::hilight_using_display_list(
	const MGColor& Hcolor,
	double span_length,	///<Line segment span length.
	int line_density	///<line density to draw a surface in wire mode.
)const{
	MGPickObject::hilight_using_display_list(Hcolor,span_length,line_density);
	const MGCurve* c=curve();
	MGPosition P=m_start_end ? c->end_point() : c->start_point();
	glColor4fv(endPointColor);
	mgGDL::MGDrawPoint(P);
}
