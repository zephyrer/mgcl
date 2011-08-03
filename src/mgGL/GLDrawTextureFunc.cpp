/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/DrawFunc.h"
#include "mg/Straight.h"
#include "mg/CParam_list.h"
#include "mg/Surface.h"
#include "Tl/TLTriangles.h"
#include "Tl/TLTexPoints.h"
#include "Tl/TLRect.h"
#include "Tl/TLDataVector.h"
#include "Tl/TLTexPlane.h"
#include "mgGL/GLDrawFunc.h"
#include "mgGL/ImageRect.h"
#include "mgGL/TextureImage2.h"
#include "mgGL/TextureImages.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

namespace mgGDL{
void draw_texture_image_set_up(
	const double back_color[4]	//back ground color data.
){
	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.);
	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_NORMALIZE);
	::glColor4dv(back_color);
	glEnable(GL_TEXTURE_2D);
}

void draw_texture_image_end_up(){
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_TEXTURE_2D);
}
}

//Draw textures fo tri. tir's texture coordinates are stored in tlTexPoints() as their
//woordl coordinate length. draw_texture_image_tri() draws texture of tri converting
//the texture coordinates to [0,1] span using the box(of the image).
//box provides the texture's (left, bottom) coordinates and the length of
//the image's [0,1] span.
void mgTLData::draw_texture_image_tri(
	const mgTLTriangle& tri,
	const MGBox& box//The box (left,bottom, right, top) of the image
)const{
	if(tri.size()<3) return;

	const MGSurface& surf=surface();
	const mgTLPoints& tlpts=tlpoints();
	const mgTLTexPoints& texpnt=tlTexPoints();
	double w=box[0].length().value();
	double h=box[1].length().value();
	MGPosition base=box.low();

	if(tri.getGeometryType()==mgTESTRIANG_FAN){
		glBegin(GL_TRIANGLE_FAN);
		//std::cout<<"start of tri:GL_TRIANGLE_FAN"<<std::endl;
	}else{
		glBegin(GL_TRIANGLE_STRIP);
		//std::cout<<"start of tri:GL_TRIANGLE_STRIP"<<std::endl;
	}
	mgTLTriangle::CIndexItr j=tri.begin(), je=tri.end(); // size_t
	for(;j!=je; ++j){
		const MGPosition& uv = tlpts[*j];
		//glNormal3dv(surf.unit_normal(uv).data());
		MGPosition st=texpnt[*j]-base;
		glTexCoord2d(st[0]/w,st[1]/h);
		glVertex3dv(surf.eval(uv).data());
		//std::cout<<"uv="<<uv<<",pos="<<surf.eval(uv);///////////
		//std::cout<<",normal="<<surf.unit_normal(uv)<<std::endl;///////////
	}
	glEnd();
}

//Trim the trangle tri by the (s,t) box of the timage.
//Function's return value is the number of vertices of the output polygon.
//When the number of the vertices is 0, st_tri was outside the timage.
int trim_triangle(
	const MGBox& box,//The box (left,bottom, right, top)
	const MGPosition st_tri[3],
	const MGPosition world_tri[3],
	std::vector<MGPosition>& st_polygon,
		//The vertices of the trimmed polygon of (s,t) coordinates will be output.
		//The maximum number of the trimmed polygon is 7.
		//The coordinates will be subtracted by (left, bottom) of the box.
	std::vector<MGPosition>& world_polygon
		//The vertices world coordinates corresponding to the trimmed polygon
		//of (s,t) will be output.
		//The maximum number of the trimmed polygon is 7.
){
	MGPosition minst=box.low(), maxst=box.high();
	double smin=minst[0], smax=maxst[0];
	double tmin=minst[1], tmax=maxst[1];
	double s0=st_tri[0][0], s1=st_tri[1][0], s2=st_tri[2][0];
	double t0=st_tri[0][1], t1=st_tri[1][1], t2=st_tri[2][1];

	if(s0<=smin && s1<=smin && s2<=smin)
		return 0;
	if(s0>=smax && s1>=smax && s2>=smax)
		return 0;
	if(t0<=tmin && t1<=tmin && t2<=tmin)
		return 0;
	if(t0>=tmax && t1>=tmax && t2>=tmax)
		return 0;

	if(s0>=smin && s1>=smin && s2>=smin
		&& s0<=smax && s1<=smax && s2<=smax
		&& t0>=tmin && t1>=tmin && t2>=tmin
		&& t0<=tmax && t1<=tmax && t2<=tmax){
		for(size_t i=0; i<3; i++){
			st_polygon.push_back(st_tri[i]);
			world_polygon.push_back(world_tri[i]);
		}
		return 3;
	}

	mgImageRect irect(box,st_tri,world_tri);
	//std::cout<<"trim_triangle ,irect="<<irect;/////
	//for(size_t ii=0; ii<3; ii++)std::cout<<std::endl<<st_tri[ii];
	//std::cout<<std::endl;
	irect.createTrimPolygon(st_polygon);
	irect.convert_to_world(st_polygon,world_polygon);
	return st_polygon.size();
}

void mgTLData::draw_texture_image_tri_with_trim(
	const MGBox& image_box,//The image box to draw.
	const mgTLTriangle& triangle,
	const MGBox& trim_box//The trimming box of the triangle.
)const{
	size_t nv=triangle.size();
	if(nv<3) return;

	double w=image_box[0].length().value();//image width.
	double h=image_box[1].length().value();//image height

	MGPosition st_tri[3];
	MGPosition world_tri[3];
	MGPosition base=image_box.low();
	/*std::cout<<" draw_texture_image_tri_with_trim::Base="
		<<base<<std::endl<<", imageBox="<<image_box<<std::endl
		<<", trim_box="<<trim_box<<std::endl;*/
	get_texcoord_world(triangle[0],st_tri[0],world_tri[0],base);
	get_texcoord_world(triangle[1],st_tri[1],world_tri[1],base);

	MGBox trim_box2=trim_box-base;
	for(size_t i=2; i<nv; i++){
		std::vector<MGPosition> st_polygon;
		std::vector<MGPosition> world_polygon;
		get_texcoord_world(triangle[i],st_tri[2],world_tri[2],base);
		int nvertices=trim_triangle(trim_box2,st_tri,world_tri,st_polygon,world_polygon);
		if(nvertices){
			glBegin(GL_POLYGON);
			for(int j=0; j<nvertices; j++){
				//std::cout<<st_polygon[j]<<",";
				glTexCoord2d(st_polygon[j][0]/w,st_polygon[j][1]/h);
				glVertex3dv(world_polygon[j].data());
			}
			//std::cout<<std::endl;
			glEnd();
		}
		if(triangle.getGeometryType()==mgTESTRIANG_STRIP){
			st_tri[0]=st_tri[1];
			world_tri[0]=world_tri[1];
		}
		st_tri[1]=st_tri[2];
		world_tri[1]=world_tri[2];
	}
}

//Draw textures within the m_texture_planes[plane_id].
void mgTLData::drawt_texture_image_part(
	const MGTextureImage2& teximage,	//texture image data.
	const double back_color[4],	//back ground color data.
	int plane_id	//indicates if only parts of texture data be drawn.
		//let i1=m_texture_planes[plane_id]->index(), and i2=m_texture_planes[plane_id+1]->index().
		//Then rects of m_rects[i1] to m_rects[i2-1] are to be drawn.
)const{
	mgGDL::draw_texture_image_set_up(back_color);
	teximage.make_texture_object();
	teximage.apply_texture_object();

	const mgTLTriangles& tris=triangles();
	mgTLTriangles::const_triIterator first=tris.begin(), last=tris.end();;
	const std::vector<mgTLTexPlane*>& texplanes=tlrects().texture_planes();
	int ntpls=texplanes.size();
	if(plane_id<ntpls){
		int i1=texplanes[plane_id]->index();
		const mgTLRect& recti=rect(i1);
		int tri_id1=rect(i1).triangle_id();
		first+=tri_id1;
		if(plane_id<ntpls-1){
			int i2=texplanes[plane_id+1]->index();
			int tri_id2=rect(i2).triangle_id();
			last=first+(tri_id2-tri_id1);
		}
	}
	for(mgTLTriangles::const_triIterator i=first; i!=last; ++i){
		const mgTLTriangle& tri=**i;
		draw_texture_image_tri(tri,teximage.box());
	}

	mgGDL::draw_texture_image_end_up();
}

//Draw textures within the m_texture_planes[plane_id].
void mgTLData::drawt_texture_image_part(
	const MGTextureImages& teximage,	//texture image data.
	const double back_color[4],	//back ground color data.
	int plane_id	//indicates if only parts of texture data be drawn.
		//let i1=m_texture_planes[plane_id]->index(), and i2=m_texture_planes[plane_id+1]->index().
		//Then rects of m_rects[i1] to m_rects[i2-1] are to be drawn.
)const{
	mgGDL::draw_texture_image_set_up(back_color);

	const mgTLTriangles& tris=triangles();
	mgTLTriangles::const_triIterator first=tris.begin(), last=tris.end();;
	const std::vector<mgTLTexPlane*>& texplanes=tlrects().texture_planes();
	int ntpls=texplanes.size();
	if(plane_id<ntpls){
		int i1=texplanes[plane_id]->index();
		const mgTLRect& recti=rect(i1);
		int tri_id1=rect(i1).triangle_id();
		first+=tri_id1;
		if(plane_id<ntpls-1){
			int i2=texplanes[plane_id+1]->index();
			int tri_id2=rect(i2).triangle_id();
			last=first+(tri_id2-tri_id1);
		}
	}
	draw_texture(teximage,first,last);

	mgGDL::draw_texture_image_end_up();
}

// OpenGL display for shading.
void mgTLData::draw_texture(
	const MGTextureImage& teximage
)const{
	const mgTLTriangles& tris=triangles();
	for(mgTLTriangles::const_triIterator i=tris.begin(); i!=tris.end(); ++i){
		draw_texture_image_tri(**i,teximage.box());
	}
}

void mgTLDataVector::draw_texture(
	const MGTextureImage2& teximage,
	const double* back_color
)const{
	mgGDL::draw_texture_image_set_up(back_color);
	teximage.make_texture_object();
	teximage.apply_texture_object();
	mgTLDataVector::const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		(*i).draw_texture(teximage);
	mgGDL::draw_texture_image_end_up();
}

void mgTLDataVector::draw_texture(
	const MGTextureImages& teximages,
	const double* back_color
)const{
	mgGDL::draw_texture_image_set_up(back_color);
	mgTLDataVector::const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		(*i).draw_texture(teximages);
	mgGDL::draw_texture_image_end_up();
}

void mgTLData::draw_texture(
	const MGTextureImages& timages
)const{
	const mgTLTriangles& tris=triangles();
	draw_texture(timages,tris.begin(),tris.end());
}

void mgTLData::draw_texture(
	const MGTextureImages& timages,
	mgTLTriangles::const_triIterator first,
	mgTLTriangles::const_triIterator last
)const{
	double totalw=timages.total_width(), totalh=timages.total_height();
	MGVector totalwh(totalw,totalh);
	double w=timages.m_1width, h=timages.m_1height;

	const mgTLTexPoints& texpnt=tlTexPoints();
	for(mgTLTriangles::const_triIterator itri=first; itri!=last; ++itri){
		const mgTLTriangle& tri=**itri;

		mgTLTriangle::CIndexItr itri_id=tri.begin(), itri_ide=tri.end();
		MGBox stbox;
		for(;itri_id!=itri_ide; ++itri_id)
			stbox.expand(texpnt[*itri_id]);

		//tri.print(std::cout,*this);/////////**********
		//std::cout<<" stbox="<<stbox<<std::endl;
		int is,it,nw,nh; //var for MGTextureImage repitition in one MGTextureImages.
		int istotal,ittotal, nwtotal, nhtotal;//var for repitition of MGTextureImages.
		MGVector baseTotal;
		int kind=timages.compute_image_range(stbox,istotal,ittotal,nwtotal,nhtotal,is,it,nw,nh,baseTotal);
		if(kind<=0)
			continue;//If outside.

		if(kind<=2){
			//The case that the whole data is inside one MGTextureImages.
			if(nw==1 && nh==1){
				//The case that the whole data is inside one MGTextureImage.
				const MGTextureImage& teximage=timages(is,it);
				teximage.apply_texture_object();
				draw_texture_image_tri(tri,teximage.box()+baseTotal);
			}else{
				//The case that the whole data overrides multiple MGTextureImage.
				for(int iw=0; iw<nw; iw++){
					int ispiw=is+iw;
					for(int ih=0; ih<nh; ih++){
						const MGTextureImage& teximage=timages(ispiw,it+ih);
						teximage.apply_texture_object();
						MGBox bx=teximage.box()+baseTotal;
						draw_texture_image_tri_with_trim(bx,tri,bx);
					}
				}
			}
		}else{
			//The case that the whole data overrides multiple MGTextureImages.
			double basetSave=baseTotal[1];
			for(int i=0; i<nwtotal; i++){
				for(int j=0; j<nhtotal; j++){
					//Loop over total width and height.
					MGBox totalBxij(baseTotal,baseTotal+totalwh);
					MGBox stbox2=stbox&totalBxij;
					MGVector base2;
					int istotal2,ittotal2, nwtotal2, nhtotal2;//var for repitition of MGTextureImages.
					kind=timages.compute_image_range(stbox2,
						istotal2,ittotal2,nwtotal2,nhtotal2,is,it,nw,nh,baseTotal);
					assert(kind!=3);
					if(kind!=0){

					//The case that the whole data overrides multiple MGTextureImage.
					for(int iw=0; iw<nw; iw++){
						int ispiw=is+iw;
						for(int ih=0; ih<nh; ih++){
							//Loop over one MTTextureImage(not MGTextureImages).
							int itpih=it+ih;
							const MGTextureImage& teximage=timages(ispiw,itpih);
							teximage.apply_texture_object();
							MGBox bx=teximage.box()+baseTotal;
							draw_texture_image_tri_with_trim(bx,tri,bx&totalBxij);
						}
					}

					}
					baseTotal(1)+=totalh;
				}
				baseTotal(0)+=totalw;
				baseTotal(1)=basetSave;
			}
		}
	}
}