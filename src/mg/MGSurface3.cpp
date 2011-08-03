/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Unit_vector.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/Transf.h"
#include "mg/LBRepEndC.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/CParam_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"
#include "topo/Face.h"

extern "C" {
#include "cskernel/Bkdnp.h"
#include "cskernel/blg4sq.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
using namespace std;

// Implementation of MGSurface.
// MGSurface is an abstract class of 3D surface.
// Surface is represented using two parameter u and v.

//*******Intersection.************

//Compute intersection points of outer boundary curves of this face 
//with face2 to compute intersections.
//Function's return value is the number of ip's obtained(appended)
//into uvuv_list, may not be equal to the enlarged size of uvuv_list.
size_t MGSurface::isect_outcurves(
	const MGFSurface& face2,
	MGPosition_list& uvuv_list,	//intersection points will be appended.
		//One member in the list is of sdim 7,
		//and the last three elements are the ip direction vector.
	size_t id1			//id of uvuv(a member of uvuv_list).
		//uvuv(id1) for this face parameter uvuv(id2) for srf or face2 parameter.
		//id2=0 if id1=2, and id2=2 if id1=0.
)const{
	size_t pnum=perimeter_num();
	if(!pnum) return 0;

	size_t id2=2; if(id1==2) id2=0;
	const MGSurface* sf2=face2.get_surface_pointer();
	double t1[4]={param_s_u(),param_s_v(),param_e_u(),param_e_v()};
	size_t i; int im1;
	double t;
	MGCurve* boundary;
	size_t numi=0;
	for(i=0; i<pnum; i++){
		im1=i-1; if(im1<0) im1=pnum-1;
		boundary=perimeter_curve(im1);//std::cout<<*boundary<<std::endl;
		MGCSisect_list cs_list=face2.isect(*boundary);
		numi+=cs_list.entries();
		while(cs_list.entries()){
			MGCSisect cs=cs_list.removeFirst();
			MGPosition uvuv(4,cs.param_surface(),id2,0);
			t=cs.param_curve();
			if(i==1 || i==3){uvuv(id1)=t; uvuv(id1+1)=t1[i];}
			else			{uvuv(id1+1)=t; uvuv(id1)=t1[i];}
			if(id1)
				uvuv_list.append(*sf2,*this,uvuv);
			else
				uvuv_list.append(*this,*sf2,uvuv);
		}
		delete boundary;
	}
	return numi;
}

//Compute intersection points of an inner parameter line of this surface and sf2.
//The intersection point is used to compute surface to surface intersection lines.
//Function's return value is at most one intersection point on uvuv_list.
//One member of uvuv_list is (u1,v1,u2,v2), where (u1,v1) is a parameter of
//this surface and (u2,v2) is a parameter of surf.
MGPosition_list MGSurface::intersectInner(
	const MGSurface& sf2	//The second surface.
) const{
	MGPosition_list uvuv_list;
	size_t nspan1u=bdim_u()-order_u()+1, nspan1v=bdim_v()-order_v()+1;
	size_t nspan2u=sf2.bdim_u()-sf2.order_u()+1,
			nspan2v=sf2.bdim_v()-sf2.order_v()+1;
	int maximum;
	if(nspan1u<nspan1v){
		if(nspan1v<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1v<nspan2v) maximum=3; else maximum=1;
		}
	}else{
		if(nspan1u<nspan2u){
			if(nspan2u<nspan2v) maximum=3; else maximum=2;
		}else{
			if(nspan1u<nspan2v) maximum=3; else maximum=0;
		}
	}
	const MGSurface *surf1, *surf2;
	int nspan;
	int isU=1; if(maximum%2) isU=0;
	if(maximum<=1){
		surf1=this; surf2=&sf2;
		if(isU) nspan=nspan1u; else nspan=nspan1v;
	}else{
		surf2=this; surf1=&sf2;
		if(isU) nspan=nspan2u; else nspan=nspan2v;
	}
	double t0,t1;
	if(isU){ t0=surf1->param_s_u(); t1=surf1->param_e_u();}
	else{ t0=surf1->param_s_v(); t1=surf1->param_e_v();}
	double tau=(t0+t1)*.5;
	double delta=(t1-t0)/double(nspan);
	int sign=-1;
	MGCSisect_list csiList;
	for(int i=0; i<nspan; i++){
		tau+=double(i*sign)*delta;
		if(tau>t0 && tau<t1){
			MGCurve* crv=surf1->parameter_curve(isU,tau);
			csiList=surf2->isect(*crv); delete crv;
			if(csiList.size()) break;
		}
		if(csiList.size()) break;
		sign*=-1;
	}
	if(!csiList.size()) return uvuv_list;
	double u1,v1;
	u1=tau; v1=csiList.front().param_curve();
	if(!isU){ u1=v1; v1=tau;}
	const MGPosition& uv2=csiList.front().param_surface();
	MGPosition uvuv(4);
	if(surf1==this){
		uvuv(0)=u1; uvuv(1)=v1; uvuv.store_at(2,uv2,0,2);
	}else{
		uvuv(2)=u1; uvuv(3)=v1; uvuv.store_at(0,uv2,0,2);
	}
	uvuv_list.append(uvuv);
	return uvuv_list;
}

//Compute the intersection lines of this surface and srf2(both are not planes).
MGSSisect_list MGFSurface::isect_with_surf(
	MGPosition_list& uvuv_list,
	//Let a member of uvuv_list be uvuv. Then, uvuv's space dimension is
	//at least 4, and the first 2 is (u,v) of this and the next 2 is (u,v) of srf2. 
	//When uvuv's space dimension is more than 4, it indicates that the uvuv
	//is used to input approximate tangent of the intersection.
	const MGFSurface& face2	//2nd surface for the intersection.
)const{
	MGSSisect_list lst(this,&face2);
	MGSSisect ssi;
	MGPosition_list::iterator uvuv_id;
	int obtained;
	size_t n;
	size_t loop_number=0, max_loop=uvuv_list.entries();
	while(n=uvuv_list.entries()){
		if(n>=2){
			const MGPosition& Ps=uvuv_list.front();
			MGVector N1=unit_normal(Ps[0],Ps[1]), N2=face2.unit_normal(Ps[2],Ps[3]);
			double saS=N1.sangle(N2);
			if(saS<.02){
				const MGPosition& Pe=uvuv_list.back();
				N1=unit_normal(Pe[0],Pe[1]); N2=face2.unit_normal(Pe[2],Pe[3]);
				double saE=N1.sangle(N2);
				if(saS<saE) uvuv_list.reverse_order();
				//If two surfaces are more parallel at the starting point than
				//at the ending point, we change the intersection line direction.
			}
		}
		MGPosition uvuvS=uvuv_list.removeFirst();
		int m1;
again:	if(obtained=isect_start(uvuvS,uvuv_list,face2,ssi,uvuv_id,m1)){
			MGPosition uvuvE;
			if(obtained==3)
				uvuvE=uvuv_list.removeAt(uvuv_id);

			if(!ssi.is_null())
				uvuv_list.removeOn(*this,face2,ssi); //Remove uvuv that is on the ssi.

			if(uvuv_list.size() && obtained==3){
				//Check if the obtained section's end point is not a terminal point
				//but only a mid point of a section(and should be neglected).
				const MGFSurface* f1;const MGFSurface* f2;
				if(m1==0){
					f1=this; f2=&face2;
				}else{
					f1=&face2; f2=this;
				}
				if(!ssi.is_null()){
				const MGCurve& iline=ssi.line();
				if(f1->uvuvE_is_a_midpoint(*f2,m1,ssi,uvuvE)){
					uvuvS=MGPosition(7,uvuvS);
					uvuvS.store_at(4,iline.eval(iline.param_s(),1));
					goto again;
						//goto again means uvuvE is a mid point and will be neglected.
				}

				//Check to the starting direction by reversing the ssi.
				ssi.negate();
				if(uvuvE_is_a_midpoint(*f2,m1,ssi,uvuvS)){
					uvuvS=MGPosition(7,uvuvE);
					uvuvS.store_at(4,iline.eval(iline.param_s(),1));
					goto again;
				}
				}
			}
			if(obtained==7){
				if(!ssi.is_null()){
				MGSSisect_list::SSiterator i=lst.find_common(ssi);
				if(i==lst.end())
					lst.append(ssi);
				else{
					MGCurve& is1=(*i).line(); int nk1=is1.bdim();
					MGCurve& is2=ssi.line(); int nk2=is1.bdim();
					if(nk1<=nk2){
						if(nk1<nk2 || is1.length()<is2.length()){
							lst.removeAt(i);
							lst.append(ssi);
						}
					}
				}
				}
			}else
				if(!ssi.is_null())
					lst.append(ssi);
		}else{
			uvuv_list.append(uvuvS);
		}
		loop_number++;
		if(loop_number>=max_loop) break;
	}
	return lst;
}

//Default intersection program of MGSurface.
//It is assumed that both this and srf2 are not a plane.
MGSSisect_list MGSurface::intersect(const MGSurface& srf2) const{
	MGSSisect_list lst(this,&srf2);
	if(!has_common(srf2)) 
		return lst;

	MGPosition_list uvuv_list;
	intersect12Boundary(srf2,uvuv_list);
	if(!uvuv_list.size())
		uvuv_list=intersectInner(srf2);

	//Compute intersection line using isect_with_surf.
	lst=isect_with_surf(uvuv_list,srf2);

	return lst;
}

//Default intersection program of MGSurface with a plane.
MGSSisect_list MGSurface::intersectPl(const MGPlane& srf2)const{
	MGSSisect_list lst(this,&srf2);
	MGPosition_list uvuv_list;
	if(!box().cutting(srf2))
		return lst;

	intersect12Boundary(srf2,uvuv_list);
	if(!uvuv_list.size())
		uvuv_list=intersectInner(srf2);
	//Compute intersection lines using isect_with_plane.
	lst=isect_with_plane(uvuv_list,srf2,srf2);
	return lst;
}

//isect_dt computes incremental values du and dv for the intersection
//computation at parameter position (u,v).
void MGFSurface::isect_dt(
	double u, double v, double& du, double& dv,
	double acuRatio	//acuracy ratio.
) const{
	const MGKnotVector& tu=knot_vector_u();
	const MGKnotVector& tv=knot_vector_v();
	size_t id,k, k_half1,k_half2;
	double alfa=isect_dt_coef(0);

	id=tu.locate(u);
	k=tu.order();
	k_half1=k/2; k_half2=k-k_half1-1;
	du=(tu(id+k_half1)-tu(id-k_half2));
//	if(k<4) k=4;//k=4;
	du=du/double(k)*alfa;
	du*=acuRatio;

	id=tv.locate(v);
	k=tv.order();
	k_half1=k/2; k_half2=k-k_half1-1;
	dv=(tv(id+k_half1)-tv(id-k_half2));
	if(k<4) k=4;//k=4;
	dv=dv/double(k)*alfa;
	dv*=acuRatio;
}

#define DELTA .05;
//isect_direction() is used by isect_startPt() to define which constant
//parameter line should be used to compute intersection, and what
//incremental value be used for the parameter.
//Function's return value is direction to get next intersection(with dt).
//When =1: u=const direction, =0: v=const, =-1: cannot get intersection.
int MGFSurface::isect_direction(
	const MGFSurface& sf2,	//Second surface for the intersection.
	size_t m1,		//id of uvuvS that indicates this surface's parameter
		//position in uvuvS. (uvuvS(m1), uvuvS(m1+1))=(u,v) of this surface.
	MGPosition& uvuvS,//start parameter (u,v) pair of this surface and sf2.
	double& du,	//Incremental value of the parameter kind of kdt will be output.
	double& dv, //Right dt will be output according to the function's output =0,1.
	double acuRatio	//acuracy ratio.
)const{
	size_t m1p1=m1+1;
	MGPosition uvuv(4),uv;
	double& u0=uvuvS(m1); double& u1=uvuv(m1);
	double& v0=uvuvS(m1p1); double& v1=uvuv(m1p1);
	double lzero=MGTolerance::line_zero();
	MGTolerance::set_line_zero(lzero*acuRatio);

	int kdt=-1;

//Normalize the starting point direction.
	int obtained, kdt2;
	size_t pnum;
	double du10, dv10;
	MGVector direction;
	if(uvuvS.sdim()>4){
	//When direction specified.
		direction.resize(3);
		direction.store_at(0,uvuvS,4,3);
	}
	if(!direction.is_null() && !direction.is_zero_vector()){
		kdt=isect_direction_with_direction(u0,v0,direction,du,dv);
		du10=du*DELTA; dv10=dv*DELTA;
		double u=u0, v=v0; if(kdt) u+=du10; else v+=dv10;
		if(!in_range(u,v)){
			kdt=!kdt;
			u=u0, v=v0; if(kdt) u+=du10; else v+=dv10;
			if(!in_range(u,v)) {kdt=-1; goto ret10;}
		}
	
		size_t m2=2; if(m1) m2=0;
		double u20=uvuvS(m2); double v20=uvuvS(m2+1);
		size_t pnum;
		if(sf2.on_a_perimeter(u20,v20,pnum)){
			size_t degu=(pnum+1)%2, degv=pnum%2;
			MGVector T=sf2.eval(u20,v20,degu,degv);if(pnum>=2) T*=-1.;
			if(!T.parallel(direction)){
				MGVector N=sf2.normal(u20,v20);
				if(N%(T*direction)<0.) {kdt=-1; goto ret10;}
			}
		}

		MGVector P00=eval(u0,v0);
		if(obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv)){
			MGVector dir1(eval(u1,v1)-P00);
			double angle1=dir1.cangle(direction), angle2;
			int kdt3=!kdt;
			MGPosition uvuv2(4);
			double u=u0, v=v0; if(kdt3) u+=du10; else v+=dv10;
			if(in_range(u,v)){
				if(isect_start_incr(sf2,uvuvS,kdt3,du10,dv10,m1,uvuv2)){
					MGVector dir2(eval(uvuv2[m1],uvuv2[m1p1])-P00);
					angle2=dir2.cangle(direction);
					if(angle1<angle2){
						kdt=kdt3;
						uvuv=uvuv2;
					}
					double du2=u1-u0; du=fabs(du); if(du2<0.) du*=-1.;
					double dv2=v1-v0; dv=fabs(dv); if(dv2<0.) dv*=-1.;
				}
			}
		}else{
			kdt=!kdt;
			obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv);
		}
		if(!obtained) {
			kdt=-1; goto ret10;
		}

		//Check if intersection is not going to opposite direction as uvS'direction.
		MGVector dir(eval(u1,v1)-P00);
		double cang=direction.cangle(dir);
		if(cang<=MGTolerance::angle_zero()) {kdt=-1; goto ret10;}
		kdt2=-1;//kdt will be defined in isect_inner_dt() below.
	}else{
		isect_dt(u0,v0,du,dv,acuRatio);
		if(on_a_perimeter(u0,v0,pnum)){
		//When on a surface perimeter.
			if(pnum==3 || pnum==1){
				kdt=1;
				if(pnum==1) du*=-1.;
			}else{
				kdt=0;
				if(pnum==2) dv*=-1.;
			}
			du10=du*DELTA; dv10=dv*DELTA;
			if(!(obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv))){
				kdt=!kdt;
				double u=u0, v=v0;
				if(kdt) u+=du10; else v+=dv10;
				if(in_range(u,v))
					obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv);
				else
					obtained=0;
				if(!obtained){
					du10*=-1.; dv10*=-1.;
					if(in_range(u,v))
						obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv);
				}
			}
		}else{
		//When direction not specified.
			kdt=1;
			du10=du*DELTA; dv10=dv*DELTA;
			if(!(obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv))){
				du10*=-1.; dv10*=-1.;
				if(!(obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv))){
					kdt=0;
					if(!(obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv))){
						du10*=-1.; dv10*=-1.;
						obtained=isect_start_incr(sf2,uvuvS,kdt,du10,dv10,m1,uvuv);
					}
				}
			}
		}
		kdt2=-1;//kdt will be defined in isect_inner_dt() below.
		if(!obtained) {
			kdt=-1; goto ret10;
		}
	}

	uv=MGPosition(2,uvuvS,0,m1);
	kdt=kdt2;
	dv=v1-v0; du=u1-u0;
	isect_inner_dt(0,uv,du,dv,kdt,acuRatio);
//	du*=10.; dv*=10.;
ret10:
	MGTolerance::set_line_zero(lzero);
	return kdt;
}

//uvuvE_is_a_midpoint() is used by isect_with_surf() to determine if
//uvuvE is a mid point at the intersectio and should be neglected.
bool MGFSurface::uvuvE_is_a_midpoint(
	const MGFSurface& f2,
	int m1,		//id of uvuvE that indicates this surface's parameter
		//position in uvuvE. (uvuvE(m1), uvuvE(m1+1))=(u,v) of this surface.
	const MGSSisect& ssi, //ssi obtained so far.
	const MGPosition& uvuvE//start parameter (u,v) pair of this surface and sf2.
)const{
	const MGSurface* sf1=get_surface_pointer();
	const MGSurface* sf2=f2.get_surface_pointer();
	const MGFace* face1=get_face_pointer();
	const MGFace* face2=f2.get_face_pointer();

	const MGCurve& iline=ssi.line();
	const MGCurve* iuvline;
	int m2;
	if(m1==0){
		m2=2;
		iuvline=&(ssi.param1());
	}else{
		m2=0;
		iuvline=&(ssi.param2());
	}

	const int m1p1=m1+1;
	MGPosition uvuv(4);
	const double& u0=uvuvE(m1); const double& u1=uvuv(m1);
	const double& v0=uvuvE(m1p1); const double& v1=uvuv(m1p1);
	double t=iline.param_e();
	MGVector direction=iline.eval(t,1), duv=iuvline->eval(t,1);

	double du,dv;
	int kdt=isect_direction_with_direction(u0,v0,direction,du,dv);
	if(duv[0]<0.) du=-fabs(du); if(duv[1]<0.) dv=-fabs(dv);

	double du10=du*DELTA; double dv10=dv*DELTA;
	double u=u0, v=v0; if(kdt) u+=du10; else v+=dv10;
	if(!sf1->in_range(u,v)){
		kdt=!kdt;
		u=u0, v=v0; if(kdt) u+=du10; else v+=dv10;
		if(!sf1->in_range(u,v))
			return false;
	}

	double u20=uvuvE(m2); double v20=uvuvE(m2+1);
	size_t pnum;
	if(sf2->on_a_perimeter(u20,v20,pnum)){
		size_t degu=(pnum+1)%2, degv=pnum%2;
		MGVector T=sf2->eval(u20,v20,degu,degv);if(pnum>=2) T*=-1.;
		if(!T.parallel(direction)){
			MGVector N=sf2->normal(u20,v20);
			if(N%(T*direction)<0.)
				return false;
		}
	}

	if(face1){
		if(!face1->in_range(u,v)){
			kdt=!kdt;
			u=u0, v=v0; if(kdt) u+=du10; else v+=dv10;
			if(!face1->in_range(u,v))
				return false;
		}
	}
	if(!isect_start_incr(f2,uvuvE,kdt,du10,dv10,m1,uvuv)){
		kdt=!kdt;
		u=u0, v=v0; if(kdt) u+=du10; else v+=dv10;
		if(!sf1->in_range(u,v))
			return false;
		if(!isect_start_incr(f2,uvuvE,kdt,du10,dv10,m1,uvuv))
			return false;;
	}

	//Check if intersection is not going to opposite direction as uvS'direction.
	MGVector dir(eval(u1,v1)-eval(u0,v0));
	if(dir.is_zero_vector())
		return false;
	double cang=direction.cangle(dir);
	if(cang<=MGTolerance::angle_zero())
		return false;

	if(face1)
		if(!face1->in_range(u1,v1))
			return false;
	if(face2)
		if(!face2->in_range(uvuv[m2], uvuv[m2+1]))
			return false;

	return true;
}

//isect_direction_with_direction() is used by isect_direction() to temporalily define
//which constant parameter line should be used to compute intersection,
//and what incremental value be used for the parameter.
//Function's return value isect_direction_with_direction() is 1 for u=const parameter
//line, and 0 for v=const parameter line.
int MGFSurface::isect_direction_with_direction(
	double u, double v,		//start parameter (u,v) of this surface.
	const MGVector& tangent,//To indicate which direction isect line
							//should march toward.
	double& du,				//Incremental value sign of the parameter kind of
	double& dv)const		//isect_direction_with_direction will be output.
{
	MGUnit_vector Su=eval(u,v,1,0), Sv=eval(u,v,0,1);
	int is_u=0;
	double err=MGTolerance::wc_zero(); err*=err*.1;
	double Sut=Su%tangent, Svt=Sv%tangent;
	if(Su%Su<err) is_u=1;
		//If length of Su is zero, use u=const parameter line.
	else if(!(Sv%Sv<err)){
		//If both length of Su and of Sv are not zero, 
		//define direction using angle between tangent and Su, Sv.
		if(fabs(Svt)<fabs(Sut)) is_u=1;
	}
	//If length of Sv is zero, is_u=0(use v=const parameter line).

	isect_dt(u,v,du,dv);
	if(Sut<0.) du*=-1.;
	if(Svt<0.) dv*=-1.;
	return is_u;
}

#define sect_id_max 5
const static size_t sect_div_id_max=sect_id_max;	//Maximum id of sect_div.
//***NOTE*** sect_div_id_max is the number isect_div_id_max() returns.
//isect_dt_coef provides coef of how fine parameter increment should be,
//given num of ip computed so far.
double MGFSurface::isect_dt_coef(
	  size_t n) const{
//sect_div[] are coefficients that determine how fine sectional line be
//approximated. Let maximum order of the two surfaces be k, then one span is
//divided into approximately sect_div[6]*2./k.
//First number of sect_div[.] is small since points that are near to
//a border line are computed fine.
//const static double sect_div[sect_id_max+1]={.15, 0.165, 0.18, 0.195, 0.21, 0.225};
//const static double sect_div[sect_id_max+1]={.15, 0.15, 0.15, 0.15, 0.15, 0.15};
const static double sect_div[sect_id_max+1]={.28, 0.28, 0.29, 0.29, 0.30, 0.30};
//const static double sect_div[sect_id_max+1]={.28, 0.28, 0.28, 0.28, 0.28, 0.28};
//const static double sect_div[sect_id_max+1]={.3, 0.3, 0.3, 0.3, 0.3, 0.3};
//const static double sect_div[sect_id_max+1]={.13, 0.17, 0.21, 0.24, 0.26, 0.28};

	size_t div_id=n; if(div_id>sect_div_id_max) div_id=sect_div_id_max;
	return sect_div[div_id];
}

//isect_div_id_max is maximum id of array of sect_div defined in
//isect_dt_coef. That is, isect_div_id_max+1 is the length of the array
//sect_div.
size_t MGFSurface::isect_div_id_max()const{
	return sect_div_id_max;
}

//"isect_inner_dt" is a dedicated function of isect_startPt,
// comutes adequate incremental parameter value(du,dv) and parameter line kind
//kdt(u=const or v=const).
void MGFSurface::isect_inner_dt(
	size_t n,	//num of i.p. obtained so far(not include uvnow).
	const MGPosition& uvnow,//intersection point obtained last(of this).
	double& du, double& dv,	//incremental length from previous to uvnow is input.
				//New du or dv will be output according to kdt's return value.
	int& kdt,	//Parameter kind used so far is input, will be output as:
				//=1:parameter line kind(u=const), =0: v=const,
				//=-1:should halt computation since incremental value is zero.
	double acuRatio	//Accurate ratio.
) const{
	double uerr=param_error_u()*acuRatio;
	double verr=param_error_v()*acuRatio;
	double abdu=fabs(du), abdv=fabs(dv);
	double ratio=acuRatio;
	if(abdu<=uerr){
		if(abdv<=verr){kdt=-1; return;}
		else kdt=0;
	}else if(abdv<=verr) kdt=1;
	else{
		MGVector dfdu=eval(uvnow,1,0), dfdv=eval(uvnow,0,1);
		double fuu=dfdu%dfdu, fuv=dfdu%dfdv, fvv=dfdv%dfdv;
		double dffu=fuu*du+fuv*dv, dffv=fuv*du+fvv*dv;
		double dffubyfvv=dffu*dffu*fvv, dffvbyfuu=dffv*dffv*fuu;
		if(kdt==-1){
			if(dffubyfvv>=dffvbyfuu) kdt=1; else kdt=0;
		}else if(dffubyfvv>=dffvbyfuu*1.8){
			//u=const and v-varying parameter line.
			if(abdu>uerr*4.) kdt=1;
		}else if(dffvbyfuu>=dffubyfvv*1.8){
			//v=const and u-varying parameter line.
			if(abdv>verr*4.) kdt=0;
		}
		if(kdt) ratio*=dffubyfvv; else ratio*=dffvbyfuu;
		ratio/=(dffubyfvv+dffvbyfuu);
	}

//Define new dt, kdt.
	const MGKnotVector* t;
	double dtold;
	if(kdt){
		t=&(knot_vector_u());
		dtold=du;
	}else{
		t=&(knot_vector_v());
		dtold=dv;
	}
	size_t k=t->order();
	size_t k_half1=k/2; size_t k_half2=k-k_half1-1;
	size_t id=(*t).locate(uvnow[(kdt+1)%2]);
	double dt=(*t)(id+k_half1)-(*t)(id-k_half2);
	if(dtold<0.) dt=-dt;
	if(k<4) k=4;
	dt/=double(k); dt*=isect_dt_coef(n)*ratio;
	if(n){
		//When this is not the 1st call of isect_inner_dt,
		//dt must not exceed twice or half of the old dt.
		double dtr1=dt*dt, dtr2=dtold*dtold;
		if(dtr1 > 2.*dtr2) dt=dtold*2.;
		else if(dtr1 < .5*dtr2) dt=dtold*.5;
	}
	if(kdt) du=dt; else dv=dt;
	return;
}

//Check if the intersection line lineb's start and end tangent vectors are accurate
//enough. If they do not have enough accuracy, isect_start_tan returns
//which end did not have the accuracy.
// 1:start, 2:end, 3:start and end. If both ends had enough accuracy, returns 0.
int isect_start_tan(
	const MGFSurface& sf1,
	const MGFSurface& sf2,
	const MGLBRep& lineb,
	MGVector* Tse[2]	//If an end had not the accuracy, accurate tangent
						//will be output. Tse[0]:start, Tse[1]:end.
){
	double t[2]={lineb.param_s(), lineb.param_e()};
	double azero=MGTolerance::angle_zero();
	double mzero=MGTolerance::mach_zero()*1.5;
	int ng=0;
	for(size_t j=0; j<2; j++){//Loop for start and end point.
		MGVector EndP=lineb.eval(t[j]);
		MGVector N1=sf1.unit_normal(EndP[0],EndP[1]),
			N2=sf2.unit_normal(EndP[2],EndP[3]);
		//double sa2=N1.sangle(N2);

		MGVector& tan=*(Tse[j]);
		tan=lineb.eval(t[j],1); //cout<<tan<<endl;////
		MGVector T, Tline(3,tan,0,4);
		MGVector N3=N1*N2;
		if(N3.len()<=mzero){
			MGVector N;
			if(N1%N2>=0.) N=N1+N2; else N=N1-N2;
			T=N*(Tline*N);
		}else{
			if(N3%Tline>0.) T=N3;
			else        T=-N3;
		}
		double sa=T.sangle(Tline);

		if(sa>=azero){
			ng+=j+1;
			T.set_unit();
			T*=Tline.len();
			tan.store_at(4,T,0,3); //cout<<tan<<endl;////////
		}
	}

	return ng;
}

//Update lineb so as to have the tangent tan for start or end according to ngtan.
void isect_start_adjustSE(
	int ngtan,	//Return value of isect_start_tan, indicates which end be
				//adjusted. =1: start tangent, =2:end tangent, =3:both tangent.
	MGNDDArray& tau,	//data point abcissa.
	MGBPointSeq& point,	//data point ordinate.
	MGLBRep& lineb,	//line b-rep obtained so far. tangent adjusted new B-rep
					//will be output.
	MGVector* tan[2]//accurate tangent data obtained by isect_start_tan.
){
	int n=point.length();
	MGENDCOND beginc=MGENDC_NO, endc=MGENDC_NO;
	if(ngtan%2){
		beginc=MGENDC_1D;
		tau(1)=tau(0);
		point.store_at(1,*(tan[0]));
	}
	if(ngtan>=2){
		endc=MGENDC_1D;
		tau(n-2)=tau(n-1);
		point.store_at(n-2,*(tan[1]));
	}
	MGBPointSeq& bp=lineb.line_bcoef(); bp.resize(n,point.sdim());
	double* ktv=lineb.knot_data();
	int irc=bp.capacity();
	double* work=new double[n*9];//cout<<lineb<<endl;
	int pointSize=point.capacity(), pointDim=point.sdim(), error;
	error=blg4sq_(beginc,endc,tau.data(),point.data(),pointSize,n,
		pointDim,irc,work,ktv,bp.data());
	//cout<<lineb<<endl;
	delete[] work;
}

//Compute the maximum difference square between the intersection line and the two surfaces,
//this and sr2. The difference evaluation is done at the data point tau.
double MGFSurface::isect_start_dif(
	const MGNDDArray& tau,	//data points
	const MGLBRep& line,		//the intersection line of this and sf2.
	const MGFSurface& sf2		//second surface.
)const{
//Check if the obtained intersection lines are on the two surfaces.	
	double devi1=-1., devi2=-1.;//get the maximum deviations in devi1,2.
	size_t n=tau.length();
	for(size_t i=1; i<n; i++){
		double t=tau(i-1)+tau(i); t*=.5;
		MGVector iP=line.eval(t);
		MGVector Qnow(3,iP,0,4), P1=eval(iP[0],iP[1]), P2=sf2.eval(iP[2],iP[3]);
		MGVector def1(P1-Qnow), def2(P2-Qnow);
		double len1=def1%def1, len2=def2%def2;
		if(devi1<len1) devi1=len1;
		if(devi2<len2) devi2=len2;
	}
	if(devi1<devi2) devi1=devi2;
	return devi1;
}

//isect_start compute one intersection line of two surfaces, this and sf2,
// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2)
// and direction in uvuv(4-6)(optionally).
//isect_start halts the computation when intersection reached to a boundary of
// this or sf2, or reached to one of the points in uvuv_list.
//The function's return value is:
// =0: Intersection was not obtained.
// !=0: Intersection was obtained as follows:
//    =1: End point is a point on a perimeter of one of the surfaces.
//    =3: End point is one of boundary points in uvuv_list.
//    =4: End point is the starting point.
//    =7: isect_start halted the computation since intersection was lost
//     during the computation.
int MGFSurface::isect_start(
	const MGPosition& uvuv_startIn, //Starting point of the intersection line.
	MGPosition_list& uvuv_list,	//isect_start will halt when ip reached one of 
		//the point in uvuv_list. isect_start does not change uvuv_list(actually
		//uvuv_list is const.)
	const MGFSurface& sf2,	//2nd surface.
	MGSSisect& ssi,			//Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id,
			//When the end point of ip was one of the points of uvuv_list,
			//uvuv_list's iterator of the point will be returned, that is, 
			//when the function's return value was 3 or 5.
			//When was not a point of uvuv_list, end() of uvuv_list will be
			//returned.
	int& m1	//id that indicates which surface was used as the main surface.
			//m1=0: this surface, m1=2: sf2.
)const{
	ssi.set_null();
	MGBPointSeq point;
	MGLBRep lineb;
	int obtained;
	double lineError=MGTolerance::line_zero();
	//MGTolerance::set_line_zero(lineError*.5);
	double err2=lineError*.5; err2*=err2;
	double acuRatio=1.;
	double errwcSave=MGTolerance::wc_zero();
	double errwc=errwcSave;//*.5;

	double deviOld;
	int accuracy_not_improved=0;
	for(size_t numRepeat=0; numRepeat<8; numRepeat++){

	MGTolerance::set_wc_zero(errwc);
	MGPosition uvuvE;
	obtained=isect_startPt(uvuv_startIn,uvuv_list,sf2,acuRatio,point,uvuv_id,m1);
	if(obtained==0 || obtained==1){ obtained=0; break;}

	int n=point.length(), pointSize=point.capacity(), pointDim=point.sdim();
	MGNDDArray tau(n);
	double tauati=0.; tau(0)=tauati;
	MGVector P1=eval(point(0,0),point(0,1)),P2=sf2.eval(point(0,2),point(0,3));
	MGVector Qnow,Qpre=(P1+P2)*.5;
	for(int i=0; i<n; i++){
		P1=eval(point(i,0),point(i,1));
		P2=sf2.eval(point(i,2),point(i,3));
		Qnow=(P1+P2)*.5;
		point.store_at(i,Qnow,4,0,3);
		tau(i)=tauati+=(Qnow-Qpre).len();
		Qpre=Qnow;
	}
	if(tau[n-1]<=errwcSave){
		MGTolerance::set_wc_zero(errwcSave);
		return obtained;
	}

	int IMLT=1;
	double ratio=4.;
	bkdnp_(&n,tau.data(),point.data(),pointSize,pointDim,IMLT,ratio);		
	//std::cout<<point<<std::endl;;
	tau.set_length(n); point.set_length(n);
	int error;
	int k=4;//Order of the intersection is always 4.
	if(k>n) k=n;
	MGKnotVector knotv(tau,k);
//	cout<<point<<endl;
	lineb=MGLBRep(tau,point,knotv,error);// std::cout<<" original="<<lineb<<endl;
	if(error) { obtained=0; break;}

//Check if the obtained intersection lines are on the two surfaces.	
	double devi=isect_start_dif(tau,lineb,sf2);//get the maximum deviations in devi.
	if(devi<=err2){
		/*if(n<4) break;
		MGVector tan1,tan2;
		MGVector* tan[2]={&tan1, &tan2};
		int ngtan=isect_start_tan(*this,sf2,lineb,tan);
		if(ngtan){
			MGLBRep lineb2(n,k,lineb.sdim());
			isect_start_adjustSE(ngtan,tau,point,lineb2,tan);
			//Adjust start and end tangent line by value tan.
			double devi2=isect_start_dif(tau,lineb2,sf2);
			//get the maximum deviations in devi.
			if(devi2<devi*4) lineb=lineb2;
		}*/
		break;
	}
	if(numRepeat>=2){
		if(devi>deviOld*.9) accuracy_not_improved++;
		if(accuracy_not_improved>=3) break;//Halt the loop if the accuracy was not improved.
	}
		
	deviOld=devi;
	acuRatio*=.3;	//Try again by raising accuracy.
	errwc*=.5;

	}

	if(obtained){
		//cout<<" before remove_knot="<<lineb<<endl;
		lineb.remove_knot(4,3);
		//cout<<" after remove_knot="<<lineb<<endl;
		MGLBRep* lineuv1=new MGLBRep(2,lineb,0,0);
		MGLBRep* lineuv2=new MGLBRep(2,lineb,0,2);
		MGLBRep* linexyz=new MGLBRep(3,lineb,0,4);
		ssi=MGSSisect(linexyz,lineuv1,lineuv2);
	}

	//MGTolerance::set_line_zero(lineError);	//Restore the error.
	MGTolerance::set_wc_zero(errwcSave);
	return obtained;
}

//"isect_start_boundary" is a dedicated function of isect_start.
//"isect_start_boundary" computes one intersection point of two surfaces,
// this surface's parameter line at uv+dt(according to kdt) and sf2,
// given previous intersetion point(uv) and incremental value dt.
// "isect_start_boundary" is used only when no intersection found at dt.
//Function's return value is 0: ip not found.
//                           2: ip found as intersection line boundary point.
int MGFSurface::isect_start_boundary(
const MGFSurface& sf2,	//2nd surface b-rep.
const MGPosition& uvuv_pre,	//Starting parameter values of ip. 
int kdt,			//kdt=true: u=const parameter line,
					//    else: v=const parameter line of this surface.
double du, double dv,//Incremental parameter length.
size_t m1,		//id of parameter of this surface in uvuv_pre or uvuv_now.
MGPosition& uvuv_now//New parameter values of ip will be output.
) const{
	int obtained;
	// i.p. will be an end of the intersection line.
	// Compute correct end point of the line by bi-section method.
	double dt;
	if(kdt) dt=du; else dt=dv;
	double left=0.,right=dt, middle=dt/2.;
	double dfdt;
	if(kdt) dfdt=(eval(uvuv_pre[m1], uvuv_pre[m1+1],0,1)).len();
	else    dfdt=(eval(uvuv_pre[m1], uvuv_pre[m1+1],1,0)).len();
	double error=MGTolerance::wc_zero();
	int loop=0;
	//dt is narrowed by half every time.
	while(dfdt*fabs(right-left)>error && loop++<32){
		if(kdt)
			obtained=isect_start_incr(sf2,uvuv_pre,kdt,middle,dv,m1,uvuv_now);
		else
			obtained=isect_start_incr(sf2,uvuv_pre,kdt,du,middle,m1,uvuv_now);
		if(obtained){
			obtained=2; left=middle;
		}else right=middle;
		middle=(left+right)/2.;
	}
	return obtained;
}

//isect_start_incr compute one intersection point of two surfaces,
//this surface's parameter line at uv1+dt(according to kdt) and sf2,
//given previous intersetion point(uv1,uv2) and incremental value dt.
//Here uv1 is uvuv_pre(m1, m1+1), and uv2 is other two values of uvuvpre.
//isect_start_incr is a dedicated function of isect_start
//and isect_start_boundary.
//Function's return value is true: if ip found,
//                           false: if ip not found.
int MGFSurface::isect_start_incr(
const MGFSurface& sf2,	//2nd surface b-rep.
const MGPosition& uvuv_pre,	//Starting parameter values of ip. 
int kdt,			//kdt=true: u=const parameter line,
					//    else: v=const parameter line of this surface.
double du, double dv,//Incremental parameter length.
size_t m1,		//id of parameter of this surface in uvuv_pre or uvuv_now.
MGPosition& uvuv_now//New parameter values of ip will be output.
) const{
	double* t;
	double u,v;	double dt;
	if(kdt){
		t=&v; dt=dv*1.5;
	}else{
		t=&u; dt=du*1.5;
	}
	size_t m2=2; if(m1==2) m2=0;
	MGPosition uv1i(2,uvuv_pre,0,m1);
	MGCurve* pline=isect_incr_pline(uv1i,kdt,du,dv,u,v);
	//cout<<endl<<"In isect_start_incr:"<<(*pline);//////////////////
	MGPosition uv2(2);
	MGPosition uv2i(2,uvuv_pre,0,m2);
	int obtained=0;
/////////////2007/01/08
	double ts=pline->param_s(), te=pline->param_e(), tt;
	if((*t-ts)<fabs(dt)){
		obtained=sf2.isect_guess(*pline,uv2i,ts,uv2,tt);
	}else if((te-*t)<fabs(dt)){
		obtained=sf2.isect_guess(*pline,uv2i,te,uv2,tt);
	}
	if(obtained){
		*t=tt;
	}else{
/////////////2007/01/08
		tt=*t;
		if(!(obtained=sf2.isect_guess(*pline,uv2i,tt,uv2,*t))){
			for(size_t incr=1; incr<=4;incr++){
				delete pline;
				pline=isect_incr_pline(uv1i,kdt,du,dv,u,v,incr);
				//cout<<"In isect_start_incr:"<<(*pline);///////
				obtained=sf2.isect_guess(*pline,uv2i,tt,uv2,*t);
				if(obtained) break;
			}
		}
	}
	if(obtained){
		uvuv_now(m1)=u; uvuv_now(m1+1)=v;
		uvuv_now(m2)=uv2[0]; uvuv_now(m2+1)=uv2[1];
	}
//	cout<<pline->eval(*t)<<","<<eval(u,v)<<endl;
	delete pline;
	return obtained;
}

//isect_startPt compute an array of parameter value pairs of this surf and sf2
//for one intersection line of the two surfaces,
// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2)
// and direction in uvuv_startI(4-6) (optionally).
//isect_startPt is a dedicated function for isect_start.
//isect_startPt halts the computation when intersection
//reached to a boundary of this or sf2, or reached to one of the points
//in uvuv_list.
//The function's return value is:
//  =0: Intersection was not obtained.
// !=0: Intersection was obtained as follows:
//  =1: End point is a point on a perimeter of one of the surfaces.
//  =3: End point is one of boundary points in uvuv_list.
//  =4: End point is the starting point.
//  =7: isect_startPt halted the computation since intersection was lost
//     during the computation.
int MGFSurface::isect_startPt(
	const MGPosition& uvuv_startIn, //Starting point of the intersection line.
	MGPosition_list& uvuv_list,	//isect_startPt will halt when ip reached one of 
		//the point in uvuv_list. isect_startPt does not change uvuv_list(actually
		//uvuv_list is const.) uvuv's space dimension is at least 4,
		//and the first 2 is (u,v) of this and the next 2 is (u,v) of sf2. 
	const MGFSurface& sf2,	//2nd surface.
	double acuRatio,//Accurate ratio, should be decreased by multiplyng .2
		//(or a number less than 1.).
	MGBPointSeq& point,		//Surface-surface intersection parameter values
		//will be returned as:point(.,0) and point(.,1) for(u,v) of this surface
		//point(.,2) and point(.,3) for(u,v) of this surface.
		//point has the dimension of(.,7).
	MGPosition_list::iterator& uvuv_id,
			//When the end point of ip was one of the points of uvuv_list,
			//uvuv_list's iterator of the point will be returned, that is, 
			//when the function's return value was 3 or 5.
			//When was not a point of uvuv_list, end() of uvuv_list will be
			//returned.
	int& m1	//id that indicates which surface was used as the main surface.
			//m1=0: this surface, m1=2: sf2.
)const{
//cout<<endl<<"isect_startPt uvuv_startIn="<<uvuv_startIn<<endl;
//cout<<",uvuv_list="<<uvuv_list;

	//Prepare work area for isect_start_boundary and isect_start_incr.
	size_t i, len1=isect_area_length(), len2=sf2.isect_area_length();
	size_t lmax=len1+len2;
	len1/=2; len2/=2;
	size_t ncd1=coef_sdim(), ncd2=sf2.coef_sdim();
	if(len1<len2) len1=len2;
	if(ncd1<ncd2) ncd1=ncd2;

	MGBox pr1=param_range(), pr2=sf2.param_range();
	double tol=MGTolerance::rc_zero();
	double dtmin1[2]=
		{pr1(1).length().value()*tol,pr1(0).length().value()*tol};
	double dtmin2[2]=
		{pr2(1).length().value()*tol,pr2(0).length().value()*tol};
	//dtmin1,2 are minimum dt of below.

	tol*=(1000.*acuRatio); if(tol > .001) tol=.001;
	double tTol[4]={
		tol*knot_vector_u().param_span(),
		tol*knot_vector_v().param_span(),
		tol*sf2.knot_vector_u().param_span(),
		tol*sf2.knot_vector_v().param_span()
	}, tTolNow[4];	//tTol is the box tolerance to 
		//terminate the marching. if a point in uvuv_list is within this
		//tolerance of current intersection point, the marching will be
		//terminated.

	MGPosition uvuv_start(4,uvuv_startIn), uvuv_now(4);
	MGPosition uv1i(2,uvuv_start,0,0);
	MGPosition uv2i(2,uvuv_start,0,2);
	MGPosition_list::iterator uvuvi=uvuv_list.begin(), uvuve=uvuv_list.end();
	uvuv_id=uvuve;//Initialize the return value of uvuv_id.
	//Get maximum dt to not pass too far.
	double dt_max1=-1., dt_max2=-1.;
	while(uvuvi!=uvuve){
		double uvlen=(uv1i-MGPosition(2,*uvuvi)).len();
		if(dt_max1<0. || dt_max1>uvlen) dt_max1=uvlen;
		uvlen=(uv2i-MGPosition(2,*uvuvi++,0,2)).len();
		if(dt_max2<0. || dt_max2>uvlen) dt_max2=uvlen;
	}
	dt_max1/=4.; dt_max2/=4.;

//Define initial parameter direction.
	int kdt,kdt2;;
		//When kdt=1: dt=u-value(isect of u=const parameter line)
		//     kdt=0: dt=v-value(isect of v=const parameter line)
	double du,dv, du3,dv3; 
	const MGFSurface* surf1; const MGFSurface* surf2;
	double* dtmin; double dtmax;
//m1,2 are used to indicate where in line, uvuv_now, and uvuv_passed
//to store of the parameter data. m1 is the id for surf1.
//*surf1 is used as parameter line and *surf2 is as surface in isectn_incr.
//(surf1,m1, kdt, dt) are one pair.

	size_t pnum;
	//Check if input point uv1i or uv2i is a degenerated point.
	MGVector S1u0=eval(uv1i.ref(0), uv1i.ref(1), 1,0);
	MGVector S1v0=eval(uv1i.ref(0), uv1i.ref(1), 0,1);
	MGVector S2u0=sf2.eval(uv2i.ref(0), uv2i.ref(1), 1,0);
	MGVector S2v0=sf2.eval(uv2i.ref(0), uv2i.ref(1), 0,1);

	int use_srf1=0;
		//=0:can use both, =1:use this surf, =2: use sf2,
		//=3: do not use on_a_perimeter()
	if(S1u0.is_zero_vector() || S1v0.is_zero_vector()){
		//uv1i is a degenerated point.
		if((S2u0.is_zero_vector() || S2v0.is_zero_vector()))
			use_srf1=3;
			//both uv1i and uv2i are degenerated points.
		else
			use_srf1=2;
	}else{
		if(S2u0.is_zero_vector() || S2v0.is_zero_vector())
			use_srf1=1;//uv2i is degenerated point.
	}

	//Check if starting point on this surface's perimeter.
	if(use_srf1==1) m1=0;
	else if (use_srf1==2) m1=2;
	else{
		const MGSBRep* thisS=dynamic_cast<const MGSBRep*>(this->get_surface_pointer());
		const MGSBRep* thatS=dynamic_cast<const MGSBRep*>(sf2.get_surface_pointer());
		if(thisS && !thatS) m1=0;
		else if(!thisS && thatS) m1=2;
		else if(use_srf1==0 && on_a_perimeter(uv1i(0),uv1i(1),pnum)) m1=0;
		else if(use_srf1==0 && sf2.on_a_perimeter(uv2i(0),uv2i(1),pnum)) m1=2;
		else m1=0;
	}

	if(m1==0){
		surf1=this; surf2=&sf2; dtmin=dtmin1; dtmax=dt_max1;
	}else{
		surf2=this; surf1=&sf2; dtmin=dtmin2; dtmax=dt_max2;
	}

//*****************Define initial direction.**********************
	MGPosition uvuv_start2(uvuv_startIn);
	kdt=surf1->isect_direction(*surf2,m1,uvuv_start2,du,dv,acuRatio);
	if(kdt==-1)
		return 0;

	double dt;
	if(kdt) dt=du; else dt=dv;
	if(dtmax>0. && dtmax<fabs(dt)){
		double save=dt; dt=dtmax; if(save<0.) dt=-dt;
	}

	point=MGBPointSeq(lmax,7);//(.,0-1) for this surface parameter (u,v),
							//(., 2-3) for sf2 surface parameter (u,v), and
							//(., 4-6) for (x,y,z) positional data. 

// Store starting point data.
//Recompute starting point by this method.
	uvuv_now=uvuv_start;
	surf1->isect_start_incr(*surf2,uvuv_start,kdt,0.,0.,m1,uvuv_now);

	point.store_at(0,uvuv_now,0,0,4);
	size_t n=1, nm1=0;
//Loop until new i.p. is on a boundary of surfaces or, i.p. reached to
// one of the points of uvuv_list.
	int obtained, try_num;
	bool kdt_changed=false;
	MGPosition uvuv_pre(uvuv_now);
	while(1){
		try_num=0;
	try_again:
		if(!(obtained=			//===First try.
			surf1->isect_start_incr(*surf2,uvuv_pre,kdt,du,dv,m1,uvuv_now))){
			du3=du*.3333; dv3=dv*.3333;
			obtained=		   //===Second try of same dt.
				surf1->isect_start_incr(*surf2,uvuv_pre,kdt,du3,dv3,m1,uvuv_now);
		}
		if(obtained){
			if(kdt_changed){
			double du4=uvuv_now[m1]-uvuv_pre[m1], dv4=uvuv_now[m1+1]-uvuv_pre[m1+1];
			if(kdt2){
				if(dv*dv4<=0.){ obtained=7;break;}
			}else{
				if(du*du4<=0.){  obtained=7;break;}
			}
			kdt_changed=false;
			}
		}else{
			kdt=!kdt; kdt_changed=true;
			if(kdt){if(fabs(du3)<=dtmin[kdt]) du3*=2.;}
			else {if(fabs(dv3)<=dtmin[kdt]) dv3*=2.;}
			//===Third try by changing kdt(u and v parameter line.)
			obtained=	   
				surf1->isect_start_incr(*surf2,uvuv_pre,kdt,du3,dv3,m1,uvuv_now);
			if(obtained){
				double du4=uvuv_now[m1]-uvuv_pre[m1], dv4=uvuv_now[m1+1]-uvuv_pre[m1+1];
				if(kdt){
					if(dv*dv4<=0.){ obtained=7;break;}
				}else{
					if(du*du4<=0.){  obtained=7;break;}
				}
				kdt_changed=false;
			}else{
				if(try_num==0 && n>=2){
					size_t m2=0; if(m1==0) m2=2;
					double du2old=uvuv_pre[m2]-point(n-2,m2);
					double dv2old=uvuv_pre[m2+1]-point(n-2,m2+1);
					kdt2=-1;
					MGPosition uv2(2,uvuv_pre,0,m2);
					double du2=du2old, dv2=dv2old;
					surf2->isect_inner_dt(n,uv2,du2,dv2,kdt2,acuRatio);
					if(kdt2!=-1){
						try_num+=1;
						//4th try by exchanging surf1 and surf2.
						const MGFSurface* stemp=surf1;surf1=surf2;surf2=stemp;
						if(m1==0) m1=2; else m1=0;
						kdt=kdt2; du=du2; dv=dv2;
						if(m1==0) dtmin=dtmin1;
						else      dtmin=dtmin2;
						kdt_changed=false;
						goto try_again;
					}
				}
				kdt=!kdt; kdt_changed=true;
					//===Final try of first kdt.
				obtained=
					surf1->isect_start_boundary(*surf2,uvuv_pre,kdt,du3,dv3,m1,uvuv_now);
				if(!obtained){
					obtained=7;break;
				}
			}
		}

		if(try_num==0){
			//Test if we passed too long surf2 distance.
			size_t m2=0; if(m1==0) m2=2;
			double du2old=uvuv_now[m2]-uvuv_pre[m2];
			double dv2old=uvuv_now[m2+1]-uvuv_pre[m2+1];
			kdt2=-1;
			MGPosition uv2(2,uvuv_now,0,m2);
			double du2=du2old, dv2=dv2old;
			surf2->isect_inner_dt(n,uv2,du2,dv2,kdt2,acuRatio);
			if(kdt2!=-1){
				if((kdt2 && fabs(du2)>2.1*fabs(du2old)) ||
				  (!kdt2 && fabs(dv2)>2.1*fabs(dv2old))){
					try_num+=1;
				//Exchange surf1 and surf2 since we passed too long surf2's ditance.
					const MGFSurface* stemp=surf1;surf1=surf2;surf2=stemp;
					if(m1==0) m1=2; else m1=0;
					kdt=kdt2; du=du2; dv=dv2;
					if(m1==0) dtmin=dtmin1;
					else          dtmin=dtmin2;
					uvuv_now=uvuv_pre;
					goto try_again;
				}
			}
		}

		for(size_t ii=0; ii<4; ii++){
			double dttemp=uvuv_pre[ii]-uvuv_now[ii];
			if(dttemp<0.) dttemp=-dttemp; tTolNow[ii]=dttemp*2.;
			if(tTolNow[ii]<tTol[ii]) tTolNow[ii]=tTol[ii];
		}
		MGBox uv_passed(uvuv_now, tTolNow);
			// uv_passed is the box of one span parameter length that passed.
			// When parameter of one of uvuv_list is included in this box, 
			// computation terminates.
			//cout<<" uv_passed="<<uv_passed<<endl;////
			//cout<<" tTolNow="<<tTolNow[0]<<","<<tTolNow[1]<<","<<tTolNow[2]<<","<<tTolNow[3]<<endl;////

		//Test if a point in uvuv_list or uvuv_start is in the last 1 span
		//parameter range(uv_passed). If so, end point of intersection point
		//is found and terminate the marching.
		if(uvuv_list.in(uv_passed,uvuv_id,4)){
			uvuv_now=MGPosition(4,*uvuv_id);
			obtained=3;
		}else if(n>3 && uv_passed>>uvuv_start){
			//cout<<" uv_passed="<<uv_passed<<endl;////
			//cout<<" uvuv_now="<<uvuv_now<<endl;////
			//cout<<" tTolNow="<<tTolNow<<endl;////
			//cout<<" uvuv_start="<<uvuv_start<<endl;////
			uvuv_now=uvuv_start;
			obtained=4;
		}

	//Store obtained data in line.
		size_t pointSize=point.capacity();
		if(n>=pointSize){
			size_t len=pointSize+lmax;
			point.reshape(len);
		}
		point.store_at(n,uvuv_now,0,0,4);
		nm1=n++; point.set_length(n);
		if(obtained!=1) break;

	//Define new surf1,surf2,dt,kdt.
		du=uvuv_now[m1]-uvuv_pre[m1]; dv=uvuv_now[m1+1]-uvuv_pre[m1+1];
		kdt2=kdt;
		MGPosition uv(2,uvuv_now,0,m1);
		uvuv_pre=uvuv_now;
//cout<<"n="<<n-1<<" kdt="<<kdt<<" du="<<du<<",dv="<<dv<<" uvuv="<<uvuv_now; ////
//cout<<" "<<eval(uvuv_now(0),uvuv_now(1))<<endl;//////////////
		surf1->isect_inner_dt(nm1,uv,du,dv,kdt2,acuRatio);
		if(kdt2==-1){ obtained=7; break;}//kdt will be used even if kdt2=-1.
		if(kdt!=kdt2){ kdt=kdt2; kdt_changed=true;}	
	}

	if(n<=1) return 0;
	if(obtained==7){
		uvuv_now=uvuv_pre;
		if(kdt_changed) kdt=!kdt;
	}
	// obtained=1: ip found as a point on a perimeter of one of the surfaces.
	//             This case occurs at on_the_perimeter() since isect_start_incr's 
	//              return value is always 1.
	//         =2: ip found by isect_start_boundary as a boundary point.
	//         =3: ip passed one of boundary points in uvuv_list.
	//         =4: ip returned to the starting point.
	//         =7: ip was lost during the computation(may be tangent).
	//cout<<"isect_startPt, endpoint="<<uvuv_now<<", obtained="<<obtained<<endl;

	MGPosition uvuv_end=uvuv_now;//Save the end point data.
//Normalize last points. The normalization is done by deleting last
//ndel points and adding (nadd+ndel) points in replace.
	size_t ndel,nadd=0;
	if(n<=isect_div_id_max()){ndel=0;}
	else if(n<isect_div_id_max()+3){ndel=n-isect_div_id_max();}
	else{ndel=3;}
	//nadd is number of points to add for the normalization.
	size_t naddby2=nadd+ndel;
	size_t nold=n; n+=nadd;
	if(n>point.capacity()) point.reshape(n);

	//Compute normalized span length of the last point (dt).
	int id_sect_div=nold-2-nadd; if(id_sect_div<0) id_sect_div=0;
	double divt,span=0., maxdivt=isect_dt_coef(id_sect_div);
	size_t j=naddby2;
	for(i=0; i<=naddby2; i++, j--){
		divt=isect_dt_coef(j); if(divt>maxdivt) divt=maxdivt;
		span=span+divt;
	}

	double dtSum=0.;
	size_t id=(kdt+1)%2+m1;
	for(i=nold-ndel-1; i<=nold-1; i++){
		double dt2=point(i,id)-point(i-1,id);
		dtSum=dtSum+fabs(dt2);
	}
	double dtsave=dtSum/span;		//dtsave is always plus.

	//Re-compute last (nadd+ndel) points in normalized spans.
	size_t inext=n-naddby2-2; j=naddby2;
	uvuv_pre(0)=point(inext,0); uvuv_pre(1)=point(inext,1);
	uvuv_pre(2)=point(inext,2); uvuv_pre(3)=point(inext,3);
	nm1=n-1;
	//cout<<"-----------------------------------------------"<<endl;
	if(kdt) du=point(inext+1,m1)-point(inext,m1);
	else dv=point(inext+1,m1+1)-point(inext,m1+1);
	du=point(inext+1,id)-point(inext,id);
	for(i=inext+1; i<nm1; i++, j--){
		divt=isect_dt_coef(j); if(divt>maxdivt) divt=maxdivt;
		double dt=dtsave*divt;//dtsave is always plus.
		if(kdt){if(du<0.) du=-dt;else du=dt;}
		else{if(dv<0.) dv=-dt;else dv=dt;}
		surf1->isect_start_incr(*surf2,uvuv_pre,kdt,du,dv,m1,uvuv_now);
		point.store_at(i,uvuv_now,0,0,4);
		//cout<<"n="<<i<<" kdt="<<kdt<<" du="<<du<<",dv="<<dv<<" uvuv="<<uvuv_now; ////
		//cout<<" "<<eval(uvuv_now(0),uvuv_now(1))<<endl;//////////////
		uvuv_pre=uvuv_now;
	}
	//Recompute End point by this method(uvuv_now includes the end point found).
	surf1->isect_start_incr(*surf2,uvuv_end,kdt,0.,0.,m1,uvuv_end);
	point.store_at(nm1,uvuv_end,0,0,4);
	point.set_length(n);
		//cout<<"n="<<i<<" kdt="<<kdt<<" uvuv="<<uvuv_end; ////
		//cout<<" "<<eval(uvuv_end(0),uvuv_end(1))<<", Obtained="<<obtained<<endl;//////////////
	return obtained;
}

//isect_startH compute one intersection line of two surfaces, this and sf2,
// given starting intersetion point uvuv((u1,v1) of this and (u2,v2) of sf2).
// isect_startH halts the computation when intersection
// reached to a boundary of this or sf2, or reached to one of the points
// in uvuv_list.
//The function's return value is:
// =0: Intersection was not obtained.
// !=0: Intersection was obtained as follows:
//    =1: End point is a point on a perimeter of one of the surfaces.
//    =3: End point is one of boundary points in uvuv_list.
//    =4: End point is the starting point.
//    =7: isect_startH halted the computation since intersection was lost
//     during the computation.
//This is a proprietary function for MGShell.
int	MGSurface::isect_startH(
	const MGPosition& uvuv_startIn, //Starting point of the intersection line.
	MGPosition_list& uvuv_list,	//isect_startH will halt when ip reached one of 
		//the point in uvuv_list. isect_startH does not change uvuv_list(actually
		//uvuv_list is const). uvuv's space dimension is at least 4,
		//and the first 2 is (u,v) of this and the next 2 is (u,v) of sf2. 
	const MGSurface& sf2,	//2nd surface.
	MGSSisect& ssi,			//Surface-surface intersection line will be output.
	MGPosition_list::iterator& uvuv_id
		//When the end point of ip was one of the points of uvuv_list, that is, 
		//when the function's return value was 3, uvuv_list's iterator
		//of the point will be returned,
		//When the end point was not a point of uvuv_list, end() of uvuv_list
		//will be returned.
) const{
	int obtained;

	const MGPlane* pl1=dynamic_cast<const MGPlane*>(this);
	if(pl1){//Case that this is a plane.
		obtained=pl1->isect_startHPL(uvuv_startIn, uvuv_list, sf2, ssi, uvuv_id);
	}else{//Case that this is not a plane.
		const MGPlane* pl2=dynamic_cast<const MGPlane*>(&sf2);
		if(pl2){//Case that sf2 is a plane.
		//Compute intersection lines using isect_startPlane.
			obtained=isect_startPlane(uvuv_startIn,uvuv_list,*pl2,ssi,uvuv_id);
		}else{//Case that both are not a plane.
		//Compute intersection line using isect_start.
			int m1;
			obtained=isect_start(uvuv_startIn,uvuv_list,sf2,ssi,uvuv_id,m1);
		}
	}
	return obtained;
}

MGSSisect_list MGSurface::isect(const MGFSurface& fsurf) const{
	MGSSisect_list list=fsurf.isect(*this);
	list.replace12();
	return list;
}

//Intersection of this plane and a Face.
MGSSisect_list MGSurface::isect(const MGFace& f)const{
	MGSSisect_list list=f.isect(*this);
	list.replace12();
	return list;
}
