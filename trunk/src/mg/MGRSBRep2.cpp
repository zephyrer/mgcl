/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position_list.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"
#include "mg/Matrix.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/SPhere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/Bluprt.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;
//***************************************************************/
//*                      All rights reserved.                   */
//***************************************************************/

// Implements MGRSBRep Class
//
// Defines Rational Surface B-Representation.
// This NURBS is of homogeneous form, i.e., B-Coefficients have
// weight included values. 
// When usual NURBS form is (xij, yij, zij, wij) ,
// MGRSBRep form is (xij*wij, yij*wij, zij*wij, wij)
//				 for i=0,..., m-1, and j=0,..., n-1.

//******** Intersection implementation *******

// Surface と Curve の交点を求める。
MGCSisect_list MGRSBRep::isect(const MGCurve& curve)const{
	return curve.isect(*this);
}
// Surface と Curve の交点を求める。
MGCSisect_list MGRSBRep::isect(const MGRLBRep& curve)const{
	return curve.isect(*this);
}
// Surface と Curve の交点を求める。
MGCSisect_list MGRSBRep::isect(const MGEllipse& curve)const{
	return curve.isect(*this);
}
// Surface と Curve の交点を求める。
MGCSisect_list MGRSBRep::isect(const MGLBRep& curve)const{
	return curve.isect(*this);
}
// Surface と Curve の交点を求める。
MGCSisect_list MGRSBRep::isect(const MGSurfCurve& curve)const{
	return curve.isect(*this);
}
// Surface と Curve の交点を求める。
MGCSisect_list MGRSBRep::isect(const MGBSumCurve& curve)const{
	return curve.isect(*this);
}

// Surface と Surface の交線を求める。
//Compute intesection of Sphere and Surface.
MGSSisect_list MGRSBRep::isect(const MGSurface& srf2) const{
	MGSSisect_list list=srf2.isect(*this);
	list.replace12();
	return list;
}
MGSSisect_list MGRSBRep::isect(const MGPlane& srf2) const{
	return intersectPl(srf2);
}
MGSSisect_list MGRSBRep::isect(const MGSphere& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGRSBRep::isect(const MGCylinder& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGRSBRep::isect(const MGSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGRSBRep::isect(const MGRSBRep& srf2) const{
	return intersect(srf2);
}
MGSSisect_list MGRSBRep::isect(const MGBSumSurf& srf2) const{
	return intersect(srf2);
}

//The following two function will be used in perps or isect
//to decide how many division of the surface along u or v direction
//should be applied before using perp_guess or isect_guess.
size_t MGRSBRep::intersect_dnum_u() const{
	return (bdim_u()+2-order_u())*order_u();
}
size_t MGRSBRep::intersect_dnum_v() const{
	return (bdim_v()+2-order_v())*order_v();
}

//"isect_incr_pline" is a dedicated function of isect_start_incr, will get
// shortest parameter line necessary to compute intersection.
MGCurve* MGRSBRep::isect_incr_pline(
	const MGPosition& uv,	//last intersection point.
	int kdt,				//Input if u=const v-parameter line or not.
							// true:u=const, false:v=const.
	double du, double dv,//Incremental parameter length.
	double& u,				//next u value will be output
	double& v,				//next v value will be output
	size_t incr			//Invremental id value of B-coef's.
)const{
	MGRLBRep* rlb=new MGRLBRep;
	m_surface.isect_incr_pline2(uv,kdt,du,dv,u,v,incr,rlb->homogeneous());
	return rlb;
}

//Return order of intersection line order of MGLBRep.
size_t MGRSBRep::isect_order() const{
	size_t ku=order_u(), kv=order_v(); 
	if(ku<kv) ku=kv;
	return ku+1;
}

//Intersection of Surface and a straight line.
MGCSisect_list MGRSBRep::isectSl(
	const MGStraight& sl,
	const MGBox& uvbox //indicates if this surface is restrictied to the parameter
					//range of uvbox. If uvbox.is_null(), no restriction.
)const{
	MGCSisect_list list(&sl,this);
	const MGBox& sbx=box();
	if(!sbx.crossing(sl))
		return list;

	MGUnit_vector SLD=sl.direction();
	MGMatrix mat; mat.to_axis(SLD,2);	//Matrix to transform SLD to be z axis.

	double u0,u1,v0,v1;
	bool uvbox_is_null=uvbox.is_null();
	if(uvbox_is_null){
		u0=param_s_u(); u1=param_e_u();
		v0=param_s_v(); v1=param_e_v();
	}else{
		u0=uvbox[0].low_point(); u1=uvbox[0].high_point();
		v0=uvbox[1].low_point(); v1=uvbox[1].high_point();
		double uspan=u1-u0, vspan=v1-v0;
		if(MGREqual_base(u0,param_s_u(),uspan) && 
			MGREqual_base(u1,param_e_u(),uspan) && 
			MGREqual_base(v0,param_s_v(),vspan) &&
			MGREqual_base(v1,param_e_v(),vspan) ) uvbox_is_null=true;
	}
	double u=(u0+u1)*.5, v=(v0+v1)*.5;
	MGVector PN=sl.eval(sl.closest(eval(u,v)));
		//PN is the nearest point of the sl from the center of this surface,
		//and is the vector to translate sl to pass through the origin.
	MGRSBRep* sf2D=static_cast<MGRSBRep*>(copy_surface());
	(*sf2D)-=PN; (*sf2D)*=mat;
	//std::cout<<(*sf2D);////////////
	sf2D->change_dimension(2);
	//sf2D is the 2D surface that is transformed from the original this
	//surface so that the straight line sl is seen as a point of origin(0.,.0.).
	if(!uvbox_is_null) sf2D->limit(uvbox);
	//std::cout<<(*sf2D);////////////

	size_t m=sf2D->bdim_u(), n=sf2D->bdim_v();
	size_t mm1=m-1, nm1=n-1;
	MGSPointSeq sp=sf2D->surface_bcoef().non_homogeneous();
	//std::cout<<(sp);////////////
	MGKnotVector& tu=sf2D->knot_vector_u();
	MGKnotVector& tv=sf2D->knot_vector_v();
	size_t kum1=tu.order()-1, kvm1=tv.order()-1;

	for(size_t i=kum1; i<m; i++){
	size_t ip1=i+1;
	while(tu(i)== tu(ip1) && ip1<m) {i=ip1; ip1++;}
	for(size_t j=kvm1; j<n; j++){
		size_t jp1=j+1;
		while(tv(j)== tv(jp1) && jp1<n) {j=jp1; jp1++;}
		MGBox bx;
		for(size_t iu=i-kum1; iu<=i; iu++)
			for(size_t ju=j-kvm1; ju<=j; ju++) bx.expand(sp(iu,ju));
		if(bx.includes_origin()){
			//Compute intersection point using isect_guess_straight.
			MGPosition uv((tu[i]+tu[ip1])*.5, (tv[j]+tv[jp1])*.5);
			double t=sl.closest(eval(uv));
			if(isect_guess_straight(sl,t,uv,t,uv)){
				if(uvbox_is_null || uvbox>>uv) list.append(sl.eval(t),t,uv);
			}
		}
	}
	}

	delete sf2D;
	return list;
}

//Test if the RSBRep is planar or not.
//Returned is 0(false) if this is not planar, 1(true) if this is planar.
int MGRSBRep::planar(
		MGPlane& plane,		//Plane that might be closest to this.
							//plane is always output even if not planar.
		double& deviation	//maximum deviation of this from the output plane.
		)const{
	double u0=param_s_u(), u1=param_e_u();
	double v0=param_s_v(), v1=param_e_v();
	MGUnit_vector N=(normal(u0,v0)+normal(u0,v1)+normal(u1,v0)+normal(u1,v1));
	MGPosition P=(eval(u0,v0)+eval(u0,v1)+eval(u1,v0)+eval(u1,v1)
					+eval((u0+u1)*.5,(v0+v1)*.5))/5.;
	plane=MGPlane(N,P);

	size_t ncdm1=sdim();
	size_t i,j,k;
//	double mt[4]; for(k=0; k<=3; k++) mt[k]=mat(k,3);
	size_t m=bdim_u(), n=bdim_v();
	deviation=0.;
	double d=P%N;//d is surface expression's d(ax+by+cz=d), of plane(P,N).
	for(i=0; i<m; i++){
	for(j=0; j<n; j++){
		double w=coef(i,j,ncdm1), x=d;
		for(k=0; k<ncdm1; k++) x-=(coef(i,j,k)/w)*N[k];
//		double x=d-coef(i,j,0)*mt[0]+coef(i,j,1)*mt[1]+coef(i,j,2)*mt[2]
//					+coef(i,j,3)*mt[3];
		if(x<0.) x=-x;
		if(deviation<x) deviation=x;
	}
	}
	return MGAZero(deviation);
}

//Test if part of the surface is planar or not within the tolerance tol.
//The part of the surface is input by the surface parameter range uvbox.
//Returned is 0(false) if this is not planar, 1(true) if planar.
int MGRSBRep::planar(
	const MGBox& uvbox,//This surface parameter range.
	double tol,	//maximum deviation allowed to regard the sub surface as a plane.
	int* divideU//Direction to subdivide will be output, if this was not planar.
				//=1: u direction, =0: v direction.
	) const{
	tol*=1.1;
	int i,j,k;

	MGBox uvb=param_range()&uvbox;
	MGUnit_vector N;
	MGPosition P;
	int direction;
	if(!flat(uvb,tol,direction,P,N)){
		if(divideU) *divideU=direction;
		return 0;
	}

	int ncdm1=sdim();
	int ncd=ncdm1+1;

	const MGInterval& urng=uvbox[0];
	double u0=urng[0].value(), u1=urng[1].value();
	const MGInterval& vrng=uvbox[1];
	double v0=vrng[0].value(), v1=vrng[1].value();
	unsigned ku=order_u(), kv=order_v();
	const MGKnotVector& tv=knot_vector_v();
	int j0=tv.locate(v0)-kv+1, j1=tv.locate(v1);
		//Coef necessary for v direction part is from j0 to j1
		//(B-rep dimension is j1-j0+1).
		//Knot vector necessary is from j0 to j1+kv.
	int nunew, nvnew=j1-j0+1;

	size_t nu, nv, sizeu, sizev;
	surface_bcoef().length(nu, nv);	surface_bcoef().capacity(sizeu, sizev);
	unsigned kmax=ku; if(kmax<kv) kmax=kv;
	double* work=new double[kmax*kmax];

	MGSPointSeq surf1(nvnew, nu, ncd); size_t irc=nvnew*nu;
	int nmax=nu; if(nmax<nvnew) nmax=nvnew;
	int nvncd=nvnew; if(nvncd<ncd) nvncd=ncd;
	MGBPointSeq temp(nmax,nvncd);
	MGKnotVector t(kmax, nmax);
	double* tpointer=t.data(); double* temppointer=temp.data();

	for(k=0; k<ncd; k++){
		bluprt_(ku, nu, knot_data_u(), coef_data(0,j0,k),
			sizeu,nvnew,u0,u1,nmax,work,&nunew,tpointer,temppointer,1);
		for(i=0; i<nunew; i++) for(j=0; j<nvnew; j++) surf1(j,i,k)=temp(i,j);
	}
	//surf1.set_length(nvnew,nunew);cout<<surf1<<endl;//////////
	const double d=P%N;//d is surface expression's d(ax+by+cz=d), of plane(P,N).
	const double* tvnew=knot_data_v()+j0;
	int nvnew2;
	double w,x;
	for(i=0; i<nunew; i++){
		bluprt_(kv,nvnew,tvnew,surf1.data(0,i,0),
			irc,ncd,v0,v1,nmax,work,&nvnew2,tpointer,temppointer,1);
//		for(size_t i1=0;i1<nvnew2;i1++){
//			for(size_t j1=0; j1<ncd; j1++) cout<<temp(i1,j1)<<",";
//			cout<<endl;
//		}
		for(j=0; j<nvnew2; j++){
			w=temp(j,ncdm1), x=d;
			for(k=0; k<ncdm1; k++) x-=(temp(j,k)/w)*N[k];
//			for(k=0; k<ncdm1; k++) cout<<(temp(j,k)/w)*N[k]<<",";
//			cout<<endl;
			if(x<0.) x=-x;
			if(x>tol) break;
		}
		if(j<nvnew2) break;
	}

	int retcode=1;
	if(i<nunew || j<nvnew2){
		if(divideU) *divideU=direction;
		retcode=0;
	}

	delete[] work;
	return retcode;
}

//Obtain 1D surface rep. of this surf which can be used for
//isect(const MGPlane& pl). This surf1D is used in isect for
//the argument of isect_startPlane, which will use surf1D to compute isect(pl).
//surf1D=0.(intersection with x=0. plane) is the intersection lines.
MGSBRep* MGRSBRep::surf1D(
	const MGPlane& pl
)const{
	MGVector C(3,pl.root_point());
	MGVector PN(3,pl.normal()); //Plane Normal.
	MGVector N(4,PN);
	N(3)=-(C%PN);
		//N is 4D vector that is normal to the 4D vector(C,1), and,
		//pl.u_deriv() and pl.v_deriv().
	MGMatrix mat; mat.to_axis(N,3);
		//Matrix to transform N to be w(weight) axis.

	size_t i,j;
	double len=N.len();
	double mt[4]; for(i=0; i<4; i++) mt[i]=mat(i,3)*len;
		//Multiplication of len is to adjust the correctness(tolerance) of
		//the intersection computation.

	//Construct a working 1D surface sf1D.
	MGSBRep* surf=new MGSBRep;
	size_t m=bdim_u(), n=bdim_v();
	MGSPointSeq& cp=surf->m_surface_bcoef;
	cp.resize(m,n,1);
	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			cp(i,j,0)=coef(i,j,0)*mt[0]+coef(i,j,1)*mt[1]
						+coef(i,j,2)*mt[2]+coef(i,j,3)*mt[3];
//		cout<<(surface_bcoef()*mat)<<cp<<endl;////
	surf->m_uknot=knot_vector_u();
	surf->m_vknot=knot_vector_v();
	return surf;
}
