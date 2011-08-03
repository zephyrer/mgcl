/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Unit_vector.h"
#include "mg/Matrix.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/PPRep.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/BSumCurve.h"
#include "mg/SurfCurve.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implements MGRLBRep Class
//
// Defines Rational Line B-Representation.
// This NURBS is of homogeneous form, i.e., B-Coefficients have
// weight included values. 
// When usual NURBS form is (xi, yi, zi, wi) ,
// MGRLBRep form is (xi*wi, yi*wi, zi*wi, wi) for i=0,..., n-1.

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
MGRLBRep& MGRLBRep::coordinate_exchange(size_t i, size_t j){
	m_line.line_bcoef().coordinate_exchange(i,j);
	update_mark();
	return *this;
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGRLBRep* MGRLBRep::clone() const{return new MGRLBRep(*this);}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGRLBRep* MGRLBRep::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2 		// Source order of this line.
)const{
	return new MGRLBRep(sdim,*this,start1,start2);
}

//Construct new curve object by copying to newed area,
//and limitting the parameter range to prange.
//Returned is newed object and must be deleted.
MGCurve* MGRLBRep::copy_limitted(const MGInterval& prange) const{
	MGInterval pr=param_range()&prange;
	return new MGRLBRep(pr.low_point(),pr.high_point(),*this);
}

//Divide this curve at the designated knot multiplicity point.
//Function's return value is the number of the curves after divided.
int MGRLBRep::divide_multi(
	MGPvector<MGCurve>& crv_list,	//divided curves will be appended.
	int multiplicity	//designates the multiplicity of the knot to divide at.
						//When multiplicity<=0, order()-1 is assumed.
						//When multiplicity>=order(), order() is assumed.
)const{
	MGPvector<MGCurve> crv_list2;
	int n=m_line.divide_multi(crv_list2,multiplicity);
	for(int i=0; i<n; i++){
		MGLBRep& lbi=*(static_cast<MGLBRep*>(crv_list2[i]));
		crv_list.push_back(new MGRLBRep(lbi,1));
	}
	return n;
}

//Modify the original NURBS by extrapolating the specified perimeter.
//The extrapolation is C2 continuous if the order >=4.
MGRLBRep& MGRLBRep::extend(
	int start,			//Flag of start or end poit of the line.
						//If start is true extend on the start point.
	double length,		//chord length to extend. 
	double dk         //Coefficient of how curvature should vary at
//    extrapolation start point. When dk=0, curvature keeps same, i.e.
//    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
//    i.e. dK/dS=-K/length at extrapolation start point.
//    (S=parameter of arc length, K=Curvature at start point)
//    That is, when dk reaches to 1 from 0, curve changes to flat.
){
	double t,tau,vlen;
	if(start){
		t=param_s();
		vlen=eval(t,1).len();
		if(MGMZero(vlen)) return *this;
		tau=t-length/vlen;
	}else{
		t=param_e();
		vlen=eval(t,1).len();
		if(MGMZero(vlen)) return *this;
		tau=t+length/vlen;
	}
	extend_with_parameter(tau,dk);
	return *this;
}

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGRLBRep::extend(
	double length,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	int se=start ? 1:0;
	extend(se,length,0.);
}

//Extrapolate the curve by the parameter value.
MGRLBRep& MGRLBRep::extend_with_parameter(
	double tau,		//The parameter value at the end of extended point.
					//When tau<param_s(), extension will be done at the starting point.
					//When tau>param_e(), extension will be done at the end point.
	double dk     //Coefficient of how curvature should vary at the connecting point.
					//See extend();
){
	MGPPRep pp;
	extrapolated_pp(tau,dk,pp);

	MGLBRep lbext(pp);//std::cout<<lbext;
	if(tau<param_s()) connect(0,0,lbext);
	else connect(0,2,lbext);
	return *this;
}

//Extracts control points.
//Fucntion's return value is 
//true if control points was obtained, false if not.
bool MGRLBRep::get_control_points(
	MGBPointSeq& cpoints	//Control points will be output.
)const{
	cpoints=non_homogeneous_bcoef();
	return true;
}

//Provide divide number of curve span for function intersect.
size_t MGRLBRep::intersect_dnum() const{
	size_t k=order();
	int nspan=bdim()+1-k;
	int km1=k-1; if(km1<=0) km1=1;
	return nspan*km1;
}

//Test if this cure is co-planar with the 2nd curve curve2.
//MGPlane expression will be out to plane if this is co-planar.
//Function's return value is true if co-planar.
bool MGRLBRep::is_coplanar(const MGCurve& curve2, MGPlane& plane)const{
	const MGStraight* sl11=dynamic_cast<const MGStraight*>(&curve2);
	if(sl11) return sl11->is_coplanar(*this,plane);

	const MGLBRep* lb11=dynamic_cast<const MGLBRep*>(&curve2);
	if(lb11) return lb11->is_coplanar(*this,plane);

	MGPosition point;
	MGStraight sl1;
	int plkind=planar(plane,sl1,point);
	if(plkind==0) return false;

	MGPosition point2;
	MGStraight sl2;
	int plkind2=-1;
	MGPlane plane2;

	MGPosition uv;
	const MGRLBRep* rlb=dynamic_cast<const MGRLBRep*>(&curve2);
	if(rlb){
		plkind2=rlb->planar(plane2,sl2,point2);
		if(plkind2==0) return false;

		if(plkind2==1){
			if(plkind==1){
				plane=MGPlane(point,point2,point);
				return true;
			}else if(plkind==2){
				plane=MGPlane(sl1,point2);
				return true;
			}else{
				return plane.on(point2,uv);
			}
		}else if(plkind2==2){
			if(plkind==1){
				plane=MGPlane(sl2,point);
				return true;
			}else if(plkind==2){
				return sl1.is_coplanar(sl2,plane);
			}else{
				return plane.on(sl2);
			}
		}else{
			if(plkind==1){
				plane=plane2;
				return plane.on(point,uv);
			}else if(plkind==2){
				plane=plane2;
				return plane.on(sl1);
			}else{
				return plane==plane2;
			}
		}
	}else{
		if(!curve2.is_planar(plane2)) return false;
		if(plkind==1){
			plane=plane2;
			return plane.on(point,uv);
		}else if(plkind==2){
			plane=plane2;
			return plane.on(sl1);
		}
		return plane2==plane;
	}
}

//Test if this cure is planar or not.
//MGPlane expression will be out to plane if this is planar.
//Function's return value is true if planar.
bool MGRLBRep::is_planar(MGPlane& plane)const{
	MGStraight line;
	MGPosition point;
	int isp=planar(plane,line,point);
	if(isp==1){//IF this is within a point.
		plane=MGPlane(mgZ_UVEC,point);
	}else if(isp==2){//IF this is within a straight.
		return line.is_planar(plane);
	}

	return isp>0;
}

//Spline と Curve の交点を求める。
//Intersection point of NURBS and curve.
MGCCisect_list MGRLBRep::isect(const MGCurve& curve2) const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection point of NURBS and straight line.
MGCCisect_list MGRLBRep::isect(const MGStraight& sline) const{
	MGCCisect_list list(this,&sline);
	if(!has_common(sline)) return list;

	MGCParam_list clist;
	double errsave=MGTolerance::wc_zero();
	MGTolerance::set_wc_zero(errsave*.3);//make error smaller.
	if(sdim()<=2 && sline.sdim()<=2) clist=isect_2D(sline);
	else {
		MGPlane plane(sline, mgORIGIN);
			//The plane that includes sline and passes through the origin.
		clist=isect_3D(plane);
	}
	MGTolerance::set_wc_zero(errsave);	//Retrieve error.

	MGPosition p; double t1,t2;
	MGCParam_list::Citerator i;
	for(i=clist.begin(); i!=clist.end(); i++){
		t1=(*i); p=eval(t1);
		if(sline.on(p,t2)) list.append(MGCCisect(p,t1,t2));
	}

	return list;
}

//Spline と Curve の交点を求める。
//Intersection point of NURBS and curve.
MGCCisect_list MGRLBRep::isect(const MGSurfCurve& curve2) const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Spline と Curve の交点を求める。
//Intersection point of NURBS and curve.
MGCCisect_list MGRLBRep::isect(const MGBSumCurve& curve2) const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection of Spline and Surface
MGCSisect_list MGRLBRep::isect(const MGSurface& surf)const{
	return surf.isect(*this);
}

MGCSisect_list MGRLBRep::isect(const MGPlane& plane) const{
	MGCSisect_list list(this,&plane);

	MGCParam_list clist=isect_3D(plane);
	double tt; MGPosition p, uvpl;
	size_t inum=clist.entries();
	for(size_t i=0; i<inum; i++){
		tt=clist.removeFirst(); p=eval(tt);
		uvpl=plane.uv(p);
		const MGCSisect isect(plane.eval(uvpl),tt,uvpl);
		list.append(isect);
	}

	return list;
}

//Intersection of Spline and Surface
MGCSisect_list MGRLBRep::isect(const MGSphere& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGRLBRep::isect(const MGCylinder& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGRLBRep::isect(const MGRSBRep& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGRLBRep::isect(const MGSBRep& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGRLBRep::isect(const MGBSumSurf& surf)const{
	return surf.intersect(*this);
}

//Compute intersection point of 1D sub NURBS of original B-rep.
//Parameter values of intersection point will be returned.
MGCParam_list MGRLBRep::intersect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
)const{
	assert(coordinate<sdim());

	MGVector N(-1., f);		//Normal of Vector(f,1.).
	return isect_nD(N,1,coordinate);
}

//Compute intersection points of 2D sub NURBS of original B-rep.
//Parameter values of this at intersection points will be returned.
//Straight line sl and this(RLBRep) will be projected to 2D plane of
//coordinate kind (coordinate, coordinate+1), then intersection will
//be computed.
MGCParam_list MGRLBRep::isect_2D(						
	const MGStraight& sl,// Straight line.
	size_t coordinate	// Coordinate kind of 2D sub space.
)const{
	assert(coordinate<sdim());

	MGVector C(3,sl.nearest_to_origin(),0,coordinate); C(2)=1.;
	MGVector D(3,sl.direction(0.),0,coordinate); D(2)=0.;
	MGVector N=C*D;		//Normal of Vector C and D.
	//I.e., normal of the plane that includes sl and passes through origin.
	return isect_nD(N,2,coordinate);
}

//Compute intersection points of 3D sub NURBS of original B-rep.
//Parameter values of this at intersection points will be returned.
//This(RLBRep) will be projected to 3D plane of coordinate kind 
//(coordinate, coordinate+1, coordinate+2), then intersection will
//be computed.
MGCParam_list MGRLBRep::isect_3D(
	const MGPlane& pl,	// Plane.
	size_t coordinate	// Coordinate kind of 3D sub space.
)const{
	assert(coordinate<sdim());

	MGVector C(3,pl.root_point(),0,coordinate);
	MGVector PN(3,pl.normal(),0,coordinate); //Plane Normal.
	MGVector N(4,PN);
	N(3)=-(C%PN);
		//N is 4D vector that is normal to C, pl.u_deriv() and pl.v_deriv().
	return isect_nD(N,3,coordinate);
}

//Compute intersection points of (n+1)-Dimensional sub NURBS
// and (n+1)-dimension plane that passes through the origin.
//(n+1)-Dimensional sub NURBS comprises from n dimension coordinate
// and weight element data.
//Parameter values of this at intersection points will be returned.
//MGVector N is the normal vector of the plane.
MGCParam_list MGRLBRep::isect_nD(						
	const MGVector& N,		// Normal of (n+1)-dimension plane.
	size_t dimension,		// Number of dimension n.
	size_t coordinate		// Coordinate kind of n-dimensional sub NURBS.
)const{
	assert(dimension<4);	//Currently dimension must be up to 3.
	double m[4]; size_t cod[4];

	size_t i,j;
	MGMatrix mat; mat.to_axis(N,dimension);
	double len=N.len();
	for(j=0; j<=dimension; j++) m[j]=mat(j,dimension)*len;
		//Multiplication of len is to adjust the correctness(tolerance) of
		//the intersection computation.
	size_t nbd=bdim(), w_id=sdim();
	size_t ncod=dimension; if(ncod>w_id) ncod=w_id;
	cod[0]=coordinate;
	for(j=1; j<ncod; j++){
		cod[j]=cod[j-1]+1; if(cod[j]>=w_id) cod[j]=0; 
	}
	for(j=ncod; j<dimension; j++) cod[j]=w_id+1;
			//To enforce retrieval of zero from m_line.line_bcoef().

	MGLBRep brep;
	MGBPointSeq& brepcf=brep.line_bcoef(); brepcf.resize(nbd,1);
	brep.knot_vector()=knot_vector();

	double min,max;
	cod[dimension]=w_id;
	const MGBPointSeq& cf=m_line.line_bcoef();
	double x=cf(0,cod[0])*m[0];
	for(j=1; j<=dimension; j++) x+=cf(0,cod[j])*m[j];
	min=max=brepcf(0,0)=x;
	for(i=1; i<nbd; i++){
		x=cf(i,cod[0])*m[0];
		for(j=1; j<=dimension; j++) x +=cf(i,cod[j])*m[j];
		brepcf(i,0)=x;
		if(x<min) min=x; if(x>max) max=x;
	}
	double error=MGTolerance::wc_zero();
	if(max<-error || min>error) return MGCParam_list();
//		cout<<brep<<endl;
	return brep.isect_1D(0.);
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> MGRLBRep::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
)const{
	size_t i, nbd=bdim();
	size_t j, w_id=sdim();
	size_t k=order();
	size_t ncod=3; if(ncod>w_id) ncod=w_id;

	const MGBPointSeq& cf=m_line.line_bcoef();
	MGLBRep* brep=new MGLBRep();
	MGBPointSeq& rcoef=brep->line_bcoef(); rcoef.resize(nbd,1);
	MGKnotVector& knotv=brep->knot_vector(); knotv.size_change(k,nbd);

	double min,max;
	const MGKnotVector& t=knot_vector();
	double rt=0.; for(j=0; j<ncod; j++) rt+=coef(0,j)*g[j];
	min=max=rt-=coef(0,w_id)*g[3];
	rcoef(0,0)=rt;
	knotv(0)=t(0);
	for(i=1; i<nbd; i++){
		rt=0.; for(j=0; j<ncod; j++) rt+=coef(i,j)*g[j];
		rcoef(i,0)=rt-=coef(i,w_id)*g[3];
		knotv(i)=t(i);
		if(min>rt) min=rt; if(max<rt) max=rt;
	}
	size_t npk=nbd+k;
	for(i=nbd; i<npk; i++) knotv(i)=t(i);
	MGInterval minmax(min,max);
	m_box=new MGBox(1,&minmax);
	return std::auto_ptr<MGCurve>(brep);
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv2's parameter
//as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list MGRLBRep::perps(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}
MGPosition_list MGRLBRep::perps(
	const MGStraight& crv2	//The second curve
)const{return perpsSl(crv2);}
