/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position_list.h"
#include "mg/Transf.h"
#include "mg/Unit_vector.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Ellipse.h"
#include "mg/RLBRep.h"
#include "mg/CParam_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
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

// MGEllipse.cc
//
// MGEllipse is a class to define an ellipse of 2D or 3D.
// Ellipse is expressed as below using parameter t:
// Point(t) = m_center + m_m * cos(t) + m_n * sin(t)

// メンバ関数

//Extrapolate this curve by an (approximate) chord length.
//The extrapolation is C2 continuous.
void MGEllipse::extend(
	double len,	//approximate chord length to extend. 
	bool start		//Flag of which point to extend, start or end point of the line.
					//If start is true extend on the start point.
){
	double t, ts=param_s(), te=param_e(), extlen;
	double ts2, te2;
	double whole_length1=length(ts, te);
	MGEllipse el2(*this); el2.unlimit();
	double whole_length2=el2.length(el2.param_s(), el2.param_e());
	double remain=whole_length2-whole_length1;
	if(remain<len)
		len=remain;

	if(start){
		unlimit_start();
		extlen=-len;
		t=ts;
		te2=te;
	}else{
		unlimit_end();
		ts2=ts;
		t=te;
		extlen=len;
	}
	double t_target=length_param(t,extlen);

	if(start)
		ts2=t_target;
	else
		te2=t_target;
	set_param(gp_to_radian(ts2), gp_to_radian(te2));
}

//Provide divide number of curve span for function intersect.
size_t MGEllipse::intersect_dnum()const{
	return size_t(1.+(m_prange[1]-m_prange[0])*4./mgPAI);
}

//Test if this cure is linear or not, that is, is straight or not.
//MGStraight expression will be out to straight if this is linear or not.
//Function's return value is true if linear.
bool MGEllipse::is_linear(MGStraight& straight)const{
	double ts=(param_s()+param_e())*.5;
	straight=MGStraight(MGUnit_vector(eval(ts,1)),MGPosition(eval(ts)));
	return false;
}

//Test if this cure is planar or not.
//MGPlane expression will be out to plane if this is planar.
//Function's return value is true if planar.
bool MGEllipse::is_planar(MGPlane& plane)const{
	plane=MGPlane(m_normal,m_center);
	return true;
}

//Test if this ellipse is whole ellipse or not.
//Function's return value is true if this is whole.
bool MGEllipse::is_whole_ellipse()const{
	double t0=gp_to_radian(param_s()), t1=gp_to_radian(param_e());
	return ((t1-t0)+MGTolerance::angle_zero())>=mgDBLPAI;
}

// Ellipse と Curve の交点を求める。
MGCCisect_list MGEllipse::isect(
	const MGCurve& curve2
)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

// Ellipse と Straight の交点を求める。
MGCCisect_list MGEllipse::isect(
	const MGStraight& st
)const{

	MGCCisect_list list(this, &st); // リストの領域を得る
	if(!has_common(st)) return list;

	if((st.direction()).orthogonal(m_normal)){
		// When st and ellipse are parallel

		MGVector vt(m_center, st.root_point());
		if(vt.is_zero_vector() | m_normal.orthogonal(vt)){
		//When ellipse and straight line lie on a same plane.
			double t[2],angle[2]; int tangen;
			MGTransf tr(m_m, m_n, m_center); // tr is transform to transform
			// ellipes to lie on xy plane and the center to be origin.
			MGPosition stpt=st.root_point()*tr;
			MGUnit_vector dir=st.direction()*tr.affine();
			double dirlen=st.direction_len();
			int n=isect2d(stpt, dir,t,angle,tangen);
			if(n){
				double p;
				for(int i=0; i<n; i++){
					if(st.in_range(p=t[i]/dirlen)){
						MGCCRELATION rel=MGCCREL_ISECT;
						if(tangen) rel=MGCCREL_TANGENT;
						MGCCisect ip(st.eval_position(p),
							radian_to_gp(angle[i]),p,rel);
						list.append(ip);
					}
				}
			}
		}
	}else{
		if(MGAZero(st.distance(m_center))) return list;
		//Compute intersection point with the plane that includes st.
		MGPlane pl1(m_normal, m_center);	//Plane that includes ellipse.
		MGPlane pl2(st, m_center);			//Plane that includes st.
		MGStraight st2;
		pl1.relation(pl2, st2);
		MGCCisect_list list2=isect(st2);
		size_t n=list2.entries(); double t;
		for(size_t i=0; i<n; i++){
			MGCCisect is=list2.removeFirst();
			MGPosition p(is.point());
			if(st.on(p,t)) list.append(MGCCisect(p,is.param1(),t,is.rel()));
		}
	}

	return list;
}

// Ellipse と Ellipse の交点を求める。
MGCCisect_list MGEllipse::isect(const MGEllipse &e2) const{
	MGCCisect_list list(this, &e2); // リストの領域を得る
	double r1=m_r, r2=e2.m_r;
	if(MGAZero(r1) || MGAZero(r2)) return list;
	if(!circle() || !e2.circle()) return intersect(e2);
	 	//When general ellipses, obtain by general intersection function.

	MGVector c1toc2=e2.m_center-m_center;
	if(c1toc2.is_zero_vector()) return list;
	if(m_normal.parallel(e2.m_normal) && c1toc2.orthogonal(m_normal)){
	//When two ellipses are on a same plane and both of them are circles.
		double d=c1toc2.len();		//d is the distance of two centers.
		if((r1+r2)>=d && (d+r1)>=r2 && (d+r2)>=r1){
		//When intersections of whole circles exist.
			const MGEllipse *e1p, *e2p; double cosval;
			if(r1>=r2){
				e1p=this; e2p=&e2;
				cosval=(d/r1-(r2*r2)/(r1*d)+r1/d)*0.5;
			}
			else{
				e2p=this; e1p=&e2;
				cosval=(d/r2-(r1*r1)/(r2*d)+r2/d)*0.5;
				c1toc2=-c1toc2;
			}
			double angle=acos(cosval);
			MGMatrix mat1;mat1.set_rotate_3D(m_normal,angle);
			MGMatrix mat2;mat2.set_rotate_3D(m_normal,-angle);
			MGVector vec1=c1toc2*mat1, vec2=c1toc2*mat2;
			MGCCisect_list list1=(*e1p).isect
				(MGStraight(MGSTRAIGHT_HALF_LIMIT,vec1,(*e1p).m_center));
			if(list1.entries()){
				MGCCisect isect=list1.first();
				double t2; MGPosition P=isect.point();
				if((*e2p).on(P,t2)){
					if(r1>=r2) list.append(P,isect.param1(),t2);
					else       list.append(P,t2,isect.param1());
				}
			}
			list1=(*e1p).isect
				(MGStraight(MGSTRAIGHT_HALF_LIMIT,vec2,(*e1p).m_center));
			if(list1.entries()){
				MGCCisect isect=list1.first();
				double t2; MGPosition P=isect.point();
				if((*e2p).on(P,t2)){
					if(r1>=r2) list.append(P,isect.param1(),t2);
					else       list.append(P,t2,isect.param1());
				}
			}
		}
		return list;
	}else return intersect(e2);
}

// Ellipse と Curve の交点を求める。
MGCCisect_list MGEllipse::isect(
	const MGSurfCurve& curve2
)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

// Ellipse と Curve の交点を求める。
MGCCisect_list MGEllipse::isect(
	const MGBSumCurve& curve2
)const{
	MGCCisect_list list=curve2.isect(*this);
	return list.replace12();
}

//Intersection of Ellipse and Surface
MGCSisect_list MGEllipse::isect(const MGSurface& surf) const{
	return surf.isect(*this);
}

//Intersection of Ellipse and a plane.
MGCSisect_list MGEllipse::isect(const MGPlane& plane) const{
	return plane.isect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGEllipse::isect(const MGSphere& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGEllipse::isect(const MGCylinder& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGEllipse::isect(const MGRSBRep& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGEllipse::isect(const MGSBRep& surf)const{
	return surf.intersect(*this);
}

//Intersection of Spline and Surface
MGCSisect_list MGEllipse::isect(const MGBSumSurf& surf)const{
	return surf.intersect(*this);
}

//Compute intersection point of 1D sub curve of original curve.
//Parameter values of intersection point will be returned.
MGCParam_list MGEllipse::intersect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
)const{
	MGCParam_list list(this);

	double a=m_m.ref(coordinate), b=m_n.ref(coordinate);
	double c=m_center.ref(coordinate);
	double ab2=a*a+b*b;
	double sq_ab2=sqrt(ab2);
	if(MGMZero(sq_ab2)) return list;

//Intersection can be obtained as the solution of
//		c+a*cos(angle)+b*sin(angle) = f	, i.e.
//		cos(angle)=(a*(f-c) +- b*sqrt(a*a+b*b-(d-c)**2))/(a*a+b*b)
//		sin(angle)=((f-c)-a*cos(angle))/b
//				  =(b*(f-c) -+ a*sqrt(a*a+b*b-(d-c)**2))/(a*a+b*b)   .

	double fmc=f-c;
	double cosa, sina, angle;
	if(MGREqual2(fabs(fmc),sq_ab2)){
		if(fmc<0.){ cosa=-a/sq_ab2; sina=-b/sq_ab2;}
		else      { cosa= a/sq_ab2; sina= b/sq_ab2;}
		angle=MGAngle(cosa,sina);
		if(in_RelativeRange_of_radian(angle))
			list.append(radian_to_gp(RelativeRange_in_radian(angle)));
	}else{
		if(fmc<(-sq_ab2) || fmc>sq_ab2) return list;

		double ab2mfmc2=sqrt(ab2-fmc*fmc);
		cosa=(a*fmc+b*ab2mfmc2)/ab2; sina=(b*fmc-a*ab2mfmc2)/ab2;
		angle=MGAngle(cosa,sina);
		if(in_RelativeRange_of_radian(angle))
			list.append(radian_to_gp(RelativeRange_in_radian(angle)));
		cosa=(a*fmc-b*ab2mfmc2)/ab2; sina=(b*fmc+a*ab2mfmc2)/ab2;
		angle=MGAngle(cosa,sina);
		if(in_RelativeRange_of_radian(angle))
			list.append(radian_to_gp(RelativeRange_in_radian(angle)));
	}
	return list;
}

// 与ポイントから楕円への垂線の足 のパラメータ値を求める。
// Function's return value is !=0 if point is obtained, 0 if not.
int MGEllipse::perp_point (
	const MGPosition &p,	// 与ポイント
	double &d,				// パラメータ値
	const double *g			// guess parameter value.
	) const {
	size_t i;
	if(ellipse_type() == MGELLIPSE_EMPTY) {d = 0.0; return 0;}

	double theta[4];
	MGPosition p2d=p*MGTransf(m_m, m_n, m_center);
	size_t nump=perp2d(p2d.ref(0), p2d.ref(1), theta);
	d=theta[0];

	if(nump){
		size_t j=0;
		if(nump>1){
			if(g){	//Compute nearest param to *g
				double tg=gp_to_radian(*g);
				double param=RelativeRange_in_radian(tg);
				double plen=param_length(param, theta[0]);
				for(i=1; i<nump; i++){
					double plent=param_length(param, theta[i]);
					if(plent<plen){ j=i; plen=plent;};
				}
			}else{  //Find minimum length point from p.
				double len=(p-eval_in_radian2(theta[0])).len();
				for(i=1; i<nump; i++){
					double lent=(p-eval_in_radian2(theta[i])).len();
					if(lent<len){ j=i; len=lent;};
				}
			}
		}
		d=theta[j];
	}
	
	d=radian_to_gp(d);
	return nump;
}
	
// 与ポイントから楕円へ下ろした垂線の足の，楕円のパラメータ値を
// すべて求める。
MGCParam_list MGEllipse::perps(
const MGPosition& point		// Given point.
	) const{
	MGPosition p2d=point*MGTransf(m_m, m_n, m_center);
	double theta[4];		// Intersection points are maixmum 4.
	size_t nump=perp2d(p2d.ref(0), p2d.ref(1), theta);
	MGCParam_list tlist(this);
	for(size_t i=0; i<nump; i++) tlist.append(radian_to_gp(theta[i]));
	return tlist;	
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv2's parameter
//as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list MGEllipse::perps(
	const MGCurve& crv2		//The second curve
)const{
	MGPosition_list list=crv2.perps(*this);
	return MGPosition_list(list,1,0);	
}

MGPosition_list MGEllipse::perps(
	const MGStraight& sl		//The second curve
)const{
	int n;
	MGPosition_list list;

	MGUnit_vector sl_dir=sl.direction(), el_normal=normal();
	if(sl_dir.parallel(el_normal)){
	//1. When sl is normal to ellipse plane.
		MGVector lvec=sl.direction();
		MGPosition lpoint=sl.root_point();
		double r=lvec%m_normal;
		if(!MGMZero(r)){
			double d=m_normal%m_center;
			double s=(d-lpoint%m_normal)/r;
			//s is parameter of sl, of intersection of ellipse plane and sl.
			if(sl.in_range(s)){
				double eparam;
				MGPosition ipoint=lpoint+lvec*s;//ipoint is the intersection.
				MGCParam_list plist=perps(ipoint);
				int m=n=plist.entries();
				while(n--){
					eparam=plist.removeFirst();
					list.append(MGPosition(eparam,s));
				}
				if(m<2){
					if(m==0){
						list.append(MGPosition(param_s(),s));
						list.append(MGPosition(param_e(),s));
					}else{
						MGPosition sp=start_point(), ep=end_point();
						MGPosition p=eval(eparam);
						double lenp=(ipoint-p).len();
						double lensp=(ipoint-sp).len();
						double lenep=(ipoint-ep).len();
						if(lenp<lensp){
							if(lensp>lenep)
								list.append(MGPosition(param_s(),s));
							else list.append(MGPosition(param_e(),s));
						}else{
							if(lensp<lenep)
								list.append(MGPosition(param_s(),s));
							else list.append(MGPosition(param_e(),s));
						}
					}
				}
			}
		}
	}else if(sl_dir.orthogonal(el_normal)){
	//2. When sl and ellipse plane are parallel.
		double h1=MGDeterminant(major_axis(),sl_dir,el_normal);
		double h2=MGDeterminant(minor_axis(),sl_dir,el_normal);
		double h=sqrt(h1*h1+h2*h2);
		double ca=h1/h, sa=h2/h;
		double t1=MGAngle(ca,sa), t2=MGAngle(-ca,-sa);
		//t1 and t2 are two angles of the ellipse whose points are
		//perpendicular to straight line sl.
		if(in_RelativeRange_of_radian(t1)){
			MGPosition p=eval_in_radian2(t1=RelativeRange_in_radian(t1));
			double s=sl.perp_param(p);
			if(sl.in_range(s)) list.append(MGPosition(radian_to_gp(t1),s));
		}
		if(in_RelativeRange_of_radian(t2)){
			MGPosition p=eval(t2=RelativeRange_in_radian(t2));
			double s=sl.perp_param(p);
			if(sl.in_range(s))
				list.append(*this,sl,MGPosition(radian_to_gp(t2),s));
		}
		//Add intersetion points.
		MGCCisect_list ilist=isect(sl);
		n=ilist.entries();
		while(n--){
			MGCCisect is=ilist.removeFirst();
			list.append(*this,sl,MGPosition(is.param1(),is.param2()));
		}
	}else{
	//3.General case.
		double ts=m_prange[0], te;
		double range=m_prange[1]-ts;
		size_t ndiv=int(1.+2.*range/mgPAI); if(ndiv>4) ndiv=4;
		double incr=range/double(ndiv);
		te=ts+incr;
		MGPosition ppair;

		//Compute by subdividing the ellipse to parts of half pai at most.
		MGEllipse el(*this, true);
		for(size_t i=0; i<ndiv; i++){
			el.m_prange[0]=ts; el.m_prange[1]=te;
			double guess=(ts+te)/2.;
			if(el.perp_guess(ts,te,sl,1.,0.,guess,0.,ppair)){
				ppair(0)=radian_to_gp(ppair[0]);
				list.append(*this,sl,ppair);
			}
			ts+=incr; te+=incr;
		}
	}
	return list;
}

// 点が楕円上にあるか調べる。楕円上であれば，そのパラメータ値を，
// そうでなくても最近傍点のパラメータ値を返す。
bool MGEllipse::on(
	const MGPosition & point,	// 指定点
	double& param				// パラメータ
	) const {
	MGTransf tr=MGTransf(m_m, m_n, m_center);
	MGPosition point2d=point*tr;
	double p=point2d.ref(0), q=point2d.ref(1);
	double a=m_m.len();
	MGVector Vb=m_n*tr.affine();
	double bx=Vb.ref(0), by=Vb.ref(1);
	double t;	//Perp point in radian.
	bool on;
	if(MGAZero(a) && MGAZero(Vb.len())){
		on=0; t=m_prange[0];
	}else{
		double d1=(p*by-q*bx); d1*=d1;
		double d2=a*a*(by*by-q*q);
		if(MGREqual2(d1, d2)){
			if(MGAZero(a)){
				if(fabs(bx)>=fabs(by)) t=asin(p/bx);
				else				t=asin(q/by);
			}else if(MGAZero(by)) t=0.;
			else t=MGAngle((p-bx/by*q)/a, q/by);
			on=in_RelativeRange_of_radian(t);
			t=RelativeRange_in_radian(t);
		}else on=0;
		if(!on){	// Compute nearest point.
			t=m_prange[0];
			double len=(point-eval_in_radian2(t)).len();
			MGCParam_list tlist=perps(point);
			MGCParam_list::Citerator i=tlist.begin(), iend=tlist.end();
			double lent;
			while(i!=tlist.end()){
				double t1=*i++;
				lent=(point-eval_position(t1)).len();
				if(lent<len){ t=gp_to_radian(t1); len=lent;};
			}
			lent=(point-end_point()).len();
			if(len>lent){
				t=m_prange[1];
				len=lent;
			}
			if(len<=MGTolerance::wc_zero()) on=1;
		}
	}
	param=radian_to_gp(t);
	return on;
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> MGEllipse::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
) const{
	MGVector G(3,g);
	MGEllipse* el1=new MGEllipse();
	MGVector one(1);
	one(0)=G%center()-g[3]; el1->m_center=one;
	el1->m_circle=0;
	if(m_gprange){
		double* rng=new double[2]; rng[0]=m_gprange[0]; rng[1]=m_gprange[1];
		el1->m_gprange=rng;
	}
	one(0)=G%major_axis(); el1->m_m=one;
	one(0)=G%minor_axis(); el1->m_n=one;
	one(0)=1.; el1->m_normal=one;
	el1->m_prange[0]=m_prange[0]; el1->m_prange[1]=m_prange[1];
	el1->m_r=(el1->m_m[0]+el1->m_n[0])*.5;
	return std::auto_ptr<MGCurve>(el1);
}

//Compute intesection points of  2D ellipse whose center is origin and
//a straight line of 2D.
int MGEllipse::isect2d(
	const MGPosition& sp,//Start point of straight line.
	const MGVector& dir, //Direction unit vector of the straight line
	double t[2],		//Parameter of intersection point of the straight.
	double angle[2],	// angles of ellipse in radian.
	int& tangen			// Return if isect is tangent point(1) or not(0).
	) const{
//Return value of isect2D is number of intersetion points: 0, 1, or 2.

	double a=m_m.len(), b=m_n.len(); //Length of x and y axis.
		//parameter arrangement value of the straight when using sp1 instead of sp.
	if(m_m.orthogonal(m_n)){//When normal ellipse.

	double dir2=dir%dir;
	double spdirbydir2=(sp%dir)/dir2;
	MGPosition sp1=sp-dir*spdirbydir2;
		//sp1 is nearest point to the origin.
		//sp1 is used instead of sp to minimize computation error.
	double sp1tosp=(sp1%dir)/dir2-spdirbydir2;

	double u=dir.ref(0), v=dir.ref(1);
	double x1=sp1.ref(0), y1=sp1.ref(1);

//Using a, b, u, v, x1, and y1, isect of the ellipse and the straight line whose 
//root point and direction are sp1 and dir can be expressed as:
//  t= (-(u*b*b*x1+v*a*a*y1)(+-)a*b*sqrt(avbu-(u*y1-v*x1)**2))/avbu
//  , where avbu=a*a*v*v+b*b*u*u.
	double ca,sa, ang;
	double avbu=a*a*v*v+b*b*u*u;
	double uyvx=(u*y1-v*x1); uyvx *= uyvx;
	double ubxvay=-(u*b*b*x1+v*a*a*y1);
	double p=uyvx/avbu, q=ubxvay/avbu;
	if(MGREqual(p, 1.)){		   // Point is tangent point.
		tangen=1;
		t[0]=q+sp1tosp;
		ca=(x1+u*q)/a; sa=(y1+v*q)/b;
		ang=MGAngle(ca,sa);
		if(in_RelativeRange_of_radian(ang)){
			angle[0]=RelativeRange_in_radian(ang); return 1;
		}else return 0;
	}else if(p>1.) return 0;
	else{					  // Points are not tangent points.
		tangen=0;
		int count=0;
		double r=a*b*sqrt((1.-p)/avbu);
		t[0]=q+r;
		ca=(x1+u*t[0])/a; sa=(y1+v*t[0])/b;
		t[0]+=sp1tosp;
		ang=MGAngle(ca,sa);
		if(in_RelativeRange_of_radian(ang)){
			angle[0]=RelativeRange_in_radian(ang); count+=1;
		}
		t[count]=q-r;
		ca=(x1+u*t[count])/a; sa=(y1+v*t[count])/b;
		t[count]+=sp1tosp;
		ang=MGAngle(ca,sa);
		if(in_RelativeRange_of_radian(ang)){
			angle[count]=RelativeRange_in_radian(ang); count+=1;
		}
		return count;
	}

	}else{//When Projected ellipse whose axises are not orthogonal.
		MGStraight sl(MGSTRAIGHT_UNLIMIT, dir, sp);
		MGStraight sl2(sl);
		double apb=a+b;
		MGBox box(MGPosition(-apb,-apb), MGPosition(apb,apb));
		sl2.limit(box);

		MGEllipse el(*this,true);//el's parameter range is in radian.
		el*=MGTransf(m_m,m_n,m_center);
		MGCCisect_list list=el.intersect(sl2);
		tangen=0; size_t i=0;
		while(list.entries()){
			MGCCisect is=list.removeFirst();
			angle[i]=is.param1(); t[i++]=sl.param(is.point());
			if(i==2) break;
		}
		return i;
	}
}

//Compute angles of elllipse points that are perpendicular to a point (p,q).
//This function is 2-D version of perp_point.
//Function's return value is number of points obtained.
size_t MGEllipse::perp2d(
	double p,			//x-coordinate of the point.
	double q,			//y-coordinate of the point.
	double theta[4])	//angles of the ellipse.
const{
	if(!m_m.orthogonal(m_n)){//When two axises are not orthogonal.
		MGPosition P(p,q);
		MGEllipse el(*this, true);
		el*=MGTransf(m_m,m_n,m_center);
		MGCParam_list list=el.MGCurve::perps(P);
		size_t i=0;
		while(list.entries()){
			theta[i++]=list.removeFirst();
			if(i==4) break;
		}
		return i;
	}

	size_t nump=0;
	int pzero=MGRZero2(p,m_r), qzero=MGRZero2(q,m_r);
	double plen,pu,qu,pangle;
	if(!pzero || !qzero){
		plen=sqrt(p*p+q*q);
		pu=p/plen, qu=q/plen;
		pangle=MGAngle(pu, qu);
	}

	if(pzero&&qzero){		//Case that (p,q) is origin.
		theta[1]=mgHALFPAI; theta[3]=mgHALFPAI+mgPAI;
		theta[0]=0.0; theta[2]=mgPAI;
		nump=4;
	} else if(pzero){		//Case that the point is on y-axis.
		theta[0]=mgHALFPAI; theta[1]=mgHALFPAI+mgPAI;
		nump=2;
	} else if(qzero){		//Case that the point is on x-axis.
		theta[0]=0.0; theta[1]=mgPAI;
		nump=2;
	} else if(m_circle){	//Case that the ellipse is a circle.
		theta[0]=pangle;
		if(pangle>mgPAI) theta[1]=pangle-mgPAI;
		else           theta[1]=pangle+mgPAI;
		nump=2;
	} else{

	// Genaral case.
	double a=m_m.len(), b=m_n.len(); //Length of x and y axis.
	double a2mb2=a*a-b*b, ap=a*p, bq=b*q, sinv, cosv;
	double sign1, sign2, ang1, ang2;
	double angs;
	double error=MGTolerance::angle_zero();
	double ang, signt;

	if(p>=0.){
		if(q>=0.) angs=0.0;
		else	  angs=mgHALFPAI+mgPAI;
	}
	else {
		if(q>=0.) angs=mgHALFPAI;
		else	  angs=mgPAI;
	}			   

	//Angle is computed using bisection algorithm. Sign of the following
	//expression is use to determine if solution can be located between
	//two angle ang1 and ang2:
	//	sign=(a*a-b*b)*cosv*sinv-a*p*sinv+b*q*cosv;
	//This sign is signed length of the vector product of the two vectors,
	//vector1:from a point P on the ellipse to the point (p,q),
	//vector2:normal vector of the ellipse on P.

	//Case1: Angle is smaller than angle of (p,q) and 
	//greater than nearest half PAI.
	ang1=angs; sinv=sin(ang1); cosv=cos(ang1);
	sign1=a2mb2*cosv*sinv-ap*sinv+bq*cosv;
	ang2=pangle; sinv=qu;  cosv=pu ;
	sign2=a2mb2*cosv*sinv-ap*sinv+bq*cosv;
	if(sign1*sign2 <= 0.0){
		while((ang2-ang1)>error){
			ang=(ang1+ang2)/2.; sinv=sin(ang); cosv=cos(ang);
			signt=a2mb2*cosv*sinv-ap*sinv+bq*cosv;
			if(sign1*signt <= 0.0) ang2=ang;
			else                   ang1=ang;
		}
		theta[nump]=(ang1+ang2)/2.;
		nump +=1;
	}

	//Case2: Angle is greater than angle of (p,q) and 
	//smaller than nearest half PAI.
	ang1=pangle; sinv=qu;  cosv=pu ;
	sign1=a2mb2*cosv*sinv-ap*sinv+bq*cosv;
	ang2=angs+mgHALFPAI; sinv=sin(ang2); cosv=cos(ang2);
	sign2=a2mb2*cosv*sinv-ap*sinv+bq*cosv;	
	if(sign1*sign2 <= 0.0){
		while((ang2-ang1)>error){
			ang=(ang1+ang2)/2.; sinv=sin(ang); cosv=cos(ang);
			signt=a2mb2*cosv*sinv-ap*sinv+bq*cosv;
			if(sign1*signt <= 0.0) ang2=ang;
			else                   ang1=ang;
		}
		theta[nump]=(ang1+ang2)/2.;
		nump +=1;
	}

	//Case3:Check all other three half pai ranges.
	for(int i=0; i<3; i++){
		angs += mgHALFPAI; angs=fmod(angs, mgDBLPAI);
		ang1=angs; sinv=sin(ang1); cosv=cos(ang1);
		sign1=a2mb2*cosv*sinv-ap*sinv+bq*cosv;	
		ang2=angs+mgHALFPAI; sinv=sin(ang2); cosv=cos(ang2);
		sign2=a2mb2*cosv*sinv-ap*sinv+bq*cosv;	
		if(sign1*sign2 <= 0.0){
			while((ang2-ang1)>error){
				ang=(ang1+ang2)/2.; sinv=sin(ang); cosv=cos(ang);
				signt=a2mb2*cosv*sinv-ap*sinv+bq*cosv;
				if(sign1*signt <= 0.0) ang2=ang;
				else                   ang1=ang;
			}
			theta[nump]=(ang1+ang2)/2.;
			nump +=1;
		}
	}
	}

	double save=RelativeRange_in_radian(theta[0]);
	size_t i=0;
	while(i<nump){
		if(!in_RelativeRange_of_radian(theta[i])){
			nump -= 1;
			for(size_t j=i; j<nump; j++) theta[j]=theta[j+1];
		}else i+=1;
	}
	if(nump==0) theta[0]=save;
	for(i=0; i<nump; i++) theta[i]=RelativeRange_in_radian(theta[i]);
	return nump;
}
