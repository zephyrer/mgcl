/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Transf.h"
#include "mg/Position.h"
#include "mg/CParam_list.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/Position_list.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Knot.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/Sphere.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/Blipp.h"
#include "cskernel/Bkdnp.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLBRep.cpp
//
// Implement MGLBRep class.

//Member Function

//Compute continuity with brep2. maximum continuity conputed is 2.
int MGLBRep::continuity(const MGLBRep& brep2, int& which, double& ratio)const{
	const double start1=param_s();
	const double end1=param_e();
	const double start2=brep2.param_s();
	const double end2=brep2.param_e();
	double* data=new double[sdim()>3?sdim():3];
	data[0]=data[1]=data[2]=0.0;
	int cn;
	MGVector deriv1(0.,0.,0.), deriv2(0.,0.,0.);
	MGVector d21(0.,0.,0.), d22(0.,0.,0.);
	double dlen1,dlen2;

	which=2; ratio=1.0;
//Determine which point of brep1 is continuous to which
//point of brep2.
	MGPosition p1[2]; MGPosition p2[2];
	p1[0]=eval(start1, 0); p1[1]=eval(end1, 0); 
	p2[0]=brep2.eval(start2, 0); p2[1]=brep2.eval(end2, 0); 

	double tratio;
// Check of C-0 continuity.
	const MGLBRep *bl1, *bl2; double t1,t2;
	if(p1[1]==p2[0]) {
		which=2;		//which=2 means end of brep1 connects
						// to start of brep2.
		bl1=this; t1=end1;
		bl2=&brep2; t2=start2;
	}else if(p1[1]==p2[1]) {
		which=3;		//which=3 means end of brep1 connects
						// to end of brep2.
		bl1=this; t1=end1;
		bl2=&brep2; t2=end2;
	}else if(p1[0]==p2[1]) {
		which=1;		//which=1 means start of brep1 connects
						// to end of brep2.
		bl1=&brep2; t1=end2;
		bl2=this; t2=start1;
	}else if(p1[0]==p2[0]) {
		which=0;		//which=0 means start of brep1 connects
						// to start of brep2.
		bl1=&brep2; t1=start2;
		bl2=this; t2=start1;
	}else {cn=-1; goto end_process;}

// Check of C-1 continuity.
	deriv1=bl1->eval(t1, 1); dlen1=deriv1.len();
	deriv2=bl2->eval(t2, 1); dlen2=deriv2.len();
	if(MGRZero(dlen1) || MGRZero(dlen2)) {cn=0; goto end_process;}

	if(which==0)		deriv1 *= -1.;
	else if(which==3)	deriv2 *= -1.;

	tratio=dlen2/dlen1;
	if(deriv1.parallel(deriv2) && (deriv1%deriv2 > 0.)){
		d21=bl1->eval(t1, 2);
		d22=bl2->eval(t2, 2);
		d22 /=(tratio*tratio);
		if(d22==d21) cn=2; else cn=1;
	} else cn=0;
	if(which<=1) ratio=dlen1/dlen2;
	else ratio=tratio;

	end_process: delete[] data;
				 return cn;
}

//Provide divide number of curve span for function intersect.
size_t MGLBRep::intersect_dnum() const{
	size_t k=order();
	int nspan=bdim()+1-k;
	int km2=k-2; if(km2<=0) km2=1;
	return nspan*km2+1;
}

//Intersection point of spline and curve.
MGCCisect_list MGLBRep::isect(const MGCurve& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//Intersection point of spline and curve.
MGCCisect_list MGLBRep::isect(const MGStraight& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//Intersection point of spline and curve.
MGCCisect_list MGLBRep::isect(const MGRLBRep& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//Intersection point of spline and curve.
MGCCisect_list MGLBRep::isect(const MGEllipse& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//Intersection point of spline and curve.
MGCCisect_list MGLBRep::isect(const MGLBRep& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//Intersection point of spline and curve.
MGCCisect_list MGLBRep::isect(const MGSurfCurve& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//Intersection point of spline and curve.
MGCCisect_list MGLBRep::isect(const MGBSumCurve& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//compute isects by splitting this curve to sub-curves that do not have
//C0 points in it.
MGCCisect_list MGLBRep::isect_by_split_to_C1(const MGCurve& curve2)const{
	MGCCisect_list list(this, &curve2);
	if(!has_common(curve2))
		return list;

	size_t k=order();
	if(k==2)
		list.append(isect_order2(curve2));
	else{
		size_t index, n=bdim(), multi_found;
		size_t orderm1=k-1, start=k;
		const MGKnotVector& t1=knot_vector();
		do{	//Compute intersections by dividing to parts of continuity>=C1.
			multi_found=t1.locate_multi(start,orderm1,index);
											//Locate C0 continuity point.
			if(start==k && index==n)
				list.append(C1isect(curve2));	//Use original.
			else{
				MGLBRep spanC1(t1[start-1],t1[index],*this);
				list.append(spanC1.C1isect(curve2));
			}
			start=index+multi_found;
		}while(index<n);
	}
	return list;
}

//isect for this LBRep that does not include C0 continuity points in it.
MGCCisect_list MGLBRep::C1isect(const MGCurve& curve2) const{
	MGCCisect_list list=curve2.isect_withC1LB(*this);
	return list.replace12();
}

//Compute intersection points of Line B-Rep and a straight line.
MGCCisect_list MGLBRep::C1isect(const MGStraight& line) const{
	MGCCisect_list list(this, &line);
	if(!has_common(line))
		return list;

	MGPosition p; double t1,t2;
	MGCParam_list clist;
	if(sdim()<=2 && line.sdim()<=2){
		MGUnit_vector direction=line.direction();
		MGVector normal(-direction.ref(1),direction.ref(0));
			//Normal to line.
		double dist=normal%line.root_point(); //Distance from the origin.
		MGMatrix mat; mat.to_axis(direction);
		MGLBRep temp=(*this)*mat;//cout<<mat<<temp<<temp.box()<<endl;
		clist=temp.isect_1D(dist,1);
//		if(clist.size()) cout<<eval(clist.front())<<endl;
	}else{
		MGPlane plane; MGStraight sl;
		MGCCisect_list cclist;
		switch(planar(plane, sl, p)){
		case 1:		//BRep is a point.
			if(line.on(p, t2))
				list.append(MGCCisect(p,param_s(),t2));
			return list;
		case 2:		//BRep is a straight line.
			cclist=sl.isect(line);
			if(cclist.entries()){
				MGCCisect isect=cclist.first(); p=isect.point();
				if(on(p, t1))
					list.append(MGCCisect(p,t1,isect.param2()));
			}
			return list;
		default:	//BRep is on a plane or a general 3D line.
			MGUnit_vector x_axis=line.direction();
			MGUnit_vector y_axis=plane.normal()*x_axis;
			MGTransf tr(x_axis, y_axis, line.root_point());
			MGLBRep lb2=(*this)*tr;
			clist=lb2.isect_1D(0.,1);
		}
	}

	MGCParam_list::Citerator i;
	for(i=clist.begin(); i!=clist.end(); i++){
		t1=(*i); p=eval(t1);
		if(line.on(p,t2)) list.append(MGCCisect(p,t1,t2));
	}

	return list;
}

//Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list MGLBRep::isect_withC1LB(const MGLBRep& curve2)const{
	MGCCisect_list list(this, &curve2);
	if(!has_common(curve2))
		return list;

	size_t k=order();
	if(k==2)
		list.append(isect_order2(curve2));
	else{
		size_t index, n=bdim(), multi_found;
		size_t orderm1=k-1, start=k;
		const MGKnotVector& t1=knot_vector();
		do{	//Compute intersections by dividing to parts of continuity>=C1.
			multi_found=t1.locate_multi(start,orderm1,index);
											//Locate C0 continuity point.
			if(start==k && index==n)
				list.append(intersect(curve2));	//Use original.
			else{
				MGLBRep spanC1(t1[start-1],t1[index],*this);
				list.append(spanC1.intersect(curve2));
			}
			start=index+multi_found;
		}while(index<n);
	}
	return list;
}

//isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list MGLBRep::isect_with_noCompoSC(const MGSurfCurve& curve2)const{
	return isect_by_split_to_C1(curve2);
}

//isect_order2 is a private function for isect, computes intersection of
//LBRep(this) of order 2, i.e. polyline B-Rep, with crv.
MGCCisect_list MGLBRep::isect_order2(const MGCurve& crv) const{
	MGCCisect_list list(this,&crv);
	for(size_t i=0; i<bdim()-1; i++){
		MGStraight sl(coef(i+1),coef(i));	//i-th line segment.
		double clen=sl.param_e();
		if(!MGMZero(clen)){
			MGCCisect_list ls=sl.isect(crv); //isect with i-th line segment.
			while(ls.entries()){
				MGCCisect is=ls.removeFirst();
				double t1=knot(i+1);
				t1=t1+(is.param1()/clen)*(knot(i+2)-t1);
				list.append(is.point(),t1,is.param2(),is.rel());
			}
		}
	}
	return list;
}

//isect with a surface.
MGCSisect_list MGLBRep::isect(const MGSurface& surf) const{
	return isect_by_split_to_C1(surf);
}

//isect with a surface.
MGCSisect_list MGLBRep::isect(const MGPlane& surf) const{
	return isect_by_split_to_C1(surf);
}

//isect with a surface.
MGCSisect_list MGLBRep::isect(const MGSphere& surf) const{
	return isect_by_split_to_C1(surf);
}

//isect with a surface.
MGCSisect_list MGLBRep::isect(const MGCylinder& surf) const{
	return isect_by_split_to_C1(surf);
}

//isect with a surface.
MGCSisect_list MGLBRep::isect(const MGSBRep& surf) const{
	return isect_by_split_to_C1(surf);
}

//isect with a surface.
MGCSisect_list MGLBRep::isect(const MGRSBRep& surf) const{
	return isect_by_split_to_C1(surf);
}

//isect with a surface.
MGCSisect_list MGLBRep::isect(const MGBSumSurf& surf) const{
	return isect_by_split_to_C1(surf);
}

//compute isects by splitting this curve to sub-curves that do not have
//C0 points in it.
MGCSisect_list MGLBRep::isect_by_split_to_C1(const MGSurface& surf) const{
	MGCSisect_list list(this, &surf);
	if(!has_common(surf))
		return list;

	size_t k=order();
	if(k==2)
		return isect_order2(surf);

	size_t index, n=bdim(), multi_found;
	size_t orderm1=k-1, start=k;
	const MGKnotVector& t1=knot_vector();
	do{	//Compute intersections by dividing to parts of continuity>=C1.
		multi_found=(knot_vector()).locate_multi(start,orderm1,index);
										//Locate C0 continuity point.
		if(start==k && index==n)
			list.append(surf.isect_withC1LB(*this));	//Use original.
		else{
			MGLBRep spanC1(t1[start-1],t1[index],*this);
			list.append(surf.isect_withC1LB(spanC1));
		}
		start=index+multi_found;
	}while(index<n);
	return list;
}

//Compute intersection points of 1D sub B-Rep of original B-rep.(BLIPP)
//Parameter values of intersection points will be returned.
//isect_1D covers this LBRep's C0 ontinuity.
MGCParam_list MGLBRep::intersect_1D(						
	double f,			// Coordinate value
	size_t coordinate	// Coordinate kind of the data f(from 0).
						// Id of m_line_bcoef.
)const{
	assert(coordinate<sdim());
	//std::cout<<(*this)<<std::endl;

	MGCParam_list tlist(this);

	const double error=MGTolerance::wc_zero()*.7;
	unsigned k=order(), n=bdim();
	if(k<2 || n<k) return tlist;
	if(k==2){

	//Case when polyline.
	size_t nm1=n-1; double f0,f1;
	for(size_t i=0; i<nm1; i++){
		f0=coef(i,coordinate); f1=coef(i+1,coordinate);
		double fmf0=f-f0, fmf1=f-f1;
		double t0=knot(i+1), t1=knot(i+2);
		if(fmf0*fmf1<=0.){
			double fdif=f1-f0;
			if(MGMZero(fdif)) tlist.append((t0+t1)*.5);
			else tlist.append(t0+(fmf0/fdif)*(t1-t0));
		}else{
			//Check if end points are within error.
			if(fabs(fmf0)<error) tlist.append(t0);
			if(fabs(fmf1)<error) tlist.append(t1);
		}
	}

	}else{

	//case when normal spline.
	double* work=new double[4*k*k+3*k];
	double* x_temp=new double[n];
	int nx; int iend;

	size_t index, orderm1=k-1, multi_found;
	size_t start=k, startId, nspan;
	double ts;
	const double* knotp; const double* coefp;
	const MGKnotVector& t=knot_vector();
	double delta=param_e()-param_s();
	do{	//Compute intersections by dividing to parts of continuity>=C1.
		multi_found=t.locate_multi(start,orderm1,index);
									//Locate C0 continuity point.
		if(start==k && index==n){	//Use original.
			nspan=n; knotp=knot_data(); coefp=m_line_bcoef.data(0,coordinate);
			ts=param_s(); ts-=delta;
			blipp_(k,nspan,knotp,coefp,f,error,ts,n,work,&nx,x_temp,&iend);
			for(int i=0; i<nx; i++) tlist.append(x_temp[i]);
			break;
		}
		startId=start-k;
		nspan=index-startId;
		if(nspan>=k && t(start-1)<t(index)){
			knotp=t.data(startId);
			coefp=m_line_bcoef.data(startId,coordinate);
			ts=t(start-1); ts-=delta;
			blipp_(k,nspan,knotp,coefp,f,error,ts,n,work,&nx,x_temp,&iend);
			for(int i=0; i<nx; i++) tlist.append(x_temp[i]);
		}
		start=index+multi_found;
	}while(index<n);

	delete[] work; delete[] x_temp;

	}
	return tlist;
}

MGCSisect_list MGLBRep::isect_order2(const MGSurface& srf) const{
	MGCSisect_list list(this,&srf);
	for(size_t i=0; i<bdim()-1; i++){
		MGStraight sl(coef(i+1),coef(i));	//i-th line segment.
		double clen=sl.param_e();
		if(!MGMZero(clen)){
			MGCSisect_list ls=sl.isect(srf); //isect with i-th line segment.
			while(ls.entries()){
				MGCSisect is=ls.removeFirst();
				double t1=knot(i+1);
				t1=t1+(is.param_curve()/clen)*(knot(i+2)-t1);
				list.append(is.point(),t1,is.param_surface(),is.rel());
			}
		}
	}
	return list;
}

//Obtain so transformed 1D curve expression of this curve that
//f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
//of oneD and xi(t) is i-th coordinate expression of this curve.
//This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> MGLBRep::oneD(
	const double g[4]			//Plane expression(a,b,c,d) where ax+by+cz=d.
) const{
	size_t i, n=bdim();
	size_t j, m=sdim(); if(m>3) m=3;
	size_t k=order();
	MGLBRep* brep=new MGLBRep();
	MGBPointSeq& rcoef=brep->line_bcoef(); rcoef.resize(n,1);
	MGKnotVector& knotv=brep->knot_vector(); knotv.size_change(k,n);

	double min,max;
	const MGKnotVector& t=knot_vector();
	double rt=0; for(j=0; j<m; j++) rt+=coef(0,j)*g[j];
	min=max=rt-=g[3];
	rcoef(0,0)=rt;
	knotv(0)=t(0);
	for(i=1; i<n; i++){
		rt=0; for(j=0; j<m; j++) rt+=coef(i,j)*g[j];
		rt-=g[3];
		rcoef(i,0)=rt;
		if(min>rt) min=rt; if(max<rt) max=rt;
		knotv(i)=t(i);
	}
	size_t npk=n+k;
	for(i=n; i<npk; i++) knotv(i)=t(i);
	MGInterval minmax(min,max);
	brep->m_box=new MGBox(1,&minmax);
	return std::auto_ptr<MGCurve>(brep);
}

// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
// すべて求める。
MGCParam_list MGLBRep::perps(
	const MGPosition& point	// 与ポイント
)const{
	MGCParam_list tlist(this);
	
	// Points normal from a point can be obtained as extremum of
	// the folowing expression G**2(square of length from the point P).
	// G=(f-p)=sum of (Pij-Pj)*Bi about i, where j is index of axis.
	// Hence the first derivative of the square is:
	// 2*(sum about j of Gj*Gj'). The extremum is obtained by solving;
	// (sum about j of Gj*Gj')=0 .

	//std::cout<<(*this);
	//1. Compute knot vector of the B-rep of (sum about j of Gj*Gj').
	size_t k=order(); size_t km1=k-1; size_t nold=bdim();
	size_t knew=k+km1; size_t nnew=nold+km1;
	MGKnotVector t(m_knot_vector);		//std::cout<<m_knot_vector;/////////////
	t.add_data(MGKnot(m_knot_vector(km1),km1), knew);//std::cout<<t;
	t.add_data(MGKnot(m_knot_vector(nold), km1),knew);//std::cout<<t;
	t.change_order(knew); t.set_bdim(nnew);	 //std::cout<<t<<std::endl;/////////////

	//2. Compute data point.
	MGNDDArray dtp; dtp.update_from_knot(t);
	//3. Compute (sum about j of Gj*Gj') for each data point.
	MGBPointSeq gbydg(nnew,1); MGVector g(sdim()), dg(sdim());
	for(size_t i=0; i<nnew; i++){
		g=eval(dtp(i))-point; dg=eval(dtp(i),1);
		gbydg(i,0)=g%dg;
	}
	//4. Obtain the B-rep and the extremum.
			//std::cout<<dtp<<gbydg;/////////
	int n=gbydg.length(), pointSize=gbydg.capacity(), pointDim=1, IMLT=1;
	double ratio=5.;
	bkdnp_(&n,dtp.data(),gbydg.data(),pointSize,pointDim,IMLT,ratio);
			//std::cout<<dtp<<gbydg;/////////
	dtp.set_length(n); gbydg.set_length(n);
	MGLBRep dglen;
	if(n==nnew){
		int error;
		dglen=MGLBRep(dtp, gbydg, t, error);
		if(error) return tlist;
	}else{
		dglen=MGLBRep(dtp,gbydg,knew);
	}
			//std::cout<<(*this)<<dglen;//////////////
	double tolold=MGTolerance::wc_zero();
	double tol=deriv_length()*tolold;
	MGTolerance::set_wc_zero(tol);
	MGCParam_list temp=dglen.isect_1D(0.);//std::cout<<"********"<<temp<<std::endl;
	MGTolerance::set_wc_zero(tolold);
	MGCParam_list::Citerator j;
	for(j=temp.begin(); j!=temp.end(); j++){
		double tau;
		if(perp_guess(1.,0.,point,*j,tau))
			tlist.append(tau);
	}

	return tlist;
}

//Compute all the perpendicular points of this curve and the second one.
//That is, if f(s) and g(t) are the points of the two curves f and g,
//then obtains points where the following conditions are satisfied:
//  fs*(f-g)=0.    gt*(g-f)=0.
//Here fs and gt are 1st derivatives at s and t of f and g.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGLBRep::perps(const MGCurve& crv2)const{
	 return perps_by_split_to_C1(crv2);
}
MGPosition_list MGLBRep::perps(const MGStraight& crv2)const{
	 return perps_by_split_to_C1(crv2);
}
MGPosition_list MGLBRep::perps(const MGRLBRep& crv2) const{
	 return perps_by_split_to_C1(crv2);
}
MGPosition_list MGLBRep::perps(const MGEllipse& crv2)const{
	 return perps_by_split_to_C1(crv2);
}
MGPosition_list MGLBRep::perps(const MGLBRep& crv2)const{
	return perps_by_split_to_C1(crv2);
}
MGPosition_list MGLBRep::perps(const MGSurfCurve& crv2)const{
	return perps_by_split_to_C1(crv2);
}
MGPosition_list MGLBRep::perps(const MGBSumCurve& crv2)const{
	return perps_by_split_to_C1(crv2);
}

//perps_by_split_to_C1 is a private function for perps, computes perpendicular
//points of LBRep(this) with crv.
MGPosition_list MGLBRep::perps_by_split_to_C1(const MGCurve& crv2)const{
	size_t k1=order();
	if(k1==2)
		return perps_order2(crv2);

	//When this order is not 2.
	MGPosition_list list;
	size_t index, nb=bdim(), multi_found;
	size_t start=k1; int orderm1=k1-1;
	const MGKnotVector& t1=knot_vector();
	do{	//Compute intersections by dividing to parts of continuity>=C1.
		multi_found=t1.locate_multi(start,orderm1,index);
											//Locate C0 continuity point.
		if(start==k1 && index==nb)
			list.append(*this, crv2, C1perps(crv2));	//Use original.
		else{
			MGLBRep spanC1(t1[start-1],t1[index],*this);
			list.append(*this, crv2, spanC1.C1perps(crv2));
		}
		start=index+multi_found;
	}while(index<nb);
	return list;
}

//Perpendicular points with C1 conitnuity LBRep lbC1.
//MGPosition P in the MGPosition_list contains this and crv's parameter
//as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list MGLBRep::perps_withC1LB(
   const MGLBRep& lbC1
)const{
	MGPosition_list list;
	size_t k=order();
	if(k==2)
		return perps_order2(lbC1);

	//When this order is not 2.
	size_t index, nb=bdim(), multi_found;
	size_t start=k; int orderm1=k-1;
	const MGKnotVector& t1=knot_vector();
	do{	//Compute intersections by dividing to parts of continuity>=C1.
		multi_found=t1.locate_multi(start,orderm1,index);
											//Locate C0 continuity point.
		if(start==k && index==nb)
			list.append(*this, lbC1, perpendiculars(lbC1));	//Use original.
		else{
			MGLBRep spanC1(t1[start-1],t1[index],*this);
			list.append(*this, lbC1, spanC1.perpendiculars(lbC1));
		}
		start=index+multi_found;
	}while(index<nb);
	return list;
}

//Perpendicular points with SurfCurve
//whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list MGLBRep::perps_with_noCompoSC(const MGSurfCurve& curve2)const{
	return perps_by_split_to_C1(curve2);
}

//Peeps of this LBRep that does not include C0 continuity points in it.
MGPosition_list MGLBRep::C1perps(const MGCurve& curve2) const{
	MGPosition_list list=curve2.perps_withC1LB(*this);
	return MGPosition_list(list,1,0);	
}

//perps_order2 is a private function for perps, computes perpendiculars of
//LBRep(this) of order 2, i.e. polyline B-Rep, with crv.
MGPosition_list MGLBRep::perps_order2(const MGCurve& crv)const{
	MGPosition_list list;
	for(size_t i=0; i<bdim()-1; i++){
		MGStraight sl(coef(i+1),coef(i));	//i-th line segment.
		double clen=sl.param_e();
		if(!MGMZero(clen)){
			MGPosition_list ls=sl.perps(crv); //isect with i-th line segment.
			while(ls.entries()){
				MGPosition P=ls.removeFirst();
				double t1=knot(i+1);
				t1=t1+(P(0)/clen)*(knot(i+2)-t1);
				list.append(*this,crv,MGPosition(t1,P(1)));
			}
		}
	}
	return list;
}

//Check if the line B-rep is planar.
//Funtion's return value is;
// 0: Not planar, nor a point, nor straight line.
// 1: B-Rep is a point.		2: B-Rep is a straight line.
// 3: B-Rep is planar.
int MGLBRep::planar(
	MGPlane& plane			//When Brep is not straight line nor a point,
							// plane is returned.
							//Even when not planar, plane nearest is returned.
	, MGStraight& line		//When Brep is a line, line is returned.
	, MGPosition& center	//Center of the B-Rep is always returned.
)const{
	return line_bcoef().planar(plane,line,center);
}

/*MGPosition_list MGLBRep::perps(
	const MGLBRep& lb2		//The second curve
)const{
	size_t k1=order(), k2=lb2.order();
	if(k1==2){
		if(k2==2){
	//When both B-Rep are order 2,
			//compute perpendiculars using relation of 2 straight lines.
			size_t i,j,lb2span=lb2.bdim()-1;
			MGStraight* sl2array=new MGStraight[lb2span];
			for(j=0; j<lb2span; j++){	//Construct lb2's line segments.
				sl2array[j]=MGStraight(lb2.coef(j+1),lb2.coef(j));
			}
			double t1,t2,len1,len2; MGCCisect is;
			MGPosition_list list;
			for(i=0; i<bdim()-1; i++){
				MGStraight sl1(coef(i+1),coef(i));	//j-th line segment of this.
				if(!MGMZero(len1=sl1.param_e())){
					for(j=0; j<lb2span; j++){
						if(!MGMZero(len2=sl2array[j].param_e())){
							MGPosition_list ll=sl1.perps(sl2array[j]);
							if(ll.entries()){
								MGPosition pp=ll.removeFirst();
								t1=knot(i+1);
								t1=t1+(pp(0)/len1)*(knot(i+2)-t1);
								t2=lb2.knot(j+1);
								t2=t2+(pp(1)/len2)*(lb2.knot(j+2)-t2);
								list.append(*this, lb2, MGPosition(t1,t2));
							}
						}
					}
				}
			}
			delete[] sl2array;
			return list;
		}else
			return perps_order2(lb2);
	//When this order is 2 and lb2's is not.
	}else if(k2==2){
	//When this order is not 2 and lb2's is 2.
		MGPosition_list list=lb2.perps_order2(*this);
		return MGPosition_list(list,0,1);
	}
	//When both orders are not 2.
	MGPosition_list list;
	size_t index, nb=bdim(), multi_found;
	size_t start=k1; int orderm1=k1-1;
	const MGKnotVector& t1=knot_vector();
	do{	//Compute intersections by dividing to parts of continuity>=C1.
		multi_found=t1.locate_multi(start,orderm1,index);
											//Locate C0 continuity point.
		if(start==k1 && index==nb)
			list.append(*this, lb2, perps_1span(lb2));	//Use original.
		else{
			MGLBRep spanC1(t1[start-1],t1[index],*this);
			list.append(*this, lb2, spanC1.perps_1span(lb2));
		}
		start=index+multi_found;
	}while(index<nb);
	return list;
}*/
