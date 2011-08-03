/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Position_list.h"
#include "mg/Transf.h"
#include "mg/Straight.h"
#include "mg/RLBRep.h"
#include "mg/CompositeCurve.h"
#include "mg/SPointSeq.h"
#include "mg/SurfCurve.h"
#include "mg/RSBRep.h"
#include "mg/Plane.h"
#include "mg/CSisect_list.h"
#include "mg/Tolerance.h"
#include "Tl/TLData.h"
#include "Tl/TLInputParam.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// Implementation of MGSurface.
// MGSurface is an abstract class of 3D surface.
// Surface is represented using two parameter u and v.

//
///////////// Constructor ///////////
//
// 初期化なしでオブジェクトを作成する。
MGSurface::MGSurface():MGGeometry(){;}

//
///////////デストラクタ /////////////
//
MGSurface::~MGSurface(){;}

///////////Operator overload./////////////

//Generate arrow data of the tangent along u and v and the normal
//at the parameter value (u,v) of the surface.
//data[0] is the origin of the u-tangent arrow, data[1] is the top of the u-tangent arrow,
//data[2], [3] are two bottoms of u-tangent arrowhead.
//data[0], [4], [5], [6] are the points of v-tangent arrow.
//data[0], [7], [8], [9] are the points of v-tangent arrow.
void MGSurface::arrow(double u,double v, MGPosition data[10])const{
	arrow(box(),u,v,data);
}

//Generate arrow data, given box. The length of the arrows are defined from 
//box.len()
void MGSurface::arrow(const MGBox& box, double u,double v, MGPosition data[10])const{
	const double arrow_length=.06//of total length of the curve
				, head_length=.3;//of arrow_length.

	data[0]=eval(u,v);
	MGVector du=eval(u,v,1,0),dv=eval(u,v,0,1);
	double len=box.len()*arrow_length;
	MGUnit_vector ndu(du), ndv(dv);
	du=ndu*len; dv=ndv*len;
	MGPosition& P=data[0];
	one_arrow(P,du,dv,data[1],data[2],data[3]);
	one_arrow(P,dv,du,data[4],data[5],data[6]);
	MGUnit_vector N(du*dv);
	one_arrow(P,N*len,du,data[7],data[8],data[9]);
}

//Return box of the parameter space of the surface.
MGBox MGSurface::box_param() const{return parameter_range();}

//Obtain ceter coordinate of the geometry.
MGPosition MGSurface::center() const{
	double u=(param_s_u()+param_e_u())*0.5;
	double v=(param_s_v()+param_e_v())*0.5;
	return eval(u,v);
}

//Obtain ceter parameter value of the geometry.
MGPosition MGSurface::center_param() const{
	double u=(param_s_u()+param_e_u())*0.5;
	double v=(param_s_v()+param_e_v())*0.5;
	return MGPosition(u,v);
}

//Compute the closest point parameter value (u,v) of this surface
//from a point.
MGPosition MGSurface::closest(const MGPosition& point) const{
	size_t i;

	MGPosition_list list=perps(point);
	size_t pnum=perimeter_num();
	MGCurve* perimtr;
	for(i=0; i<pnum; i++){
		perimtr=perimeter_curve(i);
		list.append(perimeter_uv(i,perimtr->closest(point)));
		delete perimtr;
	}

	MGPosition uv1, uv;
	double dist1,dist;
	size_t n=list.entries();
	if(n){
		uv=list.removeFirst();
		dist=(point-eval(uv)).len();
		for(size_t i=1; i<n; i++){
			uv1=list.removeFirst();
			dist1=(point-eval(uv1)).len();
			if(dist1<dist) {uv=uv1; dist=dist1;}
		}
	}
	return uv;
}

//Compute the closest point on all the perimeters of the surface.
//The point is returned as the parameter value (u,v) of this surface.
MGPosition MGSurface::closest_on_perimeter(const MGPosition& point)const{
	MGPosition uv1, uv;
	double dist1,dist=-1.;
	size_t pnum=perimeter_num();
	MGCurve* perimtr;
	for(size_t i=0; i<pnum; i++){
		perimtr=perimeter_curve(i);
		uv1=perimeter_uv(i,perimtr->closest(point));
		dist1=(point-eval(uv1)).len();
		if(dist<0.){dist=dist1; uv=uv1;}
		else if(dist1<dist){dist=dist1; uv=uv1;}
		delete perimtr;
	}
	return uv;
}

//Compute the closest point on all the perimeters of the surface from
//the straight sl.
//The point is returned as the parameter value (u,v) of this surface.
MGPosition MGSurface::closest_on_perimeter(const MGStraight& sl)const{
	MGPosition uv1, uv;
	double dist=-1.;
	size_t pnum=perimeter_num();
	MGCurve* perii;
	for(size_t i=0; i<pnum; i++){
		perii=perimeter_curve(i);
		MGPosition st=perii->closest(sl);//st[0]:param of perii, st[1]:param of sl.
		double disti=(perii->eval(st[0]) - sl.eval(st[1])).len();
		if(dist<0.){dist=disti; uv=perimeter_uv(i,st[0]);}
		else if(disti<dist){dist=disti; uv=perimeter_uv(i,st[0]);}
		delete perii;
	}
	return uv;
}

//Compute the closest point on the surface from this point.
MGPosition MGPosition::closest(const MGSurface& surf) const{
	return surf.closest(*this);
}

//Compute surface curvatures:
//value[0]=K:Gaussian curvature=k1*k2, value[1]=H:Mean curvature=(k1+k2)/2,
//value[2]=k1:minimum curvature, and value[3]=k2=maximum curvature.
//N is the unit normal vector at position (u,v).
void MGSurface::curvatures(
	const MGPosition& uv, double value[4], MGUnit_vector& N) const{
	curvatures(uv[0], uv[1], value, N);
}
void MGSurface::curvatures(
	double u, double v, double value[4], MGUnit_vector& N) const{
	double Q[6];
	fundamentals(u,v,Q,N);
	double EG=Q[0]*Q[2], EN=Q[0]*Q[5], GL=Q[2]*Q[3];
	double EGmF2=EG-Q[1]*Q[1]; double EGmF2h=2.*EGmF2;
	double mb=EN+GL-2.*Q[1]*Q[4];	//EN+GL-2FM

	value[0]=(Q[3]*Q[5]-Q[4]*Q[4])/EGmF2;//(LN-M*M)/(EG-F*F)
	value[1]=mb/EGmF2h;					//(EN+GL-2FM)/(2*(EG-F*F))

	double r=(EN-GL); r*=r;
	double s=(Q[2]*Q[4]-Q[1]*Q[5])*(Q[0]*Q[4]-Q[1]*Q[3]);//(GM-FN)*(EM-FL)
	r+=(4.*s);if(r<0.) r=0.;
	r=sqrt(r);
	if(EGmF2>0.){value[3]=(mb+r)/EGmF2h; value[2]=(mb-r)/EGmF2h;}
	else        {value[2]=(mb+r)/EGmF2h; value[3]=(mb-r)/EGmF2h;}
}

//Compute direction unit vector of the geometry.
MGUnit_vector MGSurface::direction(const MGPosition& param)const{
	return normal(param);
}

//Compute if MGSurfCurve scurve(*this, param_curve) has the same direction
//to world_curve, assuming that scurve and world_curve are the same curve.
//Function's return value is:
//1: same direction, -1:oppositie direction.
int MGSurface::equal_direction(
	const MGCurve& param_curve,	//(u,v) parameter representation curve of this.
	const MGCurve& world_curve	//world representation curve.
)const{
	MGSurfCurve pline(*this, param_curve);
	double ts=pline.param_s(), tsb=world_curve.param_s();
	MGVector dir1=pline.eval(ts,1), dir2;
	MGVector Ps=pline.eval(ts), PBs=world_curve.eval(tsb);
	if(Ps==PBs){
		dir2=world_curve.eval(tsb,1);
	}else{
		double teb=world_curve.param_e();
		MGVector PBe=world_curve.eval(teb);
		MGVector dif1(Ps,PBs), dif2(Ps,PBe);
		if(dif1%dif1<dif2%dif2){
			dir2=world_curve.eval(tsb,1);
		}else{
			dir2=world_curve.eval(teb,1,1);
		}
	}
	if(dir1%dir2>0.)
		return 1;
	else
		return -1;
}

//Evaluate all the points (ui, vj) into spoint(i,j,.),
//where ui=utau(i) for 0<=i<utau.length() and vj=vtau(j) for 0<=j<vtau.length().
void MGSurface::eval_spoint(
	const MGNDDArray&	utau,		//u方向のデータポイント
	const MGNDDArray&	vtau,		//v方向のデータポイント
	MGSPointSeq&		spoint)const//evaluated data will be output to spoint.
{
	size_t nu=utau.length(), nv=vtau.length();
	spoint.resize(nu,nv,sdim());
	for(size_t i=0; i<nu; i++){
		double u=utau[i];
		for(size_t j=0; j<nv; j++){
			spoint.store_at(i,j,eval(u,vtau[j]));
		}
	}
}

//Evaluate right continuous surface data.
//Evaluate all positional data, 1st and 2nd derivatives.
void MGSurface::eval_all(
	double u, double v,	// Parameter value of the surface.
	MGPosition& f,			// Positional data.
	MGVector&   fu,			// df(u,v)/du
	MGVector&   fv,			// df/dv
	MGVector&   fuv,		// d**2f/(du*dv)
	MGVector&   fuu,		// d**2f/(du**2)
	MGVector&   fvv			// d**2f/(dv**2)
)const{
	f=eval(u,v);
	fu=eval(u,v,1,0);
	fv=eval(u,v,0,1);
	fuv=eval(u,v,1,1);
	fuu=eval(u,v,2,0);
	fvv=eval(u,v,0,2);
}

//Evaluate right continuous surface data.
//Evaluate all positional data, 1st and 2nd derivatives.
void MGSurface::eval_all(
	const MGPosition& uv,	// Parameter value of the surface.
	MGPosition& f,			// Positional data.
	MGVector&   fu,			// df(u,v)/du
	MGVector&   fv,			// df/dv
	MGVector&   fuv,		// d**2f/(du*dv)
	MGVector&   fuu,		// d**2f/(du**2)
	MGVector&   fvv			// d**2f/(dv**2)
	) const
{
	eval_all(uv.ref(0),uv.ref(1),f,fu,fv,fuv,fuu,fvv);
}

//evaluate gap between this surface's perimeter iperi and the given curve curve.
//function's return value is the maximum gap.
double MGSurface::eval_gap(
	const MGCurve& curve, // (I/ ) curve to evaluate.
	int            iperi,   // (I/ ) perimeter number, 0: vmin, 1: umax, 2: vmax, and 3: umin. 
	MGPosition&    uv       // ( /O) the parameter of this that had the largest gap.
)const{
	int num_peri=perimeter_num();
	if(!num_peri) return 0.;
	if(iperi>=num_peri || iperi<0) return 0.;

	uv.resize(2);
	MGCurve* peri = perimeter_curve(iperi);

	int id = -1;
	switch(iperi){
	case 0:	uv(1) = param_s_v(); id = 0;break;
	case 1:	uv(0) = param_e_u(); id = 1;break;
	case 2:	uv(1) = param_e_v(); id = 0;break;
	case 3:	uv(0) = param_s_u(); id = 1;break;
	};

	MGNDDArray tau; 
	tau.update_from_knot(curve.knot_vector());

	double t0 = peri->param_s();
	double t1 = peri->param_e();
	const double ratio = (t1-t0)/(curve.param_e()-tau[0]);

	MGPosition P = curve.start_point();
	double t = peri->closest(P);
	double lenmax = (peri->eval(t) - P).len();

	size_t ntaum1=tau.length()-1;
	for(size_t i = 0; i < ntaum1; i++){
		double tmp = (tau[i] + tau[i+1]) * .5;
		P = curve.eval(tmp);
		double tg = ratio * (tmp - tau[0]) + t0;
		if(!peri->perp_guess(t0, t1, P, tg, t)){
			//std::cerr << "** perp_guess (gap) error\n";
			t = peri->closest(P);
		}
		double len = (P - peri->eval(t)).len();
		if(len > lenmax){
			lenmax = len;
			uv(id) = t;
		}
		
		P = curve.eval(tau[i+1]);
		tg = ratio * (tau[i+1] - tau[0]) + t0;
		if(!peri->perp_guess(t0, t1, P, tg, t)){
			//std::cerr << "** perp_guess (gap) error\n";
			t = peri->closest(P);
		}
		len = (P - peri->eval(t)).len();
		if(len > lenmax){
			lenmax = len;
			uv(id) = t;
		}
	}
	P = curve.end_point();
	if(!peri->perp_guess(t0, t1, P, t1, t)){
		//std::cerr << "** perp_guess  (gap) error\n";
		t = peri->closest(P);
	}
	double len = (P - peri->eval(t)).len();
	if(len > lenmax){
		lenmax = len;
		uv(id) = t;
	}
	delete peri;
	return lenmax;
}

//evaluate gap between this surface's perimeters and the given curve curve.
//evaluation is performed for the perimeter i and curve[i] for 0<=i<=4.
//function's return value is the maximum gap.
double MGSurface::eval_gap(
	const MGCurve* curve[4], // (I/ ) curves to evaluate.
	MGPosition&    uv       // ( /O) the parameter of this that had the largest gap.
)const{
	int num_peri=perimeter_num();
	double gap=-1.;
	MGPosition uv2;
	for(int i=0; i<num_peri; i++){
		double gap2=eval_gap(*curve[i],i,uv2);
		if(gap<0. || gap2>gap){
			gap=gap2; uv=uv2;
		}
	}
	return gap;
}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector MGSurface::evaluate(
	const MGPosition& t,	// Parameter value.
				//t's space dimension is geometry's manifold dimension.
	const size_t* nderiv	//Order of derivative of i-th parameter
				//in nderiv[i].
				//When nderiv=null, nderiv[i]=0 is assumed for all i.
)const{
	if(nderiv) return eval(t,nderiv[0],nderiv[1]);
	else       return eval(t);
}

//Compute 1st and 2nd fundamental quantities of the surface.
//In Q, 1st and 2nd fundamental quantities are returned as:
//Q[0]=E, Q[1]=F, Q[2]=G,
//Q[3]=L, Q[4]=M, Q[5]=N.
void MGSurface::fundamentals(
	const MGPosition&uv,	//Surface parameter value (u,v)
	double Q[6],
	MGUnit_vector& N)		//Normal vector at uv will be returned.
const{
	fundamentals(uv[0], uv[1], Q, N);
}
void MGSurface::fundamentals(
	double u, double v,		//Surface parameter value (u,v)
	double Q[6],
	MGUnit_vector& N)		//Normal vector at (u,v) will be returned.
const{
	MGPosition f;	// Positional data.
	MGVector  fu;	// df(u,v)/du
	MGVector  fv;	// df/dv
	MGVector  fuv;	// d**2f/(du*dv)
	MGVector  fuu;	// d**2f/(du**2)
	MGVector  fvv;	// d**2f/(dv**2)
	eval_all(u,v,f,fu,fv,fuv,fuu,fvv);
	Q[0]=fu%fu; Q[1]=fu%fv; Q[2]=fv%fv;
	N=unit_normal(u,v);
	Q[3]=fuu%N; Q[4]=fuv%N; Q[5]=fvv%N;
}

//Compare two parameter values. If p1 is less than p2, return true.
//Comparison is done after prjected to i-th perimeter of the surface.
bool MGSurface::less_than(
	size_t i,	//perimeter number.
	const MGPosition& p1,
	const MGPosition& p2) const{return true;}

//Transform the coordinates of boundary of this geometry so that
//new coordinate of boundary is the same coordinate as the new one of
//this geometry after negate() of this geometry is done.
//That is, boundary coordinates are of parameter world of this geometry.
void MGSurface::negate_transform(MGGeometry& boundary) const{
	MGCurve* pline=dynamic_cast<MGCurve*>(&boundary);
	if(pline) pline->coordinate_exchange(0,1);
}

double mgDistance_flat(
	const MGVector& P0, 
	const MGVector& Pm, 
	const MGVector& P1
){
	MGVector Qm=(P0+P1)*.5;
	MGVector dif=Pm-Qm;
	return dif%dif;
}

//Compute the difference of min and max of the three doubles a1, a2, and a3.
double MGCL::Max3(double a1, double a2, double a3){
	double min, max;
	if(a1>=a2){
		if(a1>=a3){ max=a1; if(a2>=a3) min=a3; else min=a2;}
		else{ max=a3; min=a2;}
	}else if(a2>=a3){ max=a2; if(a1>=a3) min=a3; else min=a1;}
	else{
		min=a1; max=a3;
	}
	return max-min;
}

//compute sample point of the surface to get the approximate plane
//of the surface within the parameter range (u0,v0) to (u1, v1).
void MGSurface::compute_sample_point(
	double u0,
	double u1,
	double v0,
	double v1,
	MGPosition Pn[9],	//9 sample points will be output.
	MGPosition& center,	//center of the sample points will be output.
	MGUnit_vector& N,	//average normal of Nn[] will be output. 
	MGVector* Nn_in		//9 normals of the surface will be output.
)const{
	MGVector NnLocal[9];
	double um=(u0+u1)*0.5, vm=(v0+v1)*0.5;

	MGVector* Nn;
	if(Nn_in) Nn=Nn_in;
	else Nn=NnLocal;

	size_t i;	//id of Pn[].
	Pn[0]=eval(u0,v0); Nn[0]=normal(u0,v0);
	Pn[1]=eval(um,v0); Nn[1]=normal(um,v0);
	Pn[2]=eval(u1,v0); Nn[2]=normal(u1,v0);
	Pn[3]=eval(u0,vm); Nn[3]=normal(u0,vm);
	Pn[4]=eval(um,vm); Nn[4]=normal(um,vm);
	Pn[5]=eval(u1,vm); Nn[5]=normal(u1,vm);
	Pn[6]=eval(u0,v1); Nn[6]=normal(u0,v1);
	Pn[7]=eval(um,v1); Nn[7]=normal(um,v1);
	Pn[8]=eval(u1,v1); Nn[8]=normal(u1,v1);
	center=Pn[0]; MGVector VN=Nn[0]; 
	for(i=1; i<9; i++){center+=Pn[i]; VN+=Nn[i];}
	center/=9.;
	N=VN;
}

//Test if the surface is flat or not within the parameter value rectangle of uvbox.
//Function's return value is:
//	true: if the surface is flat
//  false: if the surface is not falt.
//When this is not falt, the direction that indicates which direction the surface
//should be divided will be output.
//***** the flatness is tested only approximately. This is for exclusive use of
//planar().
bool MGSurface::flat(
	const MGBox& uvbox,
	double tol,		//Tolerance allowed to regart flat
					//(Allowed distance from a plane).
	int& direction,	//   1: u-direction is more non flat.
					//   0: v-direction is more non flat.
	MGPosition& P,	//Position of the flat plane will be output.
	MGUnit_vector& N//Normal of the flat plane will be output.
)const{
	const MGInterval& urng=uvbox[0];
	double u0=urng[0].value(), u1=urng[1].value();
	const MGInterval& vrng=uvbox[1];
	double v0=vrng[0].value(), v1=vrng[1].value();

	MGPosition Pn[9];
	compute_sample_point(u0,u1,v0,v1,Pn,P,N);

	double x, d=P%N;
	double dist[9];
	bool is_flat=true;
	for(int i=0; i<9; i++){
		x=dist[i]=d-Pn[i]%N;
		if(x<0.) x=-x;
		if(x>tol) is_flat=false;;
	}

	double um=(u0+u1)*0.5, vm=(v0+v1)*0.5;
	MGVector dfdu=eval(um,vm,1,0), dfdv=eval(um,vm,0,1);
	double ulen=dfdu.len()*(u1-u0), vlen=dfdv.len()*(v1-v0);
	if(ulen*5.<vlen) direction=0;
	else if(ulen>vlen*5.) direction=1;
	else{
		double udif=MGCL::Max3(dist[0],dist[1],dist[2])
					+MGCL::Max3(dist[3],dist[4],dist[5])
					+MGCL::Max3(dist[6],dist[7],dist[8]);
		double vdif=MGCL::Max3(dist[0],dist[3],dist[6])
					+MGCL::Max3(dist[1],dist[4],dist[7])
					+MGCL::Max3(dist[2],dist[5],dist[8]);
		double error=MGTolerance::wc_zero_sqr()*15.;

		double difmax;
		if(udif>=vdif){
			difmax=udif;
			direction=1;
		}else{
			difmax=vdif;
			direction=0;
		}
		if(difmax<=error){
			if(ulen>vlen) direction=1;
			else direction=0;
		}
	}
	return is_flat;
}

//Given world curve wcrv on this face, get the parameter space representation pcrv.
//Function's return value is pcrv, which is newed one. Must be deleted.
MGCurve* MGSurface::get_parameterCurve(const MGCurve& wcrv)const{
	MGPvector<MGCurve> uvs, worlds;
	int n=project(wcrv,uvs,worlds);//n=number of curves projected.
	if(n<=0){
		MGPosition uv0, uv1;
		on(wcrv.start_point(), uv0);
		on(wcrv.end_point(), uv1);
		return new MGStraight(uv1,uv0);
	}else if(n==1){
		return uvs.release(0);
	}

	MGCompositeCurve* pcrv=new MGCompositeCurve(uvs.release(0));
	for(int i=1; i<n; i++){
		MGCurve& crvi=*(uvs[i]);
		MGPosition uv0=pcrv->end_point(), uv1=crvi.start_point();
		MGPosition Pe=eval(uv0);
		MGPosition Ps=eval(uv1);
		if((Pe-Ps).len()>MGTolerance::wc_zero())
			pcrv->connect_to_end(new MGStraight(uv1,uv0));
		pcrv->connect_to_end(uvs.release(i));
	}
	return pcrv;
}

//Compute the approximate plane in the parameter range from (u0, v0) to (u1,v1)
//The plane's origin is center point of the plane when the surface is mapped
//onto the plane. The uderiv of the plane is the direction from the point(u0, v0) to (u1,v0).
void MGSurface::get_approximate_plane(
	double u0,double u1,//u range from u0 to u1.
	double v0,double v1,//v range from v0 to v1.
	MGPlane& plane,		//The plane will be output.
	double* width,	//The width and the height of the plane that include all the data
	double* height	//for the surface point to map onto the plane will be output
)const{
	MGPosition Pn[9];
	MGPosition center;
	MGUnit_vector N;
	compute_sample_point(u0,u1,v0,v1,Pn,center,N);
	
	size_t i;	//id of Pn[].
	MGPlane plane2(N,center);
	MGVector xaxis(plane2.eval(plane2.uv(Pn[2])),plane2.eval(plane2.uv(Pn[0])));
	MGVector yaxis=N*xaxis;
	plane=MGPlane(xaxis,yaxis,center);
	if(!width) return;

	MGBox uvrange;
	for(i=0; i<9; i++){
		uvrange.expand(plane.uv(Pn[i]));
	}
	double umin=fabs(uvrange[0].low_point()), umax=fabs(uvrange[0].high_point());
	if(umin<umax){
		*width=umax*2.;
	}else{
		*width=umin*2.;
	}
	double vmin=fabs(uvrange[1].low_point()), vmax=fabs(uvrange[1].high_point());
	if(vmin<vmax){
		*height=vmax*2.;
	}else{
		*height=vmin*2.;
	}
}

//Compute the approximate plane in the parameter range from Pn, Nn, center, and N.
//when the surface is within surface_tol and angle from the plane.
//The plane's origin is center point of the plane when the surface is mapped
//onto the plane.
//the uderiv is the direction from the point(u0, v0) to (u1,v0).
//Function's return value is true when the surface is within the tolerance surface_tol,
//and the surface normals are within angle from the plane's normal
//plane, width, and height are valid only when function's return value is true.
bool test_to_get_approximate_plane(
	const MGPosition Pn[9],
	const MGVector Nn[9],
	const MGPosition& center,
	const MGUnit_vector& N,
	double surface_tol,	//tolerance allowed for the deviation from the plane to the surface.
	double angle,		//angle allowed for the normal of the plane and the normals of the
						//surface.
	MGPlane& plane,		//The plane will be output.
	double& width,	//The width and the height of the plane that include all the data
	double& height	//for the surface point to map onto the plane.
){
	size_t i;	//id of Pn[].
	double x, d=center%N;
	for(i=0; i<9; i++){
		x=d-Pn[i]%N;
		if(x<0.) x=-x;
		if(x>surface_tol) return false;
		double sang=N.sangle(Nn[i]);
		if(sang>angle) return false;
	}

	MGPlane plane2(N,center);
	MGUnit_vector xaxis=MGVector(plane2.eval(plane2.uv(Pn[2])),plane2.eval(plane2.uv(Pn[0])));
	MGUnit_vector yaxis=N*xaxis;
	plane=MGPlane(xaxis,yaxis,center);
	MGBox uvrange;
	for(i=0; i<9; i++){
		uvrange.expand(plane.uv(Pn[i]));
	}
	double umin=fabs(uvrange[0].low_point()), umax=fabs(uvrange[0].high_point());
	if(umin<umax){
		width=umax*2.;
	}else{
		width=umin*2.;
	}
	double vmin=fabs(uvrange[1].low_point()), vmax=fabs(uvrange[1].high_point());
	if(vmin<vmax){
		height=vmax*2.;
	}else{
		height=vmin*2.;
	}
	return true;
}

//Compute the approximate plane in the parameter range from (u0, v0) to (u1,v1)
//when the surface is within surface_tol and angle from the plane.
//The plane's origin is center point of the plane when the surface is mapped
//onto the plane.
//the uderiv is the direction from the point(u0, v0) to (u1,v0).
//Function's return value is true when the surface is within the tolerance surface_tol,
//and the surface normals are within angle from the plane's normal
//plane, width, and height are valid only when function's return value is true.bool MGSurface::test_and_get_approximate_plane(
bool MGSurface::test_and_get_approximate_plane(
	double u0,double u1,//u range from u0 to u1.
	double v0,double v1,//v range from v0 to v1.
	double surface_tol,	//tolerance allowed for the deviation from the plane to the surface.
	double angle,		//angle allowed for the normal of the plane and the normals of the
						//surface.
	MGPlane& plane,		//The plane will be output.
	double& width,	//The width and the height of the plane that include all the data
	double& height	//for the surface point to map onto the plane.
)const{
	MGPosition Pn[9];
	MGVector Nn[9];
	MGPosition center;
	MGUnit_vector N;
	compute_sample_point(u0,u1,v0,v1,Pn,center,N,Nn);
	return test_to_get_approximate_plane(Pn,Nn,center,N,surface_tol,angle,plane,width,height);
}

//Test if (u,v) is inside the face.
//Function's return value is:
//  0:outside the face.
//  1:unknown.
//  2:inside the face, not on a boundary.
//  <0:(u,v) is on an inner boundary, and abs(return code) is the loop id.
//  4:(u,v) is on the outer boundary.
//  >=10: (u,v) is on a perimeter, (10+perimeter number) will be returned.
int MGSurface::in_range_with_on(const MGPosition& uv)const{
	if(param_range()<<uv)
		return 0;

	size_t perim;
	double u=uv[0], v=uv[1];
	if(on_a_perimeter(u,v,perim))
		return 10+perim;
	else
		return 2;
}

//Obtain i-th inner_boundary curves(world coordinates representation)
//of the FSurface. Let the output of inner_boundary(i) be wcurves and
//of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
//to pcurves[j] one to one. Number of inner_boundary can be obtained
//by the function number_of_inner_boundary().
MGPvector<MGCurve> MGSurface::inner_boundary(size_t i)const{
	return MGPvector<MGCurve>();
}

//Obtain i-th inner_boundary curves(world coordinates representation)
//of the FSurface. Let the output of inner_boundary(i) be wcurves and
//of inner_boundary_param(i) be pcurves, then wcurves[j] corresponds
//to pcurves[j] one to one. Number of inner_boundary can be obtained
//by the function number_of_inner_boundary().
MGPvector<MGCurve> MGSurface::inner_boundary_param(size_t i)const{
	return MGPvector<MGCurve>();
}

//Compute normal vector(not unit) at uv.
MGVector MGSurface::normal(const MGPosition& uv) const
{	return normal(uv.ref(0), uv.ref(1));}

//Compute normal vector(not unit) at uv.
MGVector MGSurface::normal(double u,double v) const
{	return eval(u,v,1,0)*eval(u,v,0,1);}

//Compute unit normal vector at uv.
MGUnit_vector MGSurface::unit_normal(const MGPosition& uv) const
{	return unit_normal(uv.ref(0), uv.ref(1));}

#define MAX_LOOP_NUM 100
//Compute unit normal vector at uv.
MGUnit_vector MGSurface::unit_normal(double u,double v) const{
	MGVector fu(eval(u,v,1,0)), fv(eval(u,v,0,1));
	MGVector norml=fu*fv;
	double err_sqr=MGTolerance::mach_zero();
	if((norml%norml)>=err_sqr)
		return norml;

	double rzero=MGTolerance::rc_zero();
	size_t ndu=size_t(double(intersect_dnum_u())/rzero),
			ndv=size_t(double(intersect_dnum_v())/rzero);
	size_t n=ndu;
	if(n>ndv)
		n=ndv;
	double du=param_e_u()-u, du2=param_s_u()-u;
	if(du<(-du2))
		du=du2;
	du/=double(ndu);
	double dv=param_e_v()-v, dv2=param_s_v()-v;
	if(dv<(-dv2))
		dv=dv2;
	dv/=double(ndv);
	double fulen=fu%fu, fvlen=fv%fv;
	if(fulen<=err_sqr && fvlen>err_sqr)//case that fu is zero and fv is not.
		du=0.;
	else if(fulen>err_sqr && fvlen<=err_sqr)//case that fv is zero and fu is not.
		dv=0.;

	//case that fu is parallel to fv.

	n=MAX_LOOP_NUM;
	for(size_t i=0;i<n;i++){
		u+=du; v+=dv;
		fu=eval(u,v,1,0); fv=eval(u,v,0,1);
		norml=fu*fv;
		//We expect both of df/du and df/dv are not zero at(u,v).
		if((norml%norml)>=err_sqr)
			return norml;
	}

	return norml;
}

// 点がSpline上にあるか調べる。Spline上であれば，そのパラメータ値を，
// そうでなくても最近傍点のパラメータ値を返す。
//Function's return value is true if input point is on the surface,
// and  false if the point is not on the surface.
bool MGSurface::on(
	const MGPosition &P,	// 指定点
	MGPosition& uv		// Parameter value will be returned.                    
)const{
	size_t i;

	uv=MGPosition(param_s_u(), param_s_v());
	if(perp_one(P,uv)){if(MGAZero((P-eval(uv)).len())) return 1;}
	//Now if perp point exists, it is a candidate of the nearest point.

	//Compute nearest points with four perimeters.
	double t; bool onp=false; size_t pnum=perimeter_num();
	MGPosition* uvp=new MGPosition[pnum+1];
	uvp[0]=uv;
	MGCurve* perim;
	for(i=0;i<pnum;i++){
		perim=perimeter_curve(i); onp=perim->on(P,t);
		uvp[i+1]=perimeter_uv(i,t);
		delete perim;
		if(onp){uv=uvp[i+1]; break;}
	}
	if(!onp){
	//Compute the nearest point out of uvp[.].
		double dist=(P-eval(uvp[0])).len(), dist2;
		size_t id=0;
		for(i=1;i<=pnum;i++){
			dist2=(P-eval(uvp[i])).len();
			if(dist2<dist){id=i; dist=dist2;}
		}
		uv=uvp[id];
	}
	delete[] uvp;
	return onp;
}

// 点が曲面上にあるかを調べる。曲面上にあれば，そのパラメーター値を，
// なくても最近傍点のパラメータ値を返す。
// Function's return value is >0 if the point is on the surface,
// and 0 if the point is not on the curve.
bool MGPosition::on(	const MGSurface& surf,	// Surface pointer
		MGPosition& uv	//Parameter value of the nearest point on the surface.
)const{
	return surf.on(*this, uv);
}

//Test if input (u,v) is parameter value on a perimeter of the surface.
//If u or v is on a perimeter, it will be updated to the perimeter value.
bool MGSurface::on_a_perimeter(
	double& u, double& v,		//Surface parameter (u,v)
	size_t& perim_num//if function returns true,the perimete number will be output.
		//If function returns false, the nearest perimeter number will be output.
)const{
	bool on=1;
	double u0=param_s_u(), u1=param_e_u();
	double v0=param_s_v(), v1=param_e_v();
	double uspan=u1-u0, vspan=v1-v0;
	//uspan*=.1; vspan*=.1;
	double dist, distmin;
	perim_num=3;
	if(MGREqual_base(u,u0,uspan)){
		u=u0;
		if(MGREqual_base(v,v0,vspan)){
			perim_num=0; v=v0;
		}
	}else{
		distmin=(u-u0)/uspan; if(distmin<0.) distmin=-distmin;
		if(MGREqual_base(u,u1,uspan)){
			perim_num=1; u=u1;
			if(MGREqual_base(v,v1,vspan)){
				perim_num=2; v=v1;
			}
		}else{
			dist=(u1-u)/uspan; if(dist<0.) dist=-dist;
			if(dist<distmin){ distmin=dist; perim_num=1;}
			if(MGREqual_base(v,v0,vspan)){
				perim_num=0; v=v0;
			}else{
				dist=(v-v0)/vspan; if(dist<0.) dist=-dist;
				if(dist<distmin){ distmin=dist; perim_num=0;}
				if(MGREqual_base(v,v1,vspan)){
					perim_num=2; v=v1;
				}else{
					on=0;
					dist=(v1-v)/vspan; if(dist<0.) dist=-dist;
					if(dist<distmin) perim_num=2;
				}
			}
		}
	}
	return on;
}

//Test if input x is parameter value on a perimeter of the surface.
//If x is on a perimeter, x will be updated to the perimeter value.
//Function's return value is true if on a perimeter.
bool MGSurface::on_a_perimeter2(
	int is_u,	//specify if x is u or v value, is_u!=0(true) means u value.
	double& x,	//Surface parameter (u,v)
	size_t& perim_num//if function returns true,the perimete number will be output.
)const{
	perim_num=0;
	if(!perimeter_num())
		return false;

	double u0=param_s_u(), u1=param_e_u();
	double v0=param_s_v(), v1=param_e_v();
	if(is_u){
		double ulen=u1-u0;
		if(MGREqual_base(x, u0, ulen)){
			x=u0;
			perim_num=3;
			return true;
		}else if(MGREqual_base(x, u1, ulen)){
			x=u1;
			perim_num=1;
			return true;
		}
	}else{
		double vlen=v1-v0;
		if(MGREqual_base(x, v0, vlen)){
			x=v0;
			perim_num=0;
			return true;
		}else if(MGREqual_base(x, v1, vlen)){
			x=v1;
			perim_num=2;
			return true;
		}
	}
	return false;
}

//Test if input (u,v) is on the perimeter perim_num.
//If u or v is on a perimeter, true will be returned.
bool MGSurface::on_the_perimeter(
	size_t perim_num,	//a perimete number is input.
	double u, double v	//Surface parameter (u,v)
)const{
	double span, x,y;
	if(perim_num%2){//When v perimeter(u=min or u=max perimeter).
		span=knot_vector_u().param_span();
		if(perim_num==3) x=param_s_u();
		else x=param_e_u();
		y=u;
	}else{			//When u perimeter(v=min or v=max perimeter).
		span=knot_vector_v().param_span();
		if(perim_num==0) x=param_s_v();
		else x=param_e_v();
		y=v;
	}
	return MGREqual_base(y,x,span);
}

//Test the uvcurve is on a perimeter.
//If on a perimeter, true will be returned.
bool MGSurface::on_perimeter(
	const MGCurve& uvcurve,	//curve of surface parameter (u,v)
	size_t& perim_num //if function returned true, the perimete number will be output.
)const{
	if(!perimeter_num()) return false;
	double u0=param_s_u(), u1=param_e_u();
	double error=u1-u0; error*=5.;

	const MGBox& bx=uvcurve.box();
	double crvu0=bx[0].low_point(), crvu1=bx[0].high_point();
	if(MGREqual_base(crvu0,u0,error)&&MGREqual_base(crvu1,u0,error)){
		perim_num=3;
		return true;
	}
	if(MGREqual_base(crvu0,u1,error)&&MGREqual_base(crvu1,u1,error)){
		perim_num=1;
		return true;
	}

	double v0=param_s_v(), v1=param_e_v();
	error=v1-v0; error*=5.;
	double crvv0=bx[1].low_point(), crvv1=bx[1].high_point();
	if(MGREqual_base(crvv0,v0,error)&&MGREqual_base(crvv1,v0,error)){
		perim_num=0;
		return true;
	}
	if(MGREqual_base(crvv0,v1,error)&&MGREqual_base(crvv1,v1,error)){
		perim_num=2;
		return true;
	}

	return false;
}

//Obtain outer_boundary curves(world coordinates representation) of the FSurface.
//Let the output of outer_boundary() be wcurves and of outer_boundary_param()
//be pcurves, then wcurves[i] corresponds to pcurves[i] one to one.
MGPvector<MGCurve> MGSurface::outer_boundary()const{
	MGPvector<MGCurve> perims;
	size_t n=perimeter_num();
	for(size_t i=0; i<n; i++){
		MGCurve* peri=perimeter_curve(i);
		if(i>=2) peri->negate();
		perims.push_back(peri);
	}
	return perims;
}

//Obtain boundary curves(parameter space representation) of the FSurface.
//Let the output of boundary() be wcurves and of boundary_parameter()
//be pcurves, then wcurves[i] corresponds to  pcurves[i] one to one.
MGPvector<MGCurve> MGSurface::outer_boundary_param()const{
	MGPvector<MGCurve> crvs;
	size_t n=perimeter_num();
	if(!n) return crvs;

	MGBox uvrange=param_range();
	for(size_t i=0; i<n; i++){
		MGInterval& range=uvrange[i%2];
		double t0=range.low_point(),t1=range.high_point();
		MGPosition uv0=perimeter_uv(i,t0), uv1=perimeter_uv(i,t1);
		if(i<=1) crvs.push_back(new MGStraight(uv1,uv0));
		else crvs.push_back(new MGStraight(uv0,uv1));
	}
	return crvs;
}

double MGSurface::param_error() const{
	double er2=param_e_u()-param_s_u();
	double er3=param_e_v()-param_s_v();
	double error=er2*er2+er3*er3;
	double rcz=MGTolerance::rc_zero();
	error*=(rcz*rcz);
	return sqrt(error);
}
double MGSurface::param_error_u() const{
	double er=param_e_u()-param_s_u();
	return er*MGTolerance::rc_zero();
}
double MGSurface::param_error_v() const{
	double er=param_e_v()-param_s_v();
	return er*MGTolerance::rc_zero();
}

//Compuate square of parameter span length
//from (u.min, v.min) to (u.max, v.max).
double MGSurface::param_span() const{
	MGBox prange=param_range();
	double r0=prange.ref(0).length().value();
	double r1=prange.ref(1).length().value();
	return r0*r0+r1*r1;
}

// 自身の上の指定点を表すパラメータ値を返す。
// If input point is not on the surface, return the nearest point on the
// surface.
MGPosition MGSurface::param(
	const MGPosition& point		// 指定点
)const{
	MGPosition uv;
	on(point,uv);
	return uv;
}

// Return surface's parameter value of this point.
// If this point is not on the surface, return the nearest point's parameter
// value on the surface.
MGPosition MGPosition::param(const MGSurface& srf) const{
	return srf.param(*this);
}

//Let wcurve be a world curve rep that lies on this surface, and
//pcurve is parameter (u,v) expression of wcurve. That is,
//wcurve==MGSurfCurve pline(*this,pcurve). Then, param_of_pcurve() will obtain
//the parameter tp of pcurve that represent the same point as wcurve.eval(tw).
//Let S() is this surface, fp() is pcurve, and fw() is wcurve.
//Then S(fp(tp))=fw(tw).
double MGSurface::param_of_pcurve(
	double tw,			//point parameter of wcurve to get the pcurve parameter.
	const MGCurve& wcurve,//world curve that lies on this surface.
	const MGCurve& pcurve,//This surface's parameter rep of wcurve.
	const double* guess	//guess parameter value to compute tp. When guess=null,
						//param_of_pcurve will define the guess parameter.
)const{
	MGSurfCurve pline(*this,pcurve);
	double s0=pcurve.param_s(), s1=pcurve.param_e();
	double sguess;
	if(guess){
		sguess=*guess;
	}else{
		double t0=wcurve.param_s(), t1=wcurve.param_e();
		double tspan=t1-t0, sspan=s1-s0;
		double ratio=sspan/tspan;

		//get the adequate approximate parameter of this edge.
		if(pline.eval(s0,1)%wcurve.eval(t0,1)>0.)
			sguess=s0+(tw-t0)*ratio;//if pcurve and wcurve are the same direction.
		else
			sguess=s0+(t1-tw)*ratio;
	}

	double tp=sguess;
	MGPosition P=wcurve.eval(tw);// is the point of tw.
	int pobtained=pline.perp_guess(s0,s1,P,sguess,tp);
	if(pobtained){
		MGVector dif=P-pline.eval(tp);
		if(dif%dif<=(MGTolerance::wc_zero()*2.))
			return tp;
	}
	pline.on(P,tp);
	return tp;
}

//Compute parameter value of given point. Same as param.
// 自身の上の指定点を表すパラメータ値を返す。
// If input point is not on the geometry, return the nearest point on the
// geometry.
MGPosition MGSurface::parameter(
	const MGPosition& P	//Point(指定点)
	)const{ return param(P);}

//Obtain parameter curves.
//In the case of surFSurface, parameter curve is only one. However, in the case
//of FSurface,  number of parameter curves are more than one.
MGPvector<MGCurve> MGSurface::parameter_curves(
	int is_u,		//True(!=0) if x is u-value.(i.e. obtain u=const line)
	double x)const	//parameter value. u or v-value accordint to is_u.
{
	MGPvector<MGCurve> crvs;
	crvs.push_back(parameter_curve(is_u,x));
	return crvs;
}

//Return parameter range of the geometry(パラメータ範囲を返す)
MGBox MGSurface::parameter_range() const
{	return param_range();}

// Construct perimeter (u,v) parameter position.
MGPosition MGSurface::perimeter_uv(unsigned i,double t) const
// i is perimeter number:
// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
// t is perimeter parameter line's parameter value of u or v.
{
	assert(i<4);

	MGPosition uv(2);
	switch(i){
	case 0:	  uv(0)=t;uv(1)=param_s_v(); break;
	case 1:	  uv(1)=t;uv(0)=param_e_u(); break;
	case 2:	  uv(0)=t;uv(1)=param_e_v(); break;
	default:  uv(1)=t;uv(0)=param_s_u(); break;
	}
	return uv;
}

//Return all foots of the perpendicular straight lines from P.
MGPosition_list MGSurface::perps(
	const MGPosition& P				// Point of a space(指定点)
)const{
	const MGKnotVector& tu=knot_vector_u();
	const MGKnotVector& tv=knot_vector_v();
	size_t ku=tu.order(), kv=tv.order();
	size_t k=kv; if(k<ku) k=ku; k++;
	size_t ndiv=k/2; if(!ndiv) ndiv=1;
	//ndiv=1;///****
	double dndiv=double(ndiv);
	double erroru=tu.param_error()*2., errorv=tv.param_error()*2.;
	double errorM=MGTolerance::wc_zero()*-1.;
	MGPosition P0,P1,P2,P3;
	MGVector N0,N1,N2,N3;
	MGUnit_vector PN0,PN1,PN2,PN3;

	size_t kvm1=kv-1;
	size_t bdimu=bdim_u(), bdimv=bdim_v();
	MGPosition_list list;
	for(size_t i=ku-1;i<bdimu;i++){
		double tui=tu[i];
		double uspan=(tu[i+1]-tui)/dndiv;
		if(uspan<=erroru)
			continue;

		for(size_t j=kvm1;j<bdimv;j++){
			double tvj=tv[j];
			double vspan=(tv[j+1]-tvj)/dndiv;
			if(vspan<=errorv)
				continue;

			double u=tui; 
			for(size_t ii=0;ii<ndiv;ii++){
				double un=u+uspan, u2=u+uspan*.5;
				for(size_t jj=0;jj<ndiv;jj++){
					double v=tvj+vspan*double(jj);
					double vn=v+vspan, v2=v+vspan*.5;

					P0=eval(u,v); P1=eval(un,v);
					N0=normal(u2,v);
					PN0=N0*(P1-P0);
					double d0=PN0%P-PN0%eval(u2,v);
					if(d0<errorM)
						continue;

					P2=eval(un,vn);
					N1=normal(un,v2);
					PN1=N1*(P2-P1);
					double d1=PN1%P-PN1%eval(un,v2);
					if(d1<errorM)
						continue;

					P3=eval(u,vn);
					N2=normal(u2,vn);
					PN2=N2*(P3-P2);
					double d2=PN2%P-PN2%eval(u2,vn);
					if(d2<errorM)
						continue;

					N3=normal(u,v2);
					PN3=N3*(P0-P3);
					double d3=PN3%P-PN3%eval(u,v2);
					if(d3<errorM)
						continue;

					MGPosition uvguess(u2,v2), uv;
					if(perp_guess(P,uvguess,uv))
							list.append(*this,uv);
				}
				u=un;
			}
		}
	}
	return list;
}

int MGSurface::perp_guess(
const MGPosition& uv0,
const MGPosition& uv1,//parameter range of this surface.
		//When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
		//no limit for this parameter range.
const MGCompositeCurve& crv,//MGCompositeCurve.
double t0, double t1,//parameter range of curve.
		//When t0>=t1, no limit for curve2 parameter range.
const MGPosition& tuvg,	//Guess parameter value of curve and this surface.
MGPosition& tuv	//perpendicular points' parameter values will be output.
//tuv(0): curve's parameter, (tuv(1),tuv(2)):this surface's parameter.
) const{
	if(!crv.number_of_curves()) return 0;
	double t=tuvg[0];
	size_t i=crv.find(t);
	const MGCurve& curvei=crv.curve(i);
	if(t0<t1){
		double tis=curvei.param_s(), tie=crv.curve(i).param_e();
		if(t0<tis) t0=tis;
		if(tie<t1) t1=tie;
	}
	return perp_guess_general(uv0,uv1,curvei,t0,t1,tuvg,tuv);
}

//Return the foot of the perpendicular straight line from P that is 
//nearest to point P. Computation is done from the guess parameter value.
//Function's return value is whether point is obtained(1) or not(0)
int MGSurface::perp_guess(
	const MGPosition& uv0,
	const MGPosition& uv1,		//parameter range of this surface.
					//When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
					//no limit for this parameter range.
	const MGPosition& P,		//Point(指定点)
	const MGPosition& uvguess,	// guess parameter value of surface
	MGPosition& uv				// Parameter value will be returned.                    
)const{
//Method: Let S(u,v) is the surface, then f=(S-P)**2 should be extremum. 
//If S(u,v) is the point where the shortest distance between P and S, the
//following conditions are satisfied:
//  df/du=2*Su*(S-P)=0.    df/dv=2*Sv*(S-P)=0.

// Starting with guess parameter (u,v), du and dv are
//  obtained by solving the two equations
//	E+dE(u,v)=0 and F+dF(u,d)=0, where E(u,v)=Su*(S-P), and F(u,v)=Sv*(S-P).
//	Then (u+du, v+dv) is the next guess parameter value.
//	Here dE(u,v) and dF(u,v) are total differentials of E and F.
// dE=(Suu*(S-P)+Su*Su)*du+(Suv*(S-P)+Su*Sv)*dv=-E
// dF=(Suv*(S-P)+Su*Sv)*du+(Svv*(S-P)+Sv*Sv)*dv=-F
// If we set A=Suu*(S-P)+Su*Su, B=Suv*(S-P)+Su*Sv, and C=Svv*(S-P)+Sv*Sv,
// du=(B*F-C*E)/(A*C-B*B)	dv=(B*E-A*F)/(A*C-B*B) .
//
	uv.resize(2);
	MGPosition S;
	MGVector Su,Sv,Suv,Suu,Svv,SP,dS;
	double du,dv;
	double error_sqr=MGTolerance::wc_zero_sqr(); error_sqr *=.25;

	double error_rel=MGTolerance::rc_zero();
	double u0=uv0(0), u1=uv1(0), v0=uv0(1), v1=uv1(1);
	if(u0>=u1) {u0=param_s_u(); u1=param_e_u();}
	if(v0>=v1) {v0=param_s_v(); v1=param_e_v();}
	double uspan=u1-u0, vspan=v1-v0;

	double A,B,C,E,F,ACmB2;
	int loop=0, ulow=0, uhigh=0, vlow=0, vhigh=0;

	double uold,usave, vold,vsave;
	double& u=uv(0); double& v=uv(1);
	uold=u=uvguess.ref(0), vold=v=uvguess.ref(1);
	while(loop++<16 && ulow<5 && uhigh<5 && vlow<5 && vhigh<5){
		eval_all(u,v,S,Su,Sv,Suv,Suu,Svv);
		SP=S-P; E=Su%SP; F=Sv%SP;
		A=Suu%SP+Su%Su; B=Suv%SP+Su%Sv; C=Svv%SP+Sv%Sv;
		ACmB2=A*C-B*B;
		if(MGMZero(ACmB2)) return 0;

		du=(B*F-C*E)/ACmB2; dv=(B*E-A*F)/ACmB2;
		dS=Su*du+Sv*dv;
		u+=du; v+=dv;
		if(dS%dS<=error_sqr){
			if(MGREqual_base(u0,u,uspan) || u<u0)
				u=u0;
			else if(MGREqual_base(u1,u,uspan) || u>u1)
				u=u1;
			if(MGREqual_base(v0,v,vspan) || v<v0)
				v=v0;
			else if(MGREqual_base(v1,v,vspan) || v>v1)
				v=v1;
			return true;
		}

		usave=u; vsave=v;
		if(u<u0){
			//if(uold<u0 && u<=uold) break; //If no convergence, return.
			ulow+=1; uhigh=0; u=u0;
		}
		else if(u>u1){
			//if(uold>u1 && u>=uold) break; //If no convergence, return.
			ulow=0; uhigh+=1; u=u1;
		}
		else {ulow=0; uhigh=0;}

		if(v<v0){
			//if(vold<v0 && v<=vold) break; //If no convergence, return.
			vlow+=1; vhigh=0; v=v0;
		}
		else if(v>v1){
			//if(vold>v1 && v>=vold) break; //If no convergence, return.
			vlow=0; vhigh+=1; v=v1;
		}
		else {vlow=0; vhigh=0;}

		uold=usave; vold=vsave;
	}	
	return 0;
}

//Compute perpendicular points of a curve and a surface,
//given guess starting paramter values.
int MGSurface::perp_guess(
	const MGPosition& uv0,
	const MGPosition& uv1,		//parameter range of this surface.
					//When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
					//no limit for this parameter range.
	const MGCurve& curve,		//curve.
	double t0, double t1,		//parameter range of curve.
					//When t0>=t1, no limit for curve2 parameter range.
	const MGPosition& tuvg,	//Guess parameter value of curve and this surface.
	MGPosition& tuv		//perpendicular points' parameter values
					//will be output.
	//tuv(0): curve's parameter, (tuv(1),tuv(2)):this surface's parameter.
)const{
	const MGCompositeCurve* ccrv=dynamic_cast<const MGCompositeCurve*>(&curve);
	if(ccrv) return perp_guess(uv0,uv1,*ccrv,t0,t1,tuvg,tuv);
	return perp_guess_general(uv0,uv1,curve,t0,t1,tuvg,tuv);
}

//Return the foot of the perpendicular straight line from P.
//Computation is done from the guess parameter value.
//Function's return value is whether point is obtained(true) or not(false).
bool MGSurface::perp_guess(
	const MGPosition& P,		//Point
	const MGPosition& uvguess,	// guess parameter value of the shell
	MGPosition& uv				// Parameter value will be returned.
)const{
	MGPosition uv0(param_e_u(),param_e_v()), uv1(param_s_u(), param_s_v());
	return perp_guess(uv0,uv1,P,uvguess,uv)!=0;
}

//Compute perpendicular points of a curve and the FSurface,
//given guess starting paramter values.
//Function's return value is:
//   perp_guess=true if perpendicular points obtained,
//   perp_guess=false if perpendicular points not obtained,
bool MGSurface::perp_guess(
	const MGCurve& curve,	//curve.
	const MGPosition& uvguess,	//Guess parameter value of the FSurface.
	double tguess,			//Guess parameter value of the curve.
	MGPosition& uv,			//perpendicular point's parameter values of the shell
	double& t				//will be output.
)const{
	MGPosition uv0(param_e_u(),param_e_v()), uv1(param_s_u(), param_s_v());
	double t0=curve.param_e(), t1=curve.param_s();
	MGPosition tuvg(tguess,uvguess[0],uvguess[1]);
	MGPosition tuv;
	int obtained=perp_guess_general(uv0,uv1,curve,t0,t1,tuvg,tuv);
	uv.resize(2); uv(0)=tuv[1]; uv(1)=tuv[2];
	t=tuv[0];
	return obtained!=0;
}

//Compute perpendicular points of a curve and a surface,
//given guess starting paramter values.
int MGSurface::perp_guess_general(
	const MGPosition& uv0,
	const MGPosition& uv1,		//parameter range of this surface.
					//When uv0(0)>=uv1(0) or uv0(1)>=uv1(1),
					//no limit for this parameter range.
	const MGCurve& curve,		//curve.
	double t0, double t1,		//parameter range of curve.
					//When t0>=t1, no limit for curve2 parameter range.
	const MGPosition& tuvg,	//Guess parameter value of curve and this surface.
	MGPosition& tuv	//perpendicular points' parameter values will be output.
	//tuv(0): curve's parameter, (tuv(1),tuv(2)):this surface's parameter.
)const{
//Function's return value is:
//   perp_guess=true if perpendicular points obtained,
//   perp_guess=false if perpendicular points not obtained,

//****Method****
//Let f(t) and g(u,v) are curve and surface, then h=(f-g)**2 should be
//minimized. 
//If f(t) and g(u,v) are the points where the shortest(or longest)
//distance between f and g, then the following conditions are satisfied:
//  dh/dt=2*ft*(f-g)=0.    dh/du=2*g1*(g-f)=0.    dh/dv=2*g2*(g-f)=0.
//Here ft=df/dt, g1=dg/du, g2=dg/dv.

// Starting with guess parameter (t,(u,v)), dt and (du,dv) are
//  obtained by solving the three equations
//	E+dE(t,u,v)=0 ,  F+dF(t,u,v)=0, and G+dG(t,u,v)=0, 
//  where E(t,u,v)=ft%(f-g) , F(t,u,v)=g1%(g-f), G(t,u,v)=g2%(g-f).
//	Then (t+dt,(u+du,v+dv)) is the next guess parameter value.
//	Here dE(t,u,v), dF(t,u,v), and dG(t,u,v) are total differentials of
//  E, F, and G.
//  dE= (ftt%(f-g)+ft%ft)*dt-ft%g1*du-ft%g2*dv=-E
//  dF=-ft%g1*dt+(g1%g1+(g-f)%g11)*du+(g1%g2+(g-f)%g12)*dv=-F
//  dG=-ft%g2*dt+(g1%g2+(g-f)%g12)*du+(g2%g2+(g-f)%g22)*dv=-G
//  (ft, g1, and g2 are once differential of f, g about u, ang g about v.
//  ftt, g11, g12, and g22 are twice differential of f and g.)
//  If we set 
//  A=ftt%(f-g)+ft%ft, B=-ft%g1, C=-ft%g2,
//  X=g1%g1-(f-g)%g11, Y=g1%g2-(f-g)%g12, Z=g2%g2-(f-g)%g22
//	L1=B*B-A*X, L2=B*C-A*Y,
//	L3=C*C-A*Z,
//  L4=C*X-B*Y, L5=C*Y-B*Z
//  M1=A*F-B*E, M2=A*G-C*E, M3=B*G-C*F,
//  du, dv, and dt are obtained as;
//  du=(M1*L3-L2*M2)/(L1*L3-L2*L2)	dv=(L1*M2-M1*L2)/(L1*L3-L2*L2)
//  dt=(-E-B*du-C*dv)/A.
//
	double error_sqr=(MGTolerance::wc_zero_sqr())*.25;
	//.25 is multiplied to enforce more strictness of line intersection.
	//Tolerance is made half(.25=.5*.5).
	double error_rel=MGTolerance::rc_zero();
	double terror=(t1-t0)*error_rel;
	double t0e=t0-terror, t1e=t1+terror;
	double u0=uv0(0), u1=uv1(0), v0=uv0(1), v1=uv1(1);
	if(u0>=u1) {u0=param_s_u(); u1=param_e_u();}
	if(v0>=v1) {v0=param_s_v(); v1=param_e_v();}
	double uerror=(u1-u0)*error_rel;
	double u0e=u0-uerror, u1e=u1+uerror;
	double verror=(v1-v0)*error_rel;
	double v0e=v0-verror, v1e=v1+verror;

	MGPosition f,g;
	MGVector ft,ftt,g1,g2,g11,g12,g22,fmg, df,dg1,dg2;

	double A,B,C,E,F,G,X,Y,Z;
	double A2,B2,C2;
	double L1,L2,L3,L4,L5,L1322,L1524,L2534, M1,M2,M3;
	double AL1322,AL1524,AL2534;
	int loop=0, tlow=0,thigh=0, ulow=0,uhigh=0, vlow=0, vhigh=0;
	double dfdf,dgdg;

	double dt,told,tsave, du,uold,usave, dv,vold,vsave;
	double& t=tuv(0); double& u=tuv(1); double& v=tuv(2);
	told=t=tuvg(0), uold=u=tuvg(1), vold=v=tuvg(2);
	while(loop++<16 &&
		tlow<5 && thigh<5 && ulow<5 && uhigh<5 && vlow<5 && vhigh<5){
		eval_all(u,v,g,g1,g2,g12,g11,g22);
		curve.eval_all(t,f,ft,ftt);
		fmg=f-g;
		E=ft%fmg; F=-(g1%fmg); G=-(g2%fmg);
		A=ftt%fmg+ft%ft; B=-(ft%g1); C=-(ft%g2);
		A2=A*A; B2=B*B; C2=C*C;
		X=g1%g1-fmg%g11; Y=g1%g2-fmg%g12; Z=g2%g2-fmg%g22;
		L1=B2-A*X; L2=B*C-A*Y;
		L3=C2-A*Z;
		L4=C*X-B*Y; L5=C*Y-B*Z;
		M1=A*F-B*E; M2=A*G-C*E; M3=B*G-C*F;
		L1322=L1*L3-L2*L2; L1524=L1*L5-L2*L4; L2534=L2*L5-L3*L4;
		//cout<<"f="<<f<<" g="<<g<<" ft="<<ft<<endl;//////////

		//Compute du, dv.
		AL1322=L1322; if(L1322<0.) AL1322=-L1322;
		AL1524=L1524; if(L1524<0.) AL1524=-L1524;
		AL2534=L2534; if(L2534<0.) AL2534=-L2534;
		if(AL1322>=AL1524){
			if(AL1322>=AL2534){
				if(MGMZero(L1322)) return 0;
				du=(M1*L3-L2*M2)/L1322; dv=(L1*M2-M1*L2)/L1322;
			}else{
				if(MGMZero(L2534)) return 0;
				du=(M2*L5-L3*M3)/L2534; dv=(L2*M3-M2*L4)/L2534;
			}
		}else{
			if(AL1524>=AL2534){
				if(MGMZero(L1524)) return 0;
				du=(M1*L5-L2*M3)/L1524; dv=(L1*M3-M1*L4)/L1524;
			}else{
				if(MGMZero(L2534)) return 0;
				du=(M2*L5-L3*M3)/L2534; dv=(L2*M3-M2*L4)/L2534;
			}
		}
		//Compute dt.
		if(A2>=B2){
			if(A2>=C2){
				if(MGMZero(A)) return 0;
				dt=(-E-B*du-C*dv)/A;
			}else{
				if(MGMZero(C)) return 0;
				dt=(-G-Y*du-Z*dv)/C;
			}
		}else{
			if(B2>=C2){
				if(MGMZero(B)) return 0;
				dt=(-F-X*du-Y*dv)/B;
			}else{
				if(MGMZero(C)) return 0;
				dt=(-G-Y*du-Z*dv)/C;
			}
		}

		df=ft*dt; dfdf=df%df;
		dg1=g1*du; dg2=g2*dv; dgdg=dg1%dg1+dg2%dg2;
		t+=dt; u+=du; v+=dv;	// Update t,(u,v).
		if(dfdf<=error_sqr && dgdg<=error_sqr){
			int found=curve.in_range(t) && in_range(u,v);
			if(found){
				if(t0<t1) found=(t0e<=t && t<=t1e);
				if(found){
					if(u0<u1) found=(u0e<=u && u<=u1e);
					if(found) if(v0<v1) found=(v0e<=v && v<=v1e);
				}
			}
			return found;
		}
		
		tsave=t; usave=u; vsave=v;
		if(t<t0){
			//if(told<=t0 && t<=told) return 0; //If no convergence, return.
			tlow+=1; thigh=0; t=t0;
		}else if(t>t1){
			//if(told>=t1 && t>=told) return 0; //If no convergence, return.
			tlow=0; thigh+=1; t=t1;
		}else {tlow=0; thigh=0;}

		if(u<u0){
			//if(uold<=u0 && u<=uold) break; //If no convergence, return.
			ulow+=1; uhigh=0; u=u0;
		}else if(u>u1){
			//if(uold>=u1 && u>=uold) break; //If no convergence, return.
			ulow=0; uhigh+=1; u=u1;
		}else {ulow=0; uhigh=0;}

		if(v<v0){
			//if(vold<=v0 && v<=vold) break; //If no convergence, return.
			vlow+=1; vhigh=0; v=v0;
		}else if(v>v1){
			//if(vold>=v1 && v>=vold) break; //If no convergence, return.
			vlow=0; vhigh+=1; v=v1;
		}else {vlow=0; vhigh=0;}

		told=tsave; uold=usave; vold=vsave;
	}

	return 0;
}

//指定点から最も近い、垂線の足とパラメータ値を返す。
//Return the foot of the perpendicular straight line from p that is 
//nearest to point p.
// Function's return value is whether point is obtained(1) or not(0)
int MGSurface::perp_point(
	const MGPosition& p,// 指定点(point)
	MGPosition& uv,		//Parameter value of the surface will be returned.
	const MGPosition* uvguess	// guess parameter value of surface
)const{
	MGPosition uv0(1.,1.), uv1(0.,0.);
	if(uvguess) return perp_guess(uv0,uv1,p,*uvguess,uv);
	else        return perp_one(p,uv);
}

//Compute all foot points of the perpendicular line from this point to
//a surface.
// ポイントから与曲面へ下ろした垂線の足の，曲面のパラメータ値を
// すべて求める。
MGPosition_list MGPosition::perps(
	const MGSurface& srf	//Surface
) const{
	return srf.perps(*this);
}

//Compute the parameter value of the closest point from the straight to
//this object.
//sl is the eye projection line whose direction is from yon to hither, and if
//sl had multiple intersection points, The closest point to the eye will be selected.
MGPosition MGSurface::pick_closest(const MGStraight& sl)const{
	MGCSisect_list ises=isectSl(sl);
	MGPosition uv;
	size_t n=ises.size();
	if(n){
		MGCSisect_list::CSiterator i=ises.begin(), ie=ises.end();
		double t=(*i).param_curve();
		uv=(*i).param_surface();
		for(i++; i!=ie; i++){
			double t2=(*i).param_curve();
			if(t2>t){
				t=t2; uv=(*i).param_surface();
			}
		}
	}else{
		uv=closest_on_perimeter(sl);
	}
	return uv;
}

// 入力パラメータをパラメータ範囲でまるめて返却する。
MGPosition MGSurface::range(const MGPosition& uv) const{
	double u=uv.ref(0);
	if(u<param_s_u()) u=param_s_u(); if(u>param_e_u()) u=param_e_u();
	double v=uv.ref(1);
	if(v<param_s_v()) v=param_s_v(); if(v>param_e_v()) v=param_e_v();
	return MGPosition(u,v);
}

//ノット削除関数(B表現曲線のみ)
//トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
//removal knot. line_zero tolerance is used.
void MGSurface::remove_knot(){;}

// 指定点をとおり指定方向ベクトルを持つ直線の回りを指定角度
// 回転させて自身の曲面とする。
MGSurface& MGSurface::rotate_self(
	const MGVector& v,
	double d,
	const MGPosition& p
	) {
	MGTransf t(3); t.set_rotate_3D(v, d, p);// 指定された変換を作成する。
	*this *= t;	// 自身を変換する。
	return *this;
}

//Creates a surface of revolution.
//Parameterization of the surface is:
//	u=const parameter line generates given curve(when u=0.).
//  v=const parameter line generates a circle whose center is axis.
std::auto_ptr<MGSurface> MGCL::create_revolved_surface(
	const MGCurve& curve,   // profile curve
	const MGStraight& axis, // revolution axis
	double angle            // revolution angle
	){
	assert(!axis.direction().is_zero_vector());
	//assert(!MGZero_angle(angle));
	
	if(curve.identify_type() == MGRLBREP_TID){
		const MGRLBRep& wiper = dynamic_cast<const MGRLBRep&>(curve);
		return std::auto_ptr<MGSurface>(new MGRSBRep(wiper, axis, angle));
	}
	std::auto_ptr<MGRLBRep> wiper(MGCL::convert_to_rational(curve));
	return std::auto_ptr<MGSurface>(new MGRSBRep(*wiper, axis, angle));
}
