/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Interval.h"
#include "mg/Box.h"
#include "mg/Transf.h"
#include "mg/Unit_vector.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/CCisect.h"
#include "mg/defint.h"
#include "mg/nlbit.h"
#include "mg/Plane.h"
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

// コンストラクタ
// 初期化なしで楕円を生成する。
MGEllipse::MGEllipse(): MGCurve(), m_gprange(0),m_knotV(0) {;}

//Copy constructor.
//If to_radian is true,
//the parameter range will be changed to radian values.
MGEllipse::MGEllipse(const MGEllipse& el, bool to_radian)
:MGCurve(el), m_center(el.m_center), m_m(el.m_m), m_n(el.m_n)
,m_circle(el.m_circle),m_normal(el.m_normal),m_r(el.m_r),m_knotV(0){
	m_prange[0]=el.m_prange[0]; m_prange[1]=el.m_prange[1];
	if(el.m_gprange && !to_radian){
		m_gprange=new double[3];
		m_gprange[0]=el.m_gprange[0];
		m_gprange[1]=el.m_gprange[1];
		m_gprange[2]=el.m_gprange[2];
	}else m_gprange=0;
}

// Change space dimension and ordering of axis.
MGEllipse::MGEllipse(size_t dim // new space dimension
	, const MGEllipse& ellip	// original ellipse data
	, size_t start1				// start axis position of new ellipse
	, size_t start2)			// start axis position of original ellipse
:MGCurve(ellip)
, m_center(dim,ellip.m_center,start1,start2)
, m_m(dim,ellip.m_m,start1,start2)
, m_n(dim,ellip.m_n,start1,start2),m_knotV(0){
	update_mark();
	m_prange[0]=ellip.m_prange[0]; m_prange[1]=ellip.m_prange[1];
	if(ellip.m_gprange){
		m_gprange=new double[3];
		m_gprange[0]=ellip.m_gprange[0];
		m_gprange[1]=ellip.m_gprange[1];
		m_gprange[2]=ellip.m_gprange[2];
	}else m_gprange=0;
	set_normal_r_c();
}

// 中心、２つのベクトル、パラメータ範囲から楕円を生成する。
MGEllipse::MGEllipse (
	const MGPosition & p,		// 中心点              
	const MGVector & v1,		// 第１ベクタ          
	const MGVector & v2,		// 第２ベクタ          
	const MGInterval & intrvl	// 角度の範囲(ラジアン)
):MGCurve(),m_center(p),m_m(v1), m_n(v2),m_gprange(0),m_knotV(0){
	set_normal_r_c();
	if(intrvl.finite())
		set_param(intrvl.low_point(), intrvl.high_point());
	else
		set_param(0.,mgDBLPAI);
}

//  矩形に内接する楕円を生成する
//  plane    楕円が乗る平面に平行な平面
//  corner1  楕円が内接する矩形のコーナーの座標
//  corner2  corner1 の対頂角の座標のヒント
// 
//  もし corner2 が corner1 を通る plane に平行な平面に
//  乗っていない場合、corner2 は平行な平面に面直投影され、
//  それが計算に用いられる。
MGEllipse::MGEllipse(
	const MGPlane&    plane, 
	const MGPosition& corner1, 
	MGPosition&       corner2
):MGCurve(),m_gprange(0),m_knotV(0){
	// 1. plane を法線方向に持ち上げる
	MGPlane pl(plane);
	pl += (corner1 - pl.center()).project(pl.normal());

	// 2. corner2 を今の平面に面直投影する
	corner2 = pl.eval(pl.closest(corner2));

	// 3. MGEllipse コンストラクタを用いて楕円を生成する
	MGPosition org = (corner1 + corner2) * .5;
	MGVector diag  = corner2 - corner1;
	MGVector u = diag.project(pl.u_deriv()) * .5;
	if(u.is_zero_vector()){
		return;
	}
	MGVector v = diag.project(pl.v_deriv()) * .5;
	if(v.is_zero_vector()){
		return;
	}
	// 4. operator= でメンバをセットする
	*this = MGEllipse(org, u, v, MGInterval(0., mgDBLPAI));
}

//A whole circle form center, radius, and the normal.
// 中心、法線ベクトル、半径を指定して真円を作成する。
MGEllipse::MGEllipse(
	const MGPosition & p,		// 中心点  
	double d,					// 半径    
	const MGVector & v			// 法線    
):MGCurve(),m_center(p),m_normal(v),m_r(fabs(d)),m_circle(1)
,m_gprange(0),m_knotV(0){
	m_prange[0]=0.0; m_prange[1]=mgDBLPAI;
	if(d<0.) m_normal = -m_normal;
	MGUnit_vector m,n;
	m_normal.orthonormal(m_normal, m, n);
	m_m =m_r*m; m_n =m_r*n;
}

// 中心、法線、始点と角度を指定して円弧を作成する。
MGEllipse::MGEllipse(
	const MGPosition& center,		// 中心点           
	const MGPosition& start,		// 始点             
	double d,						// Angle in radian. 
	const MGVector& v				// 法線             
):MGCurve(),m_center(center),m_normal(v),m_m(start,center),m_circle(1)
,m_gprange(0),m_knotV(0){
	if(d<0.0) {m_normal=-m_normal; d=-d;}
	if(d >= mgDBLPAI){
		m_prange[0]=0.0; m_prange[1]=mgDBLPAI;
	}else{
		m_prange[0]=0.0; m_prange[1]=d;
	}
	m_r = m_m.len();			  // 半径の設定

	// 法線と長軸から短軸を作成する。
	if (m_m.parallel(m_normal)){
		MGUnit_vector m, n;
		m_normal.orthonormal(m_normal, m, n);
		m_m =m_r*m; m_n =m_r*n;
	}else{
		m_n=m_r*((m_normal*m_m).normalize());
		m_m=m_r*((m_n*m_normal).normalize());
		size_t dim=sdim();
		m_m=MGVector(dim,m_m);
		m_n=MGVector(dim,m_n);
	}
}

// １点Pと２つのベクトルからなる２直線に接する指定半径の
// 円弧を作成する。
// 第一ベクトルとの接点を始点、他方を終点とし、
// d>0 のときこれらで分割された円弧の指定点 P側を
// d<0 のときPと反対側を作成する。
MGEllipse::MGEllipse (
	const MGPosition& p,	// 直線の始点          
	const MGVector& v1,		// 直線の方向ベクトル１
	const MGVector& v2,		// 直線の方向ベクトル２
	double d				// 半径                
):MGCurve(),m_r(fabs(d)),m_circle(1), m_gprange(0),m_knotV(0){
	MGUnit_vector m(v1); MGUnit_vector n(v2);
	MGUnit_vector v3 = (m+n)/2.0; // 2つの直線の２等分線の単位方向ベクトル
	double ang=m.angle(n);        // v1,v2のなす角度を求める。

	m_normal=n*m;
	if(m_normal.parallel(m)) {
		MGUnit_vector mm,N;
		m_normal.orthonormal(m_normal, mm, N);
		m_m =m_r*mm; m_n = m_r*N;
		m_center=p+v3;
	}
	else{
		double l1=v1.len();
		if(MGRZero(l1) || m.parallel(n)){
			m_center=p+v3;
			m_m=m-m_center;
		}
		else{
			double ang2=ang/2.0;
			double k = m_r / sin(ang2);
			m_center = p + k*v3;			// 中心点を求める
			m_m = (k*cos(ang2))*m - k*v3;	// 長軸を求める。
		}
		m_n = m_r * ((m_normal * m_m).normalize());	// 短軸を求める。
		m_n=MGVector(sdim(), m_n);
	}
	double theta=mgPAI-ang;
	if(d<0.){ m_normal=-m_normal; m_n=-m_n; theta=mgPAI+ang; }
	m_prange[0]=0.0; m_prange[1]=theta;	// パラメータ範囲の設定
}

//Arc from start point, through point, and end point.
// 始点、通過点、終点を指定して円弧を作成する。
MGEllipse::MGEllipse (
	const MGPosition& start,	//Start point:始点
	const MGPosition& through,	//Through point:通過点
	const MGPosition& end,	//End point:終点
	bool whole_circle		//true if the whole circle is to generate.
):MGCurve(),m_circle(1), m_gprange(0),m_knotV(0){
	MGVector TS(start,through);
	MGVector TE(end,through);
	MGPosition st2=(start+through)/2.0;// 始点ー通過点の中点		
	MGPosition et2=(end+through)/2.0;// 通過点ー終点の中点
	MGPosition se2=(start+end)/2.0;

	double endAngle;
	if(start == through && start == end && through == end){
		// ３点が誤差範囲内で一致する。
		m_center = (start+through+end)/3.;
		m_r = 0.0;
		m_m = MGVector(1.0,0.0);
		m_n = MGVector(0.0,1.0);
		m_normal = mgZ_UVEC;
		endAngle=mgDBLPAI;// パラメータ範囲の設定
	}else if (start == through || through == end){
		// 通過点が始点か終点に一致するときstart-endが直径の半円を作成する。
		m_center = se2;
		m_m = MGVector(start,m_center);
		MGUnit_vector N;
		if(m_center.sdim()==2){
			m_normal = mgZ_UVEC;
			N=MGVector(2,m_normal*m_m);
		}else
			MGUnit_vector(m_m).orthonormal(m_m, N, m_normal);
		m_r = m_m.len();
		m_n = m_r*N;
		endAngle=mgPAI;	// パラメータ範囲の設定
	}else if (start == end){
		// 始点と終点が一致するときstart-throughが直径の真円を作成する。
		m_center = st2;	// 始終点と通過点の中点を中心にする。
		m_m = MGVector(start,m_center);// 中心から始点へのベクトルが長軸。
		MGUnit_vector N;
		if(m_center.sdim()==2){
			m_normal = mgZ_UVEC;
			N=MGVector(2,m_normal*m_m);
		}else
			MGUnit_vector(m_m).orthonormal(m_m, N, m_normal);
		m_r = m_m.len();
		m_n = m_r*N;
		endAngle=mgDBLPAI;	// パラメータ範囲の設定
	}else if (TS.parallel(TE)){
		// ３点が直線上にのるとき 始点、終点を直径とする。半円を作成する。
		// 始点終点の中点を中心にする。
		m_center = se2;
		m_m = MGVector(start,m_center);// 中心から始点へのベクトルを長軸とする。
		MGUnit_vector N;
		if(m_center.sdim()==2){
			m_normal = mgZ_UVEC;
			N=MGVector(2,m_normal*m_m);
		}
		else MGUnit_vector(m_m).orthonormal(m_m, N, m_normal);
		m_r = m_m.len();
		m_n = m_r*N;
		endAngle=mgPAI;	// パラメータ範囲の設定
	}else{
		// ３点がばらばらで直線上に乗らないとき中心点を求める。
		// 求める点は始点ー通過点、通過点ー終点２つの線分の
		// 垂直２等分線の交点なので２つの垂直２等分線を作成する。
		MGUnit_vector N=TE*TS;
		MGStraight s1(MGSTRAIGHT_UNLIMIT,TS*N,st2);//始点ー通過点の垂直２等分線
		MGStraight s2(MGSTRAIGHT_UNLIMIT,TE*N,et2);//通過点ー終点の垂直２等分線
		MGCCisect isect;
		s1.relation(s2,isect);
		const MGPosition& C = isect.point();//２直線の交点が中心点。
		copy_ellipse_data(MGEllipse(C,start,end,N,whole_circle));
		return;
	}
	m_prange[0]=0.0;// パラメータ範囲の設定
	if(whole_circle)
		m_prange[1]=mgDBLPAI;
	else
		m_prange[1]=endAngle;
	update_mark();
}

//An arc of radius r whose start point is start, and the end point is end.
//and that is normal to the vector N.
//The circle lies on the plane whose normal is N and that passes through
//start and end.
//radius r is able to have minus value, in which case the longer part of the
//arc out of the whole circle is constructed.
//The center of the circle C is:
//C=M+MGUnit_vector(sign(r)*N*(start-end))*sqrt(r*r-d*d),
//where M=(start+end)*.5, and d is the distance between start and M.
MGEllipse::MGEllipse(
  double r,		//radius
  const MGPosition& start,
  const MGPosition& end,
  const MGVector& N,
  bool whole_circle		//true if the whole circle is to generate.
):MGCurve(),m_gprange(0),m_knotV(0){
	MGPosition M((start+end)*.5);
	double d=start.distance(M);
	if(MGAZero(d))
		whole_circle=true;
	double sign=1.; if(r<0.) sign=-1.;
	MGUnit_vector N2(N*sign);
	MGUnit_vector SE(end-start), MtoC, SE2;
	N2.orthonormal(SE,SE2,MtoC);
	double ctom=r*r-d*d;if(ctom<0.) ctom=0.;
	MGPosition C=M+MtoC*sqrt(ctom);
	double ang=mgDBLPAI;
	if(!whole_circle)
		ang=C.angle(start,end,N);
	copy_ellipse_data(MGEllipse(C,start,ang,N));
}

//Arc from center, start point, end point.
MGEllipse::MGEllipse(
	const MGPosition& center,	//center of the circle.
	const MGPosition& start,	//Start point:始点
	const MGPosition& end,		//End point:終点
	const MGVector& N,			//Normal of the plane the circle lies on.
	bool whole_circle		//true if the whole circle is to generate.
):MGCurve(),m_gprange(0),m_knotV(0){
	double ang=mgDBLPAI;
	if(start!=end){
		if(!whole_circle)
			ang=center.angle(start,end,N);
	}
	copy_ellipse_data(MGEllipse(center,start,ang,N));
}

///An arc of radius r whose start point is start, and the end point is end.
///The circle lies on the plane that the three points start, end, and reference
///lie on.
///The center of the circle C is:
///C=M+MGUnit_vector(sign(r)*N*(end-start))*sqrt(r*r-d*d),
///where M=(start+end)*.5, d is the distance between start and M,
///and N=(reference-start)*(end-start).
///That is, if r>0., the position reference indicates on which side
///the C lies against the vector (end-start), C and reference lie
///in the same half plane that the vector (end-start) divides.
///If r<0., C and reference lie in the opposite side.
///When r>0., smaller part arc of the circle is selected,
///and when r<0., larger part of the circle is selected.
MGEllipse::MGEllipse(
  double r,		//radius
  const MGPosition& start,
  const MGPosition& end,
  const MGPosition& reference,
  bool whole_circle		//true if the whole circle is to generate.
):MGCurve(),m_gprange(0),m_knotV(0){
	MGVector N=(end-start)*(reference-start);
	copy_ellipse_data(MGEllipse(r,start,end,N,whole_circle));
}

//Replace this ellipse with the arc whose start point is start,
//whose end point is end, and whose tangent at the start point is dir_s.
MGEllipse::MGEllipse(
	const MGPosition& start,
	const MGPosition& end,
	const MGVector& dir_s,
	bool whole_circle		//true if the whole circle is to generate.
):MGCurve(),m_gprange(0),m_knotV(0){
	MGVector T=end-start;
	MGVector N=dir_s*T;//normal of the arc.
	MGUnit_vector C=N*dir_s;
	double r=T.len()/2./dir_s.sangle(T);//radius of the arc.
	MGPosition center=start+r*C;//cneter of the arc.
	double ang=mgDBLPAI;
	if(!whole_circle)
		ang=angle(center, start, end, N);
	copy_ellipse_data(MGEllipse(center,start,ang,N));
}

/////////Destructor/////////
MGEllipse::~MGEllipse(){
	if(m_gprange) delete[] m_gprange;
	if(m_knotV) delete m_knotV;
}

/////////////// メンバ関数  /////////////////

//Compute if m_normal is the same as x, y, or z-axis or not.
//Return value is:
//	 0, 1, 2: m_normal is the same as x, y, or z-axis, each.
//	-1, -2, -3: m_normal is not the same but nearest to x, y, or z-axis each.
int MGEllipse::axis() const{
	double sx,sy,sz;
	sz=m_normal.sangle(mgZ_UVEC);
	if(MGZero_angle(sz)) return 2;
	else{
		sy=m_normal.sangle(mgY_UVEC);
		if(MGZero_angle(sy)) return 1;
		else{
			sx=m_normal.sangle(mgX_UVEC);
			if(MGZero_angle(sx)) return 0;
		}
	}
	sx=fabs(sx); sy=fabs(sy); sz=fabs(sz);
	if(sx<=sy){
		if(sx<=sz) return -1;
		else       return -3;
	}
	else{
		if(sy>sz)  return -3;
		else       return -2;
	}
}

// 入力パラメータ値により楕円のparameter範囲を狭めた部分の楕円の周りを
// 囲む最小のボックスを返す。
MGBox MGEllipse::box_limitted( const MGInterval& intrvl ) const{
	MGEllipse temp(*this);		// Generate temporal ellipse to force new
	temp.limit(intrvl);			// parameter range of intrvl;
	return temp.box();
}

// 楕円全体の周りを囲む最小のボックスを返す。
MGBox* MGEllipse::compute_box() const{
	if(ellipse_type()==MGELLIPSE_EMPTY)
		return new MGBox(start_point(), end_point());

	size_t dim=sdim();
	MGBox* boxt=new MGBox(dim); double ang;
	for(size_t i=0; i<dim; i++){
		double a=m_m.ref(i), b=m_n.ref(i);
		double aba=fabs(a), abb=fabs(b);
		if(aba<abb){
			if(MGMZero(abb)) ang=0.;
			else ang=atan2(b,a);
		}else{
			if(MGMZero(aba)) ang=0.;
			else ang=atan(b/a);
		}                
		if(in_RelativeRange_of_radian(ang)) (*boxt)(i)=
			MGInterval((*boxt)(i),m_center.ref(i)+a*cos(ang)+b*sin(ang));
		ang+=mgPAI;
		if(in_RelativeRange_of_radian(ang)) (*boxt)(i)=
			MGInterval((*boxt)(i),m_center.ref(i)+a*cos(ang)+b*sin(ang));
	}
	MGBox boxt2=MGBox(*boxt,start_point());
	delete boxt;
	boxt=new MGBox(boxt2,end_point());

	return boxt;
}

//Changing this object's space dimension.
MGEllipse& MGEllipse::change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2) 		// Source order of this object.
{
	m_center=MGPosition(sdim,m_center,start1,start2);
	m_m=MGVector(sdim,m_m,start1,start2);
	m_n=MGVector(sdim,m_n,start1,start2);
	set_normal_r_c();
	update_mark();
	return *this;
}

//Chane the parameter range of this ellipse to radian.
void MGEllipse::change_param_to_radian(){
	if(m_gprange){
		delete m_gprange; m_gprange=0;
		if(m_knotV){delete m_knotV; m_knotV=0;}
	}
}

//Change parameter range, be able to change the direction by providing
//t1 greater than t2.
void MGEllipse::change_range(
	double t1,		//Parameter value for the start of original. 
	double t2		//Parameter value for the end of original. 
){
	if(t1>t2){
		negate();
		double save=t1; t1=t2; t2=save;
	}

	double prlen=m_prange[1]-m_prange[0];
	if(param_range_is_in_radian()){
		if(MGREqual_base(t1,m_prange[0],prlen) && 
			MGREqual_base(t2,m_prange[1],prlen))
			return;
	}else{
		double gprlen=m_gprange[1]-m_gprange[0];
		if(MGREqual_base(t1,m_gprange[0],gprlen) && 
			MGREqual_base(t2,m_gprange[1],gprlen))
			return;
	}

	if(m_knotV){delete m_knotV; m_knotV=0;}
	assert(!MGMZero(t2-t1));
	if(!m_gprange)
		m_gprange=new double[3];
	m_gprange[0]=t1; m_gprange[1]=t2;
	m_gprange[2]=prlen/(t2-t1);
}

//Exchange ordering of the coordinates.
//Exchange coordinates (i) and (j).
MGEllipse& MGEllipse::coordinate_exchange(size_t i, size_t j){
	assert(i<sdim() && j<sdim());
	m_center.swap(i,j);
	m_m.swap(i,j); m_n.swap(i,j);
	m_normal.swap(i,j);
	update_mark();
	return *this;
}

//Construct new curve object by copying to newed area.
//User must delete this copied object by "delete".
MGEllipse* MGEllipse::clone() const{return new MGEllipse(*this);}

//copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
//When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
//Otherwise,  the new curve will be a MGLBRep.
//Returned object must be deleted.
MGCurve* MGEllipse::copy_as_nurbs() const{
	MGCurve* nurbs=new MGRLBRep(*this);
	nurbs->change_range(param_s(), param_e());
	return nurbs;
}

//Construct new curve object by changing
//the original object's space dimension.
//User must delete this copied object by "delete".
MGEllipse* MGEllipse::copy_change_dimension(
	size_t sdim,		// new space dimension
	size_t start1, 		// Destination order of new line.
	size_t start2) 		// Source order of this line.
	const{
	return new MGEllipse(sdim,*this,start1,start2);
}

//copy all the ellipse specific data into this from ellipse2.
void MGEllipse::copy_ellipse_data(const MGEllipse& ellipse2){
	update_mark();
	m_center=ellipse2.m_center;
	m_normal=ellipse2.m_normal;
	m_m=ellipse2.m_m;
	m_n=ellipse2.m_n;
	m_r=ellipse2.m_r;
	m_prange[0]=ellipse2.m_prange[0];
	m_prange[1]=ellipse2.m_prange[1];
	m_circle=ellipse2.m_circle;
	if(ellipse2.m_gprange){
		m_gprange=new double[3];
		m_gprange[0]=ellipse2.m_gprange[0];
		m_gprange[1]=ellipse2.m_gprange[1];
		m_gprange[2]=ellipse2.m_gprange[2];
	}else m_gprange=0;
	delete m_knotV; m_knotV=0;
}

//  焦点と通過点を与えて楕円を生成する
//  F0  焦点
//  F1  焦点
//  P   楕円上の点
// 
//  F0 と F1 は相異なる焦点でなければ失敗する。
//  P が2焦点を結ぶ直線上にあるとき失敗する。
bool MGEllipse::create_ellipse(
	const MGPosition& F0,	//focus 1
	const MGPosition& F1,	//focus 2
	const MGPosition& P		//a point on the ellipse that is not on the line
							//through F0 and F1.
){
	if(F0 == F1 || P.is_collinear(F0, F1)){
		return false;
	}
	// 楕円の中心は2焦点間の中点
	const MGPosition& origin = (F0 + F1) * .5;
	// 正の実数 a を 線分 F0-P と F1-P の長さの和の半分とする
	double a = (F0.distance(P) + F1.distance(P)) *.5;
	// 正の実数 L は 2焦点間の距離の半分
	double L = F0.distance(origin);

	// major axis を計算する
	const MGVector major = a * (F1 - origin) / L;
	assert(MGAEqual(major.len(), a));

	// major axis を回転することで minor axis の向きを計算する
	MGTransf trans;
	MGVector normal = MGPlane(F0, F1, P).normal();
	trans.set_rotate_3D(normal, mgHALFPAI, origin);
	MGPosition tmp(origin + major);
	tmp *= trans;

	MGVector minor = tmp - origin;
	// 長さを調整する
	minor /= minor.len();
	minor *= sqrt(fabs(a * a - L * L));

	// 楕円を生成して終了
	*this = MGEllipse(origin, major, minor, MGInterval(0., mgDBLPAI));
	return true;
}

//Compute curvilinear integral of the 1st two coordinates, i.e. 2D.
//This integral can be used to compute area sorounded by the curve.
//(線積分）を求める。
double MGEllipse::curvilinear_integral(double t1, double t2) const{
	t1=range(t1); t2=range(t2);
	t1=gp_to_radian(t1); t2=gp_to_radian(t2);
	double sval=sin(t2)-sin(t1);
	double cval=cos(t2)-cos(t1);
	double a1=m_m[1]*m_center[0]-m_m[0]*m_center[1];
	double a2=m_n[1]*m_center[0]-m_n[0]*m_center[1];
	double a3=m_m[0]*m_n[1]-m_m[1]*m_n[0];
	double val=a1*cval+a2*sval+a3*(t2-t1);
	if(!m_gprange) return val;
	else return val*m_gprange[2];
}

//Return ellipse type:楕円のタイプを返却する。
MGELLIPSE_TYPE MGEllipse::ellipse_type()const{
	if(MGREqual_base(m_prange[1],m_prange[0],mgDBLPAI)) return MGELLIPSE_EMPTY;
	if(MGREqual_base(m_prange[1]-mgDBLPAI,m_prange[0],mgDBLPAI))
		return MGELLIPSE_CLOSED;
	else return MGELLIPSE_SEGMENT;
}

// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector MGEllipse::eval(
	double t,			// Parameter value.
	size_t nderiv,		// Order of Derivative.
	int left			//Left continuous(left=true)
						//or right continuous(left=false).
)const{
	t=range(t);
	if(!m_gprange) return eval_in_radian2(t,nderiv);
	t=gp_to_radian(t);
	if(nderiv==0) return eval_in_radian2(t);
	else{
		double dpower;
		if(nderiv<=2){
			dpower=m_gprange[2];
			if(nderiv==2) dpower*=dpower;
		}else dpower=pow(m_gprange[2],double(nderiv));
		return eval_in_radian2(t,nderiv)*dpower;
	}
}

// 楕円上の与えられたパラメータ値における一次微分値を返却する。
MGVector MGEllipse::eval_deriv(double d) const{
	d=range(d);
	if(m_gprange){
		d=gp_to_radian(d);
		return (-sin(d)*m_m + cos(d)*m_n)*m_gprange[2];
	}else return (-sin(d)*m_m + cos(d)*m_n);
}

// パラメータを与え、位置、一次微分値、二次微分値を求める。
void MGEllipse::eval_all(
	double d,
	MGPosition & p,
	MGVector & v1,
	MGVector & v2
)const{
	d=range(d);
	d=gp_to_radian(d);
	double cosv=cos(d), sinv=sin(d);
	p = m_center + cosv*m_m + sinv*m_n;
	if(m_gprange){
		v1 =(-sinv*m_m + cosv*m_n); v1*=m_gprange[2];
		v2 = -m_m * cosv -m_n * sinv; v2*=(m_gprange[2]*m_gprange[2]);
	}else{
		v1 = -sinv*m_m + cosv*m_n;
		v2 = -m_m * cosv -m_n * sinv;
	}
}

// 与えられたパラメータ値に相当する楕円上の点を返却する。
MGVector MGEllipse::eval_in_radian2(
	double t,		//Parameter value in radian.
	size_t nderiv	//Order of Derivative.
 )const{
	double cosv=cos(t), sinv=sin(t);
	if(nderiv==0) return m_center + cosv*m_m + sinv*m_n;
	else{
		MGVector v;
		div_t result=div(nderiv, 4);
		switch(result.rem){
		case 0: v=  cosv*m_m+sinv*m_n; break;
		case 1: v= -sinv*m_m+cosv*m_n; break;
		case 2: v= -cosv*m_m-sinv*m_n; break;
		case 3: v=  sinv*m_m-cosv*m_n; break;
		}
		return v;
	}
}

//Evaluate ellipse data.
//Input parameter t must be in radian.
MGVector MGEllipse::eval_in_radian(
		double t,		//Parameter value in radian.
		size_t nderiv	// Order of Derivative.
)const{
	double s=radian_to_gp(t);
	s=range(s);
	t=gp_to_radian(s);
	return eval_in_radian2(t,nderiv);
}

//Evaluate ellipse data.
//Input parameter degree must be in degree.
MGVector MGEllipse::eval_in_degree(
		double degree,	//Parameter value in degree.
		size_t nderiv	// Order of Derivative.
)const{
	double t=MGCL::degree_to_radian(degree);
	return eval_in_radian(t,nderiv);
}

//Access to i-th element of knot.
//i=0, 1 and returns start or end parameter value of the ellipse.
double MGEllipse::knot(size_t i) const{
	assert(i<=2);
	if(i<=1) return param_s();
	else return param_e();
}

//Returns the knot vector of the curve.
const MGKnotVector& MGEllipse::knot_vector() const{
	if(!m_knotV){
		m_knotV=new MGKnotVector(2,2);
		(*m_knotV)(0)=(*m_knotV)(1)=param_s();
		(*m_knotV)(2)=(*m_knotV)(3)=param_e();
	}
	return *m_knotV;
}
MGKnotVector& MGEllipse::knot_vector(){
	if(!m_knotV){
		m_knotV=new MGKnotVector(2,2);
		(*m_knotV)(0)=(*m_knotV)(1)=param_s();
		(*m_knotV)(2)=(*m_knotV)(3)=param_e();
	}
	return *m_knotV;
}

// 与えられたパラメータ間の曲線に沿った代数的距離を返却する。
// パラメータ値が昇順で与えられたとき正値、降順のとき負値を
// 返す。
double MGEllipse::length(double r1, double r2) const{
	double t1=gp_to_radian(r1), t2=gp_to_radian(r2);
	if(t1>t2){ double save=t1; t1=t2; t2=save;}
	MGEllipse temp(*this, true);//temp is ellipse of parameter range in radian.
	double len;
	if(t1<temp.m_prange[0]) t1=temp.m_prange[0];
	if(t2>temp.m_prange[1]) t2=temp.m_prange[1];
	if(m_circle) len=m_r*(t2-t1);
	else{
		double eps=MGTolerance::rc_zero(); eps*=mgDBLPAI*MGEllipse::m_r;
		MGCurveLengthDrive f(this);
		len=mgDefint(f,t1,t2, eps);
	}
	if(r1>r2) len=-len;
	return len;
}

// 楕円の全体の長さを返却する。
double MGEllipse::length() const{
	double t1=param_s(), t2=param_e();
	return length(t1,t2);
}

// パラメータで示される点t_inから指定距離lenはなれた点のパラメータ
// を返す。
double MGEllipse::length_param(double t_in, double len) const{
	double ts=gp_to_radian(range(t_in));

	if(MGAZero(m_r)) return radian_to_gp(ts);//When radius=0
	double d=RelativeRange_in_radian(ts+len/m_r);
	if(m_circle || ellipse_type()==MGELLIPSE_EMPTY) //When circle.
		return radian_to_gp(d);

	double a,b;
	if(m_m.orthogonal(m_n)){ a=m_m.len(); b=m_n.len();}
	else{
		MGMatrix mat; mat.set_xy_axis(m_m, m_n);
		MGVector axis1=m_m*mat, axis2=m_n*mat;
		MGEllipse etemp
			(MGPosition(0.,0.),axis1,axis2,MGInterval(0., mgDBLPAI));
		MGBox bx=etemp.box();
		a=bx(0).length().value()/2.; b=bx(1).length().value()/2.;
	}

	MGEllipse etemp(*this,true);//Ellipse in radian range.
	double t1,t2,tsave;
	if(ellipse_type()==MGELLIPSE_SEGMENT){
		if(len>=0.) tsave=m_prange[1];
		else 		tsave=m_prange[0];
		double total=etemp.length(ts,tsave);
		if(fabs(len)>=fabs(total)){//Solution is outside range;
			return radian_to_gp(tsave); 
		}
		if(MGMZero(a)) t1=m_prange[1];else t1=RelativeRange_in_radian(ts+len/a);
		if(MGMZero(b)) t2=m_prange[1];else t2=RelativeRange_in_radian(ts+len/b);
 	}else{
		double len_round=etemp.length();
		len=fmod(len, len_round);
		if(MGMZero(a)) t1=m_prange[1]; else t1=ts+len/a;
		if(MGMZero(b)) t2=m_prange[1]; else t2=ts+len/b;
		d=ts+len/m_r;
	}
	if(t1>t2){tsave=t1; t1=t2; t2=tsave;}
	double len1=length(ts,t1)-len, len2=length(ts,t2)-len
								 , lenm=length(ts,d)-len;
	if(len1*lenm<=0.) t2=d; else t1=d;

	//Now there exists a solution between t1 and t2.
	MGCurveLenParamDrive clen(&etemp,len,ts);
	double eps=MGTolerance::rc_zero(); eps*=2.*mgPAI*MGEllipse::m_r;
	int ier;
	double x=mgNlbit(clen, t1,t2, eps, 30, ier);
	return radian_to_gp(x); 
}

// 自身の楕円に指定された範囲のlimitを付加する。
MGEllipse& MGEllipse::limit(const MGInterval& intrvl){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	double t1=param_s(), t2;
	if(intrvl.empty()) t2=t1;
	else{//現在のパラメータ範囲との積をとる。
		t2=param_e();
		if(intrvl<t1) t2=t1;
		else if(t2<intrvl) t1=t2;
		else{
			if(t1<intrvl) t1=intrvl.low_point();
			if(intrvl<t2) t2=intrvl.high_point();
		}
	}
	if(m_gprange){
		double a1=gp_to_radian(t1), a2=gp_to_radian(t2);
		m_gprange[0]=t1; m_gprange[1]=t2;
		double t2mt1=t2-t1;
		if(MGMZero(t2mt1)) m_gprange[2]=1.;
		else m_gprange[2]=(a2-a1)/t2mt1; 
		m_prange[0]=a1; m_prange[1]=a2;
	}else{
		m_prange[0]=t1; m_prange[1]=t2;
	}
	update_mark();
	return *this;
}

// 楕円の方向を反転する。方向ベクトルを逆にする。範囲があるとき
// は始終点を入れ換える。
void MGEllipse::negate(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	m_normal = -m_normal;
	//Replace m_m and m_n.
	MGVector msave=m_m, nsave=m_n;
	double t01=m_prange[0]+m_prange[1];
	double ct01=cos(t01), st01=sin(t01);
	m_m=msave*ct01+nsave*st01;
	m_n=msave*st01-nsave*ct01;
}

//Obtain parameter value if this curve is negated by "negate()".
double MGEllipse::negate_param(double t)const{
	if(m_gprange) t=(m_gprange[1]+m_gprange[0])-t;
	else{
		t=(m_prange[0]+m_prange[1])-t;
	}
	return t;
}

//Normalize parameter value t to the nearest knot if their distance is
//within tolerance. For ellipse, the knots are start and end points.
double MGEllipse::param_normalize(double t) const{
	double plen=param_span(), tnew;
	tnew=param_s(); if(MGRZero2(tnew-t,plen)) return tnew;
	tnew=param_e(); if(MGRZero2(tnew-t,plen)) return tnew;
	return t;
}

//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by the calling program, or memory leaked.
MGEllipse* MGEllipse::part(double t1, double t2, int multiple) const{
	assert(t2-t1>param_error());
	MGEllipse* el=new MGEllipse(*this);
	el->limit(MGInterval(t1,t2));
	return el;
}

//Return space dimension
size_t MGEllipse::sdim() const{
	if(m_normal.is_null()) return 0;
    double a=m_normal.ref(1)*m_normal.ref(1)+m_normal.ref(0)*m_normal.ref(0);
//sqrt(a) is sin(angle) where angle is the angle between z-axis and m_normal.
	double error=MGTolerance::angle_zero(); error *=error;
	if(a<=error) return m_center.sdim();
	else return 3;
}

// 演算子の多重定義

//Assignment.
//When the leaf object of this and crv2 are not equal, this assignment
//does nothing.
MGEllipse& MGEllipse::operator=(const MGEllipse& el2){
	if(this==&el2)
		return *this;

	MGCurve::operator =(el2);
	m_center=el2.m_center;
	m_m=el2.m_m; m_n=el2.m_n; m_normal=el2.m_normal;
	m_circle=el2.m_circle; m_r=el2.m_r;
	m_prange[0]=el2.m_prange[0]; m_prange[1]=el2.m_prange[1];
	if(el2.m_gprange){
		if(!m_gprange) m_gprange=new double[3];
		m_gprange[0]=el2.m_gprange[0];
		m_gprange[1]=el2.m_gprange[1];
		m_gprange[2]=el2.m_gprange[2];
	}else{
		if(m_gprange) delete m_gprange;
		m_gprange=0;
	}
	if(m_knotV){delete m_knotV; m_knotV=0;}
	return *this;
}
MGEllipse& MGEllipse::operator=(const MGGel& gel2){
	const MGEllipse* gel2_is_this=dynamic_cast<const MGEllipse*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

// 楕円の平行移動を行いオブジェクトを生成する。
MGEllipse MGEllipse::operator+ ( const MGVector & v) const{
	MGEllipse e(*this);
	e += v;
	return e;
}
MGEllipse operator+ (const MGVector& v, const MGEllipse& e){
	return e+v;
}

// 楕円の平行移動を行い自身の楕円とする。
MGEllipse& MGEllipse::operator+= (const MGVector& v){
	m_center += v;
	if(m_box) *m_box+=v;
	return *this;
}

// 楕円の逆方向の平行移動を行いオブジェクトを生成する。
MGEllipse MGEllipse::operator- ( const MGVector & v) const {
	MGEllipse e(*this);
	e -= v;
	return e;
}

// 楕円の逆方向の平行移動を行い自身の楕円とする。
MGEllipse& MGEllipse::operator-= (const MGVector& v) {
	m_center -= v;
	if(m_box) *m_box-=v;
	return *this;
}

//Generate scaled ellipse.
// 与えられたスケールをかけ自身の楕円とする。
MGEllipse MGEllipse::operator* (double scale)const{
	MGEllipse el(*this);
	el*=scale;
	return el;
}

//Ellipseをスケーリングしてできるオブジェクトを生成する。
//Generates an ellipse by scaling.
MGEllipse operator* (double scale, const MGEllipse& el){
	return el*scale;
}

//Update the curve by multiplying scale.
// 与えられたスケールをかけ自身の楕円とする。
MGEllipse& MGEllipse::operator*= (double scale){
	m_m*=scale;
	m_n*=scale;
	if(scale>=0.) m_r*=scale; else m_r *= -scale;
	update_mark();
	return *this;
}

// 与えられたmatrixによる変換を行いオブジェクトを生成する。
MGEllipse MGEllipse::operator* (const MGMatrix& mat) const{
	MGEllipse e= *this;
	e *= mat;
	return e;
}

// 与えられた変換による変換を行い自身の楕円とする。
MGEllipse& MGEllipse::operator*= (const MGMatrix& mat){
	m_center *= mat;
	m_m *= mat;
	m_n *= mat;
	set_normal_r_c();
	update_mark();
	return *this;
}

// 与えられた変換によるトランスフォームを行いオブジェクトを生成する。
MGEllipse MGEllipse::operator* (const MGTransf& t) const {
	MGEllipse e = *this;
	e *= t;
	return e;
}

// 与えられた変換によるトランスフォームを行い自身の楕円とする。
MGEllipse& MGEllipse::operator*= (const MGTransf& tr) {
	m_center *= tr;
	m_m *= tr.affine();
	m_n *= tr.affine();
	set_normal_r_c();
	update_mark();
	return *this;
}

// 論理演算子の多重定義
// 自身の楕円とCurveが等しいかどうか比較し判定する。
bool MGEllipse::operator==(const MGEllipse& e)const{
	if(m_center != e.m_center)
		return 0;
	if(m_m != e.m_m)
		return 0;
	if(m_n != e.m_n)
		return 0;
	double span=m_prange[1]-m_prange[0];
	if(!MGREqual2(span, e.m_prange[1]-e.m_prange[0]))
		return 0;
	if(!MGREqual_base(m_prange[1], e.m_prange[1], span))
		return 0;
	return 1;
}
bool MGEllipse::operator<(const MGEllipse& gel2)const{
	return m_r<gel2.m_r;
}
bool MGEllipse::operator==(const MGGel& gel2)const{
	const MGEllipse* gel2_is_this=dynamic_cast<const MGEllipse*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGEllipse::operator<(const MGGel& gel2)const{
	const MGEllipse* gel2_is_this=dynamic_cast<const MGEllipse*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//Normalize parameter range intrvl(in radian)of the ellipse,
//and set the parameter range in m_prange.
void MGEllipse::set_param(
	double d1,	//Input paramter range in radian.
	double d2
){
	double t0=fmod(d1,mgDBLPAI);

	double a0, a1;
	if(d1>=d2){
		a0=a1=t0;
	}else{
		if(d2-d1>=(mgDBLPAI-MGTolerance::angle_zero())){
			double d1new =t0;
			if(d1new > 0.)
				d1new=d1new-mgDBLPAI;
			a0=d1new;
			a1=d1new+mgDBLPAI;
		}else{
			double d1new =t0; double diff=d1-d1new;
			double d2new=d2-diff;
			if(d2new>mgDBLPAI){
				d1new=d1new-mgDBLPAI; d2new=d2new-mgDBLPAI;
			}
			a0=d1new;
			a1=d2new;
        }
	}

	if(m_gprange){
		double t0=radian_to_gp(a0), t1=radian_to_gp(a1);
		m_gprange[0]=t0; m_gprange[1]=t1; m_gprange[2]=(a1-a0)/(t1-t0);
	}
	m_prange[0]=a0; m_prange[1]=a1;
	update_mark();
}

//Compute m_normal, m_r, and m_circle from m_m and m_n.
void MGEllipse::set_normal_r_c(){
	double a=m_m.len(), b=m_n.len();
	m_r=sqrt((a*a+b*b)/2.);
	MGUnit_vector m,n;
	if(a>=b){
		m=m_m; m.orthonormal(m_n,n,m_normal);
	}else{
		n=m_n; n.orthonormal(m_m,m,m_normal);
		m_normal=-m_normal;
	}
	m_circle=(m_m.orthogonal(m_n) && MGREqual2(a,b));
}

//Test if input angle is in the radian parameter range.
//When m_prange[0]<=angle and angle<=m_prangle[1], return true
//(tolerance included).
bool MGEllipse::in_radian_range(double angle) const{
	double error=MGTolerance::rc_zero()*(m_prange[1]-m_prange[0]);
	double t0=m_prange[0]-error, t1=m_prange[1]+error;
	if(angle<t0) return false;
	if(angle>t1) return false;
	return true;
}

//Compute if anglei is in the parameter range.
//If in, return 1(true), if not, return 0(false).
//Radian version. Input anglei must be of radian.
bool MGEllipse::in_RelativeRange_of_radian(double anglei) const{
	double error=MGTolerance::rc_zero()*(m_prange[1]-m_prange[0]);
	double t0=m_prange[0]-error, t1=m_prange[1]+error;

	double angle=fmod(anglei,mgDBLPAI);
	bool ret_val=(t0<=angle && angle<=t1);
	if(angle<0.0) angle+=mgDBLPAI;
	else		  angle-=mgDBLPAI;
	ret_val |=(t0<=angle && angle<=t1);
	return ret_val;
}

// 入力パラメータをパラメータ範囲に丸める。
//Input p and function's return value are both in radian value.
double MGEllipse::RelativeRange_in_radian(double p) const{
	double low=m_prange[0], high=m_prange[1];
	if(low<=p && p<=high) return p;
	p = fmod(p,mgDBLPAI);
	if(low<=p && p<=high) return p;

	double pp2pai=p+mgDBLPAI, pm2pai=p-mgDBLPAI;
	MGELLIPSE_TYPE et=ellipse_type();
	if(et==MGELLIPSE_CLOSED){
		if(p>=0.0) return pm2pai;
		else       return pp2pai;
	}else if(et==MGELLIPSE_EMPTY) return (low+high)*0.5;

	double diff1,diff2;
	if(p>=0.0){
		if(low<=pm2pai && pm2pai<=high) return pm2pai;
		diff1=low-pm2pai;
		diff2=p-high;
	}else{
		if(low<=pp2pai && pp2pai<=high) return pp2pai;
		diff1=low-p;
		diff2=pp2pai-high;
	}
	if(diff1<diff2) return low;
	else            return high;
}

//Compute how distant param p1 is from param p2 in radian.
//p1 and p2 must be the values in radian.
double MGEllipse::param_length(double p1, double p2) const{
	double q1=RelativeRange_in_radian(p1), q2=RelativeRange_in_radian(p2);
	double diff;
	diff=fabs(q1-q2);if(diff>mgPAI) diff=mgDBLPAI-diff;
	return diff;
}

// 楕円からパラメータのｌｉｍｉｔを取り除く。
MGCurve& MGEllipse::unlimit(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	double a0, a1;//New parameter range in radian.
	if(m_prange[0]>=0.){
		a0=0.; a1=mgDBLPAI;
	}else if(m_prange[1]<=0.){
		a0=-mgDBLPAI; a1=0.;
	}else if(m_prange[0]>=-mgPAI && m_prange[1]<=mgPAI){
		a0=-mgPAI; a1=mgPAI;
	}else{
		a0=m_prange[0]; a1=a0+mgDBLPAI;
	}
	set_param(a0,a1);
	return *this;
}

//Unlimit parameter range of the curve to the end point direction
//(終点方向にlimitをはずす)
MGCurve& MGEllipse::unlimit_end(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	double t0=m_prange[0];//Save the original start parameter.
	set_param(t0, t0+mgDBLPAI);
	return *this;
}

//Unlimit parameter range of the curve to the start point direction
//(始点方向にlimitをはずす)
MGCurve& MGEllipse::unlimit_start(){
	if(m_knotV){delete m_knotV; m_knotV=0;}
	double t1=m_prange[1];//Save the original end parameter.
	set_param(t1-mgDBLPAI, t1);
	return *this;
}
