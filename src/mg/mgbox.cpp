/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Vector.h"
#include "mg/Position.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Straight.h"
#include "mg/Plane.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGBox.cc
// Implementation of Class MGBox

//////////// 1. Constructor ////////////

//Copy constructor.
MGBox::MGBox(const MGBox& box2):m_sdim(0), m_range(0){
	get_area(box2.sdim());
	for(size_t i=0; i<m_sdim; i++) m_range[i]=box2.m_range[i];
}

//Construct 2D Box by providing each Interval.
MGBox::MGBox(const MGInterval& xspan, const MGInterval& yspan)
:m_sdim(2), m_range(m_rData){
	m_range[0]=xspan;
	m_range[1]=yspan;
}

//Construct 3D Box by providing each Interval.
MGBox::MGBox(
	const MGInterval& xspan, const MGInterval& yspan, const MGInterval& zspan
):m_sdim(3), m_range(m_rData){
	m_range[0]=xspan;
	m_range[1]=yspan;
	m_range[2]=zspan;
}

//中心点と各辺の幅を指定しBoxを生成する。
MGBox::MGBox(const MGPosition& center, double* size)
:m_sdim(0), m_range(0){
	size_t dim=center.sdim();
	get_area(dim);
	for(size_t i=0; i<dim; i++){
		double icentr = center(i);
		assert (size[i]>=0.);
		double hsize=size[i]*0.5;
		m_range[i]=MGInterval(icentr-hsize, icentr+hsize);
	}
}

//中心点と幅を指定しBoxを生成する（すべての辺にたいして同一の値）。
//Construct Box, given center point and a size of each coordinates.
//The size is applied to all the coordinates.
MGBox::MGBox(const MGPosition& center, double size)
:m_sdim(0), m_range(0){
	assert (size>=0.);
	size_t dim=center.sdim();
	get_area(dim);
	double hsize=size*0.5;
	for(size_t i=0; i<dim; i++){
		double icentr = center(i);
		m_range[i]=MGInterval(icentr-hsize, icentr+hsize);
	}
}

//２点からBoxを生成する。
MGBox::MGBox(const MGPosition &pt1, const MGPosition &pt2)
:m_sdim(0), m_range(0){
	size_t dim1=pt1.sdim();	size_t dim2=pt2.sdim();
	size_t dim= dim1>=dim2 ? dim1:dim2;
	get_area(dim);
	for(size_t i=0; i<dim; i++){
		double t1=pt1[i], t2=pt2[i];
		if(t1<=t2) m_range[i]=MGInterval(t1,t2);
		else       m_range[i]=MGInterval(t2,t1);
	}
}

//Construct Box by providing each Interval
MGBox::MGBox(size_t dim, const MGInterval* intervl)
:m_sdim(0), m_range(0){
	assert (dim>0);
	get_area(dim);
	for(size_t i=0; i<dim; i++) m_range[i]=intervl[i];
}

//Construct Box which contains both input box and a point.
MGBox::MGBox(const MGBox& boxi, const MGPosition& point)
:m_sdim(0), m_range(0){
	size_t dim1=boxi.sdim(); size_t dim2=point.sdim();
	size_t dim= dim1>=dim2 ? dim1:dim2;
	get_area(dim);
	for(size_t i=0; i<dim; i++){
		m_range[i]=MGInterval(boxi.m_range[i],point[i]);
	}
}

//Construct Box by copying old Box, changing space dimension and
//ordering of old coordinates.
MGBox::MGBox(size_t dim, const MGBox& box, size_t start1, size_t start2)
:m_sdim(0), m_range(0){
	assert(dim>0 && start1<dim && start2<box.sdim());

	get_area(dim);
	size_t dim2=box.sdim();
	size_t dimmin= dim<=dim2 ? dim:dim2;
	size_t k,j=start2,i=start1;
	for(k=0; k<dimmin; k++){
		m_range[i++]=box.m_range[j++];
		if(i>=dim) i=0;	if(j>=dim2) j=0;
	}
	while(k++<dim){
		m_range[i++]=MGDefault::empty_interval();
		if(i>=dim) i=0;
	}
}

//////////Destructor//////////////
MGBox::~MGBox(){
	if(m_sdim>3) delete[] m_range;
}

//////////Operator Overload//////////

//Assignment.
MGBox& MGBox::operator= (const MGBox& box2){
	get_area(box2.sdim());
	for(size_t i=0; i<m_sdim; i++) m_range[i]=box2.m_range[i];
	return *this;
}

//Boxを順方向に平行移動してできるオブジェクトを生成する。 
MGBox MGBox::operator+ (const MGVector& vec) const{
	MGBox new_box = *this;
	new_box += vec;
	return new_box;
}

MGBox operator+ (const MGVector& v, const MGBox& b){
	return b+v;
}

//Boxを順方向に平行移動し自身のBoxとする。
MGBox& MGBox::operator+= (const MGVector& vec){
	for(size_t i=0; i<sdim(); i++) m_range[i] += vec[i];
	return *this;
}

//Boxを逆方向に平行移動してできるオブジェクトを生成する。
MGBox MGBox::operator- (const MGVector& vec) const{
	MGBox new_box = *this;
	new_box -= vec;
	return new_box;
}

//Boxを逆方向に平行移動し自身のBoxとする。
MGBox & MGBox::operator-= (const MGVector& vec){
	for(size_t i=0; i<sdim(); i++) m_range[i] -= vec[i];
	return *this;
}

//Boxを拡大してできるオブジェクトを生成する。
MGBox MGBox::operator* (double a) const{
	size_t spdim=sdim();
	MGBox temp(spdim);
	for(size_t i=0; i<spdim; i++) temp.m_range[i]=m_range[i]*a;
	return temp;
}

//与えられたマトリックスでBoxの変換を行いオブジェクトを生成する。
MGBox MGBox::operator* (const MGMatrix& matrix) const {
	MGBox new_box = *this;
	new_box *= matrix;
	return new_box;
}

//与えられた変換でBoxのトランスフォームを行い、
//それらを囲むBoxのオブジェクトを生成する。
MGBox MGBox::operator* (const MGTransf& transf) const {
	MGBox new_box = *this;
	new_box *= transf;
	return new_box;
}

//自身のBoxを拡大し自身のBoxとする。
MGBox& MGBox::operator*= (double scalar) {
	for(size_t i=0; i<sdim(); i++) m_range[i] *= scalar;
	return *this;
}

//与えられたマトリックスでBoxの変換を行い自身のBoxとする。
MGBox& MGBox::operator*= (const MGMatrix& matrix) {

	size_t dim1=sdim(), dim2=matrix.sdim();
	size_t dim=dim1>=dim2 ? dim1:dim2;
	if(!dim) return *this;
	if(dim>dim1) resize(dim);

	if(finite()){
		MGPosition pmin(dim), pmax(dim), p(dim);
		std::vector<MGPosition> vertices=vertex();
		size_t len=vertices.size(),i,j;
//		for(size_t k=0; k<len; k++) cout<<vertices(k);
		//Compute transformed box and minimum and maximum.
		pmin=pmax=vertices[0]*matrix;
		for(i=1; i<len; i++){
			p=vertices[i]*matrix;
			for(j=0; j<dim; j++){
				double a=p.ref(j);
				if(a<pmin.ref(j)) pmin(j)=a;
				if(a>pmax.ref(j)) pmax(j)=a;
			}
		}
		//Set transformed interval.
		for(j=0; j<dim; j++) m_range[j]=MGInterval(pmin(j), pmax(j));
	}else{
		for(size_t i=0; i<dim; i++)
			m_range[i]=MGInterval(MGINTERVAL_INFINITE);
	}
	return *this;
}

//与えられた変換でBoxのトランスフォームを行い、
//それらを囲むBoxを自身のBoxとする。
MGBox& MGBox::operator *= (const MGTransf & transf ) {
    // 自身のBoxの全辺が上下無限でも empty でもない時処理を行う
    if(finite()) {
		*this *= transf.affine();
		*this += transf.translation();
	}
	// もし、全辺上下無限の時はそのまま返却される
	return *this;
}

//Boxを縮小してできるオブジェクトを生成する。 
MGBox MGBox::operator / (double a) const{
	MGBox temp(*this);
	return temp/=a;
}

//Boxを縮小し自身のBoxとする。
MGBox& MGBox::operator /= ( double scalar ) {
	for(size_t i=0; i<sdim(); i++) m_range[i] /= scalar;
	return *this;
}

//自身のBoxと与えられたBoxを内包する最小のBoxを生成する。
MGBox MGBox::operator| (const MGBox& box2) const{
	MGBox temp(*this);
	return temp |= box2;
}

//自身のBoxと与えられたBoxを内包する最小のBoxを自身のBoxとする。 
MGBox& MGBox::operator |= (const MGBox &box2) {
	size_t dim1=sdim(); size_t dim2=box2.sdim();
	if(!dim2) return *this;
	if(!dim1) *this=box2;
	else{
		size_t dim=dim2; size_t i;
		if(dim1<dim2) {
			dim=dim1;
			resize(dim2);
			for(i=dim1; i<dim2; i++) m_range[i]=box2.m_range[i];
		}
		for(i=0; i<dim; i++) m_range[i] |= box2.m_range[i];
	}
	return *this;
}

//自身のBoxと与えられたBoxの共通部分のBoxを生成する。 
MGBox MGBox::operator& ( const MGBox& box2) const{
	MGBox temp(*this);
	return temp &= box2;
}

//自身のBoxと与えられたBoxの共通部分のBoxを自身のBoxとする。
MGBox& MGBox::operator&= ( const MGBox & box2 ) {
	size_t dim1=sdim(); size_t dim2=box2.sdim();
	size_t dim=(dim1<dim2) ? dim1:dim2;
	for(size_t i=0; i<dim; i++) m_range[i] &= box2.m_range[i];
	if(dim1>dim) resize(dim);
	return *this;
}

//自身と与えられたBoxが等しいかどうかを返却する。
//一方のBoxの全ての辺がもう一方の相当する辺に重なる場合 True を返却する。
bool MGBox::operator== (const MGBox& box2) const {
	size_t dim1=sdim(); size_t dim2=box2.sdim();
	size_t dim=dim1<dim2 ? dim1:dim2;
	bool judge = true;
	for(size_t i=0; i<dim; i++)
		if(m_range[i] != box2.m_range[i]) {judge=false; break;}
	if(judge){
		if(dim1>dim2){
			for(size_t i=dim2; i<dim1; i++)
				if(!m_range[i].empty()) {judge=false; break;}
		}
		else if(dim1<dim2){
			for(size_t i=dim1; i<dim2; i++)
				if(!box2.m_range[i].empty()) {judge=false; break;};
		}
	}
	return judge;
}

//与えられたポジションが自身のBox内に含まれているか返却する。
bool MGBox::operator >> ( const MGPosition & pt) const {
	size_t i, dim1=sdim(); if(!dim1) return false;
	size_t dim2=pt.sdim();
	// ポジションのdataがBoxの全ての辺の範囲内にある場合
	// True を返却する。
	for(i=0; i<dim1; i++) if(m_range[i] << pt[i]) return false;
	for(i=dim1; i<dim2; i++) if(!MGAZero(pt[i])) return false;
	return true;
}

//自身のBoxが与えられたBoxを囲んでいるか返却する。
//与えられたBoxが自身のBox内にある場合 True を返却する。
bool MGBox::operator >> ( const MGBox & box2 ) const {
	size_t i, dim1=sdim(); if(!dim1) return false;
	size_t dim2=box2.sdim(); if(!dim2) return true;
	size_t dim= dim1<dim2 ? dim1:dim2;
	// 与えられたBoxの全辺が自身のBoxの各辺の
	// Intervalに含まれる時、True
	for(i=0; i<dim; i++){
		if( !(m_range[i]>>box2.m_range[i]) ) return false;
	}
	//Check extra dimension of box2.
	for(i=dim1; i<dim2; i++){
		if(!(MGAZero(box2.m_range[i].high_point())&&
			MGAZero(box2.m_range[i].low_point()))) return false;
	}
	//Check extra dimension of this.
	for(i=dim2; i<dim1; i++) if(m_range[i]<<0.) return false;
	return true;
}

//Friend Function.
// Boxを拡大してできるオブジェクトを生成する。
MGBox operator*(double scale, const MGBox &box) {
	return box*scale;
}

///////////////Member Function//////////////

//Test if the line segment P01 from P0 to P1 is crossing this box or not.
//Function's return value is true if a part of P01 is included in this box,
//false if not.
bool MGBox::crossing(const MGStraight& sl)const{
	size_t sd1=sdim(), sd2=sl.sdim();
	if(!sd1) return false;
	if(!sd2) return false;

	if(sd1==1 && sd2==1){
		MGBox slbx=sl.box();
		slbx&=*this;
		return !slbx.empty();
	}
	
	if(sd1==sd2){
		if(sd1==2) return crossing2D(sl);
		for(size_t i=0; i<sd1; i++){
			MGBox bx2(2, *this, 0, i);
			MGStraight sl2(2,sl,0,i);
			if(!bx2.crossing2D(sl2)) return false;
		}
		return true;
	}

	size_t sd=sd1; if(sd1<sd2) sd=sd2;
	MGBox bxt(sd,*this);
	MGStraight slt(sd,sl);
	for(size_t i=0; i<sd; i++){
		MGBox bx2(2, bxt, 0, i);
		MGStraight sl2(2,slt,0,i);
		if(!bx2.crossing2D(sl2)) return false;
	}
	return true;
}

//Test if the straight line sl is crossing this box or not.
//Function's return value is true if a part of sl is included in this box,
//false if not. 2D version of crossing().
bool MGBox::crossing2D(const MGStraight& sl)const{
	MGSTRAIGHT_TYPE stype=sl.straight_type();
	if(stype==MGSTRAIGHT_EMPTY) return false;

	const MGInterval& ix=(*this)[0];
	const MGInterval& iy=(*this)[1];
	double bx0=ix.low_point(),bx1=ix.high_point();
	double by0=iy.low_point(),by1=iy.high_point();
	MGPosition P0;
	if(stype==MGSTRAIGHT_UNLIMIT) P0=sl.root_point();
		//Case that sl is an infinite line.
	else P0=sl.start_point();
	double x0=P0[0], y0=P0[1];
	if(bx0<=x0 && x0<=bx1 && by0<=y0 && y0<=by1) return true;

	MGVector v01;
	if(stype==MGSTRAIGHT_SEGMENT){
		MGPosition P1=sl.end_point();
		double x1=P1[0], y1=P1[1];

		if(x0<bx0 && x1<bx0) return false;
		if(x0>bx1 && x1>bx1) return false;
		if(y0<by0 && y1<by0) return false;
		if(y0>by1 && y1>by1) return false;

		if(bx0<=x1 && x1<=bx1 && by0<=y1 && y1<=by1) return true;
		if(bx0<=x0 && x0<=bx1 && bx0<=x1 && x1<=bx1) return true;
		if(by0<=y0 && y0<=by1 && by0<=y1 && y1<=by1) return true;
		v01=MGVector(P1,P0);
	}else v01=sl.direction();

	MGPosition BP0(bx0,by0), BP1(bx1,by1);
	//corner points of this box to test crossing range.
	if(x0>bx1 && by0<=y0 && y0<=by1) BP0(0)=bx1;
	if(x0<bx0 && by0<=y0 && y0<=by1) BP1(0)=bx0;
	if(y0>by1){
		if(x0>=bx0){
			BP0(1)=by1;
			if(x0>bx1) BP1(1)=by0;
		}
	}
	if(y0<by0){
		if(x0<=bx1){
			BP1(1)=by0;
			if(x0<bx0) BP0(1)=by1;
		}
	}
	MGVector v0bp0(BP0,P0), v0bp1(BP1,P0);
	if((v01*v0bp0)%(v01*v0bp1) <= 0.) return true;
	return false;
}

//Test if the plane is cutting this box or not.
//Function's return value is true: if the plane is cutting the box,
//false if all of the vertices of the box are in one side of the plane.
bool MGBox::cutting(const MGPlane& plane)const{
	std::vector<MGPosition> vertices=vertex();
	double d0=plane.distance(vertices[0]);
	double dmin=d0, dmax=d0;
	for(size_t i=1; i<8; i++){
		double di=plane.distance(vertices[i]);
		if(dmin>di) dmin=di;
		if(dmax<di) dmax=di;
	}
	double err=MGTolerance::wc_zero();
	dmin-=err; dmax+=err;
	return dmin*dmax<=0.;
}

//Compute the distance from a point to the box.
double MGBox::distance(const MGPosition& P) const{
	double dist=0.;
	size_t sd=sdim(), sd2=P.sdim(),i;
	for(i=0; i<sd; i++){
		double pi=P.ref(i);
		if(m_range[i]>pi){
			double len=(m_range[i].low_point()-pi);
			dist+=len*len;
		}else if(m_range[i]<pi){
			double len=(pi-m_range[i].high_point());
			dist+=len*len;
		}
	}
	for(i=sd; i<sd2; i++){
		double len=P.ref(i);
		dist+=len*len;
	}
	return sqrt(dist);
}

//Boxが empty かどうか返却する
bool MGBox::empty() const{
	//empty()=true when any of the intervals is empty.
	bool judge = false;
	size_t dim=sdim();
	for(size_t i=0; i<dim; i++){
		if(m_range[i].empty()){judge = true; break;}
	}
	return judge; 
}

//Expand the box by len. This box will be updated so that the center of
//the box will not be moved and the box is widened by len for each coordinate.
//That is, ref(i).high() is set to ref(i).high()+len
//and ref(i).low() is set to ref(i).low()-len for each i.
void MGBox::expand(double len){
	size_t sd=sdim();
	for(size_t i=0; i<sd; i++){
		m_range[i].low()-=len; m_range[i].high()+=len;
	}
}

//Expand the box by len[]. This box will be updated so that the center of
//the box will not be moved and the box is widened by len[i]
//for each coordinate. That is, ref(i).high() is set to ref(i).high()+len[i]
//and ref(i).low() is set to ref(i).low()-len[i] for each i.
void MGBox::expand(double* len){
	size_t sd=sdim();
	for(size_t i=0; i<sd; i++){
		m_range[i].low()-=len[i]; m_range[i].high()+=len[i];
	}

}

//Expand the box by MGTolerance::rc_zero(). That is,
//operator*=(1.+MGTolerance::rc_zero()) will be executed.
void MGBox::expand(){
	operator*=(1.+MGTolerance::rc_zero());
}

//Expand the box so that this contains the position P.
void MGBox::expand(const MGPosition& P){
	size_t dim1=sdim(); size_t dim2=P.sdim();
	size_t dim= dim1>=dim2 ? dim1:dim2;
	resize(dim);
	for(size_t i=0; i<dim; i++) m_range[i].expand(P[i]);
}

//Return true if box is finite.
bool MGBox::finite() const{
	size_t i;
	size_t dim=sdim();
	bool judge=true;
	for(i=0; i<dim; i++){
		if(!(m_range[i].finite())){judge = false; break;}
	}
	return 	judge;
}

//Get the area of m_range for the space dimension sdim.
//Result area will contain garbages.
//get_area will use m_sdim as input. m_sdim must be valid.
void MGBox::get_area(size_t dim){
	if(dim<=3){
		if(m_sdim>3) delete[] m_range;
		m_range=m_rData;
	}else{
		if(dim>m_sdim){
			if(m_sdim>3) delete[] m_range;
			m_range= new MGInterval[dim];
		}
	}
	m_sdim=dim;
}

//Boxの対角線の両端(high(), low())と中心(mid())を返却する。全ての座標値が
//最小の点が low () で、最大の点が high () で返却される。 
MGPosition MGBox::high() const {
	size_t dim=sdim();
	MGPosition high_pt(dim);

	// 全てのIntervalが少なくとも上方有限の時
	// 最大座標値が取得できる
	for(size_t i=0; i<dim; i++) high_pt(i)=m_range[i].high_point();
	return high_pt;
}

//原点が自身のBox内に含まれているか返却する。
//Test if this includes the origin(0., 0., ...).
bool MGBox::includes_origin()const{
	size_t i, dim1=sdim(); if(!dim1) return false;
	for(i=0; i<dim1; i++) if(m_range[i] << 0.0) return false;
	return true;
}

//ボックスの対角線長さを求める
//Return diagonal line length.
double MGBox::len() const{
	return (high() - low()).len();
}

MGPosition MGBox::low() const {
	size_t dim=sdim();
	MGPosition low_pt(dim);
	for(size_t i=0; i<dim; i++) low_pt(i)=m_range[i].low_point();
	return low_pt;
}

MGPosition MGBox::mid() const {
	size_t dim=sdim();
	MGPosition mid_pt(dim);
	for(size_t i=0; i<dim; i++) mid_pt(i)=m_range[i].mid_point();
	return mid_pt;
}

//Access to i-th Inteval.
const MGInterval& MGBox::ref(size_t i)const{ 
	if(i<sdim()) return m_range[i];
	else return MGDefault::empty_interval(); // Return null interval.
}

//Resize the m_range.
//If sim>current sdim(), the old data are guaranteed to hold in resized box.
void MGBox::resize(size_t dim){
	if(dim==m_sdim) return;
	if(dim<=3){
		if(m_sdim<=3){
			m_range=m_rData;
			for(size_t i=m_sdim; i<dim; i++)
				m_range[i]=MGDefault::empty_interval();
		}else{
			for(size_t i=0; i<dim; i++) m_rData[i]=m_range[i];
			delete[] m_range;
			m_range=m_rData;
		}
	}else{
		if(dim>m_sdim){
			MGInterval* rng=new MGInterval[dim];
			for(size_t i=0; i<m_sdim; i++) rng[i]=m_range[i];
			if(m_sdim>3) delete[] m_range;
			m_range=rng;
		}
	}
	m_sdim=dim;
}

//自身のBoxの最大座標値を指定された点に変更する
MGBox& MGBox::set_high(const MGPosition& high_pt){
	size_t dim1=sdim(); size_t dim2=high_pt.sdim();
	if(dim2>dim1) resize(dim2);
	for(size_t i=0; i<dim2; i++) m_range[i].set_high_point(high_pt.ref(i));
	return *this;
}

//自身のBoxの最小座標値を指定された点に変更する
MGBox& MGBox::set_low(const MGPosition &low_pt) {
	size_t dim1=sdim(); size_t dim2=low_pt.sdim();
	if(dim2>dim1) resize(dim2);
	for(size_t i=0; i<dim2; i++) m_range[i].set_low_point(low_pt.ref(i));
	return *this;
}

//Set this box as a null box.
void MGBox::set_null(){
	if(m_sdim>3) delete[] m_range;
	m_sdim=0; m_range=0;
}

// vertex computes all the vertices of the box.
std::vector<MGPosition> MGBox::vertex() const{

	size_t dimi=sdim(); 
	double a1,a2,b1,b2;

	if(!dimi) return std::vector<MGPosition>();
	else if(dimi==1){
		std::vector<MGPosition> data(2);
		double *p1=&a1, *p2=&b1;
		a1=m_range[0].low_point(); b1=m_range[0].high_point();
		data[0]=MGPosition(1,p1); data[1]=MGPosition(1,p2);
		return data;
	}
	a1=m_range[0].low_point(); b1=m_range[0].high_point();
	a2=m_range[1].low_point(); b2=m_range[1].high_point();
	std::vector<MGPosition> data(4);
	data[0]=MGPosition(a1,a2);
	data[1]=MGPosition(a1,b2);
	data[2]=MGPosition(b1,a2);
	data[3]=MGPosition(b1,b2);
	size_t dim=2, dimp1=dim+1;
	size_t len=4,len2=len*2,i;
	while(dim<dimi){
		std::vector<MGPosition> data1=data;
		data.resize(len2);
		for(i=0; i<len; i++){
			MGPosition q(dimp1,data1[i]);
			q(dim)=m_range[dim].low_point(); data[2*i]=q;
			q(dim)=m_range[dim].high_point(); data[2*i+1]=q;
		}
		dim=dimp1; dimp1=dim+1;
		len=len2; len2=len*2;
	}
	return data;
}
