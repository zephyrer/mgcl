/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Vector.h"
#include "mg/Position.h"
#include "mg/Unit_vector.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Implemetation of class Vector.

// <<< 1. Constructor >>>

//void Constructor
MGVector::MGVector(size_t dim)
:m_sdim(dim),m_element(m_data),m_length(-1.)
{
	if(dim>3) m_element=new double[dim];
	if(!dim) m_length=0.;
}

//Construct 2D vector by providing each element data.
MGVector::MGVector ( double x, double y)
:m_sdim(2),m_element(m_data),m_length(-1.) {
	m_element[0]=x; m_element[1]=y;
}

// ｘ，ｙ，ｚを指定して3Dベクトルを生成する
MGVector::MGVector(double x, double y, double z)
:m_sdim(3),m_element(m_data),m_length(-1.) {
	m_element[0]=x; m_element[1]=y;	m_element[2]=z;
}

// ｘ，ｙ，ｚ, wを指定して4Dベクトルを生成する
MGVector::MGVector(double x, double y, double z, double w)
:m_sdim(4),m_element(new double[4]),m_length(-1.) {
	m_element[0]=x; m_element[1]=y;	m_element[2]=z; m_element[3]=w;
}

// 初期値　v ですべてのエレメントを初期化してベクトルを生成する。
MGVector::MGVector(size_t dim, double v)
:m_sdim(dim),m_element(m_data),m_length(-1.) {
	if(dim>3) m_element=new double[dim];
	for(size_t i=0; i<dim; i++) m_element[i]=v;
	if(!dim) m_length=0.;
}

// double の配列でcoordinate valueを指定しベクトルを生成する。
MGVector::MGVector(size_t dim, const double* v)
:m_sdim(dim),m_element(m_data),m_length(-1.) {
	if(dim>3) m_element=new double[dim];
	for(size_t i=0; i<dim; i++) m_element[i]=v[i];
	if(!dim) m_length=0.;
}

//Construct a vector from a class MGPosition.
MGVector::MGVector(const MGPosition& point)
	:m_sdim(0),m_element(m_data)
{	*this=point.vector();}

//Construct a vector from a difference of two vectors.
MGVector::MGVector(const MGVector& dvec,	//Destination point
				   const MGVector& svec)	//Source point
:m_sdim(0),m_element(m_data),m_length(-1.) {
	size_t dimd=dvec.sdim(), dims=svec.sdim(), i;

	if(dimd==dims){
		m_sdim=dimd;
		if(dimd>3) m_element=new double[dimd];
		for(i=0; i<dimd; i++) 
			m_element[i]=dvec.m_element[i]-svec.m_element[i];
	}else if(dimd<dims){
		(*this) = -svec;
		for(i=0; i<dimd; i++) m_element[i] +=dvec.m_element[i];
	}else{
		(*this) = dvec;
		for(i=0; i<dims; i++) m_element[i] -=svec.m_element[i];
	}
	if(!sdim()) m_length=0.;
}
		 
//Construct Vector by copying old Vector, changing space dimension and
//ordering of old coordinates.
MGVector::MGVector(size_t dim, const MGVector& vec2,
				   size_t start1, size_t start2)
:m_sdim(dim),m_element(m_data),m_length(vec2.m_length) {
	if(m_sdim>3) m_element=new double[dim];

	assert (dim>0 && start1<dim && start2<vec2.sdim());

	size_t dim2=vec2.sdim(); 
	size_t dimmin= dim<dim2 ? dim:dim2; 
	size_t k,j=start2,i=start1;
	for(k=0; k<dimmin; k++) {
		m_element[i++]=vec2.m_element[j++];
		if(i>=dim) i=0; if(j>=dim2) j=0;
	}
	while(k++<dim){
		m_element[i++]=0.;
		if(i>=dim) i=0;
	}
	if(dim<dim2) m_length=-1;
}

//Construct from std::vector<double>
MGVector::MGVector(const std::vector<double>& darrays)
:m_sdim(darrays.size()),m_element(m_data),m_length(-1.){
	if(m_sdim>3) m_element=new double[m_sdim];
	for(size_t i=0; i<m_sdim; i++) m_element[i]=darrays[i];
	if(!m_sdim) m_length=0.;
}

/*
//Construct from std::valarray<double>
MGVector::MGVector(const std::valarray<double>& darrays)
:m_sdim(darrays.size()),m_element(m_data),m_length(-1.){
	if(m_sdim>3) m_element=new double[m_sdim];
	for(size_t i=0; i<m_sdim; i++) m_element[i]=darrays[i];
	if(!m_sdim) m_length=0.;
}*/

// Copy Constructor
MGVector::MGVector(const MGVector& v)
:m_sdim(v.m_sdim),m_length(v.m_length){
	if(m_sdim>3) m_element=new double[m_sdim];
	else m_element=m_data;
	for(size_t i=0; i<m_sdim; i++) m_element[i]=v.m_element[i];
}

//Destructor
//

//
//Operator Overload
//

//Assignment
MGVector& MGVector::operator =(const MGVector& v){
	if(this==&v) return *this;
	m_length=v.m_length;
	size_t i,dim=v.m_sdim;
	if(dim<=m_sdim){
		if(dim<=3 && m_sdim>3){
			delete[] m_element; m_element=m_data;
		}
		for(i=0; i<dim; i++) m_element[i]=v.m_element[i];
	}else{
		if(dim>3){
			if(m_sdim>3) delete[] m_element;
			m_element=new double[dim];
		}
		for(i=0; i<dim; i++) m_element[i]=v.m_element[i];
	}
	m_sdim=dim;
	return *this;
}

//Access to i-th Inteval.
double& MGVector::operator()(size_t i) { 
	assert(i<sdim());
	m_length=-1.;
	return m_element[i];
}

//Update vector data by array of double.
MGVector& MGVector::operator=(const double* data){
	for(size_t i=0; i<m_sdim; i++) m_element[i]=data[i];
	m_length=-1.;
	return *this;
}

// ふたつのベクトルの加算
//Addition of two vectors.
MGVector operator+(const MGVector& vec1,const MGVector& vec2){
	size_t dim=vec1.m_sdim;
	size_t dimmin=dim;
	if(dim<vec2.m_sdim)
		dim=vec2.m_sdim;
	else
		dimmin=vec2.m_sdim;

	MGVector temp(dim);
	size_t i;
	for(i=0; i<dimmin; i++)
		temp.m_element[i]=vec1.m_element[i]+vec2.m_element[i];
	for(; i<dim; i++)
		temp.m_element[i]=vec1.ref(i)-vec2.ref(i);
	return temp;
}

// 自身のベクトルに与えられたベクトルを加算して自身のベクトルとする
MGVector& MGVector::operator+= (const MGVector &vec2) {
	size_t dim1=sdim(); size_t dim2=vec2.sdim();
	size_t i;
	if(dim1>=dim2){
		for(i=0; i<dim2; i++) m_element[i] += vec2.m_element[i];
		m_length=-1.;
	}else{
		for(i=0; i<dim1; i++) m_element[i] += vec2.m_element[i];
		resize(dim2);
		for(; i<dim2; i++) m_element[i] = vec2.m_element[i];
	}
	return *this;
}

// 単項マイナス。自身のベクトルを反転し、オブジェクトを生成
MGVector MGVector::operator- () const{
	MGVector temp(m_sdim);
	for(size_t i=0; i<m_sdim; i++) temp.m_element[i] = -m_element[i];
	temp.m_length=m_length;
	return temp;
}

// ベクトルの減算
//Subtraction of two vectors.
MGVector operator-(const MGVector& vec1,const MGVector& vec2){
	size_t dim=vec1.m_sdim;
	size_t dimmin=dim;
	if(dim<vec2.m_sdim)
		dim=vec2.m_sdim;
	else
		dimmin=vec2.m_sdim;

	MGVector temp(dim);
	size_t i;
	for(i=0; i<dimmin; i++)
		temp.m_element[i]=vec1.m_element[i]-vec2.m_element[i];
	for(; i<dim; i++)
		temp.m_element[i]=vec1.ref(i)-vec2.ref(i);
	return temp;
}

// 自身のベクトルと与えられたベクトルの減算を行い自身のベクトルとする
MGVector& MGVector::operator-= (const MGVector &vec2) {
	size_t dim1=sdim(); size_t dim2=vec2.sdim();
	size_t i;
	if(dim1>=dim2){
		for(i=0; i<dim2; i++) m_element[i] -= vec2.m_element[i];
		m_length=-1.;
	}else{
		for(i=0; i<dim1; i++) m_element[i] -= vec2.m_element[i];
		resize(dim2);
		for(; i<dim2; i++) m_element[i] = -vec2.m_element[i];
	}
	return *this;
}

// スカラーの乗算を行いオブジェクトを生成
//Scalar multiplication.
MGVector operator*(const MGVector& vec1,double scale){
	size_t sd=vec1.m_sdim;
	MGVector new_vec(sd);
	for(size_t i=0; i<sd; i++)
		new_vec.m_element[i]=vec1.m_element[i]*scale;
	new_vec.m_length=fabs(scale)*vec1.m_length;
	return new_vec;
}

// スカラーの乗算を行い自身のベクトルとする
MGVector& MGVector::operator*= (double scale){
	for(size_t i=0; i<m_sdim; i++) m_element[i] *= scale;
	if(m_length>=0.) m_length*=fabs(scale);
	return *this;
}

//ベクタの外積
//vector product of two vectors.
MGVector operator*(const MGVector& vec1,const MGVector& vec2){
	MGVector v(vec1); 
	return v*=vec2;
}

// ベクトルの外積を行い自身のベクトルとする
MGVector& MGVector::operator*= (const MGVector &vec2){
	double d0,d1,d2;
	size_t dim=sdim(), dim2=vec2.sdim();
	if(dim<dim2) dim=dim2;
	if(dim<=3){
		d0=vec2.ref(2)*ref(1) - vec2.ref(1)*ref(2);
		d1=vec2.ref(0)*ref(2) - vec2.ref(2)*ref(0);
	    d2=vec2.ref(1)*ref(0) - vec2.ref(0)*ref(1);
		m_element[0]=d0; m_element[1]=d1; m_element[2]=d2;
		m_sdim=3; m_length=-1;
		return *this;
	} else {
		//When dim>3, find three id of this and vec2 that construct
		//3D vectors whose inner products are zero,
		//and whose vector products has maximum vector length.
		size_t i,j,k, i1,j1,k1;
		double v10,v11,v12, v20,v21,v22;
		double d0s,d1s,d2s;
		double dmaxt,dmaxs=-1.;
		for(i=0; i<dim-3; i++){
			v10=ref(i); v20=vec2.ref(i);
			for(j=i+1; j<dim-2; j++){
				v11=ref(j); v21=vec2.ref(j);
				d2=v21*v10-v20*v11;
				for(k=j+1; k<dim-1; k++){
					v12=ref(k); v22=vec2.ref(k);
					d0=v22*v11-v21*v12;
					d1=v20*v12-v22*v10;
					dmaxt=d0*d0+d1*d1+d2*d2;
					if(dmaxt>dmaxs){
						dmaxs=dmaxt;
						i1=i; j1=j; k1=k;
						d0s=d0; d1s=d1; d2s=d2;
					}
				}
			}
		}
		double sa=sangle(vec2);
		double v1v2stheta=len()*vec2.len()*sa;
		double vlen=sqrt(dmaxs);
		if(m_sdim>3) delete[] m_element;
		m_sdim=k1+1; m_length=-1.;
		if(m_sdim<=3) m_element=m_data; else m_element=new double[m_sdim];
		for(i=0; i<m_sdim; i++) m_element[i]=0.;
		if(MGMZero(vlen)){
			m_element[i1]=d0s; m_element[j1]=d1s; m_element[k1]=d2s;
		}else{
			double one_vlen=v1v2stheta/vlen;
			m_element[i1]=d0s*one_vlen;
			m_element[j1]=d1s*one_vlen;
			m_element[k1]=d2s*one_vlen;
		}
		return *this;
	}
}

// トランスフォームを行いオブジェクトを生成

// ベクトルの内積
//Inner product of two vectors.
double operator%(const MGVector& vec1,const MGVector& vec2){
	size_t dim=vec1.sdim(); size_t dim2=vec2.sdim();
	if(dim>dim2) dim=dim2;
	double product=0.;
	for(size_t i=0; i<dim; i++)
		product+=vec1.m_element[i]*vec2.m_element[i];
	return product;
}

// スカラー除算を行いオブジェクトを生成
//Scalar division.
MGVector operator/(const MGVector& vec1,double scale){
	size_t sd=vec1.m_sdim;
	MGVector new_vec(sd);
	for(size_t i=0; i<sd; i++)
		new_vec.m_element[i] = vec1.m_element[i]/scale;
	if(vec1.m_length>0.)
		new_vec.m_length=vec1.m_length/fabs(scale);
	return new_vec;
}

// スカラーの除算を行い自身のベクトルとする
MGVector& MGVector::operator/= (double scalar){
	for(size_t i=0; i<m_sdim; i++) m_element[i] /= scalar;
	if(m_length>0.) m_length/=fabs(scalar);
	return *this;
}

// 与えられたベクトルの成分の値を比較し、同じであれば TRUE を返却
//Test if two vectors are equal.
bool operator==(const MGVector& v1,const MGVector& v2){
	// ベクトルの差分を取得し、そのベクトルが０ベクトルの時等しい
	double len1=v1.len(), len2=v2.len();
	double dif=(v1 - v2 ).len();
	if(len1>=len2){
		if(len1<=MGTolerance::mach_zero()) return 1;
		else return MGRZero2(dif,len1);
	}else{
		if(len2<=MGTolerance::mach_zero()) return 1;
		else return MGRZero2(dif,len2);
	}
}

//Member Function
//

// 自身のベクトルと与えられたベクトルのなす角度を Radian で返却
//Compute angle in radian of two vectors.
// 0<= angle <pai.
double MGVector::angle(const MGVector& vec2) const{
	double angle;
	double ca=cangle(vec2);
	if(sdim()>3 || vec2.sdim()>3) return acos(ca); 
	double sa=sangle(vec2); // Note that sa >=0 always holds.
	if(ca>=sa)                angle=asin(sa);          // 0<=  <=pai/4.
	else if(ca>=0. && sa>=ca) angle=acos(ca);          // pai/4 <=  <=pai/2.
	else if(sa>=-ca)          angle=acos(ca);          // pai/2<=  <=3*pai/4.
	else                      angle=mgPAI-asin(sa);   // 3*pai/4<= <=pai.
	return angle;
}

//Compute angle in radian measured from this to v2 around the normal N.
//The angle's range is 0<= angle <2*pai.
//Although N is assumed to be parallel to N2=(*this)*v2, N may not perpendicular
//to v1 and v2, in which case, the projected N to N2 is used to measure the angle.
//v1.angle2pai(v2,N)+v2.angle2pai(v1,N)=2*pai always holds.
double MGVector::angle2pai(const MGVector& v2, const MGVector& N)const{
	double ca=cangle(v2);
	double sa=sangle(v2); // Note that sa >=0 always holds.
	MGVector v1v2=(*this)*v2;
	if(v1v2%N<0.) sa*=-1.;
	return MGAngle(ca,sa);
}

// 自身のベクトルと与えられたベクトルのなす角度を cosΘ で返却する
// 自身か与えられたベクトルが零ベクトルの時は、cosΘは 1.0 とする
double MGVector::cangle ( const MGVector & vec2 ) const {
	double cos_theta;
	double ll=len()*vec2.len();
	if(MGMZero(ll)) cos_theta=1.0;
	else{
		cos_theta = (*this)%vec2;
		cos_theta/=ll;
	}
	return cos_theta;
}
	
//Clear all the elements by the value init.
MGVector& MGVector::clear(double init){
	for(size_t i=0; i<m_sdim; i++) m_element[i]=init;
	m_length=-1.;
	return *this;
}

// Generate a vector by interpolating two vectors. Input scalar is a ratio
// and when zero, output vector is a copy of the own vector.
MGVector MGVector::interpolate(double t2, const MGVector& vec2) const{
	double t1=1.0-t2; 
	size_t dim1=sdim(); size_t dim2=vec2.sdim(); size_t i;
	size_t dim = dim1<dim2 ? dim2:dim1;
	MGVector temp(dim);
	if(dim1==dim2){
		for(i=0; i<dim; i++)
			temp.m_element[i]=t1*m_element[i]+t2*vec2.m_element[i];
	}else if(dim1<dim2){
		for(i=0; i<dim1; i++)
			temp.m_element[i]=t1*m_element[i]+t2*vec2.m_element[i];
		for(i=dim1; i<dim2; i++)
			temp.m_element[i]=t2*vec2.m_element[i];
	}else{
		for(i=0; i<dim2; i++)
			temp.m_element[i]=t1*m_element[i]+t2*vec2.m_element[i];
		for(i=dim2; i<dim1; i++)
			temp.m_element[i]=t1*m_element[i];
	}
	return temp;
}

// Generate a vector by interpolating two vectors by rotation.
//Input scalar t is a ratio and when zero, output vector is a copy of *this.
// New vector vnew=a*(*this)+b*vec2, where
// a=sin(theta2)/sin(theta), b=sin(theta1)/sin(theta). Here,
// theta=angle of *this and vec2. theta1=t*theta, theta2=theta-theta1.
// theta may be zero.
//When ratio is not null, ratio[0]=a and ratio[1]=b will be returned.
MGVector MGVector::interpolate_by_rotate(
	double t, const MGVector& vec2,
	double* ratio
)const{
	double theta=angle(vec2);
	double stheta=sin(theta);
	double a[2];
	if(!ratio){
		ratio=a;
	}
	if(MGMZero(stheta)){
		ratio[1]=t; ratio[0]=1.-t;
	} else{
		double theta1=t*theta;
		double theta2=theta-theta1;
		ratio[0]=sin(theta2)/stheta;
		ratio[1]=sin(theta1)/stheta;
	}
	size_t sd=sdim(), sd2=vec2.sdim();
	if(sd2>sd)
		sd=sd2;
	MGVector vec(sd);
	for(size_t i=0; i<sd; i++)
		vec.m_element[i]=ratio[0]*ref(i)+ratio[1]*vec2.ref(i);
	return vec;
}

//Test if this, v2, and v3 are on a single straight line.
//Function's return value is true if the three points are on a straight,
//false if not.
bool MGVector::is_collinear(
	const MGVector& v2,
	const MGVector& v3
)const{
	return (v2 - *this).parallel(v3-*this);
}

// ベクトルの長さを返却する 
double MGVector::len() const {
	if(m_length < 0.){
		double a=0., b;
		for(size_t i=0; i<m_sdim; i++){
			b=m_element[i];
			a+=(b*b);
		}
		m_length=sqrt(a);
	}
	return m_length;
}

// 自身のベクトルを単位ベクトル化しオブジェクトを生成する
MGUnit_vector MGVector::normalize() const{
	return MGUnit_vector(*this);
}

// 自身のベクトルと与えられたベクトルが垂直かどうか返却する
// 垂直の時、True(1) を返却
bool MGVector::orthogonal(const MGVector &vec2) const{
	// *this と vec2 のcosΘを取得し、π/2 の時 垂直 （True(1)）
	return MGRight_angle(cangle(vec2));
}

//Compute the vector that is orthogonal to vec2 and is closest to this.
//"closest" means that the angle of the two vectors is minimum and
//the two vector length are equal.
MGVector MGVector::orthogonize(const MGVector& vec2)const{
	MGVector v212=vec2*(*this)*vec2;
	return MGUnit_vector(v212)*len();
}

// 自身のベクトルと与えられたベクトルが平行かどうか返却する
// 平行の時、True(1) を返却 
bool MGVector::parallel(const MGVector &vec2) const{
	// *this と vec2 のsinΘを取得し、零の時 平行
	return MGZero_angle(sangle(vec2));
}

//Reference to i-th element.
//double MGVector::ref(size_t i) const{ 
//	if(i<sdim()) return m_element[i];
//	else         return 0.;
//}

// 自身のベクトルをベクトル(v2)に射影したベクトルを求める。
// v2 が 零ベクトルのとき(*this)が返る。
//Project this onto the vector v2.
MGVector MGVector::project(const MGVector& v2) const{
	double v2ip = v2 % v2;
	if(MGMZero(v2ip)) return *this;//If v2 is zero vector.
	return ((*this%v2)/v2ip)*v2;
}

//Resize the vector, that is , change the space dimension.
//When this is enlarged, the extra space will contain garbages.
void MGVector::resize(size_t new_sdim){
	if(m_sdim==new_sdim) return;
	if(new_sdim<=3){
		if(m_sdim>3){
			for(size_t i=0; i<new_sdim; i++) m_data[i]=m_element[i];
			delete[] m_element; m_element=m_data;
		}
	}else{
		if(m_sdim<new_sdim){
			double* data=new double[new_sdim];
			for(size_t i=0; i<m_sdim; i++) data[i]=m_element[i];
			if(m_sdim>3) delete[] m_element;
			m_element=data;
		}
	}
	m_sdim=new_sdim; m_length=-1.;
}

// 自身のベクトルと与えられたベクトルのなす角度の sinΘ の絶対値を返却する
double MGVector::sangle(const MGVector &vec2) const{
	double sin_theta;
	// 自身か与えられたベクトルが零ベクトルの場合 sinΘを0.0 にする
	double len12=len()*vec2.len();
    if(MGMZero(len12)) sin_theta = 0.0;
    else{
		if(sdim()<=3 && vec2.sdim()<=3){
			sin_theta = ((*this)*vec2).len()/len12;
		} else {
			sin_theta=sin(angle(vec2));
		}
	}
	return sin_theta;
}

//Set this as a null vector.
void MGVector::set_null(){
	if(m_sdim>3) delete[] m_element;
	m_element=m_data;
	m_sdim=0;
	m_length=-1;
}

//Change this to a unit vector.
void MGVector::set_unit(){
	size_t dim=sdim();
	if(!dim){(*this)=MGUnit_vector(); return;}
	double length=len();
	// 自身のベクトルが零ベクトルの時はデフォルトベクトルを生成
	size_t i;
	if(MGMZero(length)){
		for(i=0; i<dim-1; i++) m_element[i]=0.;
		m_element[dim-1]=1.;
	}else {
	// 零ベクトル以外は与えられたベクトルの成分を長さで割る
		for(i=0; i<dim; i++) m_element[i]=m_element[i]/length;
	}
	m_length=1.;
}

//Store vec2 data into *this.
void MGVector::store_at(
	size_t i,				//Displacement of *this.
	const MGVector& vec2,	//Vector 2.
	size_t j)				//Displacement of vec2.
{
	size_t len=vec2.sdim();
	store_at(i,vec2,j,len);
}

//Store vec2 data into *this.
void MGVector::store_at(
	size_t i,				//Displacement of *this.
	const MGVector& vec2,	//Vector 2.
	size_t j,				//Displacement of vec2.
	size_t len)				//Length to store 
{
	size_t n1=sdim(), n2=vec2.sdim();
	if(len>n1) len=n1;
	for(size_t n=0; n<len; n++){
		if(i>=n1) i=0; if(j>=n2) j=0;
		m_element[i++]=vec2.ref(j++);
	}
	m_length=-1.;
}

//swap the coordinates.
//swap coordinates (i) and (j).
void MGVector::swap(size_t i, size_t j){
	assert(i<sdim() && j<sdim());
	double x=m_element[i];
	m_element[i]=m_element[j]; m_element[j]=x;
}

// ベクトルが単位ベクトルかどうか返却する。単位ベクトルの時、True(1) を返却
bool MGVector::is_unit_vector() const {
	double length=len();
	// 長さが Tolerance を考慮した 1.0 の時 true;
	return MGREqual(length, 1.);
}

// 自身のベクトルが零ベクトルかどうか返却する
bool MGVector::is_zero_vector() const{ return MGAZero(len()); }

// ３つのベクトルから求められる行列の行列式の値を返却する
double MGDeterminant(const MGVector& v1, const MGVector& v2,
                     const MGVector& v3 ) 
{
    return v1.ref(0)*v2.ref(1)*v3.ref(2) - v3.ref(0)*v2.ref(1)*v1.ref(2)
		+  v3.ref(0)*v1.ref(1)*v2.ref(2) - v1.ref(0)*v3.ref(1)*v2.ref(2)
		+  v2.ref(0)*v3.ref(1)*v1.ref(2) - v2.ref(0)*v1.ref(1)*v3.ref(2);
}

// ベクトルのスカラーの乗算を行いオブジェクトを生成
MGVector operator *(double scal, const MGVector& vec) {
	return vec*scal;
}

// V1をベクトル(v2)に射影したベクトルを求める。
// v2 が 零ベクトルのときV1が返る。
MGVector project(const MGVector& V1, const MGVector& V2){
	return V1.project(V2);
}
