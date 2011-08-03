/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/NDDArray.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/SPointSeq.h"
#include "mg/Vector.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Straight.h"
#include "mg/Plane.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/bkdtpg.h"
#include "cskernel/Bvavpl.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implementation of Class MGSPointSeq

//Constructor
MGSPointSeq::MGSPointSeq(size_t sizeu, size_t sizev, size_t dim)
:m_capacityu(sizeu),m_capacityv(sizev),m_sdim(dim)
,m_lengthu(sizeu),m_lengthv(sizev){
	size_t len=sizeu*sizev*dim;
	if(len) m_spoint=new double[len];
	else{
		m_capacityu=m_capacityv=m_lengthu=m_lengthv=0;
		m_spoint=0;
	}
}

MGSPointSeq::MGSPointSeq(
		size_t dim,				// New Space Dimension.
		const MGSPointSeq& old,	// Origianl SpointSeq.
		size_t start1,			// Destination start order to store.
		size_t start2			// Source start order to retrieve.
		)
	//Construct a MGSPointSeq by copying original MGSPointSeq.
	//Can change the order of coordinates 
	: m_capacityu(old.m_lengthu),m_capacityv(old.m_lengthv)
	, m_sdim(dim)
	, m_lengthu(old.m_lengthu), m_lengthv(old.m_lengthv){
	assert(start1<dim && start2<old.m_sdim);

	size_t len=m_lengthu*m_lengthv*m_sdim;
	if(len){
		m_spoint=new double[len];
		size_t dim2=old.m_sdim; 
		size_t dimmin= dim<dim2 ? dim:dim2; 
		size_t k, k2=start2, k1=start1;
		for(k=0; k<dimmin; k++) {
			for(size_t i=0; i<m_lengthu; i++)
				for(size_t j=0; j<m_lengthv; j++) (*this)(i,j,k1)=old(i,j,k2);
			k1 +=1; if(k1>=dim) k1=0; 
			k2 +=1; if(k2>=dim2) k2=0;
		}
		while(k++<dim){
			for(size_t i=0; i<m_lengthu; i++)
				for(size_t j=0; j<m_lengthv; j++) (*this)(i,j,k1)=0.0;
			k1 +=1; if(k1>=dim) k1=0; 
		}
	}else m_spoint=0;
}

//Copy constructor.
MGSPointSeq::MGSPointSeq(const MGSPointSeq& rhs)
: m_capacityu(rhs.m_capacityu),m_capacityv(rhs.m_capacityv)
, m_sdim(rhs.m_sdim)
, m_lengthu(rhs.m_lengthu), m_lengthv(rhs.m_lengthv)
{
	m_spoint=new double[m_capacityu*m_capacityv*m_sdim];
	for(size_t k=0; k<m_sdim; k++)
		for(size_t i=0; i<m_lengthu; i++)
			for(size_t j=0; j<m_lengthv; j++)
				(*this)(i,j,k)=rhs.ref(i,j,k);
}

//Destructor

//Member Function

//compute an average plane of the point sequence.
//Function's return value is:
// 0: Number of points included is zero.
// 1: Point seq is a point.		2: Point seq is on a line.
// 3: Plane is output.
int MGSPointSeq::average_plane(
	MGPosition& center		// center of point seq will be output.
	, MGPlane& plane		// Plane will be output, when average_plane=3.
	, MGStraight& line		// Straight line will be output            =2.
	, double& deviation)	// Maximum deviation from point, line
	 const{					// , or plane will be output.
	const double error=MGTolerance::wc_zero();	//error allowed.
	size_t nu=length_u(), nv=length_v(), ncd=sdim(), ipu=capacity_u(), ipv=capacity_v();
	int kplane=0; double g[4], fcenter[3];
	if(ncd){
		bvavpl_(error,nu,nv,data(),ncd,ipu,ipv,&kplane,g,fcenter,&deviation);
		MGVector direction(3,g);
		center=MGPosition(ncd, fcenter);
		switch (kplane){
		case 1:	break;	//When Point.
		case 2:			// Straight line
			line=MGStraight(MGSTRAIGHT_UNLIMIT, direction, center);
			break;
		case 3:			//Plane will be output.
			plane=MGPlane(direction, g[3]);
			break;
		default:
			break;
		}
	}
	return kplane;
}

//Compute minimum box sorrounding the all the points.
MGBox MGSPointSeq::box()const{
	size_t dim=sdim(); MGInterval* intrvl=new MGInterval[dim];
	size_t lenu=length_u(), lenv=length_v();
	for(size_t k=0; k<dim; k++){
		double max,min; max=min=ref(0,0,k);
		for(size_t j=0; j<lenv; j++){
			for(size_t i=0; i<lenu; i++){
				double dijk=ref(i,j,k);
				max= max>=dijk ? max:dijk;
				min= min<=dijk ? min:dijk;
			}
		}
		intrvl[k]=MGInterval(min,max);
	}
	MGBox ubox;
	if(dim) ubox=MGBox(dim, intrvl);
	delete[] intrvl;
	return ubox;
}

//Compute minimum box sorrounding the all the points.
MGBox* MGSPointSeq::compute_box()const{
	size_t dim=sdim(); MGInterval* intrvl=new MGInterval[dim];
	size_t lenu=length_u(), lenv=length_v();
	for(size_t k=0; k<dim; k++){
		double max,min; max=min=ref(0,0,k);
		for(size_t j=0; j<lenv; j++){
			for(size_t i=0; i<lenu; i++){
				double dijk=ref(i,j,k);
				max= max>=dijk ? max:dijk;
				min= min<=dijk ? min:dijk;
			}
		}
		intrvl[k]=MGInterval(min,max);
	}
	MGBox* ubox;
	if(dim) ubox=new MGBox(dim, intrvl);
	else ubox=new MGBox();
	delete[] intrvl;
	return ubox;
}

//Returns a pointer to the area.
const double* MGSPointSeq::data(size_t i, size_t j, size_t k) const
{return &m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];}

//Returns a pointer to the area.
double* MGSPointSeq::data(size_t i, size_t j, size_t k)
{return &m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];}

//Transformation for rational(MGRLBRep) Control Polygon.
//When rational control polygon, coordinates are of homogeneous form,
//i.e., the last space dimension id is for weights,
//and other elements include weight multiplied.
MGSPointSeq& MGSPointSeq::homogeneous_transform(double scale)
//1. Scaling.
{
	--m_sdim; operator*=(scale); ++m_sdim;
	return *this;
}
//2. Vector addition.
MGSPointSeq& MGSPointSeq::homogeneous_transform(const MGVector& vec)
{
	size_t i,j,r;
	size_t dim1=sdim()-1, dim2=vec.sdim();
	size_t dim= dim1>=dim2 ? dim1:dim2;
	double* v=new double[dim];
	for(j=0; j<dim; j++) v[j]=vec.ref(j);
	size_t lenu=length_u(),lenv=length_v();
	double weight;
	if(dim>dim1){
		MGSPointSeq bnew(lenu,lenv, dim+1);
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			weight=ref(i,j,dim1);
			for(r=0; r<dim1; r++) bnew(i,j,r)=ref(i,j,r)+v[r]*weight;
			for(r=dim1; r<dim; r++) bnew(i,j,r)=v[r]*weight;
			bnew(i,j,dim)=weight;
		}
		}
		*this=bnew;
	}
	else{
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			weight=ref(i,j,dim1);
			for(r=0; r<dim2; r++) (*this)(i,j,r)+=v[r]*weight;
		}
		}
	}
	delete[] v;
	return *this;
}
//3. Matrix multiplication.
MGSPointSeq& MGSPointSeq::homogeneous_transform(const MGMatrix& mat){
	size_t dim1=sdim()-1; size_t dim2=mat.sdim();
	size_t dim= dim1>=dim2 ? dim1:dim2;
	size_t i,j,r1,r2;	double a;
	size_t lenu=length_u(),lenv=length_v();
	if(dim>dim1){
		MGSPointSeq bnew(lenu,lenv, dim+1);
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(r1=0; r1<dim; r1++){
				a=0.;
				for(r2=0; r2<dim1; r2++) a+=ref(i,j,r2)*mat.ref(r2,r1);
				bnew(i,j,r1)=a;
			}
			bnew(i,j,dim)=ref(i,j,dim1);
		}
		}
		*this=bnew;
	}else{
		double* v=new double[dim];
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(r1=0; r1<dim; r1++){
				a=0.;
				for(r2=0; r2<dim; r2++) a+=ref(i,j,r2)*mat.ref(r2,r1);
				v[r1]=a;
			}
			for(r1=0; r1<dim; r1++) (*this)(i,j,r1)=v[r1];
		}
		}
		delete[] v;
	}

	return *this;
}
//4. Transformation multiplication.
MGSPointSeq& MGSPointSeq::homogeneous_transform(const MGTransf& tr)
{
	homogeneous_transform(tr.affine());
	homogeneous_transform(tr.translation());
	return *this;
}
//Generate data point abscissa from data point ordinates SPointSeq.
//SPointSeq(*this) must be homogeneous, i.e. if SPointSeq is positional
//data, all of them must be positional data, should not include
//derivative data.
void MGSPointSeq::make_data_point(MGNDDArray& utau, MGNDDArray& vtau) const
{
	size_t i,j,k, lenu,lenv, ipu,ipv,ipuv, mid; 
	length(lenu, lenv); const size_t ncd=sdim();
	if((ncd*lenu*lenv)==0){
		utau=vtau=MGNDDArray();
		return;
	}
	double d1,d2, d;
 	capacity(ipu, ipv); ipuv=ipu*ipv;
	MGNDDArray tau1(lenu), tau2(lenu), tau3(lenu);

	//Compute u data points.
	mid=lenv/2;
	bkdtpg_(data(0,0,0), lenu, ncd, ipuv, &tau1(0));
	bkdtpg_(data(0,mid,0), lenu, ncd, ipuv, &tau2(0));
	bkdtpg_(data(0,lenv-1,0), lenu, ncd, ipuv, &tau3(0));
	utau.reshape(lenu);
	for(i=0; i<lenu; i++) utau(i)=(tau1(i)+tau2(i)+tau3(i))/3.;
	utau.set_length(lenu);

	//Compute v data points and knot vector.
	vtau.reshape(lenv);
	mid=lenu/2;
	size_t id[3]={0, mid, lenu-1};
	vtau(0)=0.;
	for(j=1; j<lenv; j++){
		d2=0.0;
		for(size_t m=0; m<3; m++){
			d1=0.;
			for(k=0; k<ncd; k++){
				d = ref(id[m],j,k)-ref(id[m],j-1,k);
				d1+=d*d;
			}
			d2 += sqrt(d1);
		}
		vtau(j)=vtau(j-1)+d2/3.;
	}
	vtau.set_length(lenv);
}

//Compute non_homogeneous coordonate data without w coordinate element,
//assumed that this is homogeneous coordinate data,
//i.e., maximum space dimension element is w(weight) coordinate.
//Result data does not include weight elements.
MGSPointSeq MGSPointSeq::non_homogeneous() const{
	size_t i,j,r, sd=sdim()-1, m=length_u(), n=length_v();
	MGSPointSeq cp(m,n,sd);
	for(i=0; i<m; i++){
	for(j=0; j<n; j++){
		double weight=ref(i,j,sd);
		for(r=0; r<sd; r++) cp(i,j,r)=ref(i,j,r)/weight;
	}
	}
	return cp;
}

double MGSPointSeq::ref(size_t i, size_t j,size_t k) const{
//Return (i,j,k)-th element data.
// When k>=sdim(), return 0.0  .
	assert(i<length_u() && j<length_v());
	if(k>=sdim()) return 0.;
	else          return m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];
}

void MGSPointSeq::reshape(size_t sizeu, size_t sizev
						 , size_t startu, size_t startv)
//Change size. Change of sdim not allowed.
{
	assert(startu<=sizeu && startv<=sizev);

	size_t nlenu=m_lengthu+startu; size_t nlenv=m_lengthv+startv;
	if(sizeu<nlenu) nlenu=sizeu; if(sizev<nlenv) nlenv=sizev;
	size_t nu=nlenu-startu, nv=nlenv-startv;
	//nu,nv are the numbers of data to move.
	
	//Reshape of m_spoint.
	size_t i,j,k;
	size_t sizeuv=sizeu*sizev;
	double* data=new double[sizeuv*m_sdim];
	for(k=0; k<m_sdim; k++){
		size_t sizeuvk=sizeuv*k;
		for(j=0; j<nv; j++){
			size_t js=startu+(j+startv)*sizeu+sizeuvk;
			for(i=0; i<nu; i++) data[i+js]=ref(i,j,k);
		}
	}
	delete[] m_spoint;
	m_spoint=data;
	m_capacityu=sizeu; m_capacityv=sizev; m_lengthu=nlenu; m_lengthv=nlenv;
}
	
//Change the size of the array to
//m_lengthu=m_capacityu=lenu, m_lengthv=m_capacityv=lenv, m_sdim=dim.
void MGSPointSeq::resize(size_t lenu, size_t lenv,  size_t dim){
	if(m_spoint) delete[] m_spoint;
	m_spoint=new double[lenu*lenv*dim];
	m_capacityu=m_lengthu=lenu;
	m_capacityv=m_lengthv=lenv;
	m_sdim=dim;
}

//Reverse the ordering of points.
void MGSPointSeq::reverse(
		int is_u)				//if true, u-drection. if not, v.
{
	size_t i,i2,j,j2,k; double save;
	size_t ncd=sdim(), lenu=length_u(), lenv=length_v();
	size_t half;
	if(is_u){
		half=lenu/2;
		for(k=0; k<ncd; k++){
			for(j=0; j<lenv; j++){
				i2=lenu-1;
				for(i=0; i<half; i++){
					save=ref(i2,j,k);
					(*this)(i2--,j,k)=ref(i,j,k);
					(*this)(i,j,k)=save;
				}
			}
		}
	}else{
		half=lenv/2;
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				j2=lenv-1;
				for(j=0; j<half; j++){
					save=ref(i,j2,k);
					(*this)(i,j2--,k)=ref(i,j,k);
					(*this)(i,j,k)=save;
				}
			}
		}
	}
}

//Set the length of effective data.
void MGSPointSeq::set_length(size_t lengthu,
					size_t lengthv )	
{	assert(lengthu<=m_capacityu && lengthv<=m_capacityv);
	m_lengthu=lengthu; m_lengthv=lengthv;
}

//Set this as a null.
void MGSPointSeq::set_null(){
	if(m_spoint) delete[] m_spoint;
	m_spoint=0;
	m_lengthu=m_lengthv=m_sdim=m_capacityu=m_capacityv=0;
}

//Store vector data vector(from+r) to this(i,j,to+r) for 0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_at(
	size_t i, size_t j,	//id of this which indicates the placement.
	const MGVector& vector,	//Input vector.
	size_t to,		//Indicates to where of this in the space dimension id.
	size_t from)	//Indicates from where of vector in the space dimension.
{
	assert(i<m_lengthu && j<m_lengthv);
	size_t sd1=sdim(), sd2=vector.sdim();
	size_t len=sd1; if(len>sd2) len=sd2;
	for(size_t r=0; r<len; r++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*(to++)]=vector.ref(from++);
	}
}

//Store vector data vector(from+r) to this(i,j,to+r) for 0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_at(
	size_t i, size_t j,	//id of this which indicates the placement.
	const MGVector& vector,	//Input vector.
	size_t to,		//Indicates to where of this in the space dimension id.
	size_t from,	//Indicates from where of vector in the space dimension.
	size_t len)		//Length of data to store.
{
	assert(i<m_lengthu && j<m_lengthv);
	size_t sd1=sdim(), sd2=vector.sdim();
	if(len>sd1) len=sd1;
	for(size_t r=0; r<len; r++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*(to++)]=vector.ref(from++);
	}
}

//Store data[r] to this(i,j,to+r)  0<=r<sdim().
//When (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_at(
	size_t i, size_t j,		//id of this which indicates the placement.
	const double* data,		//Input data.
	size_t to)		//Indicates to where of this in the space dimension id.
{
	assert(i<m_lengthu && j<m_lengthv);
	size_t sd1=sdim();
	for(size_t r=0; r<m_sdim; r++){
		if(to>=sd1) to=0;
		m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*(to++)]=data[r];
	}
}

//Store BPointSeq along u at v's id j. That is, 
//data bp(i,from+r) to this(i,j,to+r) for i=0,..., length_u()-1, and  0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_BP_along_u_at(
	size_t j,	//id of this which indicates the placement.
	const MGBPointSeq& bp,	//Input vector.
	size_t to,	//Indicates to where of this in the space dimension id.
	size_t from)	//Indicates from where of vector in the space dimension.
{
	assert(j<m_lengthv);
	size_t sd1=sdim(), sd2=bp.sdim();
	size_t len=sd1; if(len>sd2) len=sd2;
	size_t szubyj=m_capacityu*j, szubyszv=m_capacityu*m_capacityv;
	for(size_t r=0; r<len; r++, to++, from++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		size_t juvto=szubyj+szubyszv*to;
		for(size_t i=0; i<m_lengthu; i++) m_spoint[i+juvto]=bp(i,from);
	}
}

//Store BPointSeq along v at u's id i. That is, 
//data bp(j,from+r) to this(i,j,to+r) for j=0,..., length_v()-1, and  0<=r<sdim().
//When (form+r) or (to+r) reached to maximum space dimension id, next id
//becomes 0(form the start).
void MGSPointSeq::store_BP_along_v_at(
	size_t i,	//id of this which indicates the placement.
	const MGBPointSeq& bp,	//Input vector.
	size_t to,	//Indicates to where of this in the space dimension id.
	size_t from)	//Indicates from where of vector in the space dimension.
{
	assert(i<m_lengthu);
	size_t sd1=sdim(), sd2=bp.sdim();
	size_t len=sd1; if(len>sd2) len=sd2;
	size_t szubyszv=m_capacityu*m_capacityv;
	for(size_t r=0; r<len; r++, to++, from++){
		if(to>=sd1) to=0; if(from>=sd2) from=0;
		size_t iuvto=i+szubyszv*to;
		for(size_t j=0; j<m_lengthv; j++) m_spoint[iuvto+m_capacityu*j]=bp(j,from);
	}
}

//Operator Definition

double& MGSPointSeq::operator()(size_t i, size_t j, size_t k) 
//Access to (i,j)th element
{
	assert(i<m_capacityu && j<m_capacityv && k<m_sdim);
	return m_spoint[i+m_capacityu*j+m_capacityu*m_capacityv*k];
}

MGVector MGSPointSeq::operator()(size_t i, size_t j) const
//Extract (i,j,k) elements for 0<=k<sdim() as a vector.
{
	MGVector v(m_sdim);
	for(size_t k=0; k<m_sdim; k++) v(k)=ref(i,j,k);
	return v;
}

// 曲線の平行移動を行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator + ( const MGVector& vec) const{
	MGSPointSeq bnew(*this);
	return bnew+=vec;
}

// 与ベクトルだけ曲線を平行移動して自身とする。
MGSPointSeq& MGSPointSeq::operator+= (const MGVector& vec){
	size_t i,j,k; size_t dim1=sdim(), dim2=vec.sdim();
	size_t dim= dim1>=dim2 ? dim1:dim2;
	double* v=new double[dim];
	for(k=0; k<dim; k++) v[k]=vec.ref(k);
	size_t lenu=length_u(), lenv=length_v();
	if(dim>dim1){
		MGSPointSeq bnew(lenu,lenv, dim);
		for(k=0; k<dim; k++)
			for(i=0; i<lenu; i++) 
				for(j=0; j<lenv; j++) bnew(i,j,k)=ref(i,j,k)+v[k];
		bnew.set_length(lenu,lenv);
		*this=bnew;
	}else{
		for(k=0; k<dim; k++)
			for(i=0; i<lenu; i++)  
				for(j=0; j<lenv; j++) (*this)(i,j,k)+=v[k];
	}
	delete[] v;
	return *this;
}

// 曲線の逆方向に平行移動を行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator- (const MGVector& vec) const{
	MGVector v=-vec;
	return (*this)+v;
}

// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
MGSPointSeq& MGSPointSeq::operator-= (const MGVector& vec){
	MGVector v=-vec;
	return *this += v;
}

// 与えられたスケーリングで曲線の変換を行いオブジェクトを生成する。
//Scaling.
MGSPointSeq MGSPointSeq::operator* (double scale) const{
	MGSPointSeq sp(*this);
	sp*=scale;
	return sp;
}

// 与えられたスケーリングで曲線の変換を行いオブジェクトを生成する。
//Scaling.
MGSPointSeq operator* (double scale, const MGSPointSeq& sp){
	return sp*scale;
}

// 与えられたスケーリングで曲線の変換を行い自身の曲線とする。
//Scaling.
MGSPointSeq& MGSPointSeq::operator*= (double scale){
	for(size_t k=0; k<sdim(); k++)
		for(size_t i=0; i<length_u(); i++)  
			for(size_t j=0; j<length_v(); j++) (*this)(i,j,k)*=scale;
	return *this;
}

// 与えられた変換で曲線の変換を行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator* (const MGMatrix& mat) const{
	MGSPointSeq bnew(*this);
	return bnew*=mat;
}

// 与えられた変換で曲線の変換を行い自身の曲線とする。
MGSPointSeq& MGSPointSeq::operator*= (const MGMatrix& mat){
	size_t dim1=sdim(); size_t dim2=mat.sdim();
	size_t dim= dim1>=dim2 ? dim1:dim2;
	size_t i,j,k,k2;	double a;
	size_t lenu=length_u(), lenv=length_v();
	if(dim>dim1){
		MGSPointSeq bnew(lenu, lenv, dim);
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(k=0; k<dim; k++){
				a=0.;
				for(k2=0; k2<dim; k2++) a=a+ref(i,j,k2)*mat.ref(k2,k);
				bnew(i,j,k)=a;
			}
		}
		}
		bnew.set_length(lenu,lenv);
		*this=bnew;
	}else{
		double* v=new double[dim];
		for(i=0; i<lenu; i++){
		for(j=0; j<lenv; j++){
			for(k=0; k<dim; k++){
				a=0.;
				for(k2=0; k2<dim; k2++) a=a+ref(i,j,k2)*mat.ref(k2,k);
				v[k]=a;
			}
			for(k=0; k<dim; k++) (*this)(i,j,k)=v[k];
		}
		}
		delete[] v;
	}

	return *this;
}

// 与えられた変換で曲線のトランスフォームを行いオブジェクトを生成する。
MGSPointSeq MGSPointSeq::operator* (const MGTransf& tr) const{
	return (*this)*tr.affine()+tr.translation();
}

// 与えられた変換で曲線のトランスフォームを行い自身とする。
MGSPointSeq& MGSPointSeq::operator*= (const MGTransf& tr){
	return ((*this)*=tr.affine())+=tr.translation();
}

//Assignment
MGSPointSeq& MGSPointSeq::operator=(const MGSPointSeq& sp2){
	size_t lenu=sp2.length_u(), lenv=sp2.length_v(), dim=sp2.sdim();
	size_t len=capacity_u()*capacity_v()*sdim();
	if(len<lenu*lenv*dim) resize(lenu,lenv,dim);
	else{
		m_sdim=dim;
		m_lengthu=m_capacityu=lenu; m_lengthv=lenv;
		m_capacityv=len/(dim*lenu);
	}
	for(size_t m=0; m<m_sdim; m++)
		for(size_t i=0; i<lenu; i++)
			for(size_t j=0; j<lenv; j++)
				(*this)(i,j,m)=sp2(i,j,m);
	return *this;
}

//Compare two SPointSeq if they are equal.
bool MGSPointSeq::operator== (const MGSPointSeq& spoint) const{
	size_t lenu=length_u(), lenv=length_v();
	if(lenu!=spoint.length_u() || lenv!=spoint.length_v()) return 0;
	if(lenu<=0 && lenv<=0) return 1;

	size_t dim=sdim(), dim2=spoint.sdim();
	if(dim<dim2) dim=dim2;
	size_t i,j,k; double a,b; 
	double error=MGTolerance::wc_zero(); error=error*error;
	for(j=0; j<lenv; j++){
	for(i=0; i<lenu; i++){
		a=0.;
		for(k=0; k<dim; k++){
			b=ref(i,j,k)-spoint.ref(i,j,k); a += b*b;
		}
		if(a>error) return 0;
	}
	}
	return 1;
}

//Add and subtract operation of two MGSPointSeq.
MGSPointSeq& MGSPointSeq::operator+= (const MGSPointSeq& sp2){
	size_t sd=sdim(), sd2=sp2.sdim();
	if(sd2>sd) sd=sd2;
	size_t nu1=length_u(), nu2=sp2.length_u();
	size_t nu=nu1;	if(nu2>nu) nu=nu2;
	size_t nv1=length_v(), nv2=sp2.length_v();
	size_t nv=nv1;	if(nv2>nv) nv=nv2;
	size_t i,j;
	const MGVector zero(0.,0.,0.);
	if(sd>sdim()){
		MGSPointSeq sp(nu,nv,sd);
		for(i=0; i<nu1; i++)
			for(j=0; j<nv1; j++)
				sp.store_at(i,j,(*this)(i,j));//copy
		for(i=nu1; i<nu; i++)
			for(j=nv1; j<nv; j++)
				sp.store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++)
			for(j=0; j<nv2; j++)
				sp.store_at(i,j,sp(i,j)+sp2(i,j));
		(*this)=sp;
	}else{
		reshape(nu,nv);
		for(i=nu1; i<nu; i++)
			for(j=nv1; j<nv; j++)
				store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++)
			for(j=0; j<nv2; j++)
				store_at(i,j,(*this)(i,j)+sp2(i,j));
	}
	return *this;
}
MGSPointSeq MGSPointSeq::operator+ (const MGSPointSeq& sp2)const{
	MGSPointSeq sp(*this);
	return sp+=sp2;
}
MGSPointSeq MGSPointSeq::operator- (const MGSPointSeq& sp2) const{
	MGSPointSeq sp(*this);
	return sp-=sp2;
}
MGSPointSeq& MGSPointSeq::operator-= (const MGSPointSeq& sp2){
	size_t sd=sdim(), sd2=sp2.sdim();
	if(sd2>sd) sd=sd2;
	size_t nu1=length_u(), nu2=sp2.length_u();
	size_t nu=nu1;	if(nu2>nu) nu=nu2;
	size_t nv1=length_v(), nv2=sp2.length_v();
	size_t nv=nv1;	if(nv2>nv) nv=nv2;
	size_t i,j;
	const MGVector zero(0.,0.,0.);
	if(sd>sdim()){
		MGSPointSeq sp(nu,nv,sd);
		for(i=0; i<nu1; i++) for(j=0; j<nv1; j++) sp.store_at(i,j,(*this)(i,j));//copy
		for(i=nu1; i<nu; i++) for(j=nv1; j<nv; j++) sp.store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++) for(j=0; j<nv2; j++) sp.store_at(i,j,sp(i,j)-sp2(i,j));
		(*this)=sp;
	}else{
		reshape(nu,nv);
		for(i=nu1; i<nu; i++) for(j=nv1; j<nv; j++) store_at(i,j,zero);//Clear
		for(i=0; i<nu2; i++) for(j=0; j<nv2; j++) store_at(i,j,(*this)(i,j)-sp2(i,j));
	}
	return *this;
}
