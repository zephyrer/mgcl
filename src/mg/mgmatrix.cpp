/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Vector.h"
#include "mg/Unit_vector.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGMatrix.cc
// Implementatation of MGMatrix
//

//
// Constructor.

//void constructor:void コンストラクタ
MGMatrix::MGMatrix(size_t sdim)
:m_sdim(sdim){
	if(sdim) m_matrix=new double[sdim*sdim];
	else m_matrix=0;
}

//  Construct 2D matrix from two 2D vectors.
MGMatrix::MGMatrix(const MGVector& vec1, const MGVector& vec2)
	:m_sdim(2), m_matrix(new double[4]){
	(*this)(0,0)=vec1.ref(0); (*this)(0,1)=vec1.ref(1);
	(*this)(1,0)=vec2.ref(0); (*this)(1,1)=vec2.ref(1);
}

//  ３つの与えられた行ベクトルから3D Matrixを生成する。
MGMatrix::MGMatrix(
	const MGVector& v1, const MGVector& v2, const MGVector& v3)
	:m_sdim(3), m_matrix(new double[9]){
	(*this)(0,0)=v1.ref(0); (*this)(0,1)=v1.ref(1); (*this)(0,2)=v1.ref(2);
	(*this)(1,0)=v2.ref(0); (*this)(1,1)=v2.ref(1);	(*this)(1,2)=v2.ref(2);
	(*this)(2,0)=v3.ref(0); (*this)(2,1)=v3.ref(1);	(*this)(2,2)=v3.ref(2);
}

//  各軸方向で等しい Scaling のためのMatrixを生成。
MGMatrix::MGMatrix(size_t dim, double scale):m_sdim(dim){
	if(dim){
		m_matrix=new double[dim*dim];
		for(size_t i=0; i<dim; i++){
			for(size_t j=0; j<dim; j++){
				if(j==i) (*this)(i,j)=scale;
				else     (*this)(i,j)=0.;
			}
		}
	}else m_matrix=0;
}

//Construct dim dimension matrix from the array of double values.
MGMatrix::MGMatrix(
	size_t dim,			//dimension of the matrix.
	const double* values,//array of values[dim*dim].
	bool column_wise	//If column_wise=true, (*this)(i,j)=values[i+dim*j],
						//else(row_wise), (*this)(i,j)=values[j+dim*i], for 0<=i,j<=dim-1.
):m_sdim(dim), m_matrix(new double[dim*dim]){	
	if(column_wise){
		for(size_t i=0; i<dim; i++){
			for(size_t j=0; j<dim; j++){
				(*this)(i,j)=values[i+dim*j];
			}
		}
	}else{
		for(size_t i=0; i<dim; i++){
			for(size_t j=0; j<dim; j++){
				(*this)(i,j)=values[j+dim*i];
			}
		}
	}
}

//  与えられたVector回りに指定の角度回転させるMatrixを作成する。
//　Rotation matrix around vec. Space dimension of vec can be
//  any number, i.e. can be more than 3.
MGMatrix::MGMatrix(const MGVector& vec, double angle)
:m_sdim(0), m_matrix(0){
	to_axis(vec,0);
    // this = To transform vec to x-axis.

	MGMatrix temp; temp.set_rotate_2D(angle);
	MGMatrix m2=MGMatrix(3,temp,1,0);
	//m2= Rotation about x-axis.

	MGMatrix m3; m3.from_axis(vec,0);
	// m3 = To transform x-axis to vec.

	(*this)*=m2;
	(*this)*=m3;
}

//Construct a matrix to rotate in (axis1, axis2) plane.
MGMatrix::MGMatrix(size_t dim,	//Space dimension(can be more than 3).
	size_t axis1,		//axis number 1
	size_t axis2,		//axis number 2
	double cosv, double sinv) //cosv=cos(angle), sinv=sin(angel).
						//That is cosv*cosv+sinv*sinv must be 1.
:m_sdim(0), m_matrix(0){
	assert(axis1<dim && axis2<dim);
	//assert(MGREqual(1.,cosv*cosv+sinv*sinv));

	*this=MGMatrix(dim,1.);
	(*this)(axis1,axis1)=cosv; (*this)(axis1,axis2)=sinv;
	(*this)(axis2,axis1)=-sinv; (*this)(axis2,axis2)=cosv;
}

//Construct Matrix by copying old Matrix, changing space dimension and
//ordering of old coordinates.
MGMatrix::MGMatrix(size_t dim, const MGMatrix& mat2,
	size_t start1, size_t start2):m_sdim(dim){
	assert (dim>0 && start1<dim && start2<mat2.sdim());

	if(!dim){
		m_matrix=0;
		return;
	}

	m_matrix=new double[dim*dim];
	set_scale(1.0);
	size_t dim2=mat2.sdim(); 
	size_t i1,i2,j1,j2;
	size_t dimmin= dim<dim2 ? dim:dim2 ;
	i1=start1; i2=start2;
	for(size_t i=0; i<dimmin; i++){
		j1=start1; j2=start2;
		for(size_t j=0; j<dimmin; j++){
			(*this)(i1,j1)=mat2(i2,j2);
			j1 +=1; if(j1>=dim)  j1=0;
			j2 +=1; if(j2>=dim2) j2=0;
		}
		i1 +=1; if(i1>=dim)  i1=0;
		i2 +=1; if(i2>=dim2) i2=0;
	}
}

//Copy constructor.
MGMatrix::MGMatrix(const MGMatrix& mat)
:m_sdim(mat.m_sdim){
	size_t len=m_sdim*m_sdim;
	m_matrix=new double[len];
	for(size_t i=0; i<len; i++) m_matrix[i]=mat.m_matrix[i];
}

//////////// Destructor.////////

//
// メンバ関数
//

//Convert this transf matrix to OpenGL Matrix.
void MGMatrix::convert_to_glMatrix(
	double glMat[16]	//OpenGL Matrix will be output.
)const{
	for(int i=0; i<3; i++){
		int i4=i*4;
		for(int j=0; j<3; j++){
			glMat[i4+j]=(*this)(i,j);
		}
		glMat[i4+3]=0.;
	}
	glMat[12]=glMat[13]=glMat[14]=0.; glMat[15]=1.;
}

//  行列式の値を返却。
double MGMatrix::determinant() const{
	size_t dim=sdim();
	double value;
	switch (dim){
	case 0: value=0.; break;
	case 1: value=m_matrix[0]; break;
	case 2: value=m_matrix[0]*m_matrix[3]-m_matrix[1]*m_matrix[2]; break;
	case 3: value=
				(m_matrix[0]*m_matrix[4]*m_matrix[8]
				-m_matrix[0]*m_matrix[7]*m_matrix[5])+
				(m_matrix[3]*m_matrix[7]*m_matrix[2]
				-m_matrix[3]*m_matrix[1]*m_matrix[8])+
				(m_matrix[6]*m_matrix[1]*m_matrix[5]
				-m_matrix[6]*m_matrix[4]*m_matrix[2]); break;
	default: value=0.; double sign=1.;
			MGMatrix mat(dim-1);
			for(size_t m=0; m<dim; m++){
				for(size_t i=1; i<dim; i++){
					size_t im1=i-1;
					for(size_t j=0; j<dim; j++){
						if(j<m) mat(im1,j)=ref(i,j);
						else if(j>m) mat(im1,j-1)=ref(i,j);
					}
				}
				value+=mat.determinant()*sign*ref(0,m);
				sign=-sign;
			}; break;
	}
	return value;
}

//Construct a matrix to transform one of the axises to a vector 'uvec',
// and replace own matrix with it. Inverse matrix of to_axis.
// axis can be any number(can be more than 2).
MGMatrix& MGMatrix::from_axis(
	const MGUnit_vector& uvec,	// Unit vector to be an axis.
	size_t axis)				// Axis kind 0:x, 1:y, 2:z, ...
{
//This program produce matrix by inversing to_axis matrix,
//is based on to_axis.
	size_t dim=uvec.sdim(); if(axis>=dim) dim=axis+1;
	size_t i,j;
	size_t axis1=axis, axis2;
	MGUnit_vector unit=uvec;
	double len, len1, len2, max_len, cosv, sinv;
	MGMatrix mat(dim,1.), mat_inverse(dim,1.), m1,m2;
	//Compute matrix mat_inverse to rotate for axis1 element to be zero
	//from axis+1 to axis-1. from_axis is the inverse of mat_inverse.
	for(i=0; i<dim-1; i++){
		axis1++; if(axis1>=dim) axis1=0;
		//Get maximum length element id in axis2.
		axis2=0; max_len=fabs(unit.ref(0));
		for(j=1; j<dim; j++){
			len=fabs(unit.ref(j));
			if(len>max_len){ axis2=j; max_len=len;}
		}
		if(axis1==axis2) axis2=axis;

		len1=unit.ref(axis1); len2=unit.ref(axis2);
		len=sqrt(len1*len1+len2*len2);
		cosv=len1/len; sinv=len2/len;
		m1=MGMatrix(dim,axis1,axis2,sinv,cosv); mat_inverse*=m1;
		//Above mat_inverse is to rotate for axis1 element to be zero.
		m2=MGMatrix(dim,axis1,axis2,sinv,-cosv);
		mat=m2*mat;
		unit=uvec*mat_inverse;
	}
	return *this=mat;
}

// Multiply matrix m1 and m2.
MGMatrix MGMatrix::multiply(const MGMatrix& m2) const{
	size_t dim=sdim(), dim2=m2.sdim(),i,j,k;
	if(dim<dim2) dim=dim2;
	double a;
	MGMatrix m(dim);
	for(i=0; i<dim; i++){
		for(j=0; j<dim; j++){
			a=0.;
			for(k=0; k<dim; k++) a += ref(i,k)*m2.ref(k,j);
			m(i,j)=a;
		}
	}
	return m;
}

//Reference to (i,j)-th element of the matarix.
double MGMatrix::ref(size_t i, size_t j) const
{
	if(i<m_sdim && j<m_sdim) return m_matrix[i+j*m_sdim];
	else if(i != j)          return 0.;
	else                     return 1.;
}

//Resize, i.e., change, the space dimension.
void MGMatrix::resize(size_t nsdim){
	if(m_matrix) delete[] m_matrix;
	m_matrix=new double[nsdim*nsdim];
	m_sdim=nsdim;
}

//Construct a mirror reflection matrix about a plane whose normal
//is vec and that passes through the origin, then
//replace own matrix with it.
MGMatrix& MGMatrix::reflection(const MGVector& vec)
{
	to_axis(vec,0);
    // this = To transform vec to x-axis.

	for(size_t i=0; i<sdim(); i++) (*this)(i,0) *= -1.;
	//Reflection about x-axis.

	MGMatrix m2; m2.from_axis(vec,0);
	// m2 = To transform x-axis to vec.

	return *this *=m2;
}

//Obtain the scaling factor of this matrix.
double MGMatrix::scale()const{
	double a=0.;
	size_t m=sdim();
	if(!m) return a;
	for(size_t i=0; i<m; i++){
		for(size_t j=0; j<m; j++){
			double b=(*this)(i,j);
			a+=b*b;
		}
	}
	return sqrt(a);
}

// Construct 2D space Matrix to transform for 'unit' to be x-coordimate, and
// replace own Matrix.
MGMatrix& MGMatrix::set_x_axis
	(const MGUnit_vector& unit) //unit vector to be x-coordinate
{
	assert(unit.sdim()==2);

	if(m_sdim!=2) resize(2);
    double x=unit(0); double y=unit(1);
	(*this)(0,0)=x; (*this)(0,1)=-y;
	(*this)(1,0)=y; (*this)(1,1)=x;
	return *this;
}

//  原点を通り、指定Vectorに関して鏡面変換する 2D Matrixを作成し，
//  既存のMatrixと入れ換える。
MGMatrix& MGMatrix::set_reflect_2D(const MGVector& vec1)
{
     const MGUnit_vector uvec1(MGVector(2,vec1));
	 double x=uvec1(0); double y=uvec1(1); double twoxy=2.*x*y;
	 double x2my2=x*x-y*y;
	if(m_sdim!=2) resize(2);

	 (*this)(0,0)=x2my2; (*this)(0,1)=twoxy;
	 (*this)(1,0)=twoxy; (*this)(1,1)=-x2my2;

     return *this;
}

//（原点を基点とする）Vector V0 を V1に変換するMatrixを作成し、
//  自身のMatrixと入れ換える。
//Construct the matrix to rotate and scale that transform vector V0 to V1,
//and replace this matrix with the matrix.
//Space dimension of V0 and V1 can be any number greater than 1.
MGMatrix& MGMatrix::set_rotate(const MGVector& V0, const MGVector& V1){
	MGVector V2=V0*V1;
	double angle=V0.angle(V1);
	MGMatrix m(V2, angle);
	double l0=V0.len(), l1=V1.len();
	double s;
	if(MGMZero(l0)) s=1.;
	else s=l1/l0;
	m*=s;
	*this=m;
	return *this;
}

//  原点の回りに指定の角度回転させる2D Matrixを作成し,
//  既存の Matrixと入れ換える。
MGMatrix& MGMatrix::set_rotate_2D(double angle)
{
	double cosval=cos(angle); double sinval=sin(angle);
	if(m_sdim!=2) resize(2);
	(*this)(0,0)= cosval; (*this)(0,1)=sinval;
	(*this)(1,0)=-sinval; (*this)(1,1)=cosval;
     return *this;
}

//Rotation 2D matrix around origin by an angle.
//The angle is given by cval as cos(angle) and sval as sin(angle).
MGMatrix& MGMatrix::set_rotate_2D(double cval, double sval)
{
	if(m_sdim!=2) resize(2);
	(*this)(0,0)= cval; (*this)(0,1)=sval;
	(*this)(1,0)=-sval; (*this)(1,1)=cval;
     return *this;
}

//Construct a 3D matrix to transform a vector to be one of the axises,
// and replace own matrix.
MGMatrix& MGMatrix::set_axis(
	const MGUnit_vector& uvec,	//Unit vector to be an axis.
	size_t axis)				// Axis number 0:x, 1:y, 2:z.
{
	assert(axis<3);

	size_t i,j,k;
	i=axis; if(axis>=3) i=2;
	j=i+1; if(j>=3) j=0; k=j+1; if(k>=3) k=0;
	MGMatrix m1(3,1.), m2(3,1.); double a,d;
	
	double d2[3];
	for(size_t n=0; n<3; n++){
		a=uvec.ref(n); d2[n]=a*a;
	}

	double dij=d2[i]+d2[j]; double dki=d2[k]+d2[i];
	if(dij>=dki){
		d=sqrt(dij);
		// 1. Rotate around k-axis.
		m1(i,i)=uvec.ref(i)/d; m1(i,j)=-uvec.ref(j)/d;
		m1(j,i)=-m1(i,j);      m1(j,j)=m1(i,i);
		// 2. Rotate around j-axis.
		m2(i,i)=d;             m2(i,k)=-uvec.ref(k);
		m2(k,i)=uvec.ref(k);   m2(k,k)=m2(i,i);
	}else{
		d=sqrt(dki);
		// 1. Rotate around j-axis.
		m1(i,i)=uvec.ref(i)/d; m1(i,k)=-uvec.ref(k)/d;
		m1(k,i)=-m1(i,k);      m1(k,k)=m1(i,i);
		// 2. Rotate around k-axis.
		m2(i,i)=d;             m2(i,j)=-uvec.ref(j);
		m2(j,i)=uvec.ref(j);   m2(j,j)=m2(i,i);
	}
	(*this)=m1*m2;
	return *this;
}

//Construct a matrix to transform a unit vector on an axis
// to a vector 'uvec', and replace own matrix with it.
//Inverse matrix of set_axis.
MGMatrix& MGMatrix::set_vector(
	const MGUnit_vector& uvec,	// Unit vector to be an axis.
	size_t axis)				// Axis number 0:x, 1:y, 2:z.
{
	assert(uvec.sdim()<=3 && axis<3);

	size_t i,j,k;
	i=axis; j=i+1; if(j>=3) j=0; k=j+1; if(k>=3) k=0;
	MGMatrix m1(3,1.), m2(3,1.); double a,d;
	
	double d2[3];
	for(size_t n=0; n<3; n++){
		a=uvec(n); d2[n]=a*a;
	}

	double dij=d2[i]+d2[j]; double dki=d2[k]+d2[i];
	if(dij>=dki){
		d=sqrt(dij);
		// 1. Rotate around j-axis.
		m2(i,i)=d;        m2(i,k)=uvec(k);
		m2(k,i)=-uvec(k);  m2(k,k)=m2(i,i);
		// 2. Rotate around k-axis.
		m1(i,i)=uvec(i)/d; m1(i,j)=uvec(j)/d;
		m1(j,i)=-m1(i,j);  m1(j,j)=m1(i,i);
	}else{
		d=sqrt(dki);
		// 1. Rotate around k-axis.
		m2(i,i)=d;         m2(i,j)=uvec(j);
		m2(j,i)=-uvec(j);   m2(j,j)=m2(i,i);
		// 2. Rotate around j-axis.
		m1(i,i)=uvec(i)/d; m1(i,k)=uvec(k)/d;
		m1(k,i)=-m1(i,k);  m1(k,k)=m1(i,i);
	}
	(*this)=m2*m1;
	return *this;
}

//  与えられた２つの単位ベクトルを各々、X軸、Y軸にするよう原点の周りに
//  回転させる 3D Matrix を生成し，自身のMatrixと入れ替える。もし、
//  ２つ目の単位ベクトルが１つ目の単位ベクトルと直交しない場合は、２つのベ
//  クトルのかわりに両ベクトルを含む平面内で直交するよう変換したベクトルを
//  使用する。
MGMatrix& MGMatrix::set_xy_axis
	(const MGUnit_vector& uvecx,	//Unit vector 1 for x axis.
     const MGUnit_vector& uvecy)	//Unit vector 2 for y axis.
{
	assert(uvecx.sdim()<=3 && uvecy.sdim()<=3);

	double a,b;
	MGMatrix mx(3,1.), mat;

	mat.set_axis(uvecx, 0);
	// mat=the matrix to transform uvecx as x-axis.

	MGVector vecy=uvecy*mat;
	a=vecy[1]; b=vecy[2]; double dvecy=sqrt(a*a+b*b);
	if(dvecy>MGTolerance::mach_zero()){
		mx(1,1)=a/dvecy;  mx(1,2)=-b/dvecy;
		mx(2,1)=-mx(1,2); mx(2,2)=mx(1,1);
		//mx= the matrix to rotate around x-axis for uvecy*mat to be y-axis.
		(*this)=mat*mx;
	}else (*this)=mat;

	return *this;
}

//  X軸、Y軸を各々、与えられた２つの単位ベクトルにするよう原点の周りに
//  回転させる 3D Matrix を生成し，自身のMatrixと入れ替える。もし、
//  ２つ目の単位ベクトルが１つ目の単位ベクトルと直交しない場合は、２つのベ
//  クトルのかわりに両ベクトルを含む平面内で直交するよう変換したベクトルを
//  使用する。
// This is the inverse matrix of set_xy_axis().
MGMatrix& MGMatrix::set_xy_vector
	(const MGUnit_vector& uvecx,	//Unit vector 1 for x axis.
     const MGUnit_vector& uvecy)	//Unit vector 2 for y axis.
{
	assert(uvecx.sdim()<=3 && uvecy.sdim()<=3);

	MGVector ay(3);
	MGMatrix mat1; mat1.set_vector(uvecx, 0);
	// mat1=the matrix to transform x-axis as uvecx.
	ay(0)=mat1(1,0); ay(1)=mat1(1,1); ay(2)=mat1(1,2);
	//ay is the transformed vector of y-axis by matrix mat1.

	MGUnit_vector vy=(uvecx*uvecy)*uvecx;//vy is normalized uvecy.
	double cval=ay.cangle(vy); double sval=ay.sangle(vy);
	if((ay*vy)%uvecx < 0.) sval=-sval;
	MGMatrix mat2; mat2.set_rotate_3D(uvecx, cval, sval);
	// mat2 is the matrix to transform ay to be vy(normalized uvecy).

	return *this=mat1*mat2;
}

//  原点を通り、指定ベクトルに垂直な平面に関して鏡面変換する 3D Matrix
//  を作成し，既存のMatrixと入れ換える。
MGMatrix& MGMatrix::set_reflect_3D(const MGVector& vec1)
{
	set_axis(vec1,0);
    // this = To transform vec1 to x-axis.

	(*this)(0,0) *= -1.;(*this)(1,0) *= -1.;(*this)(2,0) *= -1.;
	//Reflection about x-axis.

	MGMatrix m2; m2.set_vector(vec1,0);
	// m2 = To transform x-axis to vec1.

	return (*this)*=m2;
}

//  与えられたベクトル回りに指定の角度回転させる3D Matrixを作成し既存の
//  Matrixと入れ換える。
MGMatrix& MGMatrix::set_rotate_3D(
					const MGVector& vec, //Rotate Axis vector
					double angle	)	 //Angle
{
	set_axis(vec,0);
    // this = To transform vec1 to x-axis.

	MGMatrix temp; temp.set_rotate_2D(angle);
	MGMatrix m2=MGMatrix(3,temp,1,0);
	//m2= Rotation about x-axis.

	MGMatrix m3; m3.set_vector(vec,0);
	// m3 = To transform x-axis to vec1.

	(*this)*=m2;
	return (*this)*=m3;
}

//3D rotation matrix around vec.
//The angle is given by cval as cos(angle) and sval as sin(angle).
MGMatrix& MGMatrix::set_rotate_3D
	(const MGVector& vec,		//Rotate Axis vector
	double cval, double sval)	//Angle in cos() and sin() 
{
	set_axis(vec,0);
    // this = To transform vec1 to x-axis.

	MGMatrix temp; temp.set_rotate_2D(cval, sval);
	MGMatrix m2=MGMatrix(3,temp,1,0);
	//m2= Rotation about x-axis.

	MGMatrix m3; m3.set_vector(vec,0);
	// m3 = To transform x-axis to vec1.

	(*this)*=m2;
	(*this)*=m3;
	return *this;
}

//  各軸方向で等しい Scaling のためのMatrixを生成し、既存の
//  Matrixと入れ換える。Not change space dimension.
MGMatrix& MGMatrix::set_scale(double scale){
	size_t dim=sdim();
	for(size_t i=0; i<dim; i++){
		for(size_t j=0; j<dim; j++){
			if(j==i) (*this)(i,j)=scale;
			else     (*this)(i,j)=0.;
		}
	}
	return *this;
}

//  各軸方向で異なる Scaling のためのMatrixを生成し、既存の
//  Matrixと入れ換える。
MGMatrix& MGMatrix::set_diff_scale(double* scale){
	size_t dim=sdim();
	for(size_t i=0; i<dim; i++){
		for(size_t j=0; j<dim; j++){
			if(j==i) (*this)(i,j)=scale[i];
			else     (*this)(i,j)=0.;
		}
	}
	return *this;
}

//Set up this matrix from OpenGL matrix.
MGMatrix& MGMatrix::set_glMatrix(const double glMat[16]){
	resize(3);
	for(int i=0; i<3; i++){
		int i4=i*4;
		for(int j=0; j<3; j++){
			(*this)(i,j)=glMat[i4+j];
		}
	}
	return *this;
}

//Set this as a null matrix.
void MGMatrix::set_null(){
	if(m_matrix) delete[] m_matrix;
	m_matrix=0;
	m_sdim=0;
}

//Construct a matrix to transform a vector 'uvec' to be one of the axises,
// and replace own matrix. Inverse matrix of from_axis.
// axis can be any number(can be more than 2).
MGMatrix& MGMatrix::to_axis(
	const MGUnit_vector& uvec,	// Unit vector to be an axis.
	size_t axis)				// Axis kind 0:x, 1:y, 2:z, ...
{
	size_t dim=uvec.sdim(); if(axis>=dim) dim=axis+1;
	size_t i,j;
	MGUnit_vector unit=uvec;
	double len, len1, len2, max_len, cosv, sinv;
	MGMatrix mat(dim,1.);
	size_t axis1=axis, axis2;
	//Compute matrix to rotate for axis1 element to be zero
	//from axis+1 to axis-1.
	for(i=0; i<dim-1; i++){
		axis1++; if(axis1>=dim) axis1=0;
		//Get maximum length element id in axis2.
		axis2=0; max_len=fabs(unit.ref(0));
		for(j=1; j<dim; j++){
			len=fabs(unit.ref(j));
			if(len>max_len){ axis2=j; max_len=len;}
		}
		if(axis1==axis2) axis2=axis;

		len1=unit.ref(axis1); len2=unit.ref(axis2);
		len=sqrt(len1*len1+len2*len2);
		cosv=len1/len; sinv=len2/len;
		mat*=MGMatrix(dim,axis1,axis2,sinv,cosv);
		//Above MGMatrix(...) is to rotate for axis1 element to be zero.
		unit=uvec*mat;
//		cout<<" axis from i="<<i<<":"<<endl<<mat<<unit<<endl;/////
	}
	return *this=mat;
}

//  転置Matrixを生成。
MGMatrix MGMatrix::transpose() const{
	size_t dim=sdim();
    MGMatrix mat1(dim);
	for(size_t i=0; i<dim; i++){
		for(size_t j=0; j<dim; j++) mat1(j,i)=(*this)(i,j);
	}
    return mat1;
}

//
// 演算子の多重定義
//

// Assignment
MGMatrix& MGMatrix::operator=(const MGMatrix& mat){
	size_t dim=mat.m_sdim;
	if(m_sdim<dim) resize(dim); else m_sdim=dim;
	size_t len=m_sdim*m_sdim;
	for(size_t i=0; i<len; i++) m_matrix[i]=mat.m_matrix[i];
	return *this;
}

//自身のMatrixと与えられたscaleの乗算を行いオブジェクトを生成。
//Scaling of the matrix.
MGMatrix MGMatrix::operator* (double a) const{
	MGMatrix mat(*this);
	mat *= a;
	return mat;
}

 //自身のMatrixと与えられたscaleの乗算を行いオブジェクトを生成。
 //Scaling of the matrix.
MGMatrix operator* (double scale, const MGMatrix& mat){
	return mat*scale;
}

//自身のMatrixと与えられたscaleの乗算し自身のMatrixとする。
//Scaling of the matrix.
MGMatrix& MGMatrix::operator*= (double a){
	size_t len=sdim()*sdim();
	for(size_t i=0; i<len; i++) m_matrix[i] *=a;
	return *this;
}
	
// マトリックスによるベクトルの変換を行いオブジェクトを生成
//Matrix transformation of the vector.
MGVector operator* (const MGVector& v, const MGMatrix& m){
	MGVector vec2(v);
	vec2 *= m;
	return vec2;
}

//  マトリックスによるベクトルの変換を行い自身のベクトルとする。
MGVector& operator*= (MGVector& v,const MGMatrix& mat1){
	size_t dim=v.sdim();
	size_t i,j;	double a;
	const MGVector temp(v);	//Save the vector v.
	for(i=0; i<dim; i++){
		a=0.;
		for(j=0; j<dim; j++) a=a+temp[j]*mat1.ref(j,i);
		v(i)=a;
	}
	return v;
}

// マトリックスによるベクトルの変換を行い自身のベクトルとする
//Update own vector by matrix transformation.
//The result is unit of transformed vector.
MGUnit_vector& operator*= (MGUnit_vector& v,const MGMatrix& mat){
	return v=MGVector(v)*mat;
}

//  自身のMatrixと与えられたMatrixの乗算を行いオブジェクトを
//  生成。
MGMatrix MGMatrix::operator* (const MGMatrix& mat2) const{
	return multiply(mat2);
}

//  自身のMatrixと与えられたMatrixを乗算し自身のMatrix
//  とする。
MGMatrix& MGMatrix::operator*= (const MGMatrix&  mat2)
{
	return *this=multiply(mat2);
}

//自身のMatrixと与えられたTransfの乗算を行いオブジェクトを生成。
//Matrix and Transf multiplication.
MGTransf MGMatrix::operator* (const MGTransf& tr) const
{
	return MGTransf((*this)*tr.affine(), tr.translation());
}

//  Boolean 演算

//  自身のMatrixと与えられたMatrixが等しいかどうか
//  比較を行う。
bool MGMatrix::operator== (const MGMatrix& mat2) const{
	bool equal=true;
	size_t dim=sdim(); size_t dim2=mat2.sdim();
	dim= dim<=dim2 ? dim2:dim;
	if(dim<=0) return equal;

	MGVector v1(dim), v2(dim); size_t i,j;
	for(i=0; i<dim; i++){
		for(j=0; j<dim; j++){
			v1(j)=ref(i,j);
			v2(j)=mat2.ref(i,j);
		}
		if(v1 != v2) {equal=false; break;}
	}

	return equal;
}

//  自身のMatrixと与えられたMatrixが等しいかどうか
//  比較を行う。
bool MGMatrix::operator!= (const MGMatrix& mat2) const{
	return !((*this)==mat2);
}
