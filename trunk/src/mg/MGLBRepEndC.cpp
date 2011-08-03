/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/NDDArray.h"
#include "mg/BPointSeq.h"
#include "mg/LBRepEndC.h"
#include "mg/Curve.h"

extern "C"{
#include "cskernel/Bvltan.h"
#include "cskernel/bvutan.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLBRepEndC.cc
//
// Implements MGLBRepEndC class.

//Member Data
//	MGENDCOND m_cond;	//Type of end condition.
//	MGVector  m_1deriv;	//1st derivative stored 
//						//when m_cond=MGENDC_1D or MGENDC_12D
//	MGVector  m_2deriv;	//2nd derivative stored

//Constructor

MGLBRepEndC::MGLBRepEndC(	//Construct of m_cond=MGENDC_1D or MGENDC_2D.
	MGENDCOND cond,			//Type of end condition
							//(MGENDC_1D or MGENDC_2D)
	const MGVector& deriv)	//Derivative inf according to cond
	: m_cond(cond) {
	switch (cond){
	case MGENDC_1D:
		m_1deriv=deriv; break;
	case MGENDC_2D:
		m_2deriv=deriv; break;
	default:
		m_cond=MGENDC_NO; break;
	}
}

MGLBRepEndC::MGLBRepEndC(	//Construct of m_cond=MGENDC_12D
	const MGVector& first_deriv,	//1st derivative
	const MGVector& second_deriv)	//2nd derivative
	: m_cond(MGENDC_12D)
	, m_1deriv(first_deriv), m_2deriv(second_deriv) {;}

MGLBRepEndC::MGLBRepEndC(		//BVLTAN
	int start,					//Indicates start(start==true) condition
								// or end.
	const MGNDDArray& tau,		//Data point abscissa
	const MGBPointSeq& points,	//Point seq data
	int &error)					//Error flag.
	:m_cond(MGENDC_1D)
{
	assert(tau.length()==points.length() && tau.length()>1);

	const size_t np=points.length();
	const size_t ip=points.capacity();
	const size_t ncd=points.sdim();
	int ise=2; if(start) ise=1;
	double* work12=new double[8*ncd+81]; 
	double* tangen=new double[ncd];
	bvltan_(np,points.data(),tau.data(),ip,ncd,ise,
			work12,work12+70,work12+81,tangen,&error);
	if(error==1){
		error=0;
		m_1deriv=MGVector(ncd,tangen);
	}
	delete[] work12; delete[] tangen;
}

MGLBRepEndC::MGLBRepEndC(		//BVUTAN
	int start,			//Indicates start(start==true)
						// condition or end.
	const MGBPointSeq& points)	//Point seq data
	:m_cond(MGENDC_1D)
{
	assert(points.length()>1);

	const size_t np=points.length();
	const size_t ip=points.capacity();
	const size_t ncd=points.sdim();
	int ise=2; if(start) ise=1;
	double* work=new double[8*ncd+81];
	double* tangen=new double[ncd];
	int error;
	bvutan_(ise,np,points.data(),ip,ncd,work,tangen,&error);
	m_1deriv=MGVector(ncd,tangen);
	delete[] work; delete[] tangen;
}

// Given MGCurve, construct the curve's end condition.
MGLBRepEndC::MGLBRepEndC(
	int start,		//Indicates start(start==true) condition or end.
	MGENDCOND cond,	//Type of end condition(MGENDC_1D, MGENDC_2D, or MGENDC_12D)
	const MGCurve& curve//Curve
):m_cond(cond){
	double tau;
	if(start) tau=curve.param_s();
	else tau=curve.param_e();

	if(cond==MGENDC_1D || cond==MGENDC_12D){
		m_1deriv=curve.eval(tau,1);
	}
	if(cond==MGENDC_2D || cond==MGENDC_12D){
		m_2deriv=curve.eval(tau,2);
	}
}

//Destructor
//	~MGLBRepEndC();	  We use default destructor.

//Member Function

//Initialize the instance. Will be set to the same as constructed by the void 
//constructor.
void MGLBRepEndC::initialize(){
	m_cond=MGENDC_NO;
	m_1deriv.set_null();
	m_2deriv.set_null();
}

void MGLBRepEndC::set_1st(const MGVector& first_deriv)
//Set 1st deriv and change condition type to MGENDC_1D or MGENDC_12D.
{
	switch (m_cond){
	case MGENDC_NO:
		m_cond=MGENDC_1D; m_1deriv=first_deriv; break;
	case MGENDC_2D:
		m_cond=MGENDC_12D; m_1deriv=first_deriv; break;
	default:
		m_1deriv=first_deriv; break;
	}
}

void MGLBRepEndC::set_2nd(const MGVector& second_deriv)
//Set 2nd deriv and change condition type to MGENDC_2D or MGENDC_12D.
{
	switch (m_cond){
	case MGENDC_NO:
		m_cond=MGENDC_2D; m_2deriv=second_deriv; break;
	case MGENDC_1D:
		m_cond=MGENDC_12D; m_2deriv=second_deriv; break;
	default:
		m_2deriv=second_deriv; break;
	}
}

//Operator overload.
