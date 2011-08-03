/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Declaration for class MGIgesPD104.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgiges/IgesGsec.h"
#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesPD104.h"
#include "mg/Tolerance.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

using namespace MGIges;
//!	@brief MGIgesPD104 is the class for Iges parameter data type 104(conic arc).

// Constructors.

//! Constructs an object of class MGIgesPD104.
MGIgesPD104::MGIgesPD104(MGIgesDirectoryEntry* DEpointer)
:MGIgesPD(CONIC_ARC,DEpointer){
}

//Construct PD100, supplying 2D coordinate data in each array.
MGIgesPD104::MGIgesPD104(
	const double coef[6],//Coefficients of the conic equation.
	double Zt,//Z coordinate of (x,y) plane of the above equation.
	const double start[2], const double terminate[2]
		//the start and terminate point coordinates of the conic arc in (x,y) plane.
):MGIgesPD(CONIC_ARC), m_Zt(Zt), m_X1(start[0]), m_Y1(start[1]),
m_X2(terminate[0]), m_Y2(terminate[1]){
	for(size_t i=0; i<6; i++)
		m_coef[i]=coef[i];
}

//Read in parameter data from string stream data.
void MGIgesPD104::read_in(
	char pDelimeter,
	std::istringstream& pdstream
){
	for(size_t i=0; i<6; i++){//Read A,B,C,D,E,F
		get_real(pDelimeter,pdstream,m_coef[i]);
	}
	get_real(pDelimeter,pdstream,m_Zt);
	get_real(pDelimeter,pdstream,m_X1);
	get_real(pDelimeter,pdstream,m_Y1);
	get_real(pDelimeter,pdstream,m_X2);
	get_real(pDelimeter,pdstream,m_Y2);
}

//Write out this PD as MGIgesParamLine's(into plines).
//Except for string data, one integer or double data is output
//into one MGIgesParamLine, not striding over more than one line.
//Only when string data is output(to Holleris string), the data
//may stride over more than one lines.
//plines[i] for 0<=i<plines.size() are valid.
void MGIgesPD104::write_out_into_string(
	const MGIgesGSec& gsec,	//Input gsec to input delimeter_param and delimeter_record;
	MGPvector<std::string>& plines //output plines.
)const{
	for(size_t i=0; i<6; i++){//Write A,B,C,D,E,F
		put_real(m_coef[i],gsec,plines);
	}
	put_real(m_Zt,gsec,plines);
	put_real(m_X1,gsec,plines);
	put_real(m_Y1,gsec,plines);
	put_real(m_X2,gsec,plines);
	put_real(m_Y2,gsec,plines);
}

//Convert de to MGObject(a newed object). de must be of type 104(conic arc).
//Output MGObject is a MGEllipse, or MGRLBRep(when the arc is parabola or hyperbola).
MGObject* MGIgesIfstream::convert_conic_arc(
	const MGIgesDirectoryEntry& de
)const{
	double rzero=MGTolerance::rc_zero();
	double zero=rzero*rzero;
		//This is because A,B,C are square of zero.

	int fnum=de.FormNumber();
	const std::auto_ptr<MGIgesPD>& pd=de.paramData();
	const MGIgesPD104* pd104=static_cast<const MGIgesPD104*>(pd.get());
	const double* coef=pd104->m_coef;
	double maxCoef=coef[0]; if(maxCoef<0.) maxCoef*=-1.;
	for(size_t i=1; i<=2; i++){
		double absCoef=coef[i];
		if(absCoef<0.)
			absCoef*=-1.;
		if(maxCoef<absCoef)
			maxCoef=absCoef;
	}
	if(maxCoef<zero)
		return 0;

	MGObject* obj=0;
	double A=coef[0]/maxCoef, B=coef[1]/maxCoef, C=coef[2]/maxCoef;
	double D=coef[3]/maxCoef ,E=coef[4]/maxCoef, F=coef[5]/maxCoef;
	double Zt=pd104->m_Zt;
	double X1=pd104->m_X1, X2=pd104->m_X2, Y1=pd104->m_Y1, Y2=pd104->m_Y2;

	double Q2=A*C-B*B/4.;
	double Xm=(X1+X2)*.5, Ym=(Y1+Y2)*.5;
	MGPosition P0(X1, Y1, Zt);//Start point.
	MGPosition P1(0.,0.,Zt);//Mid point.
	MGPosition P2(X2, Y2, Zt);//Terminate point.
	MGVector T0(1.,1., 0.), T2(1.,1.,0.);//Tangent vector at start and terminate point.
	double zero2=zero*zero;
	if(Q2<zero2){
		if(Q2>(-zero2)){
		//Parabola
			if(!MGRZero(E)){
			//Parabola of Y axis
				double AoverE=A/E;
				P1(0)=Xm; P1(1)=-AoverE*Xm*Xm;
				if(X1<X2){
					T0(1)=-2.*AoverE*X1;
					T2(1)=-2.*AoverE*X2;
				}else{
					T0(0)=-1.; T0(1)=2.*AoverE*X1;
					T2(0)=-1.; T2(1)=2.*AoverE*X2;
				}
			}else if(!MGRZero(D)){
			//Parabola of X axis
				double CoverD=C/D;
				P1(0)=-CoverD*Ym*Ym; P1(1)=Ym;
				if(Y1<Y2){
					T0(1)=-2.*CoverD*Y1;
					T2(1)=-2.*CoverD*Y2;
				}else{
					T0(0)=2.*CoverD*Y1; T0(1)=-1.;
					T2(0)=2.*CoverD*Y2; T2(1)=-1.;
				}
			}
		}else{
		//Hyperbola
			if(!MGRZero(A) && !MGRZero(C)){
				double FoverA=F/A, FoverC=F/C;
				if(FoverA<0. && FoverC>0.){
					double a=sqrt(-FoverA), b=sqrt(FoverC);
					double t1=atan2(Y1,b), t2=atan2(Y2,b);
					double ct1=cos(t1), ct2=cos(t2);
					double ct1ct1=ct1*ct1, ct2ct2=ct2*ct2;
					double tant1=tan(t1), tant2=tan(t2);
					double tm=(t1+t2)*.5;
					P1(0)=a/cos(tm); P1(1)=b*tan(tm);
					if(t1<t2){
						T0(0)=a*tant1/ct1; T0(1)=b/ct1ct1;
						T2(0)=a*tant2/ct2; T2(1)=b/ct2ct2;
					}else{
						T0(0)=-a*tant1/ct1; T0(1)=-b/ct1ct1;
						T2(0)=-a*tant2/ct2; T2(1)=-b/ct2ct2;
					}
				}else{
					double a=sqrt(FoverA), b=sqrt(-FoverC);
					double t1=atan2(X1,a), t2=atan2(X2,a);
					double ct1=cos(t1), ct2=cos(t2);
					double ct1ct1=ct1*ct1, ct2ct2=ct2*ct2;
					double tant1=tan(t1), tant2=tan(t2);
					double tm=(t1+t2)*.5;
					P1(0)=a*tan(tm); P1(1)=b/cos(tm);
					if(t1<t2){
						T0(0)=a/ct1ct1; T0(1)=b*tant1/ct1;
						T2(0)=a/ct2ct2; T2(1)=b*tant2/ct2;
					}else{
						T0(0)=-a/ct1ct1; T0(1)=-b*tant1/ct1;
						T2(0)=-a/ct2ct2; T2(1)=-b*tant2/ct2;
					}
				}
			}
		}
		obj=new MGRLBRep(P0,T0,P1,P2,T2);//std::cout<<(*obj)<<std::endl;//////////********
	}else{
	//Ellipse
		double FoverA=F/A, FoverC=F/C;
		double a=sqrt(-FoverA), b=sqrt(-FoverC);
		double t1=MGAngle(X1/a, Y1/b), t2=MGAngle(X2/a, Y2/b);
		double azero=mgDBLPAI*rzero;
		if(-azero<t1 && t1<azero)
			t1=0.;
		if(-azero<t2 && t2<azero)
			t2=mgDBLPAI;
		double t2mt1=t2-t1;
		if(-azero<t2mt1 && t2mt1<azero)
			t2=t1+mgDBLPAI;
		else if(t2<t1)
			t2+=mgDBLPAI;
		MGInterval rng(t1,t2);
		obj=new MGEllipse(MGPosition(0.,0.,Zt),mgX_UVEC*a, mgY_UVEC*b,rng);
	}
	return obj;
}
