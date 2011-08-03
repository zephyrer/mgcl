/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/nlbit.h"
#include "mg/CurveParameter.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"
#include "mg/CCisect_list.h"
#include "mg/CCisect.h"
#include "mg/CParam_list.h"
#include "mg/Tolerance.h"

//#if defined(_DEBUG)
//#define new DEBUG_NEW
//#undef THIS_FILE
//static char THIS_FILE[] = __FILE__;
//#endif

using namespace std;

// Implementation of MGEllipse.

//Utility class for MGEllipse constructor of curve1, 2, and tangent point on crv1.
class MGEllipseTTR: public MGCurveParameter{

	const MGCurve& m_curve1;		//curve1.
	const MGCurve& m_curve2;		//curve2.
	double m_r;	//radius of the target arc.
	double m_rabs;//absolute value of m_r.
	double m_t2guess;// initial guess parameter of m_curve2.
	const MGUnit_vector& m_normal;//The normal of the plane where the arc lies on.
	mutable double m_t2obtained;
		//When solution found, m_t2obtained will be the tangent point parameter of curve2.
	mutable MGPosition m_center;
		//When solution found, m_center will be the center of the arc.
	double m_r2;	//square of m_r.

public:

MGEllipseTTR(
	const MGCurve& crv1,//curve1.
	double t2guess,		//guess parameter value of crv2.
	const MGCurve& crv2,//curve2.
	double r,	//radius of the target arc.
	const MGUnit_vector& normal//The normal of the plane where the arc lies on.
):MGCurveParameter(crv1.param_range(),MGTolerance::wc_zero_sqr())
,m_curve1(crv1),m_t2guess(t2guess),m_curve2(crv2),m_r(r),m_rabs(r),m_normal(normal){
	if(m_rabs<0.) m_rabs*=-1.;
	m_r2=r*r;
};

//return the difference of the two radii. That is,
//let m_center be the point normally away by r from curve's point P1(of parameter t),
//and P2 is the foot point of perpendicular line from m_center to curve2.
//Then this operator returns a*a-b*b, where a=the distance from m_center to
//P1, and b=distance from P2 to m_center.
double operator()(
	double t	//is a parameter of m_curv2;
)const{
	MGVector tan1=m_curve1.eval(t,1);
	MGPosition P1=m_curve1.eval(t);
	MGUnit_vector N1=m_normal*tan1;
	if((N1%(m_curve2.eval(m_t2guess)-P1))<0.){
		N1.negate();
	}
	m_center=P1+N1*m_rabs;
	if(!m_curve2.perp_guess(1.,0.,m_center,m_t2guess,m_t2obtained)){
		m_t2obtained=m_curve2.closest(m_center);
	}
	MGVector diff=m_curve2.eval(m_t2obtained)-m_center;
	return diff%diff-m_r2;
};

MGPosition& center(){return m_center;};
double curve2_param(){return m_t2obtained;};

};

//�ړI�F
//		��{���Q�{�Ɣ��a����R�[�i�[�q���쐬����
//		�����p�����[�^�̓R�[�i�[�q���쐬���������ɐݒ肵�A��{���͓��ꕽ�ʏ�ɂ��邱��
//�߂�l�F
//		0:			����I��
//		-1:			���͒l���s��
//		-2:			��{���̋ȗ����a��蔼�aR���傫��
//		-3:			��{�����m�̌�_�����܂�Ȃ�
//Start point of the generated ellipse is either t1 side or t2,
//which depends on normal direction. The ellipse direction is defined from normal,
//If the arc length from t1 to t2(arond normal) is smaller than from t2 to t1(around normal),
//the start point is on t1 side. If the arc length from t1 to t2 is longer than from t2 to t1,
//the start point is on t2 side.
MGEllipse::MGEllipse(
	const MGCurve&			crv1,	//I:��{��1
	const MGCurve&			crv2,	//I:��{��2
	const MGUnit_vector&	normal,	//I:��{���̃m�[�}���x�N�g��
	double					r,//I:�R�[�i�[�q�̔��a
	double&					t1,		//I:��{���P�̏����p�����[�^(�R�[�i�[�q���쐬���鑤���w��)
									//O:�R�[�i�[�q�Ɛڂ����{���P�̃p�����[�^�l
	double&					t2,		//I:��{���Q�̏����p�����[�^(�R�[�i�[�q���쐬���鑤���w��)
									//O:�R�[�i�[�q�Ɛڂ����{���Q�̃p�����[�^�l
	int&					rc)		//O:���^�[���R�[�h
:MGCurve(crv1),m_gprange(0),m_knotV(0){
	MGEllipseTTR ttr(crv1,t2,crv2,r,normal);
	ttr.set_delta(crv1.param_span()/crv1.divide_number());
	rc=ttr.getCurveParameter(t1);
	if(rc)
        return;
	
	//Construct arc data.
	MGPosition P1=crv1.eval(t1);
	t2=ttr.curve2_param();
	MGPosition P2=crv2.eval(t2);
	MGPosition& Center=ttr.center();
	MGVector vec1(P1 - Center), vec2(P2 - Center);
	double dAng = vec1.angle(vec2);
	MGVector V12(vec2 - vec1), B1(normal*vec1);
	if((V12%B1)>0.){
		copy_ellipse_data(MGEllipse(Center, P1, dAng, normal));//start is P1.
	}else{
		copy_ellipse_data(MGEllipse(Center, P2, dAng, normal));//start is P2.
	}
	rc=0;//Normal return;
}

//Utility class for MGEllipse constructor of curve1, 2, and tangent point on crv1.
class MGEllipseTTP: public MGCurveParameter{
	const MGUnit_vector& m_normal;//The normal of the plane where the arc lies on.
	const MGUnit_vector& m_N1;	//curve1's normal at the point P1.
	const MGPosition& m_P1;		//curve1's point which is the end of the arc.
	const MGCurve& m_crv2;		//curve2.

public:
MGEllipseTTP(
	const MGUnit_vector& normal,//The normal of the plane where the arc lies on.
	const MGUnit_vector& N1,	//curve1's normal at the point P1.
	const MGPosition& P1,		//curve1's point which is the end of the arc.
	const MGCurve& crv2			//curve2.
):MGCurveParameter(crv2.param_range(),MGTolerance::wc_zero_sqr())
,m_normal(normal), m_N1(N1), m_P1(P1), m_crv2(crv2){;};

//return the difference of the two radii. That is,
//let Q be the intersection of the straight from m_P1 to m_N1 and from P2 to normal at P2
//of the curve2. Then this operator returns a*a-b*b, where a=signed distance from
//P1 to Q, and b=signed distance from P2 to Q.
double operator()(
	double t	//is a parameter of m_curv2;
)const{
	MGUnit_vector N2=normal2(t);
	MGVector p12=getP12(t);
	double p12n1=p12%m_N1, p12n2=p12%N2, n12=m_N1%N2;
	return (p12n1*p12n1-p12n2*p12n2)/(1.-n12*n12);
};

MGUnit_vector normal2(double t2)const{
	return m_normal*m_crv2.eval(t2,1);
};
MGVector getP12(double t2)const{
	return m_crv2.eval(t2)-m_P1;
};
double radius(double t2)const{
	MGUnit_vector N2=normal2(t2);
	double n12=m_N1%N2;
	MGVector p12=getP12(t2);
	return (p12%m_N1-(p12%N2)*n12)/(1.-n12*n12);
};

};

//Construct an arc from 2 curves and a tangent point on the curve 1.
//The 2 curves must be planar and must lie on the sampe plane.
//return code rc:
//	=0: normally constructed.
//	=1:no arcs found that are tangent to crv1 and 2.
//	=-2: error detected in mgNlbit
//       (This case does not occur normally, maybe some bugs are included).
//Start point of the generated ellipse is either t1 side or t2,
//which depends on normal direction. The ellipse direction is defined from normal,
//If the arc length from t1 to t2(arond normal) is smaller than from t2 to t1(around normal),
//the start point is on t1 side. If the arc length from t1 to t2 is longer than from t2 to t1,
//the start point is on t2 side.
MGEllipse::MGEllipse(
	const MGCurve&			crv1,	//curve1
	const MGCurve&			crv2,	//curve2
	const MGUnit_vector&	normal,	//normal of the plane 2 curvs lie on
	double					t1,		//tangent point's parameter of crv1.
	double&					t2,		//guess parameter value of crv2 is input initially,
							//correct tangent point's parameter will be output when rc=0.
	int&					rc		//return code.
):MGCurve(crv1),m_gprange(0), m_knotV(0){
	MGVector tan1=crv1.eval(t1,1);
	MGPosition P1=crv1.eval(t1);
	MGUnit_vector N1=normal*tan1;
	MGEllipseTTP radius_diff(normal,N1,P1,crv2);
	double delta=crv2.param_span()/crv2.divide_number();
	radius_diff.set_delta(delta);
	rc=radius_diff.getCurveParameter(t2);
	if(rc)
        return;
	
	//Construct arc data.
	MGPosition P2=crv2.eval(t2);
	MGPosition Center=P1+N1*radius_diff.radius(t2);
	MGVector vec1(P1 - Center), vec2(P2 - Center);
	double dAng = vec1.angle(vec2);
	MGVector V12(vec2 - vec1), B1(normal*vec1);
	if((V12%B1)>0.){
		copy_ellipse_data(MGEllipse(Center, P1, dAng, normal));
	}else{
		copy_ellipse_data(MGEllipse(Center, P2, dAng, normal));
	}
	rc=0;//Normal return;
}

//�ړI�F
//		��{���R�{����R�[�i�[�q���쐬����
//		crv1, crv2, crv3�����ꂼ��n�_�A�ʉߓ_�A�I�_�ƂȂ�悤�ɂ���B
//		�����p�����[�^�̓R�[�i�[�q���ڂ���ߕӂɐݒ肵�A��{���͓��ꕽ�ʏ�ɂ��邱��
//�߂�l�F
//		0:			����I��
//		-1:			�A���������������Ȃ�����
//		-3:			�������Ȃ�����(�������[�v�ɂ���������)
MGEllipse::MGEllipse(
	const MGCurve& crv1,	//I:��{��1(�n�_)
	const MGCurve& crv2,	//I:��{��2(�ʉߓ_)
	const MGCurve& crv3,	//I:��{��3(�I�_)
	const MGUnit_vector& normal,//I:��{���̃m�[�}���x�N�g��
	double& t1,		//I:��{��1�̏����p�����[�^(�R�[�i�[�q�̎n�_�ߕ�)
		 			//O:�R�[�i�[�q�̎n�_�ɂ������{���P�̃p�����[�^�l
	double& t2,		//I:��{��2�̏����p�����[�^(�R�[�i�[�q�̒ʉߓ_�ߕ�)
					//O:�R�[�i�[�q�̒ʉߓ_�ɂ������{���Q�̃p�����[�^�l
	double& t3,		//I:��{��3�̏����p�����[�^(�R�[�i�[�q�̏I�_�ߕ�)
					//O:�R�[�i�[�q�̏I�_�ɂ������{���Q�̃p�����[�^�l
	int& 	rc		//O:���^�[���R�[�h
):MGCurve(crv1),m_gprange(0), m_knotV(0){
	int siNit=20;
	double d1, d2, d3, r;
	double error=MGTolerance::wc_zero_sqr();
	MGPosition P1,P2,P3;
	rc=-3;
	while(siNit-->0){
		P1=crv1.eval(t1), P2=crv2.eval(t2), P3=crv3.eval(t3);
		MGPosition posM((P1+P2+P3)/3.);	//�d�S
		MGVector tan1=crv1.eval(t1,1), tan2=crv2.eval(t2,1), tan3=crv3.eval(t3,1);
		MGUnit_vector T1(tan1), T2(tan2), T3(tan3);
		MGUnit_vector N1(normal*T1), N2(normal*T2), N3(normal*T3);
		//�����x�N�g���́A�R�_�̏d�S�����ɂȂ�悤�ɂ���
		if(((posM-P1)%N1) < 0.0)
			N1.negate();
		if(((posM-P2)%N2) < 0.0)
			N2.negate();
		if(((posM-P3)%N3) < 0.0)
			N3.negate();

		//�ȉ��̎���d1,d2,d3��0�ƂȂ�r�����߂�
		//P2 + d2*T2 + r*N2 = P1 + d1*T1 + r*N1 = P3 + d3*T3 + r*N3
		double	a00 = T1%N2, a01 = (N1-N2)%N2,
				a10 = T1%N3, a11 = (N1-N3)%N3,
				b0 = (P2-P1)%N2, b1 = (P3-P1)%N3;
		double det = a00*a11-a10*a01;
		if(!MGMZero(det)){
		    d1 =  (b0*a11-b1*a01)/det;
			r =  (a00*b1-a10*b0)/det; 
		}else{
			rc = -1;
			break;
		}
		MGPosition center(P1+r*N1);
		d2 = (center-P2)%T2;
		d3 = (center-P3)%T3;
		if((d1*d1+d2*d2+d3*d3)<error){
			//���������̂ŃR�[�i�[�q�𐶐�����
			rc=0;
			break;
		}
		t1=crv1.param_round_into_range(t1+d1/tan1.len());
		t2=crv2.param_round_into_range(t2+d2/tan2.len());
		t3=crv3.param_round_into_range(t3+d3/tan3.len());
	}
	copy_ellipse_data(MGEllipse(P1, P2, P3));
}

//Utility class for MGEllipse constructor of curve, radius of the arc, and
//end point of the arc(TPR).
class MGEllipseTPR: public MGCurveParameter{

	const MGCurve& m_curve;	//curve.
	double m_r;	//radius of the target arc.
	const MGPosition& m_P2;	//end point of the arc.
	const MGUnit_vector& m_normal;//The normal of the plane where the arc lies on.
	mutable MGPosition m_center;
		//When solution found, m_center will be the center of the arc.
	double m_r2;	//square of m_r.

public:

MGEllipseTPR(
	const MGCurve& crv,//curve.
	double r,	//radius of the target arc.
	const MGPosition& P2,
	const MGUnit_vector& normal//The normal of the plane where the arc lies on.
):MGCurveParameter(crv.param_range(),MGTolerance::wc_zero_sqr())
,m_curve(crv),m_P2(P2),m_r(fabs(r)),m_normal(normal),m_r2(r*r){
};

//return the difference of the two radii. That is,
//let m_center be the point normally away by r from curve's point P1(of parameter t),
//Then this operator returns a*a-b*b, where a=the distance from m_center to
//P1, and b=distance from m_P2 to m_center.
double operator()(
	double t	//is a parameter of m_curv2;
)const{
	MGVector tan1=m_curve.eval(t,1);
	MGPosition P1=m_curve.eval(t);
	MGUnit_vector N1=m_normal*tan1;
	MGVector p12=m_P2-P1;
	if(N1%p12<0.){
		N1.negate();
	}
	m_center=P1+N1*m_r;
	MGVector cp2=m_P2-m_center;
	return cp2%cp2-m_r2;
};

MGPosition& center(){return m_center;};

};

//�ړI�F
//		��{���Ɖ~�ʒ[�_�Ɣ��a����R�[�i�[�q���쐬����
//		�����p�����[�^�̓R�[�i�[�q���쐬���������ɐݒ肵�A��{���Ɖ~�ʒ[�_�͓��ꕽ�ʏ�ɂ��邱��
//�߂�l�F
//		0:			����I��
//		-1:			���͒l���s��
//		-2:			���a�l������������
//		-3:			���a�l���傫������
//		-4:			��{�����I�t�Z�b�g�����J�[�u�Ɖ~�ʒ[�_����쐬�����~�̌�_�����܂�Ȃ�����
MGEllipse::MGEllipse (
	const MGCurve&			crv,				//I:��{��
	const MGPosition&		P2,				//I:�~�ʒ[�_
	const MGUnit_vector&	normal,				//I:��{���̃m�[�}���x�N�g��
	double					r,			//I:�R�[�i�[�q�̔��a
	double&					t,				//I:��{���̏����p�����[�^(�R�[�i�[�q���쐬���鑤���w��)
												//O:�R�[�i�[�q�Ɛڂ����{���̃p�����[�^�l
	int&					rc					//O:���^�[���R�[�h
):MGCurve(crv),m_gprange(0), m_knotV(0){
	rc = 0;

	MGEllipseTPR tpr(crv,r,P2,normal);
	tpr.set_delta(crv.param_span()/crv.divide_number());
	rc=tpr.getCurveParameter(t);
	if(rc)
        return;
	
	//Construct arc data.
	MGPosition P1=crv.eval(t);
	MGPosition& Center=tpr.center();
	MGVector vec1(P1 - Center), vec2(P2 - Center);
	double dAng = vec1.angle(vec2);
	MGVector V12(vec2 - vec1), B1(normal*vec1);
	if((V12%B1)>0.){
		copy_ellipse_data(MGEllipse(Center, P1, dAng, normal));
	}else{
		copy_ellipse_data(MGEllipse(Center, P2, dAng, normal));
	}
	rc=0;//Normal return;
}

//�ړI�F
//��{���Ɗ�{����q�~�܂�Ɖ~�ʒ[�_����R�[�i�[�q���쐬����
//�����p�����[�^�̓R�[�i�[�q���쐬���������ɐݒ肵�A��{���Ɖ~�ʒ[�_�͓��ꕽ�ʏ�ɂ��邱��
MGEllipse::MGEllipse (
	const MGCurve&			crv,	//I:��{��
	const MGPosition&		P2,		//I:�~�ʒ[�_
	const MGUnit_vector&	normal,	//I:��{���̃m�[�}���x�N�g��
	double t				//I:R�~�܂�_�̃p�����[�^
):MGCurve(crv),m_gprange(0), m_knotV(0){
	MGPosition P1=crv.eval(t);
	MGPosition Q=(P1+P2)*.5;
	MGVector V=Q-P1;

	MGVector T=crv.eval(t,1);
	MGUnit_vector N=normal*T;//Normal at t.
	if(N%V<0.)
		N.negate();//The arc is always on the side of the vector V.
	if(T%V<0.)
		T.negate();
	double vlen=V.len();
	double r=vlen*vlen/(N%V);//radius of the arc.
	MGPosition C=P1+N*r;//center.
	double angle=acos(vlen/r);//angle of the arc.
	copy_ellipse_data(MGEllipse(C, P1, mgPAI-2.*angle, T*N));
}
