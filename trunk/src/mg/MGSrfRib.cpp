/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Pvector.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/SPointSeq.h"
#include "mg/Surface.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

extern "C"{
#include "cskernel/Blgi2d.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//�����J�[�u�̎�ނœ����m�b�g�x�N�g���������Ă���X�v���C���Ȑ��񂩂ǂ������ׂ�
bool isSameKnotVector(const MGPvector<MGCurve>& pvCurves){
	std::vector<MGCurve*>::const_iterator citer = pvCurves.begin();
	MGCurve& curve0=**citer;
	long CrvId =curve0.type();
	//�X�v���C���J�[�u�ȊO���܂܂�Ă�����false��Ԃ�
	if(CrvId != MGCURVE_SPLINE && CrvId != MGCURVE_RSPLINE)
		return false;

	std::vector<MGCurve*>::const_iterator cend = pvCurves.end();
	const MGKnotVector& t0 = curve0.knot_vector();//1st knot vector.
	for(citer++; citer!=cend; citer++){
		MGCurve& curvei=**citer;
		//�J�[�u�̎�ނ�����Ă�����A�m�b�g�x�N�g��������Ă�����false��Ԃ�
		if(curvei.type()!=CrvId || curvei.knot_vector()!=t0)
			return false;
	}
	return true;
}

//���u�Ȑ��񂩂�ʂ��쐬����
//���u�Ȑ��̑S�Ẵm�b�g�������X�v���C���̎�MGSBRep��MGRSBRep���ԋp�����
//����ȊO�̏ꍇ�́A�Ȑ���LBRep�ōč\�����Ėʂ��쐬����̂�MGSBRep���ԋp�����
//�쐬����ʂ̃m�b�g�x�N�g���̓��u�Ȑ��̌�����u,���u�������v�Ƃ���
//Let v0=start parameter value, v1=terminate parameter value along v, then
//v=v0 const parameter line is curves[0], and v=v1 const parameter line is
//curves[n-1], where n=curves.size(). n must be greater or equal to 2.
//When n==2, the surface is a ruled surface(that is, order_u() is 2).
std::auto_ptr<MGSurface> MGCL::createSurfaceFromRibs(
	const MGPvector<MGCurve>& curves,//���u�Ȑ���
	bool direction_adjustment	//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	//���u��̃m�b�g�x�N�g�������낦���Ă��邩���ׂ�(B�X�v���C���ȊO��false���Ԃ�)
	if(isSameKnotVector(curves)){
		const MGCurve& crv0 = *(curves.front());
		if(crv0.type() == MGCURVE_RSPLINE){//Rational or non Rational�𒲂ׂ�
			size_t nRibNum = curves.size();
			std::vector<const MGRLBRep*> vecPtrRLBReps(nRibNum);
			for(size_t i=0; i<nRibNum; i++)
				vecPtrRLBReps[i]=static_cast<const MGRLBRep*>(curves[i]);
			return std::auto_ptr<MGSurface>(
				new MGRSBRep(vecPtrRLBReps,direction_adjustment));
		}
	}

	//�m�b�g�x�N�g�����قȂ������u�̏ꍇ�Ȑ������r���h����
	return std::auto_ptr<MGSurface>(new MGSBRep(curves,direction_adjustment));
}

MGDECL std::auto_ptr<MGSurface> MGCL::createSurfaceFromRibs(
	const std::vector<const MGCurve*>& curves,//���u�Ȑ���
	bool direction_adjustment	//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	size_t n=curves.size();
	MGPvector<MGCurve> curves2(n);
	for(size_t i=0; i<n; i++)
		curves2[i]=const_cast<MGCurve*>(curves[i]);
	std::auto_ptr<MGSurface> srf=
		MGCL::createSurfaceFromRibs(curves2,direction_adjustment);
	curves2.release_all();
	return srf;
}

//���u�Ȑ���̗��[�_�A���_�̋����𕽋ς���data points tau���쐬����
void createRibDataPoints(const MGPvector<MGLBRep>& curves, MGNDDArray& tau){
	int nBdim = curves.size();
	tau.resize(nBdim);
	double taui=0.0;
	int nBdimM1=nBdim-1; 
	double error=MGTolerance::wc_zero()*10.;
	for(int i=0; i<nBdimM1; i++){
		tau[i]=taui;
		const MGCurve& curvei=*(curves[i]); const MGCurve& curveiP1=*(curves[i+1]);
		double dStartDist = curvei.start_point().distance(curveiP1.start_point());
		double dEndDist = curvei.end_point().distance(curveiP1.end_point());
		double dCenterDist = curvei.center().distance(curveiP1.center());
		double dist=(dStartDist+dEndDist+dCenterDist)/3.0;
		if(dist<error)
			dist=error;
		taui += dist;
	}
	tau[nBdimM1]=taui;
}

// Creates a ruled surface.
std::auto_ptr<MGSurface> MGCL::create_ruled_surface(
	const MGCurve& cross1,  // a curve as Edge No.0.
	const MGCurve& cross2,   // another curve as Edge No.2.
	bool direction_adjustment	//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	std::vector<const MGCurve*> curves(2);
	curves[0]=&cross1; curves[1]=&cross2;
	return MGCL::createSurfaceFromRibs(curves,direction_adjustment);
}

///���u�Ȑ��񂩂�ʂ��쐬����
///�쐬����ʂ̃m�b�g�x�N�g���̓��u�Ȑ��̌�����u,���u�������v�Ƃ���
///This constructor only generates MGSBRep even if curves are MGRLBRep of the same
///knot configuration. To avoid this, use createSurfaceFromRibs() that generates
///MGRSBRep when curves are MGRLBRep of the same knot configuration.
///Let v0=start parameter value, v1=terminate parameter value along v, then
///v=v0 const parameter line is curves[0], and v=v1 const parameter line is
///curves[n-1], where n=curves.size(). n must be greater or equal to 2.
///When n==2, the surface is a ruled surface(that is, this->order_u() is 2).
///When direction_adjustment=true, curves[i] directions are adjusted to line up
///to the same direciton as the curves[i-1]'s.
///When direction_adjustment=false, no curves of curves[i] are negated and
///the directions of curves[i] are the direcitons of ribs.
MGSBRep::MGSBRep(const MGPvector<MGCurve>& curves,bool direction_adjustment)
:MGSurface(){
	size_t i=0, nRibNum = curves.size();
	if(nRibNum < 2)
		return;

	MGPvector<MGLBRep> lbs=rebuild_knot(curves);	//���r���h���s
	MGLBRep& lbs0=*(lbs[0]);
	if(direction_adjustment){
		double tmid=(lbs0.param_s()+lbs0.param_e())*.5;//param of mid point.
		MGVector dir1=lbs0.eval(tmid,1);
		for(i=1; i<nRibNum; i++){
			MGLBRep& lbsi=*(lbs[i]);
			MGVector dir2=lbsi.eval(tmid,1);
			if(dir1%dir2<0.){//If directions are opposite.
				lbsi.negate();
				dir2.negate();
			}
			dir1=dir2;
		}
	}
	m_uknot = lbs0.knot_vector();
	
	MGNDDArray vtau;
	createRibDataPoints(lbs,vtau);
	size_t lenu = m_uknot.bdim(), lenv = vtau.length();
	size_t nvOrder = 4;//Default order is 4.
	if(lenv < nvOrder)
		nvOrder = lenv;
	m_vknot=MGKnotVector(vtau,nvOrder);

	//�ʍ쐬����
	size_t len=lenu;
	if(lenv > len)
		len = lenv;	//len��lenu,lenv�̑傫�����̒l������
	size_t lenuv = lenu * lenv;
	size_t nOrder = m_uknot.order(), nSdim = lbs[0]->sdim();
	if(nOrder<nvOrder)
		nOrder=nvOrder;
	m_surface_bcoef = MGSPointSeq(lenu, lenv, nSdim);
	double *q = new double[len * (2 * nOrder - 1)];
	double *wk = new double[len * 2];
	double *wk2 = new double[lenuv * nSdim];

	//�Ȑ����猳�ƂȂ�R���g���[���|�C���g���쐬����
	for(i=0; i<lenu; i++){
		for(size_t j=0; j<lenv; j++){
			for(size_t k=0; k<nSdim; k++){
				wk2[j+lenv*i+lenuv*k] = lbs[j]->coef(i,k);
			}
		}
	}
	int blg4spError = 2;
	for(i=0; i<nSdim; i++){
		blgi2d_(&blg4spError,vtau.data(),wk2+lenuv*i,m_vknot.data(),
			nvOrder,lenv,lenu,lenv,lenu,wk,q,&m_surface_bcoef(0,0,i));
	}
	delete[] q;
	delete[] wk;
	delete[] wk2;
}
MGSBRep::MGSBRep(const std::vector<const MGCurve*>& curves,bool direction_adjustment){
	size_t n=curves.size();
	MGPvector<MGCurve> curves2(n);
	for(size_t i=0; i<n; i++) curves2[i]=const_cast<MGCurve*>(curves[i]);
	operator=(MGSBRep(curves2,direction_adjustment));
	curves2.release_all();
}

//���u�Ȑ��񂩂�ʂ��쐬����
//�쐬����ʂ̃m�b�g�x�N�g���̓��u�Ȑ��̌�����u,���u�������v�Ƃ���
//This constructor only generates MGSBRep even if curves are MGRLBRep of the same
//knot configuration. To avoid this, use createSurfaceFromRibs() that generates
//MGRSBRep when curves are MGRLBRep of the same knot configuration.
//Let v0=start parameter value, v1=terminate parameter value along v, then
//v=v0 const parameter line is curves[0], and v=v1 const parameter line is
//curves[n-1], where n=curves.size(). n must be greater or equal to 2.
//When n==2, the surface is a ruled surface(that is, this->order_u() is 2).
MGRSBRep::MGRSBRep(
	const std::vector<const MGRLBRep*>& vecPtrRibRLBReps,
	bool direction_adjustment//=true, curves[.] direction are adjusted to line
								//to the same direction.
){
	size_t i=0, nRibNum=vecPtrRibRLBReps.size();
	std::vector<const MGCurve*> vecPtrRibLBReps(nRibNum);
	for(i=0; i<nRibNum; i++){
		vecPtrRibLBReps[i] = &vecPtrRibRLBReps[i]->homogeneous();
	}
	*this = MGRSBRep(MGSBRep(vecPtrRibLBReps),direction_adjustment);
}
