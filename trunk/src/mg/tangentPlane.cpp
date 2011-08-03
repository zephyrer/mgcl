/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// tangentCurve.cpp : MGCurve��MGSurface����Tangent Plane�����֐�
//
#include "MGCLStdAfx.h"
#include "mg/LBRep.h"
#include "mg/Surface.h"
#include "mg/Position.h"
#include "mg/SurfCurve.h"
#include "mg/Interval.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MGLBRep TP_from_world_curve(
	const MGSurface& srf,	//�ڑ���
	const MGCurve& crv,		//�ڑ��Ȑ�
	size_t order,			//�쐬����J�[�u�̃I�[�_
	int& error)	//�G���[�R�[�h  0:OK, !=0 error code when constructing LBRep.
{
	if(order < 2)return MGLBRep();
	const MGKnotVector& tempKnotVector = crv.knot_vector();
	MGKnotVector knotVector(tempKnotVector, order);	//�w��I�[�_�[�̃m�b�g�x�N�g���ɍ��ւ���
	size_t i;
	for(i = tempKnotVector.order() - 1; i < tempKnotVector.bdim(); i++){
		double	tpara = 0.0,					//�e���|����
				spara = tempKnotVector(i),		//�X�p���̎n�_
				epara = tempKnotVector(i + 1);	//�X�p���̏I�_
		if(epara - spara < crv.param_error())continue;	//�}���`�m�b�g�̂Ƃ��̏���(RLBRep�̂�)
		//1�X�p���̕����������肷��
		MGInterval interval(spara, epara);
		int ndiv = crv.calc_div_num(interval);			//�I�t�Z�b�g�Ȑ��̕�����
		double shortspan = (epara - spara) / ndiv;
		tpara = spara;
		for (int j = 0; j < ndiv; j++){knotVector.add_data(tpara); tpara += shortspan;}
	}

	// on���g����� perps_guess���g�����ق��������̂�
	// ��ɑ��_�����߂Ă���
	MGPosition guess_prm(2), on_prm(2);
	MGPosition crv_pt(crv.start_point()), srf_tan(srf.sdim());
	MGPosition uv_min(srf.param_s_u(), srf.param_s_v()), 
			uv_max(srf.param_e_u(), srf.param_e_v());
	srf.on(crv_pt, guess_prm);
	//����_�𐶐�����
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);
	int len = dataPoint.length();
	MGBPointSeq bp1(len, srf.sdim());
	for(int j = 0; j < len; j++){
		srf.perp_guess(uv_min, uv_max, 
			crv.eval(dataPoint(j)),guess_prm, on_prm);
		MGPosition pos = srf.unit_normal(on_prm);
		guess_prm = on_prm;
		bp1.store_at(j, pos);
	}

	//���x�\���̋Ȑ��𐶐�����
	MGLBRep brep;
	brep = MGLBRep(dataPoint, bp1, knotVector, error);
	//�]����Knot���폜
	MGTolerance::push();
	// 2.0* sin��/2 �� �� �Ƃ���B
	MGTolerance::set_line_zero(MGTolerance::angle_zero());
	brep.remove_knot();
	MGTolerance::pop();
	return brep;
}

MGLBRep TP_from_parameter_curve(
	const MGSurface& srf,	//�ڑ���
	const MGCurve& pcrv,	//�ڑ��Ȑ�(�ڑ��ʏ�̃p�����[�^�J�[�u)
	const MGCurve& wcrv,	//�ڑ��Ȑ�(���E���W��̃J�[�u)
	size_t order,			//�쐬����J�[�u�̃I�[�_
	int& error)	//�G���[�R�[�h  0:OK, !=0 error code when constructing LBRep.
{
	// WorldCurve����\���ׂ���KnotVector���쐬����B
	if(order < 2)return MGLBRep();
	const MGKnotVector& tempKnotVector = wcrv.knot_vector();
	MGKnotVector knotVector(tempKnotVector, order);	//�w��I�[�_�[�̃m�b�g�x�N�g���ɍ��ւ���
	size_t i;
	for(i = tempKnotVector.order() - 1; i < tempKnotVector.bdim(); i++){
		double	tpara = 0.0,					//�e���|����
				spara = tempKnotVector(i),		//�X�p���̎n�_
				epara = tempKnotVector(i + 1);	//�X�p���̏I�_
		if(epara - spara < wcrv.param_error())continue;	//�}���`�m�b�g�̂Ƃ��̏���(RLBRep�̂�)
		//1�X�p���̕����������肷��
		MGInterval interval(spara, epara);
		int ndiv = wcrv.calc_div_num(interval);			//�I�t�Z�b�g�Ȑ��̕��������p
		double shortspan = (epara - spara) / ndiv;
		tpara = spara;
		for (int j = 0; j < ndiv; j++){knotVector.add_data(tpara); tpara += shortspan;}
	}

	//����_�𐶐�����
	// SurfCurve(src,pcrv)��on()��
	// wcrv�̈���p�����[�^t�ɑΉ�����pcrv�̃p�����[�^�lpt������
	// pt����Surface��uv�l�����Normal�����߂�
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);
	int len = dataPoint.length();
	MGBPointSeq bp1(len, srf.sdim());
	MGSurfCurve scrv(srf, pcrv);
	for(int j = 0; j < len; j++){
		double pt;
		scrv.on(wcrv.eval(dataPoint(j)),pt);
		MGPosition pos = srf.unit_normal(pcrv.eval(pt));
		bp1.store_at(j, pos);
	}

	// �璷��LBRep�����
	MGLBRep brep;
	brep = MGLBRep(dataPoint, bp1, knotVector, error);

	//�]���ȃm�b�g�𗎂Ƃ�
	MGTolerance::push();
	// 2.0* sin��/2 �� �� �Ƃ���B
	MGTolerance::set_line_zero(MGTolerance::angle_zero());
	brep.remove_knot();
	MGTolerance::pop();
	return brep;
}
