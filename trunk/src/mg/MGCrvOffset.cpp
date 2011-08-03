/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Knot.h"
#include "mg/Curve.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Straight.h"
#include "mg/Pvector.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of curve offset.

#define POW_LOW 0.34	//�ŏ��ɃI�t�Z�b�g������|�C���g�������߂�W��(���x�d��)
#define POW_HIGH 0.50   //�ŏ��ɃI�t�Z�b�g������|�C���g�������߂�W��(���x�d��)
#define NUM_DIV 50      //1�X�p���̕�������������z�����Ƃ�POW_LOW���g�p����悤�ɂ���

//����������2�{��B�\���Ȑ���ڑ�����(������ނ̂Ƃ�)
MGLBRep* join2LBRep(const MGLBRep& crv1, const MGLBRep& crv2){
	int which, cont;
	double ratio;
	cont = crv1.continuity(crv2, which, ratio);
	if(cont < 0 || which == 0 || which == 3)
		return NULL;

	return new MGLBRep(crv1, cont, which, crv2);
}

//���I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B
//�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B
//�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
MGPvector<MGCurve> MGCurve::offset(
	double ofs_value,			//�I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	int error;
	MGBPointSeq bp1(2, 1);					//ofs_value���̒����𐶐�����
	bp1.store_at(0, &ofs_value);
	bp1.store_at(1, &ofs_value);
	MGLBRep ofs_value_lb(bp1, error, 2);	//�I�[�_�[�Q�̒��������
	if(error){
		return MGPvector<MGCurve>();
	}
	ofs_value_lb.change_range(param_s(), param_e());
	return offset(ofs_value_lb, norm_vector);	//�σI�t�Z�b�g���g�p����
}

//�σI�t�Z�b�g�֐�
//�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B
//�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B
//�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
MGPvector<MGCurve> MGCurve::offset(
	const MGLBRep& ofs_value_lb,	//��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector		//�@���x�N�g��
)const{
	MGPvector<MGCurve> ofs_crvl;
	if(norm_vector == mgNULL_VEC){
		offset_proc(ofs_value_lb, ofs_crvl);
	}else{
		offset_norm_proc(ofs_value_lb, norm_vector, ofs_crvl);
	}
	return ofs_crvl;
}

//���I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B
//�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B
//�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
MGPvector<MGCurve> MGStraight::offset(
	double ofs_value,			//�I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	MGPvector<MGCurve> ofs_crvl;
	MGStraight *st1 = clone();
	if(norm_vector == mgNULL_VEC){
		*st1 += MGUnit_vector() * ofs_value;
	}else{
		MGUnit_vector tangent = st1->direction(param_s()), ofs_dir;
		ofs_dir = norm_vector * tangent;
		*st1 += ofs_dir * ofs_value;
	}
	ofs_crvl.push_back(st1);
	return ofs_crvl;
}

//�@���x�N�g�����w�肳��Ă���Ƃ��̏���
int MGCurve::offset_norm_proc(
	const MGLBRep& ofs_value_lb,	//�I�t�Z�b�g��
	const MGVector& norm_vector,	//�@���x�N�g��
	MGPvector<MGCurve>& ofs_crvl	//�I�t�Z�b�g�J�[�u���X�g
)const{
	//�Ȑ���܂ꂪ��������imultiple knot�j�ŕ�������
	MGPvector<MGCurve> crv_list, tmp_list;
	int nDiv = divide_multi_ofs(ofs_value_lb, crv_list);

	ofs_crvl.clear();
	for(int i = 0;i < nDiv; i++){
		MGLBRep tmp_brep;
		if(!crv_list[i]->offset_norm_c2_proc(ofs_value_lb, norm_vector, tmp_brep))return false;
		tmp_list.push_back(tmp_brep.clone());
	}
	int num = join(tmp_list, ofs_crvl);
	return num;
}

//�@���x�N�g�����w�肳��Ă��Ȃ��Ƃ��̏���
int MGCurve::offset_proc(
	const MGLBRep& ofs_value_lb,		//�I�t�Z�b�g��
	MGPvector<MGCurve>& ofs_crvl	//�I�t�Z�b�g�J�[�u���X�g
)const{
	//�Ȑ���܂�imultiple knot�j�ŕ�������
	MGPvector<MGCurve> crv_list, tmp_list;
	int nDiv = divide_multi_ofs(ofs_value_lb, crv_list);

	ofs_crvl.clear();
	MGUnit_vector T, B, preN;	//preN : �O��̃m�[�}���x�N�g��
	int freverse = 0;			//�������t�ɂ��Ă���Ƃ����t���O
	double curvature, torsion;
	Frenet_frame(param_s(), T, preN, B, curvature, torsion);	//�J�[�u�n�_�̃m�[�}�����ŏ��ɗ^����

	for(int i = 0;i < nDiv; i++){
		MGLBRep tmp_brep;
		if(!crv_list[i]->offset_c2_proc(ofs_value_lb, tmp_brep, preN, freverse))return false;
		tmp_list.push_back(tmp_brep.clone());
	}
	int num = join(tmp_list, ofs_crvl);
	return num;
}

//C2�A���Ȑ��̈��I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B
//�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B
//�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
MGLBRep MGCurve::offset_c2(
	double ofs_value,				//�I�t�Z�b�g��
	const MGVector& norm_vector		//�@���x�N�g��
)const{
	int error;
	MGBPointSeq bp1(2, 1);					//ofs_value���̒����𐶐�����
	bp1.store_at(0, &ofs_value);
	bp1.store_at(1, &ofs_value);
	MGLBRep ofs_value_lb(bp1, error, 2);	//�I�[�_�[�Q�̒��������
	if(error)return MGLBRep();
	ofs_value_lb.change_range(param_s(), param_e());
	return offset_c2(ofs_value_lb, norm_vector);	//�σI�t�Z�b�g���g�p����
}

//C2�A���Ȑ��̉σI�t�Z�b�g�֐�
//�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B
//�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B
//�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
MGLBRep MGCurve::offset_c2(
	const MGLBRep& ofs_value_lb,	//��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector		//�@���x�N�g��
)const{
	int rc = 0;
	MGLBRep ofs_brep;
	if(norm_vector == mgNULL_VEC){
		int freverse = 0;			//�������t�ɂ��Ă���Ƃ����t���O
		MGUnit_vector T, B, preN;	//preN : �O��̃m�[�}���x�N�g��
		double curvature, torsion;
		Frenet_frame(param_s(), T, preN, B, curvature, torsion);	//�J�[�u�n�_�̃m�[�}�����ŏ��ɗ^����
		rc = offset_c2_proc(ofs_value_lb, ofs_brep, preN, freverse);
	}else{
		rc = offset_norm_c2_proc(ofs_value_lb, norm_vector, ofs_brep);
	}
	if(!rc)return MGLBRep();
	ofs_brep.remove_knot();
	return ofs_brep;
}

//�@���x�N�g�����w�肳��Ă���C2�A���Ȑ��̃I�t�Z�b�g
int MGCurve::offset_norm_c2_proc(
	const MGLBRep& ofs_value_lb,//�I�t�Z�b�g��
	const MGVector& norm_vector,//�@���x�N�g��
	MGLBRep& ofs_brep			//�I�t�Z�b�g�J�[�u
)const{
	MGKnotVector knotVector = offset_make_knotvector(ofs_value_lb);	//�\�����������m�b�g�x�N�g�������߂�
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);

	//����_�𐶐�����
	int len = dataPoint.length();
	MGBPointSeq bp1(len, sdim());
	for(int i = 0; i < len; i++){
		MGUnit_vector T, N, B, ofs_dir;
		double curvature, torsion, ofs_value = ofs_value_lb.eval_position(dataPoint(i)).ref(0);
		Frenet_frame(dataPoint(i), T, N, B, curvature, torsion);
		ofs_dir = norm_vector * T;
		double cosine = N % ofs_dir;
		int fneg = (cosine < 0);	//�R�T�C���̐����t���O
		if(!MGMZero(curvature))		//�[�����`�F�b�N
			if((ofs_value * (-2. * fneg + 1)) > ((1. / curvature) / fabs(cosine)))return false;
		MGPosition pos;
		pos = eval(dataPoint(i));
		pos += ofs_dir * ofs_value;
		bp1.store_at(i, pos);
	}
	//�m�b�g�x�N�g�������߂ăI�t�Z�b�g�Ȑ��𐶐�����
	int error;
	ofs_brep = MGLBRep(dataPoint, bp1, knotVector, error);	//���x�\���̋Ȑ��𐶐�����
	if(error)return false;
	return true;
}

//�@���x�N�g�����w�肳��Ă��Ȃ�C2�A���Ȑ��̃I�t�Z�b�g
int MGCurve::offset_c2_proc(
	const MGLBRep& ofs_value_lb,	//�I�t�Z�b�g��
	MGLBRep& ofs_brep,			//�I�t�Z�b�g�J�[�u
	MGUnit_vector& preN,		//�O��̃m�[�}���x�N�g��
	int& freverse			//�������t�ɂ��Ă���Ƃ����t���O
)const{
	MGKnotVector knotVector = offset_make_knotvector(ofs_value_lb);	//�\�����������m�b�g�x�N�g�������߂�
	MGNDDArray dataPoint;
	dataPoint.update_from_knot(knotVector);

	//����_�𐶐�����
	int len = dataPoint.length();
	MGBPointSeq bp1(len, sdim());
	for(int i = 0; i < len; i++){
		MGUnit_vector T, N, B;
		double curvature, torsion, ofs_value = ofs_value_lb.eval_position(dataPoint(i)).ref(0);
		Frenet_frame(dataPoint(i), T, N, B, curvature, torsion);
		if(MGAZero(curvature))N = preN;		//�ȗ����������Ƃ��m�[�}���͑O�̂��g��
		//�m�[�}�����S�����������ɂȂ�悤�ɂ���(180�x�ȏ�J���Ȃ�)
		if((preN % N) < 0)freverse = !(freverse);
		int fdir = -2 * freverse + 1;	//freverse��false�̂Ƃ�1, true�̂Ƃ�-1
		if(!MGMZero(curvature))	//�[�����`�F�b�N
			if(ofs_value * fdir > (1. / curvature))return false;	//�ȗ����a���傫���I�t�Z�b�g�͔F�߂Ȃ�
		preN = N;	//�m�[�}��������Ă���
		MGPosition pos;
		pos = eval(dataPoint(i));
		pos += fdir * N * ofs_value;
		bp1.store_at(i, pos);
	}

	//�m�b�g�x�N�g�������߂ăI�t�Z�b�g�Ȑ��𐶐�����
	int error;
	ofs_brep = MGLBRep(dataPoint, bp1, knotVector, error);	//���x�\���̋Ȑ��𐶐�����
	if(error)return false;
	return true;
}

//�Ȑ���܂�ŕ�������(�I�t�Z�b�g�ʂ������Ȑ��̐܂���Ă���������)
int MGCurve::divide_multi_ofs(
	const MGLBRep& ofs_value_lb,		//�I�t�Z�b�g�ʂ������Ȑ�
	MGPvector<MGCurve>& crv_list		//���������Ȑ����X�g
)const{
	crv_list.clear();
	//�I�t�Z�b�g�ʂ������Ȑ��ƌ��̋Ȑ��̃p�����[�^�����W���Ⴄ�ƃG���[
	if(param_range() != ofs_value_lb.param_range())return 0;
	size_t	start_index = ofs_value_lb.order() - 1, index = 0, count = 0, multi = 0, vbdim = ofs_value_lb.bdim();
	MGKnotVector value_knot = ofs_value_lb.knot_vector();
	do{
		if(ofs_value_lb.order() == 2){	//�I�[�_�[2�̂Ƃ��̏���
			index = start_index + 1; multi = 1;
		}else{
			multi = value_knot.locate_multi(start_index, 2, index);
		}
		MGCurve *temp_crv = part(value_knot(start_index), value_knot(index));
		count += temp_crv->divide_multi(crv_list);
		delete temp_crv;
		start_index = index + multi - 1;
	}while(index != vbdim && value_knot(start_index) < param_e());	//���d�x��������Ȃ�������I���
	return count;
}

//1�X�p���̕����������߂�
int MGCurve::offset_div_num(
	const MGInterval& interval	//�����������߂�p�����[�^�͈�
)const{
	//�Q��order()�ŕ������čő�2�������l�����߂�
	int ord = order();
	if(ord<=2) ord = 4;		//Ellipse, Straight�̂Ƃ�
	int ord2=ord*2;
	double max_deriv = 0.0, tpara = interval.low_point();
	double shortspan = (interval.high_point() - tpara)/ord2;
	double oneSpanLength = interval.length().value();
	double oneSpanLength2=oneSpanLength*oneSpanLength;
	for(int j = 0; j <=ord2; j++){
		//�p�����[�^�͈�(0,1)��2�������l�ɒ������߂Ƀp�����[�^�͈͂�2���������
		double deriv = eval(tpara, 2).len() * oneSpanLength2;
		if(deriv > max_deriv) max_deriv = deriv;
		tpara += shortspan;
	}

	double onelzero=1./MGTolerance::line_zero();
	//�������� (1 / tol)^POW * sqrt(max_deriv / 8)�ōŏ��l�� order()�ł���
	double sqrt_max_deriv8=sqrt(max_deriv / 8.);
	int ndiv = int(pow(onelzero, POW_HIGH) * sqrt_max_deriv8);
	if(ndiv > NUM_DIV)
		ndiv = int(pow(onelzero, POW_LOW) * sqrt_max_deriv8);
	if(ndiv < ord2) ndiv = ord2;
	else if(ndiv > NUM_DIV*2) ndiv=NUM_DIV*2;
	return ndiv;
}

//����������B�\���Ȑ����X�g��ڑ�����(LBRep�̂�)�Bjoin_crvl�ɐڑ������Ȑ����X�g������B
//�߂�l�́A�����̋Ȑ����X�g�̌������Ⴄ�Ƃ��A����B�\�����m�łȂ������Ƃ�false���Ԃ�B
int join(
	MGPvector<MGCurve>& crvl,
	MGPvector<MGCurve>& join_crvl
){
	int num = crvl.size(), rc = 0;
	MGCurve *pre_pcrv = 0, *cur_pcrv, *next_pcrv;
	if(!num)
		return false;

	if(num == 1){	//�Ȑ����P�{�̂Ƃ��̏���
		cur_pcrv = crvl[0]->clone();
		cur_pcrv->remove_knot();
		join_crvl.push_back(cur_pcrv);
		return 1;
	}

	cur_pcrv = crvl[0]->clone();
	for(int i = 0; i < num - 1; i++){
		next_pcrv = crvl[i+1];
		if(!cur_pcrv || !next_pcrv)
			return false;

		MGLBRep* lb1=dynamic_cast<MGLBRep*>(cur_pcrv);
		MGLBRep* lb2=dynamic_cast<MGLBRep*>(next_pcrv);
		if(!lb1 || !lb2)
			return false;

		pre_pcrv = join2LBRep(*lb1,*lb2);
		if(!pre_pcrv){
			rc++;
			cur_pcrv->remove_knot();
			join_crvl.push_back(cur_pcrv);
			delete pre_pcrv;
			cur_pcrv = next_pcrv->clone();
			continue;
		}
		delete cur_pcrv;
		cur_pcrv = pre_pcrv;
	}
	cur_pcrv->remove_knot();
	join_crvl.push_back(cur_pcrv);
	rc++;
	return rc;
}

/*//����������2�{��B�\���Ȑ���ڑ�����(������ނ̂Ƃ�)
MGCurve* MGCurve::join(const MGCurve& crv1)const{
	return NULL;
}*/

/*//����������2�{��B�\���Ȑ���ڑ�����(������ނ̂Ƃ�)
MGCurve* MGRLBRep::join(const MGCurve& crv1) const{
	const MGRLBRep *prlbrep1 = dynamic_cast<const MGRLBRep*>(&crv1);
	if(!prlbrep1)return NULL;
	int which, cont;
	double ratio;
	cont = this->continuity(*prlbrep1, which, ratio);
	if(cont < 0 || which == 0 || which == 3)return NULL;
	return new MGRLBRep(*this, cont, which, *prlbrep1);
}*/

//�Ȑ����I�t�Z�b�g����̂ɏ\�����������m�b�g�x�N�g����ԋp����
//�I�t�Z�b�g�ʋȐ����l���ɓ���A�������̑������ɂ��킹�Ă���B
MGKnotVector MGCurve::offset_make_knotvector(
	const MGLBRep& ofs_value_lb
)const{
	//���ƂȂ�m�b�g�x�N�g���𐶐�����
	const MGKnotVector& tempKnotVector = knot_vector();
	MGKnotVector knotVector(tempKnotVector, 4);	//�I�[�_�[4�̃m�b�g�x�N�g���ɍ��ւ���
	for(size_t i=tempKnotVector.order()-1; i<tempKnotVector.bdim(); i++){
		double	tpara = 0.0,					//�e���|����
				spara = tempKnotVector(i),		//�X�p���̎n�_
				epara = tempKnotVector(i + 1);	//�X�p���̏I�_
		if(epara - spara < param_error())
			continue;	//�}���`�m�b�g�̂Ƃ��̏���(RLBRep�̂�)

		//1�X�p���̕����������肷��
		MGInterval interval(spara, epara);
		size_t ndiv = offset_div_num(interval),					//�I�t�Z�b�g�Ȑ��̕�����
			tmp_ndiv = ofs_value_lb.offset_div_num(interval);	//�I�t�Z�b�g�ʋȐ��̕�����
		if(tmp_ndiv>ndiv)
			ndiv=tmp_ndiv;						//�������̑�������p����
		double shortspan = (epara - spara) / ndiv;
		tpara = spara;
		for(size_t j=0; j<ndiv; j++){
			knotVector.add_data(tpara);
			tpara += shortspan;
		}
	}
	return knotVector;
}
