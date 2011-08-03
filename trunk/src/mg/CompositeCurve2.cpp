/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/CCisect_list.h"
#include "mg/CompositeCurve.h"
#include "mg/LBRep.h"
#include "mg/Plane.h"
#include "mg/SurfCurve.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//Implements MGCompositeCurve Class.
//MGCompositeCurve is a composite of other leaf curves.
//Assumedly they are connected as C0 continuity. However, MGCompositeCurve
//does not check their continuity, but only put first or last as the user says
// (in connect_to_end or in connect_to_start), except two curves connecting
//are both MGLBRep or MGRLBRep. When the two connecting curves are both
//MGLBRep or MGRLBRep and they have more than C0 continuity,
//they are changed to one MGLBRep orMGRLBRep representation.
//Parameter ranges of the member curves are always continuous, from param_s() of the
//1st curve to param_e() of the last form MGCompositeCurve's paramter range.

//�֐����F
//		common
//�ړI�F
//		�^����ꂽ�Ȑ��Ǝ��g�̌�_�������͋��ʕ��������邩�ǂ������ׂ�B
//�����F
//		const MGCurve&			crv2,		(I/ )	�^������Ȑ�
//		std::vector<double>&	vec_param	( /O)	���ʕ����̃p�����[�^�͈�
//		MGCCisect_list&			isect		( /O)	��_
//				 4n�̔z��ŁAt(4*i+0),t(4*i+1)�����g�̃p�����[�^�͈�(t(4*i+0) < t(4*i+1))�A
//							 t(4*i+2),t(4*i+3)���^�Ȑ��̃p�����[�^�͈�(f(t(4*i+0))=f(t(4*i+2))
//�߂�l�F
//		3:��_�����ʕ��������܂���
//		2:��_�݂̂����܂���
//		1:���ʕ����݂̂����܂���
//		0:��_�����ʕ������Ȃ�����
//		-1:���ʃG�b�W�̎����v�Z�G���[
//		-2:���ʃG�b�W���S�ȏ㋁�܂���(�̂��Ă��Ȃ��ƌ��Ȃ�)
//�ǋL�F
//	�Ȑ������ʂ��ǂ����̌덷�ɂ�line_zero()�A���p�����[�^�͈͂̎����v�Z��
//	�덷�ɂ́A�p�����[�^�͈�*rc_zero()���g�p����
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan,
	MGCCisect_list& isect
)const{
	int obtainedIsect=0, obtainedCommon=0;
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		MGCCisect_list isect2;
		int obtained2=(**i).common(crv2,vecComSpan2,isect2);
		if(obtained2<0)
			return obtained2;
		int ncom=vecComSpan2.size();
		if(ncom){
			obtainedCommon=1;
			for(int j=0; j<ncom; j++)
				vecComSpan.push_back(vecComSpan2[j]);
		}
		int nisect=isect2.size();
		if(nisect){
			obtainedIsect=2;
			isect.append(isect2);
		}
	}
	return obtainedCommon+obtainedIsect;
}

//�֐����F
//		common
//�ړI�F
//		�^����ꂽ�Ȑ��Ǝ��g�̋��ʕ��������邩�ǂ������ׂ�B
//�����F
//		const MGCurve&			crv2,		(I/ )	�^������Ȑ�
//		std::vector<double>&	vec_param	( /O)	���ʕ����̃p�����[�^�͈�
//				 4n�̔z��ŁAt(4*i+0),t(4*i+1)�����g�̃p�����[�^�͈�(t(4*i+0) < t(4*i+1))�A
//							 t(4*i+2),t(4*i+3)���^�Ȑ��̃p�����[�^�͈�(f(t(4*i+0))=f(t(4*i+2))
//�߂�l�F
//		���ʕ����̐�:	���ʕ��������܂���
//		0:				���ʕ������Ȃ�����
//		-1:				���ʃG�b�W�̎����v�Z�G���[
//		-2:				���ʃG�b�W���S�ȏ㋁�܂���(�̂��Ă��Ȃ��ƌ��Ȃ�)
//�ǋL�F
//	�Ȑ������ʂ��ǂ����̌덷�ɂ�line_zero()���A�p�����[�^�͈͂̎����v�Z�̌덷�ɂ́A
//  �p�����[�^�͈�*rc_zero()���g�p����
int MGCompositeCurve::common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++){
		std::vector<double> vecComSpan2;
		int obtained2=(**i).common(crv2,vecComSpan2);
		if(obtained2<0)
			return obtained2;
		int ncom=vecComSpan2.size();
		for(int j=0; j<ncom; j++)
			vecComSpan.push_back(vecComSpan2[j]);

	}
	return vecComSpan.size()/4;
}

//���I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
//costant offset curve. if the norm_vector is given, the positive offset direction decide
//to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGPvector<MGCurve> MGCompositeCurve::offset(
	double ofs_value,			//�I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{	
return MGPvector<MGCurve>(static_cast<MGCurve*>(0));
	}

//�σI�t�Z�b�g�֐�
//�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
//valuable offset curve. if the norm_vector is given, the positive offset direction decide
//to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGPvector<MGCurve> MGCompositeCurve::offset(
	const MGLBRep& ofs_value_lb,					//��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	return MGPvector<MGCurve>(static_cast<MGCurve*>(0));
}

//C2�A���Ȑ��̈��I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
//costant offset curve of C2 continuous curve. if the norm_vector is given, the positive offset direction
//decide to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGLBRep MGCompositeCurve::offset_c2(
	double ofs_value,			//�I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	return MGLBRep();
}

//C2�A���Ȑ��̉σI�t�Z�b�g�֐�
//�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
//�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
//�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
//�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
//valuable offset curveof C2 continuous curve. if the norm_vector is given, the positive offset direction
//decide to left hand side from ahead, or the direction to center of curvature at start parameter.
//the offset value is less than radius of curvature. line_zero() is used.
MGLBRep MGCompositeCurve::offset_c2(
	const MGLBRep& ofs_value_lb,					//��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector	//�@���x�N�g��
)const{
	return MGLBRep();
}

//Return sweep surface from crv
//Returned is a newed MGSurface, must be deleted.
MGSurface* MGCompositeCurve::sweep(
	const MGUnit_vector& uvec,			//Sweep Direction.
	double start_dist,					//distance to start edge.
	double end_dist			//distance to end edge.
)const{
	return new MGPlane();
}

//Obtain the projected curve of this onto the surface srf.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
int MGCompositeCurve::project_onto_surface(
	const MGFSurface& srf,
	MGPvector<MGCurve>& vec_crv_uv,	
		//Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	const_iterator i=begin(), ie=end();
	for(; i!=ie; i++)
		srf.project(**i,vec_crv_uv,vec_crv,vec);
	return vec_crv.size();
}
