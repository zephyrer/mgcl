/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Straight.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SPointSeq.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of surface offset.

#define POW_LOW 0.34   //�ŏ��ɃI�t�Z�b�g������|�C���g�������߂�W��(���x�d��)
#define POW_HIGH 0.50  //�ŏ��ɃI�t�Z�b�g������|�C���g�������߂�W��(���x�d��)
#define NUM_DIV 25     //1�p�b�`�̕�������������z�����Ƃ�POW_LOW���g�p����悤�ɂ���

//���I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�m�[�}�������𐳂Ƃ���B�ȗ����a���傫���I�t�Z�b�g�͍s��Ȃ�
//�߂�l�́A�I�t�Z�b�g�Ȗʃ��X�g���ԋp�����B�G���[�̂Ƃ�����0�̋Ȗʃ��X�g���Ԃ�B
//�g�������X��line_zero()���g�p���Ă���B
MGPvector<MGSurface> MGSurface::offset(
	double ofs_value,			//�I�t�Z�b�g��
	int& error) const			//�G���[�R�[�h 0:���� -2:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
{
	//���u�����̐܂�̕����ŕ�������
	MGPvector<MGSurface> vecSrf, ofs_srfl, null_srf_vec;
	int div = divide_multi_knot(vecSrf);
	for(int i = 0; i < div; i++){
		std::auto_ptr<MGSurface> pofs_srf = vecSrf[i]->offset_c1(ofs_value, error);	//�I�t�Z�b�g���s��
		if(error < 0)return null_srf_vec;
		ofs_srfl.push_back(pofs_srf.release());
	}
	return ofs_srfl;
}

//Offset.
//distance is plus value if the direction is toward normal vector of the
//FSurface. Minus if against the normal vector.
//�G���[�R�[�h 0:���� -1:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
int MGSurface::offset_fs(double distance, MGPvector<MGFSurface>& vecOfsFSurface)const{
	int error;
	MGPvector<MGSurface> srfs=offset(distance,error);
	if(error) return error;
	size_t n=srfs.size();
	for(size_t i=0; i<n; i++)
		vecOfsFSurface.push_back(srfs.release(i));
	return 0;
}

//C1�A���Ȗʂ̈��I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�m�[�}�������𐳂Ƃ���B�ȗ����a���傫���I�t�Z�b�g�͍s��Ȃ��B
//�߂�l�́A�I�t�Z�b�g�����Ȗʂ̃I�[�g�|�C���^�ł���B�G���[�̂Ƃ��k�����Ԃ�B
//�g�������X��line_zero()���g�p���Ă���B
std::auto_ptr<MGSurface> MGSurface::offset_c1(
	double ofs_value,		//�I�t�Z�b�g��
	int& error) const		//�G���[�R�[�h 0:���� -1:�ʂɂ��ꂪ����
 					// -2:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
{
	//C0�A���̂Ƃ��G���[
	std::auto_ptr<MGSurface> pofs_srf;
	const MGKnotVector &knotu = knot_vector_u(), &knotv = knot_vector_v();	//SBRep, RSBRep�݂̂�����
	size_t startu = knotu.order() - 1, startv = knotv.order() - 1, index = 0;
	if(type() == MGSURFACE_SPLINE){	//�m�����V���i���̂Ƃ���������̃`�F�b�N���s��
		if(knotu.order() != 2){
			if(knotu.locate_multi(startu, knotu.order() - 1, index)){error = -1; return pofs_srf;}
		}else if(knotu.bdim() > 2){error = -1; return pofs_srf;}
		if(knotv.order() != 2){
			if(knotv.locate_multi(startv, knotv.order() - 1, index)){error = -1; return pofs_srf;}
		}else if(knotv.bdim() > 2){error = -1; return pofs_srf;}
	}

	//�I�t�Z�b�g���ȗ����a�ȓ����ǂ������ׂ�
	if(!offset_check_curva(ofs_value)){error = -2; return pofs_srf;}

	//�ʂ𕪊�����m�b�g�x�N�g���ƃf�[�^�|�C���g���v�Z����
	MGKnotVector knotVec_u, knotVec_v;
	MGNDDArray data_point_u, data_point_v;
	offset_calc_knot_vec(knotVec_u, knotVec_v, data_point_u, data_point_v);

	//�I�t�Z�b�g�ʂ��쐬����
	size_t	ulen = knotVec_u.bdim(), vlen = knotVec_v.bdim();
	MGSPointSeq sb1(ulen, vlen, sdim());
	for(size_t i = 0; i < ulen; i++){
		for(size_t j = 0; j < vlen; j++){
			double data_u = data_point_u(i), data_v = data_point_v(j);
			MGPosition pos, ofs_pos;
			pos = eval(data_u, data_v, 0, 0);
			MGVector N(unit_normal(data_u, data_v));
			ofs_pos = pos + (N* ofs_value);
			sb1.store_at(i, j, ofs_pos);
		}
	}
	int err = 0;
	MGSBRep ofs_sbrep = MGSBRep(data_point_u, data_point_v, sb1, knotVec_u, knotVec_v, error);
	if(err){error = -3; return pofs_srf;}
	ofs_sbrep.remove_knot();
	std::auto_ptr<MGSurface> tmp_surf = std::auto_ptr<MGSurface>(new MGSBRep(ofs_sbrep));
	pofs_srf = tmp_surf;
	error = 0;
	return pofs_srf;
}

//C1�A���Ȗʂ̈��I�t�Z�b�g�֐�
//�I�t�Z�b�g�����́A�m�[�}�������𐳂Ƃ���B�ȗ����a���傫���I�t�Z�b�g�͍s��Ȃ�
//�߂�l�́A�I�t�Z�b�g�����Ȗʂ̃I�[�g�|�C���^�ł���B�G���[�̂Ƃ��̓k�����Ԃ�B
//�g�������X��line_zero()���g�p���Ă���B
std::auto_ptr<MGSurface> MGPlane::offset_c1(
	double ofs_value,				//�I�t�Z�b�g��
	int& error) const				//�G���[�R�[�h 0:���� -1:�ʂɂ��ꂪ����
									// -2:�ȗ����a�ȏ�̃I�t�Z�b�g�s�� -3:�ʐ����R���X�g���N�^�G���[
{
	MGPlane ofs_plane = *this + (this->normal() * ofs_value);
	std::auto_ptr<MGSurface> pofs_srf(new MGPlane(ofs_plane));
	error = 0;
	return pofs_srf;
}

//C1�A���Ȗʂ̈��I�t�Z�b�g�֐��Ŏg�p����m�b�g�x�N�g���ƃf�[�^�|�C���g�����߂�֐�
void MGSurface::offset_calc_knot_vec(
	MGKnotVector& knot_vec_u,			//���܂���u�m�b�g�x�N�g��
	MGKnotVector& knot_vec_v,			//���܂���v�m�b�g�x�N�g��
	MGNDDArray& data_point_u,			//���܂���u�f�[�^�|�C���g
	MGNDDArray& data_point_v) const		//���܂���v�f�[�^�|�C���g
{
	int num_div = offset_div_num();		//2�������l�ɉ����ĕ����������߂�
	int ord_u = order_u(), ord_v = order_v();
	int knotNumU = (bdim_u() - ord_u + 1) * num_div, knotNumV = (bdim_v() - ord_v + 1) * num_div;
	if(ord_u < 4)ord_u = 4; if(ord_v < 4)ord_v = 4;
	knot_vec_u = MGKnotVector(knot_vector_u(), ord_u);
	knot_vec_v = MGKnotVector(knot_vector_v(), ord_v);
	knot_vec_u.change_knot_number(knotNumU);
	knot_vec_v.change_knot_number(knotNumV);
	data_point_u.update_from_knot(knot_vec_u);
	data_point_v.update_from_knot(knot_vec_v);
}

//�ȗ����a�ƃI�t�Z�b�g�l�̔�����s��
//error�̂Ƃ�false���ԋp�����
int MGSurface::offset_check_curva(
	double ofs_value) const		//�I�t�Z�b�g��
{
	offset_check_curva_one(ofs_value);
	MGSurface *temp_srf = this->copy_surface();
	temp_srf->exchange_uv();

	int rtn = temp_srf->offset_check_curva_one(-ofs_value);	//uv����ւ����̂�normal���t�ɂȂ�
	delete temp_srf;
	return rtn;
}

//�ȗ����a�ƃI�t�Z�b�g�l�̔�����s��(u=const�̃p�����[�^�Ȑ����g�p����)
//error�̂Ƃ�false���ԋp�����
int MGSurface::offset_check_curva_one(
	double ofs_value) const		//�I�t�Z�b�g��
{
	//u�����X�p����2*order�ŕ�������p�����[�^�Ȑ��̋ȗ����a���I�t�Z�b�g�l�����傫�����ǂ������ׂ�
	size_t i, j, k, l;
	for(i = order_u() - 1; i < bdim_u(); i++){
		double uspan = (knot_u(i + 1) - knot_u(i)) / (2. * order_u());
		size_t ndiv_u = 2 * order_u();
		if(i == bdim_u() - 1)ndiv_u++;
		for(k = 0; k < ndiv_u; k++){
			double uparam = knot_u(i) + (uspan * k);
			MGCurve* pcrv = parameter_curve(1, uparam);	//u = const�̃p�����[�^�Ȑ�
			for(j = order_v() - 1; j < bdim_v(); j++){
				double vspan = (pcrv->knot(j + 1) - pcrv->knot(j)) / (2. * pcrv->order());
				size_t ndiv_v = 2 * pcrv->order();
				if(j == bdim_v() - 1)ndiv_v++;
				for(l = 0; l < ndiv_v; l++){
					double vparam = pcrv->knot(j) + (vspan * l), curva, torsion, radius;
					MGUnit_vector T, N, B, norm;
					pcrv->Frenet_frame(vparam, T, N, B, curva, torsion);
					norm = unit_normal(uparam, vparam);
					if((N % norm) * ofs_value <= 0 || MGMZero(curva))continue;
					radius = 1. / curva;
					if((radius * radius) < (ofs_value * ofs_value)){delete pcrv; return false;}
				}
			}
			delete pcrv;
		}
	}
	return true;
}

//�I�t�Z�b�g����T���v���|�C���g��1�p�b�`���Ƃ̕����������߂�
//�S�Ẵp�b�`���̕������ōő�̒l��Ԃ�
int MGSurface::offset_div_num() const{
	int max_div = 0;
	int bdu=bdim_u(), bdv=bdim_v();
	for(int i = order_u() - 1; i < bdu; i++){
		MGInterval itvl_u(knot_u(i), knot_u(i + 1));
		for(int j = order_v() - 1; j < bdv; j++){
			MGInterval itvl_v(knot_v(j), knot_v(j + 1));
			int div = offset_div_num_one(MGBox(itvl_u, itvl_v));
			if(div > max_div)max_div = div;
		}
	}
	return max_div;
}

//�I�t�Z�b�g����T���v���|�C���g��1�p�b�`���Ƃ̕����������߂�
//�S�Ẵp�b�`���̕������ōő�̒l��Ԃ�
int MGPlane::offset_div_num() const
{
	return 1;
}

//�I�t�Z�b�g����T���v���|�C���g��1�p�b�`���Ƃ̕����������߂�
//������n = sqrt(1 / tol) * sqrt((M1 + 2 *M2 + M3) / 8)�͈ȏ�̎��ŋ��܂�
int MGSurface::offset_div_num_one(
	const MGBox& param_range) const{
	int orderu=order_u(), orderv=order_v();

	//calculate maximum second derivative
	int div_u = 2 * orderu, div_v = 2 * orderv;
	double span_u = fabs(param_range.high().ref(0) - param_range.low().ref(0)) / double(div_u);
	double span_v = fabs(param_range.high().ref(1) - param_range.low().ref(1)) / double(div_v);
	double start_u = param_range.low().ref(0), start_v = param_range.low().ref(1);
	double M1 = 0.0, M2 = 0.0, M3 = 0.0;	//M1: Maximum Second derivative of Suu, M2: Suv, M3: Svv
	for(int i = 0; i <= div_u; i++){
		for(int j = 0; j <= div_v; j++){
			double param_u = start_u + (span_u * i), param_v = start_v + (span_v * j);
			double d1 = eval(param_u, param_v, 2, 0).len() * span_u * span_u;	//Second derivative of Suu
			double d2 = eval(param_u, param_v, 2, 2).len() * span_v * span_v;	//Second derivative of Suv
			double d3 = eval(param_u, param_v, 0, 2).len() * span_u * span_v;	//Second derivative of Svv
			if(d1 > M1)M1 = d1;
			if(d2 > M2)M2 = d2;
			if(d3 > M3)M3 = d3;
		}
	}
	int n = int(pow((1. / MGTolerance::line_zero()), POW_HIGH) * sqrt((M1 + (2. * M2) + M3) / 8.));
	if(n > NUM_DIV)n = int(pow((1. / MGTolerance::line_zero()), POW_LOW) * sqrt((M1 + (2. * M2) + M3) / 8.));
	if(n < 2*orderu) n=2*orderu;
	if(n < 2*orderv) n=2*orderv;
	return n;
}

//u�܂���v�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
//�߂�l�́A��������ԋp����
int MGSBRep::divide_multi_knot(
	MGPvector<MGSurface>& srfl) const		//���������Ȗʃ��X�g
{
	//���u�����̐܂�̕����ŕ�������
	MGPvector<MGSBRep> srf_u;
	int udiv = divide_multi_knot_u(srf_u), num = 0;
	for(int i = 0; i < udiv; i++){
		//����v�����̐܂�ŕ�������
		MGPvector<MGSurface> srf_v;
		num += srf_u[i]->divide_multi_knot_v(srfl);
	}
	return num;
}

//u�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
//�߂�l�́A��������ԋp����
int MGSBRep::divide_multi_knot_u(
	MGPvector<MGSBRep>& srfl) const		//���������Ȗʃ��X�g
{
	size_t start_index = 0, index = 0, count = 0, multi = 0, bdim = 0, order = 0;
	MGInterval intv = param_range().ref(1);
	const MGKnotVector &knot_vector = knot_vector_u();
	start_index = order_u() - 1;
	bdim = bdim_u();
	order = order_u();
	do{		//u,v�����ɐ܂�(C0�A��)������Ƃ��ʂ𕪊�����
		if(order == 2){		//�I�[�_�[2�̂Ƃ��̏���
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		MGBox uv_range = MGBox(MGInterval(knot_u(start_index), knot_u(index)), intv);
		srfl.push_back(new MGSBRep(uv_range, *this));
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//���d�x��������Ȃ�������I���
	return count;
}

//v�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
//�߂�l�́A��������ԋp����
int MGSBRep::divide_multi_knot_v(
	MGPvector<MGSurface>& srfl) const		//���������Ȗʃ��X�g
{
	size_t start_index = 0, index = 0, count = 0, multi = 0, bdim = 0, order = 0;
	MGInterval intv = param_range().ref(0);
	const MGKnotVector &knot_vector = knot_vector_v();
	start_index = order_v() - 1;
	bdim = bdim_v();
	order = order_v();
	do{		//u,v�����ɐ܂�(C0�A��)������Ƃ��ʂ𕪊�����
		if(order == 2){		//�I�[�_�[2�̂Ƃ��̏���
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		MGBox uv_range = MGBox(intv, MGInterval(knot_v(start_index), knot_v(index)));
		srfl.push_back(new MGSBRep(uv_range, *this));
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//���d�x��������Ȃ�������I���
	return count;
}

//u�܂���v�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
//�߂�l�́A��������ԋp����
int MGRSBRep::divide_multi_knot(
	MGPvector<MGSurface>& srfl) const		//���������Ȗʃ��X�g
{
	//���u�����̐܂�̕����ŕ�������
	MGPvector<MGRSBRep> srf_u;
	int udiv = divide_multi_knot_u(srf_u), num = 0;
	for(int i = 0; i < udiv; i++){
		//����v�����̐܂�ŕ�������
		MGPvector<MGSurface> srf_v;
		num += srf_u[i]->divide_multi_knot_v(srfl);
	}
	return num;
}

//u�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
//�߂�l�́A��������ԋp����
int MGRSBRep::divide_multi_knot_u(
	MGPvector<MGRSBRep>& srfl) const		//���������Ȗʃ��X�g
{
	MGInterval intv = param_range().ref(1);
	const MGKnotVector &knot_vector = knot_vector_u();
	size_t start_index = order_u() - 1, index = 0, count = 0, multi = 0, bdim = bdim_u(), order = order_u();
/*	size_t first_index = start_index;	//�擪�̃C���f�b�N�X�A�r���Ő܂�Ă����炻����擪�ɂ���
	//U�����̍\���������o��
	MGPvector<MGCurve> vecCrvKousei;
	size_t tmpIndex = order_v() - 1, tmpBdim = bdim_v(), tmpOrder = order_v();
	for(; tmpIndex <= tmpBdim; tmpIndex++){
		double dV = knot_v(tmpIndex);
		vecCrvKousei.push_back(parameter_curve(true, dV));
	}*/
	do{		//u�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
		if(order == 2){		//�I�[�_�[2�̂Ƃ��̏���
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		//�}���`�m�b�g�Ő܂ꂪ�Ȃ�������p�X
/*		if(index != bdim){
			bool bG1 = true;	//G1�A�����ǂ��� G1�s�A���̎�false
			double dU = knot_u(index);
			//�S�\������dU�ɂ�����G1�s�A���łȂ����ǂ������ׂ�
			MGPvector<MGCurve>::const_iterator iter = vecCrvKousei.begin();
			for(; iter != vecCrvKousei.end(); iter++){
				MGVector vecRightContinuous = (*iter)->eval(dU,1,0), vecLeftContinuous = (*iter)->eval(dU,1,1);
				double dAngle = vecRightContinuous.angle(vecLeftContinuous);
				if(dAngle > MGTolerance::angle_zero()){bG1 = false; break;}	//����Ă����ꍇ�̏���
			}
			if(bG1){start_index = index + multi - 1; continue;}	//G1�A���Ȃ̂ŕ������Ȃ�
		}
		MGBox uv_range = MGBox(MGInterval(knot_u(first_index), knot_u(index)), intv);*/
		MGBox uv_range = MGBox(MGInterval(knot_u(start_index), knot_u(index)), intv);
		srfl.push_back(new MGRSBRep(uv_range, *this));
		//first_index = start_index = index + multi - 1;
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//���d�x��������Ȃ�������I���
	return count;
}

//v�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
//�߂�l�́A��������ԋp����
int MGRSBRep::divide_multi_knot_v(
	MGPvector<MGSurface>& srfl) const		//���������Ȗʃ��X�g
{
	MGInterval intv = param_range().ref(0);
	const MGKnotVector &knot_vector = knot_vector_v();
	size_t start_index = order_v() - 1, index = 0, count = 0, multi = 0, bdim = bdim_v(), order = order_v();
/*	size_t first_index = start_index;	//�擪�̃C���f�b�N�X�A�r���Ő܂�Ă����炻����擪�ɂ���
	//V�����̍\���������o��
	MGPvector<MGCurve> vecCrvKousei;
	size_t tmpIndex = order_u() - 1, tmpBdim = bdim_u(), tmpOrder = order_u();
	for(; tmpIndex < tmpBdim; tmpIndex++){
		double dU = knot_u(tmpIndex);
		vecCrvKousei.push_back(parameter_curve(false, dU));
	}*/
	do{		//u�����ɐ܂�(�}���`�m�b�g)������Ƃ��ʂ𕪊�����
		if(order == 2){		//�I�[�_�[2�̂Ƃ��̏���
			index = start_index + 1; multi = 1;
		}else{
			multi = knot_vector.locate_multi(start_index, order - 1, index);
		}
		//�}���`�m�b�g�Ő܂ꂪ�Ȃ�������p�X
/*		if(index != bdim){
			bool bG1 = true;	//G1�A�����ǂ��� G1�s�A���̎�false
			double dV = knot_v(index);
			//�S�\������dV�ɂ�����G1�s�A���łȂ����ǂ������ׂ�
			MGPvector<MGCurve>::const_iterator iter = vecCrvKousei.begin();
			for(; iter != vecCrvKousei.end(); iter++){
				MGVector vecRightContinuous = (*iter)->eval(dV,1,0), vecLeftContinuous = (*iter)->eval(dV,1,1);
				double dAngle = vecRightContinuous.angle(vecLeftContinuous);
				if(dAngle > MGTolerance::angle_zero()){bG1 = false; break;}	//����Ă����ꍇ�̏���
			}
			if(bG1){start_index = index + multi - 1; continue;}	//G1�A���Ȃ̂ŕ������Ȃ�
		}
		MGBox uv_range = MGBox(intv, MGInterval(knot_v(first_index), knot_v(index)));*/
		MGBox uv_range = MGBox(intv, MGInterval(knot_v(start_index), knot_v(index)));
		srfl.push_back(new MGRSBRep(uv_range, *this));
		//first_index = start_index = index + multi - 1;
		start_index = index + multi - 1;
		count++;
	}while(index != bdim);	//���d�x��������Ȃ�������I���
	return count;
}

//u�܂���v�����ɐ܂�(C0�A��)������Ƃ��ʂ𕪊�����
//�߂�l�́A��������ԋp����
int MGSurface::divide_multi_knot(
	MGPvector<MGSurface>& srfl) const		//���������Ȗʃ��X�g
{
	srfl.push_back(copy_surface());
	return 1;
}

