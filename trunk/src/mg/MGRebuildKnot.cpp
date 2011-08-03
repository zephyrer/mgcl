/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/LBRep.h"
#include "mg/Tolerance.h"
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implementation of knot rebuild.

//�Ȑ���̃m�b�g�x�N�g�����č\�z����
//���͂��ꂽ�����Ȑ�����w��I�[�_�[�ōč\�z����B�g�������X��line_zero()���g�p���Ă���B
//�I�[�_�[���w�肳��Ă��Ȃ��Ƃ��Ȑ���̂����ōł��傫���I�[�_�[���g�p����B���̂Ƃ��A
//Ellipse, Straight�̃I�[�_�[��4�Ƃ��čl����B
//�p�����[�^�͈͂�1�������l�̑傫�����P�ɂȂ�悤�ɂ����Ƃ��̒����̕��ς��g�p���Ă���B
//�߂�l�͍č\�z��̋Ȑ��񂪕ԋp�����B�G���[�̂Ƃ��k�����ԋp�����B
MGPvector<MGLBRep> rebuild_knot(
	const std::vector<const MGCurve*>& crvl,	//���͋Ȑ���
	size_t order ,			//�w��I�[�_�[
	MGLBRep**		tp	//�ڑ���	input and output.
		//if tp[i] for crvl[i] was not null, converted new tp will be output.
){
	//�g�p����I�[�_�[�Ƌ�Ԏ��������߂�
	int nCrv = crvl.size(), maxOrder = 0, maxSdim = 0, i = 0;
	int idMaxOrder=0;
	MGPvector<MGLBRep> rtnBrepList(nCrv);
	double allParam = 0.0;
	for(i = 0; i < nCrv; i++){
		int ord = crvl[i]->order(), sdim = crvl[i]->sdim();
		if(sdim < 1){
			rtnBrepList.clear();
			return rtnBrepList;	//��Ԏ���0�̓G���[�Ƃ���
		}
		if(ord > maxOrder){
			maxOrder = ord; idMaxOrder=i;
		}
		if(sdim > maxSdim)
			maxSdim = sdim;
	}
	if(!order)
		order = maxOrder;		//�I�[�_�[���w�肳��Ă��Ȃ��Ƃ��̏���
	if(order <= 2)
		order = 4;	//Ellipse, Straight�̂Ƃ��I�[�_�[4�ɂ���

	//LBRep�łȂ����̂́ALBRep�ɕϊ����ă��X�g�ɓ����
	//�S�ăm�b�g�x�N�g����������LBRep�̂Ƃ����r���h���Ȃ��ł��̂܂܂̋Ȑ���ԋp����
	bool fsame = true;
	for(i = 0; i < nCrv; i++){
		const MGLBRep *ptempBRep = dynamic_cast<const MGLBRep*>(crvl[i]);
		if(ptempBRep){
			rtnBrepList.reset(i, new MGLBRep(*ptempBRep));
			if(i && fsame){
				// KnotVector�����������ǂ����̃`�F�b�N
				const MGKnotVector& kv0 = rtnBrepList[i-1]->knot_vector();
				const MGKnotVector& kv1 = rtnBrepList[i]->knot_vector();
				if (kv0.param_s() != kv1.param_s() ||
					kv0.param_e() != kv1.param_e() ||
					kv0 != kv1) fsame = false;  
//				if(rtnBrepList[i]->knot_vector() != rtnBrepList[i - 1]->knot_vector())fsame = false;
			}
		}else{
			rtnBrepList.reset(i, new MGLBRep(*crvl[i], order));
			fsame = false;
		}
		//1�������l�̑傫�����P�ɂȂ�悤�ȃp�����[�^�͈͂����߂�
		double t0=rtnBrepList[i]->param_s(), t1=rtnBrepList[i]->param_e();
		double t2=(t0+t1)*.5;
		double mag=rtnBrepList[i]->eval(t0,1).len();
		mag+=rtnBrepList[i]->eval(t1,1).len();
		mag+=rtnBrepList[i]->eval(t2,1).len();
		mag /= 3.;
		allParam += rtnBrepList[i]->param_span() * mag;
	}
	if(fsame && tp==0) return rtnBrepList;
	double avgSpan = allParam / nCrv;

	//�p�����[�^�͈͂����킹���m�b�g�x�N�g�����쐬����B�S�Ẵm�b�g�x�N�g���𑫂����킹��B
	for(i = 0; i< nCrv; i++){
		rtnBrepList[i]->change_range(0.0, avgSpan);
		if(tp){if(tp[i]) tp[i]->change_range(0.0, avgSpan);}
	}

	MGNDDArray tau_temp;
	tau_temp.update_from_knot(rtnBrepList[idMaxOrder]->knot_vector());
	if(tau_temp.length()<order)
		tau_temp.change_number(order);
	MGKnotVector mixKnotVector(tau_temp,order);
	for(i = 1; i < nCrv; i++){
		MGKnotVector& ti=rtnBrepList[i]->knot_vector();
		size_t orderti=ti.order();
		size_t orderi= order<=orderti ? order: orderti;
		mixKnotVector.mix_knot_vector(MGKnotVector(tau_temp.update_from_knot(ti),orderi));
		//if(tp){
		//	if(tp[i]) mixKnotVector.mix_knot_vector(MGKnotVector(tp[i]->knot_vector(), order));
		//}
	}
	//cout<<mixKnotVector<<endl;///////////
	//�m�b�g�x�N�g�����\���ɑ��₷�Boffset_div_num()���g�p����B
	MGKnotVector tempKnotVector = mixKnotVector;
	int bd=tempKnotVector.bdim();
	for(i=order-1; i<bd; i++){
		int maxNumDiv = 0, j = 0;
		double	spara = tempKnotVector(i),
				epara = tempKnotVector(i + 1);
		if(MGRZero((epara - spara) / avgSpan))continue;	//�}���`�m�b�g�̂Ƃ��̏���
		for(; j < nCrv; j++){	//�e�X�p���̍ő啪���������߂�
			int ndiv = rtnBrepList[j]->offset_div_num(MGInterval(spara, epara));
			if(ndiv > maxNumDiv)maxNumDiv = ndiv;
		}
		double span = (epara - spara) / maxNumDiv, tParam = tempKnotVector(i);
		for(j = 0; j < maxNumDiv - 1; j++){
			tParam += span;
			mixKnotVector.add_data(tParam, order - 1);	//�m�b�g�x�N�g���𑝂₷
		}
	}
	//cout<<mixKnotVector<<endl;///////////

	//���₵���m�b�g�ŋȐ����č쐬����
	MGNDDArray tau;
	tau.update_from_knot(mixKnotVector);
	int n = tau.length();
	double tplen=0.;
	for(i = 0; i < nCrv; i++){
		double tpleni=0.;
		MGBPointSeq bp(n, maxSdim);
		for(int j = 0; j < n; j++){
			bp.store_at(j, rtnBrepList[i]->eval(tau(j)));
		}
		int error;
		rtnBrepList.reset(i, new MGLBRep(tau, bp, mixKnotVector, error));
		if(error){rtnBrepList.clear(); return rtnBrepList;}
		if(tp){
			if(tp[i]){
				MGBPointSeq bptp(n, maxSdim);
				for(int j = 0; j < n; j++){
					MGVector tpatj=tp[i]->eval(tau(j));
					bptp.store_at(j, tpatj);
					tpleni+=tpatj.len();
				}
				tp[i]=new MGLBRep(tau, bptp, mixKnotVector, error);
				tpleni/=double(n);
			}
		}
		tplen+=tpleni;
	}
	tplen/=double(nCrv);

	//�Ȑ���̃m�b�g�x�N�g���𐸓x��ێ����č폜���s��
	::remove_knot_curves(rtnBrepList,tp,tplen);
	//cout<<rtnBrepList[0]->knot_vector()<<endl;//////////////
	return rtnBrepList;
}

MGPvector<MGLBRep> rebuild_knot(
	const MGPvector<MGCurve>& brepl,	//���͋Ȑ���
	size_t order,		//�w��I�[�_�[
	MGLBRep**	tp	//�ڑ���
){
	size_t n=brepl.size();
	std::vector<const MGCurve*> curves(n);
	for(size_t i=0; i<n; i++) curves[i]=brepl[i];
	return rebuild_knot(curves,order,tp);
}

//Approximate the curve by MGLBRep.
//The parameter of the original will not be changed.
std::auto_ptr<MGLBRep> MGCurve::rebuild(
	size_t order,		//order of the new MGLBRep, must be >=4.
	double tol		//tolerance allowed for the approximation
)const{
	const MGKnotVector& t=knot_vector();
	size_t k=t.order(), n=t.bdim(), orderm1=order-1;
	size_t nnew=n+order-k;
	MGKnotVector t2(order,nnew);
	size_t l1,l2;
	double params=t[k-1], perame=t[n];
	for(l2=0; l2<order;l2++){
		t2(l2)=params;
		t2(nnew+l2)=perame;
	}
	for(l2=order,l1=k; l2<nnew; l2++,l1++)
		t2(l2)=t[l1];
	//cout<<t2<<endl;///////////

	//�m�b�g�x�N�g�����\���ɑ��₷�Boffset_div_num()���g�p����B
	double pspan=t.param_span();
	for(size_t i=k-1; i<n; i++){
		double spara=t(i), epara=t(i+1);
		double diff=epara-spara;
		if(MGRZero2(diff,pspan))
			continue;	//�}���`�m�b�g�̂Ƃ��̏���
		int ndiv = offset_div_num(MGInterval(spara, epara));
		double span = diff/ndiv;
		for(int j=0; j<ndiv-1; j++){
			spara += span;
			t2.add_data(spara, orderm1);	//�m�b�g�x�N�g���𑝂₷
		}
	}
	//cout<<t2<<endl;///////////

	//���₵���m�b�g�ŋȐ����č쐬����
	MGNDDArray tau; tau.update_from_knot(t2);
	MGBPointSeq bp(tau.length(), sdim());
	eval_line(tau,bp);
	int error;
	std::auto_ptr<MGLBRep> lbrep(new MGLBRep(tau, bp, t2, error));
	if(error){
		return lbrep;
	}

	//�Ȑ���̃m�b�g�x�N�g���𐸓x��ێ����č폜���s��
	double tol_old=MGTolerance::set_line_zero(tol);
	lbrep->remove_knot();
	//cout<<lbrep->knot_vector()<<endl;//////////////
	MGTolerance::set_line_zero(tol_old);
	return lbrep;
}
