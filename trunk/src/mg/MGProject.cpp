/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Bisection.h"
#include "mg/Box.h"
#include "mg/Knot.h"
#include "mg/Unit_vector.h"
#include "mg/CParam_list.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"
#include "mg/SurfCurve.h"
#include "mg/FSurface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Pvector.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;
#define EXTEND_COEF 2.0		//�x�N�g�����e�̃X�C�[�v�ʂ̒��������߂�̂Ɏg�p

// Implementation of MGFSurface projection.

//����������2�{��B�\���Ȑ���ڑ�����(������ނ̂Ƃ�)
MGLBRep* join2LBRep(const MGLBRep& crv1, const MGLBRep& crv2);

//�����������T��(xyzuv)��B�\���Ȑ����X�g��ڑ����Axyz������line_zero�Ń����[�u�m�b�g����
//join_crvl_uv,join_crv_xyz�ɐڑ������Ȑ����X�g������B
//Function's return value is the number of output curves.
void prjJoin(
	MGPvector<MGLBRep>& crvl,
	MGPvector<MGCurve>& join_crvl_uv,
	MGPvector<MGCurve>& join_crvl_xyz
){
	int num = crvl.size();
	if(!num)
		return;

	if(num==1){	//�Ȑ����P�{�̂Ƃ��̏���
		const MGLBRep& lbrep=*(crvl[0]);
		join_crvl_uv.push_back(new MGLBRep(2,lbrep,0,3));
		join_crvl_xyz.push_back(new MGLBRep(3,lbrep,0,0));
		return;
	}

	MGLBRep *cur_pcrv;
	MGLBRep *next_pcrv;
	cur_pcrv = crvl.release(0);
	for(int i=1; i<num; i++){
		next_pcrv = crvl.release(i);
		MGLBRep *pre_pcrv = join2LBRep(*cur_pcrv,*next_pcrv);
		if(pre_pcrv){//If successfully joined.
			delete cur_pcrv;
			delete next_pcrv;
			cur_pcrv = pre_pcrv;
		}else{
			join_crvl_uv.push_back(new MGLBRep(2,(*cur_pcrv),0,3));
			join_crvl_xyz.push_back(new MGLBRep(3,(*cur_pcrv),0,0));
			delete cur_pcrv;
			cur_pcrv = next_pcrv;
		}
	}
	join_crvl_uv.push_back(new MGLBRep(2,*cur_pcrv,0,3));
	join_crvl_xyz.push_back(new MGLBRep(3,*cur_pcrv,0,0));
	delete cur_pcrv;
}

//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is two kind of curves:
//one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
//the parameter space of the surfaces(vec_crv_uv).
//vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
int MGFSurface::project(
	const MGCurve& crv,
	MGPvector<MGCurve>& vec_crv_uv,	
		//Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	const MGCompositeCurve* ccrv=dynamic_cast<const MGCompositeCurve*>(&crv);
	if(ccrv){
		ccrv->project_onto_surface(*this,vec_crv_uv,vec_crv,vec);
		return vec_crv.size();
	}
	//������
	assert(crv.sdim() > 1);
	std::auto_ptr<MGCurve> crvRemoveKnot(crv.clone());
	crvRemoveKnot->remove_knot();

	//�x�N�g�����e���ǂ������ׂ�
	if(!vec.is_null())	//�x�N�g�����e
		return projVector(*crvRemoveKnot, vec_crv_uv, vec_crv, vec);

	const MGPlane* pln=dynamic_cast<const MGPlane*>(this);
	if(pln)
		//Plane�Ńx�N�g�����w�肳��Ă��Ȃ��Ƃ���plane�̃m�[�}���Ńx�N�g�����e����
		return projVector(*crvRemoveKnot, vec_crv_uv, vec_crv, pln->normal());

	MGPvector<MGCurve> divided_curves;
	//�}���`�m�b�g�ɂȂ�Ƃ��͐؂��ē��e����
	int n=crvRemoveKnot->divide_multi(divided_curves);
	MGPvector<MGLBRep> crv_xyzuv_vector;//Projected curve of(x,y,z,u,v).
	for(int i=0; i<n; i++){
		project1span(*(divided_curves[i]),crv_xyzuv_vector);
	}

	//���e�Ȑ���̏璷�m�b�g���폜���Đڑ�����
	prjJoin(crv_xyzuv_vector, vec_crv_uv, vec_crv);
	return vec_crv.size();
}

//Obtain the projected curve of a curve onto the surface.
//project1span does not divide the curve at the c0 continuity, which is treated by
//project(). See projec().
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project1span' is curves of space dimension 5, which are (x,y,z,u,v).
//Here (x,y,z) is the world coordinates, and (u,v) is parameter of this surface.
void MGFSurface::project1span(
	const MGCurve& crv,	//The target curve to project.
	MGPvector<MGLBRep>& crv_xyzuv_vector//Projected curve of(x,y,z,u,v) will be appended.						
)const{
	//�p�����[�^�͈͂����߂�
	std::deque<MGPosition> ranges;//ranges[2*i] to ranges[2*i+1] is the projection range.
			//Let ranges[.]=tuv, then tuv[0] is the parameter of the curve, and (tuv[1], tuv[2])
			//is the parameter of this surface (u,v).

	double sparam = crv.param_s();
	int counter=0, rc;
	do{
		MGPosition tuv[2];
		rc = prj2GetParamRange(crv,counter,tuv,counter);
		if(rc==1 || rc==-1){
			ranges.push_back(tuv[0]);
			ranges.push_back(tuv[1]);
		}
	}while(rc < 0);	//crv1�̓r���ł���ΌJ��Ԃ�

	//���e���s��
	while(ranges.size()>=2){
		//ranges will be pop_front by prj2OneCurve.
		MGLBRep* projectedCurve;
		prj2OneCurve(crv,ranges,projectedCurve);

		//���e�Ȑ���̏璷�m�b�g���폜���Đڑ�����
		if(projectedCurve)
			crv_xyzuv_vector.push_back(projectedCurve);
	}
}

//Obtain the projected curve of a curve onto the surface.
//The direction of the projection is along the vector vec if the vec is not
//NULL, and normal to the surface if the vec is NULL.
//Output of 'project' is general world coordinate curves('vec_crv')
int MGFSurface::project(
	const MGCurve& crv,			//given curve.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec			//projection vector.
		//if vec = NULL then calculate perpendicular project.
)const{
	MGPvector<MGCurve> vec_crv_uv;
	return project(crv,vec_crv_uv,vec_crv,vec);
}

//�ʒ��ɓ��e�����_��ԋp����
//�߂�l�́A��_�܂��͖ʒ��_�����܂����Ƃ���1�A���܂�Ȃ������Ƃ���0��ԋp����
int MGFSurface::project_normal(
	const MGPosition& pos,
	const MGPosition& uv_guess,	//���ʃp�����[�^
	MGPosition& uv
)const{
	double dPreTol = MGTolerance::set_wc_zero(MGTolerance::wc_zero()*0.5);
	int rc = perp_point(pos, uv, &uv_guess);
	MGTolerance::set_wc_zero(dPreTol);
	return rc;
}

//�x�N�g�����e�́A�J�[�u��܂�ŕ������čs���A��Őڑ�����
int MGFSurface::projVector(
	const MGCurve& crv,
	MGPvector<MGCurve>& vec_crv_uv,
		//Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		//Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const{
	MGPvector<MGCurve> mult_crvl;
	MGPvector<MGLBRep> temp_vec_crv;
	//�}���`�m�b�g�ɂȂ�Ƃ��͐؂��ē��e���čēx�ڑ�����
	int num_mult = crv.divide_multi(mult_crvl);
	for(int i=0; i<num_mult; i++)
		projVectorProc(*(mult_crvl[i]), temp_vec_crv, vec);
	prjJoin(temp_vec_crv, vec_crv_uv, vec_crv);
	return vec_crv.size();
}

//�x�N�g�����e�́Acrv��vec�����ɃX�C�[�v�����ʂƌ��̖ʂƂ̌���ŋ��߂�
int MGFSurface::projVectorProc(
	const MGCurve& crv,
	MGPvector<MGLBRep>& crv_xyzuv_vector,//Projected curve of(x,y,z,u,v) will be appended.						
	const MGVector& vec
)const{
	MGBox box(this->get_box() | crv.box());
	double extLeng = box.len() * EXTEND_COEF;	//���S�̈�
	const double srfError=param_error();	
	//�����炭plane�̎��̏�����ʂɏ����Ă��Ȃ��Ƃ��܂��s���Ȃ��Ǝv��
	//���x�I�ɖ��Ȃ���΂��̏������R�����g�ɂ��Ă����v����2001/11/26
	if(get_surface_pointer()->type() == MGSURFACE_PLANE){
		//�Ȑ��̒��S����x�N�g�������̒����𐶐�����
		MGStraight centerStraight(MGSTRAIGHT_UNLIMIT,vec,crv.center());
		//�����ƖʂƂ̌�_�����߂�
		MGCSisect_list isectList = isect(centerStraight);
		if(isectList.size() > 0){
			std::list<MGCSisect>::const_iterator iterIsect = isectList.begin();
			MGPosition posIsect = iterIsect->point();
			//�������钷���𒼐�����+�Ȑ������ɂ���
			extLeng = (posIsect - crv.center()).len() + crv.length();
			extLeng *= EXTEND_COEF;	//���S�̈�
		}
	}
	if(MGAZero(extLeng / EXTEND_COEF))extLeng = 1.0;		//��΃[���ȉ��̏ꍇ�̏���
	std::auto_ptr<MGSurface> apCutSrf(crv.sweep(vec, extLeng, -extLeng));
	if(!apCutSrf.get())
		return 0;
	double errsave=MGTolerance::set_wc_zero(MGTolerance::line_zero());
	MGSSisect_list isectl = isect(*apCutSrf);//cout<<isectl<<endl;/////////////********
	MGTolerance::set_wc_zero(errsave);
	int i;
	MGSSisect_list::SSiterator ssiter;
	for(ssiter = isectl.begin(), i = 0; ssiter != isectl.end(); ssiter++, i++){
		//�J�[�u�̌��������킹�Ă���MGPvector��push����
		std::auto_ptr<MGCurve> apCrv2d(ssiter->release_param1());
		std::auto_ptr<MGCurve> apCrv3d(ssiter->release_line());
		MGCurve &crvOth2d = ssiter->param2();
		assert(apCrv2d.get() && apCrv3d.get());
		//u�p�����[�^�������������n�_�Ƃ���(���J�[�u�̌����Ɠ����ɂ���)
		double spt = (crvOth2d.start_point())(0), ept = (crvOth2d.end_point())(0);
		if(spt > ept){
			apCrv2d->negate();
			apCrv3d->negate();
		}
		//cout<<*apCrv2d<<*apCrv3d<<endl;
		double	param_len = apCrv2d->box().len(),	//�V���[�g�A�[�N�̃`�F�b�N������
				world_len = apCrv3d->box().len();
		if(param_len < srfError || world_len < MGTolerance::wc_zero())continue;
		const MGLBRep *pCrv2d = dynamic_cast<const MGLBRep*>(apCrv2d.get());
		const MGLBRep *pCrv3d = dynamic_cast<const MGLBRep*>(apCrv3d.get());
		assert(pCrv2d && pCrv3d);
		size_t bdim = pCrv2d->bdim();
		MGBPointSeq bpXYZUV(5,pCrv3d->line_bcoef());
		const MGBPointSeq &bpUV = pCrv2d->line_bcoef();
		for(size_t i2=0; i2<bdim; i2++)
			bpXYZUV.store_at(i2,bpUV(i2),3,0);
		//cout<<bpXYZUV<<endl;
		crv_xyzuv_vector.push_back(new MGLBRep(pCrv2d->knot_vector(), bpXYZUV));
	}
	return i;
}

//�g�������X�ɉ������ו������f�[�^�|�C���g�����߂�
void MGFSurface::prjGetDataPoints(
	const MGCurve& curve,
	const MGPosition& tuv0,//parameter range of curve to get the data point,
	const MGPosition& tuv1,//From tuv0 to tuv1.
	MGNDDArray& tau		//data points will be output.
)const{
	double t0=tuv0[0], t1=tuv1[0];
	MGInterval cspan(t0, t1);
	int n1=curve.offset_div_num(cspan);

	const MGSurface* srf=get_surface_pointer();

	MGCurve* pcrv1=srf->parameter_curve(1,(tuv0[1]+tuv1[1])*.5);
	MGInterval cspan2(tuv0[1], tuv1[1]);
	int n2=pcrv1->offset_div_num(cspan2);
	delete pcrv1;

	MGCurve* pcrv2=srf->parameter_curve(0,(tuv0[2]+tuv1[2])*.5);
	MGInterval cspan3(tuv0[2], tuv1[2]);
	int n3=pcrv2->offset_div_num(cspan3);
	delete pcrv2;

	int n=n1;
	if(n<n2)
		n=n2;
	if(n<n3)
		n=n3;

	tau.resize(n+1);
	double delta=(t1-t0)/n;
	double t=t0;
	for(int i=0; i<n; i++, t+=delta){
		tau(i)=t;
	}
	tau(n)=t1;
}

//Get world position data at (u, v) and store the coordinates
//and the parameter in bp.
void MGFSurface::prj_store_bp(
	int i,
	double u, double v,
	MGBPointSeq& bp
)const{
	MGPosition Psrf=eval(u, v);	//�n�_�����߂�
	bp.store_at(i,Psrf,0,0,3);
	bp(i,3)=u; bp(i,4)=v;
}

//���e�����s���ĂT�����Ȑ����擾����
//curve�ɂ�����f�[�^�|�C���g�Ɠ��e�x�N�g�����瓊�e�_����쐬��
//�_�񂩂�T����(xyz,uv)�̋Ȑ�����쐬���ĕԋp����
void MGFSurface::prj2OneCurve(
	const MGCurve& curve,	//target curve to prject.
	std::deque<MGPosition>& ranges,//start(ranges[0]) and end(ranges[1]) point parameter
							//of the curve and the face.
							//On return ranges will be so updated that processed rages[0] to [1]
							//are pop_front().
	MGLBRep*& crvProjected	//newed object will be returend when obtained.
							//When not obtained, null will be returned.
)const{
	crvProjected=0;
	MGPosition tuv0=ranges.front(); ranges.pop_front();
	MGPosition tuv1=ranges.front(); ranges.pop_front();

	MGNDDArray tau;//data point to get the projected points.
				//tau[0] is tuv[0][0], and tau[n-1]=tuv[1][0].
	prjGetDataPoints(curve,tuv0,tuv1,tau);

	int ntau = tau.length();
	int ntaum1=ntau-1;

	MGBPointSeq bp_add(ntau, 5);	//xyz,uv�����킹���a�\��[x,y,z,u,v]
	prj_store_bp(0,tuv0[1],tuv0[2],bp_add);//�n�_�����߂�

	const MGSurface* srf=get_surface_pointer();
	MGPosition uv, uv_guess(tuv0[1],tuv0[2]);
	int i=1;
	for(; i<ntaum1 ; i++){
		MGPosition Pcrv = curve.eval_position(tau(i));
		project_normal(Pcrv, uv_guess, uv);
		//	MGPosition tuv2=proj2GetParamIter(curve,tuv1,tau[i],tuv1[0]);
		//	ranges.push_front(tuv1);//to process in the next call of prjOneCurve.
		//	ranges.push_front(tuv2);//Next process's start .
		//	MGPosition tuv(tau[i-1],uv_guess[0],uv_guess[1]);
		//	tuv1=proj2GetParamIter(curve,tuv,tau[i-1],tau[i]);//New end.
		//	break;
		
		size_t periNum;
		srf->on_a_perimeter(uv(0), uv(1), periNum);
			//�ӂɂ̂��Ă����炻�̃p�����[�^���g�p����(uv��update����Ă���j
		uv_guess = uv;		//���ʓ_���X�V����
		prj_store_bp(i,uv[0],uv[1],bp_add);//�n�_�����߂�
	}
	prj_store_bp(i++,tuv1[1],tuv1[2],bp_add);//End point�����߂�
	bp_add.set_length(i);

	//�Ȑ����쐬���ăx�N�g���ɑ}������
	MGBox bxAddXYZ(3,bp_add.box());	//xyz�̃{�b�N�X�����߂�
	//�V���[�g�A�[�N�̃`�F�b�N������
	if(bxAddXYZ.len() <= MGTolerance::wc_zero())
		return;

	MGBPointSeq tmpBpXYZ(3,bp_add);
	MGNDDArray tauXYZ(tmpBpXYZ);
	crvProjected = new MGLBRep(tauXYZ, bp_add, 4, 4.);
	crvProjected->remove_knot(0,3);	//�ŏ��̂R����(XYZ)��line_zero�Ńm�b�g�폜
		//cout<<(*crvProjected);
}

//Normalize the parameter of the input range[.][0], i.e. if the value is alomost equalt to
//a knot value, round the value into the knot value. This is usefull especially for
//start or end parameter value.
int normalize_check_isect(
	const MGCurve& curve,
	int retval,
	MGPosition range[2]
){
	//�p�����[�^�l�𐳋K������
	MGPosition& tuv0=range[0];
	MGPosition& tuv1=range[1];
	double ts=tuv0(0)=curve.param_normalize(tuv0[0]);
	double te=tuv1(0)=curve.param_normalize(tuv1[0]);
	if((te-ts)<=curve.param_error()){
		if(retval==1)
			return 0;	//���e�͈͂��덷�͈͂̏ꍇ��_�ƂȂ�̂�0��Ԃ�
		else
			return -2;
	}
	return retval;
}

//�X�^�[�g�p�����[�^��^�����e�\�ȃp�����[�^�͈͂�1�����擾����B
//�߂�l��	1:�͈͂����܂���(this�̍Ō�܂�)
//			0:�͈͂����܂�Ȃ�����(this�̍Ō�܂�)
//			-1:�͈͂����܂���(this�̓r���܂ŁA�܂��ʒ��͈͂����邩������Ȃ�)
//			-2;�͈͂����܂�Ȃ��������Athis�̓r���ŁA�܂��ʒ��͈͂����邩������Ȃ�
int MGFSurface::prj2GetParamRange(
	const MGCurve& curve,
	int start_counter,	//input start counter of the curve parameter incrementation.
	MGPosition range[2],
		//range[0][0]=start point parameter of curve,
		//range[0][1-2]=the corresponding surface(this surface's) parameter (u,v)
		//of the start point of the range.
		//Regarding to range[1] : the same.
	int& next_counter	//Updated curve parameter incremental counter will be output.
)const{
	//������
	//std::cout<<curve<<std::endl;
	int ndiv = curve.intersect_dnum()+2;

	//curve�𕪊����A���e�\�����ׂ�
	const double delta=curve.param_span()/ndiv;	//�`�F�b�N�|�C���g�̃p�����[�^�X�p��
	const double tend=curve.param_e();

	int t_is_on_surface;
	double t=curve.param_s()+delta*double(start_counter), t_pre;
	MGPosition uv;
	int i=start_counter;//counter for t's incremental number.
	if(t_is_on_surface=perp_one(curve.eval(t),uv)){
		range[0]=MGPosition(t,uv[0],uv[1]);
	}else{
		//Find the 1st t that is on surface.
		while(!t_is_on_surface && i<ndiv){
			t_pre=t;
			t+=delta; i++;
			t_is_on_surface=perp_one(curve.eval(t),uv);
		}
		if(t_is_on_surface){
			range[0]=proj2GetParamIter(curve,MGPosition(t,uv[0],uv[1]),t_pre,t);
		}else{
			if(perp_one(curve.eval(tend),uv)){
				range[0]=proj2GetParamIter(curve,MGPosition(tend,uv[0],uv[1]),tend-delta,tend);
				range[1]=MGPosition(tend,uv[0],uv[1]);
				return normalize_check_isect(curve,1,range);
			}
			return 0;
		}
	}

	//Here t_is_on_surface=true always holds.
	MGPosition uv_save;
	while(t_is_on_surface && t<tend && i<ndiv){
		uv_save=uv;
		t_pre=t;
		t+=delta; i++;
		t_is_on_surface=perp_point(curve.eval(t),uv,&uv_save);
	}

	int retval;
	if(!t_is_on_surface){
		range[1]=proj2GetParamIter(curve,MGPosition(t_pre,uv_save[0],uv_save[1]),t_pre,t);
		retval=-1;
	}else{
		uv_save=uv;
		if(perp_point(curve.eval(tend),uv,&uv_save)){
			range[1]=MGPosition(tend,uv[0],uv[1]);
		}else{
			range[1]=proj2GetParamIter(curve,MGPosition(t,uv[0],uv[1]),t,tend);
		}
		retval=1;
	}

	//�I������
	next_counter=i;		//����̌����J�n�p�����[�^��Ԃ�
	return normalize_check_isect(curve,retval,range);
}

class MGPrjBisect{
public:

MGPrjBisect(
	double ts,	//parameter range from ts to te.
	double te
);

//Virtual Destructor
virtual ~MGPrjBisect(){;};

//compare with the previous function value(the initial value is set
//by set_initial_t) and replace t with the previous one if necessary.
//The function's return value is the new parameter value.
virtual void compare_replace(
	double t	//parameter value to compare at.
)=0;

int solve(
	double tolerance//The tolerance to halt the bisection iteration.
);

protected:
	double m_ts, m_te;
		//the curve's parameter range from m_ts to m_te, or from m_te to m_ts.
		//note that m_ts<m_te does not always hold.
};

//A virtual super class to solve non-linear equations by bicection methos.
MGPrjBisect::MGPrjBisect(
	double ts,	//parameter range from ts to te.
	double te
):m_ts(ts), m_te(te){}

//Compute the fn(t)'s parameter value that is the maxima 
//Function's return value will be the solution obtained.
int MGPrjBisect::solve(
	double tolerance//The tolerance to halt the bisection iteration.
){
	assert(tolerance>0.);
	int nrepition=0;
	double span=m_te-m_ts;
	tolerance+=MGTolerance::mach_zero()*2.;
	while(fabs(span)>=tolerance){//m_ts > m_te may happend.
		span*=.5; nrepition++;
	    double t=m_ts+span;
		compare_replace(t);
	}
	return nrepition;
}

class MGProjBoundaryParam: public MGPrjBisect{
public:
	MGProjBoundaryParam(
		const MGFSurface& surf,
		const MGCurve& curve,
		const MGPosition& tuv,
		double ts, double te
	):MGPrjBisect(ts,te),m_surf(surf),m_curve(curve),m_uv(tuv[1],tuv[2]){
		if(tuv[0]==ts){
			m_te=ts;
			m_ts=te;
		}
	//m_te holds the curve parameter value whose point has
	//a perpendicular point ont the surface.
	}

	//compare with the previous function value(the initial value is set
	//by set_initial_t) and replace t with the previous one if necessary.
	//The function's return value is the new parameter value.
	virtual void compare_replace(
		double t	//parameter value to compare at.
	){
		MGPosition uv;
		MGPosition P=m_curve.eval(t);
		if(m_surf.project_normal(P,m_uv, uv)){
			m_uv=uv;
			m_te=t;
		}else{
			m_ts=t;
		}
	}

	MGPosition get_result(){return MGPosition(m_te,m_uv[0],m_uv[1]);};
	const MGFSurface& m_surf;
	const MGCurve& m_curve;
	MGPosition m_uv;
};

//curve��̓��e�\�ȃp�����[�^�Ɠ��e�s�\�ȃp�����[�^��^���āA
//�Ԃɂ��铊�e�\�ȋ��E�p�����[�^�l�����߂ĕԋp����
//�`�F�b�N�|�C���g�̈ړ��� < (�p�����[�^�͈�*rc_zero())�ɂȂ�ΏI���B
//Function's return value is MGPosition tuv, where
//tuv[0]=start point parameter of curve,
//(tuv[1], tuv[2])=the corresponding surface(this surface's) parameter (u,v)
//of the start point of the range.
MGPosition MGFSurface::proj2GetParamIter(
	const MGCurve& curve,//target curve to project.
	const MGPosition& tuv,//tuv[0] is the curve's parameter value that has a perp point onto the
				//surface. and (tuv[1], tuv[2]) is the surface's parameter of the perp point.
				//tuv[0] is eithere ts or te.
	double ts,	//(ts, te) is curve's parameter range that indicates
	double te	//between ts and te there must be the boundary point.
				//That is one of the following situations occurs:
				//(1) there are no perp points at ts and there is a perp point at te,
				//(2) there is a perp point at ts and there are no perp points at te,
)const{
	//������
	const double crvError = curve.param_error()*2.;//, srfError = get_surface_pointer()->param_error();
	MGProjBoundaryParam solver(*this,curve,tuv,ts,te);
	solver.solve(crvError);
	return solver.get_result();
}
