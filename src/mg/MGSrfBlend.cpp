/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Pvector.h"
#include "mg/Tolerance.h"
#include "mg/Box.h"
#include "mg/KnotArray.h"
#include "mg/SPointSeq.h"
#include "mg/Straight.h"
#include "mg/CCisect_list.h"
#include "mg/LBRepEndC.h"
#include "mg/SBRepEndC.h"
#include "mg/SBRepTP.h"
#include "mg/Coons.h"
#include "mg/SBRep.h"

extern "C" {
#include "cskernel/Bvstan.h"
}
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
//using namespace std;

//Compute 2nd derivative.
MGVector get_2nd(const MGVector& first,const MGVector& second_guess,double curvature){
	if(MGRZero(second_guess.len())) return second_guess;
	MGUnit_vector e2(second_guess);
	double len=first.len()*curvature;
	return e2*len;
}

//ノットをあわせた後の境界線、ブレンド関数、データポイントから点列を求める。
bool bilinear_spoint_proc(
	const MGPvector<MGLBRep>& brepl,//ノットベクトルをあわせた後の境界線
	const MGCurve& blendCrvU,	//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve& blendCrvV,	//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGNDDArray& utau,		//u方向のデータポイント
	const MGNDDArray& vtau,		//v方向のデータポイント
	MGSPointSeq& spoint)		//点列
{
	spoint = MGSPointSeq(utau.length(), vtau.length(), brepl[0]->sdim());

	//4隅の点は境界線の端点の中間とする
	MGPosition	p00 = (brepl[0]->start_point() + brepl[3]->start_point()) / 2.0,
				p01 = (brepl[2]->start_point() + brepl[3]->end_point()) / 2.0,
				p10 = (brepl[0]->end_point() + brepl[1]->start_point()) / 2.0,
				p11 = (brepl[1]->end_point() + brepl[2]->end_point()) / 2.0;

	//境界線上の点を求める
	size_t i, nu=utau.length(), nv=vtau.length();
	size_t num1=nu-1, nvm1=nv-1;
	for(i=0; i<nu; i++){
		spoint.store_at(i, 0, brepl[0]->eval(utau(i)));
		spoint.store_at(i, nvm1, brepl[2]->eval(utau(i)));
	}
	for(i=0; i<nv; i++){
		spoint.store_at(0, i, brepl[3]->eval(vtau(i)));
		spoint.store_at(num1, i, brepl[1]->eval(vtau(i)));
	}

	//内部点を求める
	double spanu = brepl[0]->param_span(), spanv = brepl[1]->param_span();
	for(i=1; i<num1; i++){
		double blendu = (blendCrvU.eval(utau(i) / spanu))(0);
		double onembu=1.-blendu;
		MGPosition	pu0(spoint(i, 0)),
					pu1(spoint(i, nvm1));
		for(size_t j=1; j<nvm1; j++){
			double blendv = (blendCrvV.eval(vtau(j) / spanv))(0);
			double onembv=1.-blendv;
			MGPosition	p0v(spoint(0, j)),
						p1v(spoint(num1, j));
			spoint.store_at(i,j,onembu*p0v+blendu*p1v+onembv*pu0+blendv*pu1
				- (onembv*(onembu*p00+blendu*p10)+blendv*(onembu*p01+blendu*p11)));
		}
	}
	return true;
}

//境界線、ブレンド関数、接続面を与え、対辺が同じノットベクトルのスプライン曲線
//になるように再作成して、点列、ノットベクトル、データポイントを求め、境界線の
//パラメータにあわせてあった接続面をリビルドした後のパラメータに変更する。
//境界線はC1連続であり、vmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminから
//umaxの方向にumin,umaxの向きをvminからvmaxの方向になっているものとする。
//境界線のノットベクトルをあわせるときの誤差はline_zero()を使用している。
//		戻り値：
//		0:			正常終了
//		-1:			境界線がC1連続でなかった
//		-2:			接続面のパラメータ範囲が境界線と違った
//		-3:			境界線のノットベクトルをあわせられなかった
//		-4:			接続面のパラメータ範囲を変更できなかった
//		-5:			点列が求まらなかった
int bilinear_spoint(
	const MGCurve&			blendCrvU,		//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve&			blendCrvV,		//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGPvector<MGLBRep>& perimeters,
	const MGSBRepTP& tp,	//接続面(パラメータ範囲は境界線と同じ)
	MGSPointSeq&			spoint,			//点列
	MGNDDArray&				utau,			//u方向のデータポイント
	MGNDDArray&				vtau			//v方向のデータポイント
){
	//端末ベクトル、データポイントを求める。
	MGENDCOND condk[4]={MGENDC_NO,MGENDC_NO,MGENDC_NO,MGENDC_NO};
	for(size_t i=0; i<4; i++) if(tp.specified(i)) condk[i] = MGENDC_1D;

	//点列を求める
	MGKnotVector& knotu = perimeters[0]->knot_vector();
	MGKnotVector& knotv = perimeters[1]->knot_vector();
	utau = MGNDDArray(condk[3], condk[1], knotu);//cout<<"utau="<<utau;
	vtau = MGNDDArray(condk[0], condk[2], knotv);//cout<<"vtau="<<vtau<<endl;
	if(!bilinear_spoint_proc(perimeters, blendCrvU, blendCrvV, utau, vtau, spoint))return -5;
	return 0;
}

//境界線、接続面、点列データポイント、接ベクトルから接続条件(MGSBRepEndC)を求める
bool bilinear_endc(
	MGSPointSeq&		spoint,		//点列
	const MGSBRepTP&		tp,			//接続面(パラメータ範囲は境界線と同じ)
	const MGPvector<MGLBRep>& perimeters,
	const MGCurve&			blendCrvU,		//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve&			blendCrvV,		//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGNDDArray&		utau,		//u方向のデータポイント
	const MGNDDArray&		vtau,		//v方向のデータポイント
	MGSBRepEndC&			endc)	//端末条件
{
	if(utau.length() != spoint.length_u() || vtau.length() != spoint.length_v() || spoint.sdim() != 3)return false;
	const MGLBRep& peri0=*(perimeters[0]);
	const MGLBRep& peri1=*(perimeters[1]);
	const MGLBRep& peri2=*(perimeters[2]);
	const MGLBRep& peri3=*(perimeters[3]);

	MGENDCOND condk[4]={MGENDC_NO,MGENDC_NO,MGENDC_NO,MGENDC_NO};
	const size_t dim=3; int m;
	size_t i,p;
	size_t lenu=utau.length(), lenv=vtau.length();
	size_t lenum1=lenu-1, lenvm1=lenv-1;
	double ur[2]={utau(0),utau(lenum1)}, vr[2]={vtau(0),vtau(lenvm1)};
	const int ipse[2]={1,2};

	size_t ktp, ntp, n;
	const double* ttp; const double* rtp;
	int itanse[2]={1,1}; double tanse[6];
	size_t irtp, ip1, ip2;
	double work[105], tangen[3];

	size_t tps[4]={0,0,0,0};
	for(i=0; i<4; i++) if(tp.specified(i)) tps[i]=1;

	if(tps[0] || tps[2]){

	const MGKnotVector& knotv = peri1.knot_vector();
	MGNDDArray vtau2(condk[0], condk[2], knotv);//cout<<"vtau="<<vtau<<endl<<"vtau2="<<vtau2<<endl;
	MGSPointSeq spoint2v;
	bilinear_spoint_proc(perimeters, blendCrvU, blendCrvV, utau, vtau2, spoint2v);
	spoint2v.capacity(ip1, ip2);
	n=vtau2.length();
	double spanu = peri0.param_span();
	for(m=0; m<2; m++){
		//Process of perimeter num 0 and 2(v=min and max line)
		size_t perimeter=m*2;
		if(tps[perimeter]){
			const MGLBRep& tpm = tp.TP(perimeter);
			ktp=tpm.order(); ntp=tpm.bdim();
			ttp=tpm.knot_data(); rtp=tpm.coef_data();
			irtp=tpm.line_bcoef().capacity();
		//We have to generate first and 2nd derivative data for this perimeter.
			double v=vr[m];
			MGVector deri1[2]={peri3.eval(v,1),peri1.eval(v,1)};

			MGBPointSeq first(lenu,dim);
			for(p=0; p<dim; p++) tanse[p]=first(0,p)=deri1[0][p];
			for(p=0; p<dim; p++) tanse[p+3]=first(lenum1,p)=deri1[1][p];
			for(i=1; i<lenum1; i++){
				double tptau=utau(i);
				bvstan_(ur,ktp,ntp,ttp,rtp,tptau,n,vtau2.data(),
					spoint2v.data(i,0,0),ipse[m],itanse,tanse,irtp,
					ip1,ip2,work,tangen);
				MGVector deri1ati(dim,tangen);
				first.store_at(i,deri1ati);
			}
			endc.set_1st(perimeter,first);
		}
	}

	}

	size_t ntemp=lenv-tps[0]-tps[2];
	if((tps[0] || tps[2]) && ntemp>=2){

	MGVector deri2sv[2]={peri3.eval(vr[0],2),peri1.eval(vr[0],2)};
	MGVector deri2ev[2]={peri3.eval(vr[1],2),peri1.eval(vr[1],2)};
	double crvtr0s=peri3.curvature(vr[0]), crvtr1s=peri1.curvature(vr[0]);
	double crvtr0e=peri3.curvature(vr[1]), crvtr1e=peri1.curvature(vr[1]);
	double spanu = peri0.param_span();
	MGNDDArray taut(ntemp);
	taut(0)=vtau(0);
	for(size_t j=1 ; j<ntemp-1; j++) taut(j)=vtau(j+tps[0]);
	taut(ntemp-1)=vtau(lenv-1);//cout<<taut<<endl;
	for(i=1; i<lenum1; i++){
		double tptau=utau(i);
		MGLBRepEndC endc0, endc2;
		MGBPointSeq bpt(ntemp,3);
		bpt.store_at(0,spoint(i,0));
		double blendu = (blendCrvU.eval(tptau / spanu))(0);
		if(tps[0]){
			MGVector deri1=MGVector(3,endc.first(0)(i));
			endc0.set_1st(deri1);
			MGVector deri2ati0=deri2sv[0].interpolate_by_rotate(blendu,deri2sv[1]);
			deri2ati0=get_2nd(deri1,deri2ati0,crvtr0s+(crvtr1s-crvtr0s)*blendu);
			endc0.set_2nd(deri2ati0);
		}
		for(size_t j=1 ; j<ntemp-1; j++) bpt.store_at(j,spoint(i,j+tps[0]));
		if(tps[2]){
			MGVector deri1=MGVector(3,endc.first(2)(i));
			endc2.set_1st(deri1);
			MGVector deri2ati0=deri2ev[0].interpolate_by_rotate(blendu,deri2ev[1]);
			deri2ati0=get_2nd(deri1,deri2ati0,crvtr0e+(crvtr1e-crvtr0e)*blendu);
			endc2.set_2nd(deri2ati0);
		}
		bpt.store_at(ntemp-1,spoint(i,lenv-1));
		size_t k=6;
		size_t new_bdim=ntemp+tps[0]*2+tps[2]*2;
		if(new_bdim<k) k=new_bdim;
		int error;
		MGLBRep lbt(k,endc0,endc2,taut,bpt,error);
		if(tps[0]) spoint.store_at(i,1,lbt.eval(vtau(1)));
		if(tps[2]) spoint.store_at(i,lenv-2,lbt.eval(vtau(lenv-2)));
	}

	}

	if(tps[3] || tps[1]){

	const MGKnotVector& knotu = peri0.knot_vector();
	MGNDDArray utau2(condk[3], condk[1], knotu);//cout<<"utau="<<utau<<endl<<"utau2="<<utau2<<endl;
	MGSPointSeq spoint2u;
	bilinear_spoint_proc(perimeters, blendCrvU, blendCrvV, utau2, vtau, spoint2u);
	size_t psizeu, psizev;
	spoint2u.capacity(psizeu, psizev);
	n=utau2.length();
	ip1=1; ip2=psizeu*psizev;
	const size_t iv[2]={3,1};
	double spanv = peri1.param_span();
	for(m=0; m<2; m++){
		//Process of perimeter num 3 and 1(u=min and max line)
		size_t perimeter=iv[m];
		if(tps[perimeter]){
			const MGLBRep& tpm = tp.TP(perimeter);
			ktp=tpm.order(); ntp=tpm.bdim();
			ttp=tpm.knot_data(); rtp=tpm.coef_data();
			irtp=tpm.line_bcoef().capacity();
		//We have to generate first and 2nd derivative data for this perimeter.
			double u=ur[m];
			MGVector deri1[2]={peri0.eval(u,1),peri2.eval(u,1)};

			MGBPointSeq first(lenv,dim);
			for(p=0; p<dim; p++) tanse[p]=first(0,p)=deri1[0][p];
			for(p=0; p<dim; p++) tanse[p+3]=first(lenvm1,p)=deri1[1][p];
			for(i=1; i<lenvm1; i++){
				double tptau=vtau(i);
				bvstan_(vr,ktp,ntp,ttp,rtp,tptau,n,utau2.data(),
					spoint2u.data(0,i,0),ipse[m],itanse,tanse,irtp,
					ip1,ip2,work,tangen);
				MGVector deri1ati(dim,tangen);
				first.store_at(i,deri1ati);
			}
			endc.set_1st(perimeter,first);
		}
	}

	}

	ntemp=lenu-tps[3]-tps[1];
	if((tps[3] || tps[1]) && ntemp>=2){

	MGVector deri2su[2]={peri0.eval(ur[0],2),peri2.eval(ur[0],2)};
	MGVector deri2eu[2]={peri0.eval(ur[1],2),peri2.eval(ur[1],2)};
	double crvtr0s=peri0.curvature(ur[0]), crvtr1s=peri2.curvature(ur[0]);
	double crvtr0e=peri0.curvature(ur[1]), crvtr1e=peri2.curvature(ur[1]);

	double spanv = peri1.param_span();
	MGNDDArray taut(ntemp);
	taut(0)=utau(0);
	for(size_t j=1 ; j<ntemp-1; j++) taut(j)=utau(j+tps[3]);
	taut(ntemp-1)=utau(lenu-1);//cout<<taut<<endl;
	for(i=1; i<lenvm1; i++){
		double tptau=vtau(i);
		MGLBRepEndC endc0, endc2;
		MGBPointSeq bpt(ntemp,3);
		bpt.store_at(0,spoint(0,i));
		double blendv = (blendCrvV.eval(tptau / spanv))(0);
		if(tps[3]){
			MGVector deri2ati0=deri2su[0].interpolate_by_rotate(blendv,deri2su[1]);
			endc0.set_1st(MGVector(3,endc.first(3)(i)));
			endc0.set_2nd(deri2ati0);
		}
		for(size_t j=1 ; j<ntemp-1; j++) bpt.store_at(j,spoint(j+tps[3],i));
		if(tps[1]){
			MGVector deri2ati0=deri2eu[0].interpolate_by_rotate(blendv,deri2eu[1]);
			endc2.set_1st(MGVector(3,endc.first(1)(i)));
			endc2.set_2nd(deri2ati0);
		}
		bpt.store_at(ntemp-1,spoint(lenu-1,i));
		size_t k=6;
		size_t new_bdim=ntemp+tps[3]*2+tps[1]*2;
		if(new_bdim<k) k=new_bdim;
		int error;
		MGLBRep lbt(k,endc0,endc2,taut,bpt,error);
		if(tps[3]) spoint.store_at(1,i,lbt.eval(utau(1)));
		if(tps[1]) spoint.store_at(lenu-2,i,lbt.eval(utau(lenu-2)));
	}

	}
	return true;
}

double get_deri_coef(double t0, double t1,double alpha, double t){
	if(t<=.5) return (2.*t*(1.-alpha)+2.*alpha*t0-(t0+t1))/(t0-t1);
	return (2.*t*(1.-alpha)+2.*alpha*t1-(t0+t1))/(t1-t0);
}

void get_derivatives(
	bool along_u,			//true if for perimeter 0 and 2.
							//false if for perimeter 3 and 1.
	const MGSBRep& surf,	//Original surface.
	const MGSBRepTP& tp,	//接続面(パラメータ範囲は境界線と同じ)
	MGPvector<MGLBRep>& derivatives,
				//array of derivatives[4], must have size 4.
	const double* alpha		//derivative magnitude coefficients for
		//perimeter i in alpha[i].
		//=1. is normal, if<1., curvature will be made large, and
		//if>1. curvature will be made small.
){
	size_t sd=surf.sdim();
	int error;
	const MGKnotVector& tu=surf.knot_vector_u();
	const MGKnotVector& tv=surf.knot_vector_v();

	if(along_u){

	size_t nu=surf.bdim_u();
	MGNDDArray utau; utau.update_from_knot(tu);
	MGBPointSeq deris_u(nu,sd);
	double v[2]={tv.param_s(),tv.param_e()};
	double u0=tu.param_s(), u1=tu.param_e();

	for(size_t m=0; m<2; m++){//for perimeter 0 and 2.
		size_t perim=m*2;
		double vm=v[m];
		for(size_t i=0; i<nu; i++){
			double utaui=utau[i];
			MGVector deri=surf.eval(utaui,vm,0,1);//cout<<deri;
			if(tp.specified(perim)){
				MGVector N=tp.TP(perim).eval(utaui);//Normal of tangent plane.
				double vlen=deri.len();
				deri=MGUnit_vector(N*deri*N)*vlen;//cout<<deri<<endl;
			}
			if(alpha) deri*=get_deri_coef(u0,u1,alpha[perim],utaui);
			deris_u.store_at(i,deri);
		}
		derivatives.reset(perim,new MGLBRep(utau,deris_u,tu,error));
		//cout<<"driv at perim="<<perim<<","<<*(derivatives[perim])<<endl;
	}

	}else{

	size_t nv=surf.bdim_v();
	MGNDDArray vtau; vtau.update_from_knot(tv);
	MGBPointSeq deris_v(nv,sd);
	double u[2]={tu.param_s(),tu.param_e()};
	double v0=tv.param_s(), v1=tv.param_e();
	size_t peri[2]={3,1};

	for(size_t m=0; m<2; m++){//for perimeter 3 and 1.
		size_t perim=peri[m];
		double um=u[m];
		for(size_t j=0; j<nv; j++){
			double vtauj=vtau[j];
			MGVector deri=surf.eval(um,vtauj,1,0);//cout<<deri;
			if(tp.specified(perim)){
				MGVector N=tp.TP(perim).eval(vtauj);
					//Normal of tangent plane at surf(um,vtau[j]).
				double vlen=deri.len();
				deri=MGUnit_vector(N*deri*N)*vlen;//cout<<deri<<endl;
			}
			if(alpha) deri*=get_deri_coef(v0,v1,alpha[perim], vtauj);
			deris_v.store_at(j,deri);
		}
		derivatives.reset(perim,new MGLBRep(vtau,deris_v,tv,error));
		//cout<<"driv at perim="<<perim<<","<<*(derivatives[perim])<<endl;
	}
	
	}
}

//境界線、ブレンド関数、接続面を与え、対辺が同じノットベクトルのスプライン曲線
//になるように再作成して、点列、ノットベクトル、データポイントを求め、境界線の
//パラメータにあわせてあった接続面をリビルドした後のパラメータに変更する。
//境界線はC1連続であり、vmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminから
//umaxの方向にumin,umaxの向きをvminからvmaxの方向になっているものとする。
//境界線のノットベクトルをあわせるときの誤差はline_zero()を使用している。
//		戻り値：
//		0:			正常終了
//		-1:			境界線がC1連続でなかった
//		-2:			接続面のパラメータ範囲が境界線と違った
//		-3:			境界線のノットベクトルをあわせられなかった
//		-4:			接続面のパラメータ範囲を変更できなかった
int rebuild_curve(
	const MGCurve*			edge_crvl[4],	//境界線リスト(vmin,umax,vmax,uminの順)
	MGSBRepTP&				tp,				//接続面(パラメータ範囲は境界線と同じ)
	MGPvector<MGLBRep>& perimeters
){
	//境界線がC1連続か、接続面のパラメータ範囲が境界線と同じかどうか調べる
	for(int i=0; i<4; i++){
		//if(!edge_crvl[i]->cn_continuity(1)) return -1;	//C1連続性のチェック
		if(!tp.specified(i)) continue;
		if(edge_crvl[i]->param_range() != tp.TP(i).param_range())return -2;
	}

	//Make the tolerances smaller.
	double lzero=MGTolerance::line_zero();//save line zero;
	MGTolerance::set_line_zero(lzero*.3);
	double azero=MGTolerance::angle_zero();//save line zero;
	MGTolerance::set_angle_zero(azero*.3);

	//リビルドを行う
	int k=4;	//オーダーは４とする 
	std::vector<const MGCurve*> temp_crvl(2);
	perimeters.resize(4);
	temp_crvl[0] = edge_crvl[0];	temp_crvl[1] = edge_crvl[2];
	MGLBRep* tplb[2]={0,0};
	if(tp.specified(0)) tplb[0]=&(tp.TP(0));
	if(tp.specified(2)) tplb[1]=&(tp.TP(2));
	MGPvector<MGLBRep> temp_brepl = rebuild_knot(temp_crvl, k, tplb);
	if(!temp_brepl.size()){MGTolerance::set_line_zero(lzero);return -3;}
	tp.set_TP(0,std::auto_ptr<MGLBRep>(tplb[0]));
	tp.set_TP(2,std::auto_ptr<MGLBRep>(tplb[1]));
	perimeters.reset(0, temp_brepl.release(0));
	perimeters.reset(2, temp_brepl.release(1));

	temp_crvl[0] = edge_crvl[1];	temp_crvl[1] = edge_crvl[3];
	tplb[0]=tplb[1]=0;
	if(tp.specified(1)) tplb[0]=&(tp.TP(1));
	if(tp.specified(3)) tplb[1]=&(tp.TP(3));
	MGPvector<MGLBRep> temp_brepl2 = rebuild_knot(temp_crvl, k, tplb);
	if(!temp_brepl2.size()){MGTolerance::set_line_zero(lzero);return -3;}
	tp.set_TP(1,std::auto_ptr<MGLBRep>(tplb[0]));
	tp.set_TP(3,std::auto_ptr<MGLBRep>(tplb[1]));
	perimeters.reset(1, temp_brepl2.release(0));
	perimeters.reset(3, temp_brepl2.release(1));

	//Restore the original tolerance.
	MGTolerance::set_line_zero(lzero);
	MGTolerance::set_angle_zero(azero);

	MGLBRep& l0=*(perimeters[0]);
	MGLBRep& l1=*(perimeters[1]);
	MGLBRep& l2=*(perimeters[2]);
	MGLBRep& l3=*(perimeters[3]);

	//make the knots fine.
	int error;
	MGKnotVector tu(l0.knot_vector()), tv(l1.knot_vector());
	tu.change_knot_number(tu.bdim()*3);
	tv.change_knot_number(tv.bdim()*3);
	MGCurve* crv2[4];
	for(int i=0; i<4; i++)
		crv2[i]=edge_crvl[i]->clone();
	double u0=tu.param_s(), u1=tu.param_e();
	double v0=tv.param_s(), v1=tv.param_e();
	crv2[0]->change_range(u0,u1);
	crv2[1]->change_range(v0,v1);
	crv2[2]->change_range(u0,u1);
	crv2[3]->change_range(v0,v1);
	l0=MGLBRep(*(crv2[0]),tu);
	l1=MGLBRep(*(crv2[1]),tv);
	l2=MGLBRep(*(crv2[2]),tu);
	l3=MGLBRep(*(crv2[3]),tv);
	for(int i=0; i<4; i++){
		delete crv2[i];
		if(!tp.specified(i))continue;
		MGKnotVector* t=&tu; if(i%2) t=&tv;
		MGLBRep& tpi=tp.TP(i);
		tpi=MGLBRep(tpi,*t,error);
	}

	//Adjust the corner points.
	MGPosition P0=(l0.start_point()+l3.start_point())*.5;
	MGPosition P1=(l0.end_point()+l1.start_point())*.5;
	MGPosition P2=(l1.end_point()+l2.end_point())*.5;
	MGPosition P3=(l3.end_point()+l2.start_point())*.5;

	double fixp[2];
	//P0
	fixp[0]=l0.param_e(); l0.move(2,l0.param_s(),P0,fixp);
	fixp[0]=l3.param_e(); l3.move(2,l3.param_s(),P0,fixp);
	//P1
	fixp[0]=l0.param_s(); l0.move(2,l0.param_e(),P1,fixp);
	fixp[0]=l1.param_e(); l1.move(2,l1.param_s(),P1,fixp);
	//P2
	fixp[0]=l1.param_s(); l1.move(2,l1.param_e(),P2,fixp);
	fixp[0]=l2.param_s(); l2.move(2,l2.param_e(),P2,fixp);
	//P3
	fixp[0]=l2.param_e(); l2.move(2,l2.param_s(),P3,fixp);
	fixp[0]=l3.param_s(); l3.move(2,l3.param_e(),P3,fixp);

	return 0;
}

int rebuild_curve(
	MGPvector<MGLBRep>& edges,
	MGSBRepTP&				tp,				//接続面(パラメータ範囲は境界線と同じ)
	MGPvector<MGLBRep>& perimeters
){
	const MGCurve*	edge_crvl[4];
	for(size_t i=0; i<4; i++)
		edge_crvl[i]=edges[i];
	return rebuild_curve(edge_crvl,tp,perimeters);
}

//4本の境界線、ブレンド関数、接続面を与えて面を生成する。
//境界線はvmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminからumaxの方向に
//umin,umaxの向きをvminからvmaxの方向になっているものとする。境界線のノットベクトル
//をあわせるときの誤差はline_zero()を使用している。ブレンド曲線はパラメータ範囲0,1
//で値域も0,1である。接続面(MGSBRepTP)のパラメータ範囲は各境界線と同じとする。
//		エラーコード：
//		0:			正常終了
//		-1:			境界線がC1連続でなかった
//		-2:			接続面のパラメータ範囲が境界線と違った
//		-3:			境界線のノットベクトルをあわせられなかった
//		-4:			接続面のパラメータ範囲を変更できなかった
//		-5:			点列が求まらなかった
//		-6:			端末条件が求まらなかった
//		-7:			面が生成できなかった
MGSBRep::MGSBRep(
	const MGCurve*					edge_crvl[4],	//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGCurve&					blendCrvU,		//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve&					blendCrvV,		//空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGSBRepTP&				tp,				//接続面(パラメータ範囲は境界線と同じ)
	int&							error):MGSurface()	//エラーコード
{
	//点列、ノットベクトル、データポイントを求める
	//for(size_t iii=0;iii<4; iii++) cout<<(*edge_crvl[iii]);

	MGSBRepTP tempTP(tp);
	MGSPointSeq spoint;
	MGNDDArray utau, vtau;
	MGPvector<MGLBRep> perimeters;

	//1. Adjust parameter range of the tp and edge_crvl.
	error=rebuild_curve(edge_crvl,tempTP,perimeters);
	if(error)return;

	//2. Build spoint data.
	error=bilinear_spoint(blendCrvU,blendCrvV,perimeters,tempTP,spoint,utau,vtau);
	if(error<0)return;

	//3. 接続条件(MGSBRepEndC)を求める
	MGSBRepEndC endc;
	if(!bilinear_endc(spoint, tempTP,perimeters, blendCrvU, blendCrvV,
		utau, vtau, endc)){error = -6; return;}

	//4. 点列、データポイント、ノットベクトル、接続条件より面を生成する
	MGKnotVector& knotu=perimeters[0]->knot_vector();
	MGKnotVector& knotv=perimeters[1]->knot_vector();
	//cout<<utau<<knotu<<vtau<<knotv<<endl;//////////////////////////////////
	*this = MGSBRep(endc, utau, vtau, spoint, knotu, knotv, error);
	if(error)error = -7;
}

//Easy to use version of the above constructor.
//When blendCCrvU,V were null, straight line from 0. to 1. will be used.
MGSBRep::MGSBRep(
	MGPvector<MGLBRep>& edges,	//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	int&			error	,	//エラーコードが出力される
	const MGCurve*	blendCrvU,//空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve*	blendCrvV //空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
){
	MGPosition zero(1); zero(0)=0.;
	MGPosition one(1); one(0)=1.;
	MGStraight zero_one(one,zero);

	const MGCurve* blendu;
	const MGCurve* blendv;
	if(blendCrvU) blendu=blendCrvU; else blendu=&zero_one;
	if(blendCrvV) blendv=blendCrvV; else blendv=&zero_one;

	MGSBRepTP tp;
	MGSPointSeq spoint;
	MGNDDArray utau, vtau;
	MGPvector<MGLBRep> perimeters;

	//1. Adjust parameter range of the edge_crvl.
	error=rebuild_curve(edges,tp,perimeters);
	if(error)return;

	//2. Build spoint data.
	error=bilinear_spoint(*blendu,*blendv,perimeters,tp,spoint,utau,vtau);
	if(error<0) return;

	//3. 点列、データポイント、ノットベクトル、接続条件より面を生成する
	MGKnotVector& knotu=perimeters[0]->knot_vector();
	MGKnotVector& knotv=perimeters[1]->knot_vector();
	//cout<<utau<<knotu<<vtau<<knotv<<endl;//////////////////////////////////
	MGSBRepEndC endc;
	*this = MGSBRep(endc, utau, vtau, spoint, knotu, knotv, error);
	if(error) error = -7;
}

MGSBRep::MGSBRep(
	const MGCurve*	edge_crvl[4],	//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGSBRepTP& tp,			//接続面(パラメータ範囲は境界線と同じ)
	int& error,	//エラーコード
	const double* alpha	//derivative magnitude coefficients for perimeter i in alpha[i].
		//=1. is normal, if<1., curvature will be made large, and
		//if>1. curvature will be made small.
		//If alpha=null, all of the coefs are assumed to be 1.
):MGSurface(){
//	cout<<*(edge_crvl[0])<<endl;
	//cout<<MGTolerance::instance();
	MGSBRepTP tp2(tp);
	MGPvector<MGLBRep> perimeters;
	if(rebuild_curve(edge_crvl,tp2,perimeters)) return;
	//cout<<tp2<<endl;

	//Save the original parameter range.
	MGKnotVector& tu=perimeters[0]->knot_vector();
	double u0=tu.param_s(), u1=tu.param_e();
	MGKnotVector& tv=perimeters[1]->knot_vector();
	double v0=tv.param_s(), v1=tv.param_e();

	//Change parameter range to (0., 1.).
	size_t i;
	for(i=0; i<4; i++){
		perimeters[i]->change_range(0., 1.);
		if(tp2.specified(i)) tp2.TP(i).change_range(0., 1.);
	}
	//cout<<tp2<<endl;

	MGPvector<MGLBRep> derivatives(4);
	bool along_u=true;

//Method 1.
	MGSBRep ruled0(along_u,perimeters);//cout<<ruled0;
	MGSBRep ruled1(!along_u,perimeters);//cout<<ruled1;
//	get_derivatives(along_u,ruled0,tp2,derivatives,alpha);
//	get_derivatives(!along_u,ruled1,tp2,derivatives,alpha);

//Method 2.
	MGSPointSeq sp(ruled0.surface_bcoef()+ruled1.surface_bcoef());
	sp*=.5;
	MGSBRep ruled01(sp,ruled0.knot_vector_u(), ruled0.knot_vector_v());//cout<<ruled01;
	get_derivatives(along_u,ruled01,tp2,derivatives,alpha);
	get_derivatives(!along_u,ruled01,tp2,derivatives,alpha);

//Method 3.
/*	MGSBRep bilinear(perimeters,error);
	bilinear.change_range(1,0.,1.);	bilinear.change_range(0,0.,1.);
	get_derivatives(along_u,bilinear,tp2,derivatives,alpha);
	get_derivatives(!along_u,bilinear,tp2,derivatives,alpha);
*/
	double taumax[4];
/*	double cosmax[4];
	tp2.get_perimeters_max_cos(derivatives,taumax,cosmax);
	for(i=0; i<4; i++)
		cout<<endl<<"i="<<i<<":taumax="<<taumax[i]<<", cosmax="<<cosmax[i];cout<<endl;
*/
	//Save the old knot vector.
	MGKnotVector tu2(tu), tv2(tv);
	//Construct a coons' patch.
	MGCoons coons(perimeters,derivatives);

	MGNDDArray utau(3,tu2.bdim()-2,tu2);
	tv2.change_knot_number(tv2.bdim()*3);
	MGNDDArray vtau(3,tv2.bdim()-2,tv2);
	MGSPointSeq spoint(utau.length(),vtau.length(),3);
	coons.eval(utau,vtau,spoint);
	MGSBRepEndC endc(utau,vtau,coons);
	MGSBRep* surf=new MGSBRep(endc,utau,vtau,spoint,tu2,tv2,error);
	if(error){ error = -7; return ;}

	double sinmax[4];
	double error_max=tp2.get_perimeters_max_sin(*surf,taumax,sinmax);
	error_max*=1.005;
	double azero=MGTolerance::angle_zero();
	if(error_max<azero) error_max=azero;

	double line0=MGTolerance::line_zero();
	double line02=line0*2.;
	MGSBRep* surf2=0;
	//Remove knot by the minimum line zero that satisfy the angle continuity.
	for(size_t q=0; q<4; q++){
		delete surf2;
		surf2=new MGSBRep(*surf);
		line02*=.5;
		MGTolerance::set_line_zero(line02);
		surf2->remove_knot();
		double max2=tp2.get_perimeters_max_sin(*surf2,taumax,sinmax);
		if(max2<=error_max) break;
	}
	delete surf;
	surf=surf2;

	MGTolerance::set_line_zero(line0);
	surf->change_range(1,u0,u1);
	surf->change_range(0,v0,v1);
	*this=*surf;
	delete surf;
}
