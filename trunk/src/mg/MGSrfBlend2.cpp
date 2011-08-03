/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Pvector.h"
#include "mg/SPointSeq.h"
#include "mg/LBRepEndC.h"
#include "mg/SBRepEndC.h"
#include "mg/SBRepTP.h"
#include "mg/SBRep.h"
#include "mg/BSumSurf.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/blg4sp2.h"
}

using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

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
);

void get_all_derivatives(
 	const MGPvector<MGLBRep>& perimeters,//perimeters.
	const MGSBRep& surf,	//Original surface.
	const MGSBRepTP& tp,	//接続面(パラメータ範囲は境界線と同じ)
	MGPvector<MGLBRep>& derivatives
				//array of derivatives[4], must have size 4.
){
	size_t sd=surf.sdim();
	int error;
	double azero=MGTolerance::angle_zero()*.05;
	const MGKnotVector& tu=surf.knot_vector_u();
	const MGKnotVector& tv=surf.knot_vector_v();

	size_t m;
///////////perimeter 0 and 2.
	size_t nu=surf.bdim_u();size_t num1=nu-1;
	MGNDDArray utau; utau.update_from_knot(tu);
	MGBPointSeq deris_u(nu,sd);
	MGBPointSeq dpu(nu,1);
	double v[2]={tv.param_s(),tv.param_e()};

	for(m=0; m<2; m++){//for perimeter 0 and 2.
		size_t perim=m*2;
		if(tp.specified(perim)){
			double vm=v[m];
			deris_u.store_at(0,perimeters[3]->eval(vm,1));
			for(size_t i=1; i<num1; i++){
				double utaui=utau[i];
				MGVector deri=surf.eval(utaui,vm,0,1);//cout<<deri;
				MGVector N=tp.TP(perim).eval(utaui);//Normal of tangent plane.
				double vlen=deri.len();
				deri=MGUnit_vector(N*deri*N)*vlen;//cout<<deri<<endl;
				deris_u.store_at(i,deri);
				dpu(i,0)=vlen;
			}
			deris_u.store_at(num1,perimeters[1]->eval(vm,1));
			dpu(0,0)*=.2;dpu(num1,0)*=.2;//Make precise at the ends.
			dpu(1,0)*=.4;dpu(nu-2,0)*=.4;//Make precise at the ends.
			derivatives.reset(perim,new MGLBRep(utau,deris_u,tu,error));//***NON Smoothing version///
			//cout<<"driv at perim="<<perim<<","<<*(derivatives[perim])<<endl;
		}
	}

////////// perimeter 1 and 3.
	size_t nv=surf.bdim_v();size_t nvm1=nv-1;
	MGNDDArray vtau; vtau.update_from_knot(tv);
	MGBPointSeq deris_v(nv,sd);
	double u[2]={tu.param_s(),tu.param_e()};
	size_t peri[2]={3,1};

	for(m=0; m<2; m++){//for perimeter 3 and 1.
		size_t perim=peri[m];
		if(tp.specified(perim)){
			double um=u[m];
			deris_v.store_at(0,perimeters[0]->eval(um,1));
			for(size_t j=1; j<nvm1; j++){
				double vtauj=vtau[j];
				MGVector deri=surf.eval(um,vtauj,1,0);//cout<<deri;
				MGVector N=tp.TP(perim).eval(vtauj);
					//Normal of tangent plane at surf(um,vtau[j]).
				double vlen=deri.len();
				deri=MGUnit_vector(N*deri*N)*vlen;//cout<<deri<<endl;
				deris_v.store_at(j,deri);
			}
			deris_v.store_at(nvm1,perimeters[2]->eval(um,1));
			derivatives.reset(perim,new MGLBRep(vtau,deris_v,tv,error));
			//cout<<"driv at perim="<<perim<<","<<*(derivatives[perim])<<endl;
		}
	}

	size_t i,l;
//*****Twist adjustment. Here dudv and dvdu are only averaged, may be able to improve.
//Adjust the twist vector at the corners.
	MGVector twist[4];
	for(i=0; i<4; i++){
		size_t im1=(i+3)%4;
		if(tp.specified(im1) && tp.specified(i)){
			double tim1=derivatives[im1]->param_e(),ti=derivatives[i]->param_s();
			if(im1>=2) tim1=derivatives[im1]->param_s();
			if(i>=2) ti=derivatives[i]->param_e();
			MGVector deriim1=derivatives[im1]->eval(tim1,1), derii=derivatives[i]->eval(ti,1);
			twist[i]=(deriim1+derii)*.5;
			//cout<<endl<<"****in get_all_derivatives,corner twist at "<<i<<"::"<<endl;
			//cout<<"i-1="<<deriim1<<",i="<<derii<<",twist="<<twist[i]<<endl;
		}
	}
	MGLBRepEndC ECS, ECE;
	MGNDDArray utau2(3,tu.bdim()-2,tu);
	MGNDDArray vtau2(3,tv.bdim()-2,tv);
	size_t ntau2,nutau2=utau2.length(), nvtau2=vtau2.length();
	MGBPointSeq dyu(nutau2,1), dyv(nvtau2,1);
	for(l=0; l<nutau2; l++) dyu(l,0)=1.;
	dyu(0,0)=.01;dyu(nutau2-1,0)=.01;
	for(l=0; l<nvtau2; l++) dyv(l,0)=1.;
	dyv(0,0)=.01;dyv(nvtau2-1,0)=.01;
	MGNDDArray* tau2;
	const MGKnotVector* t;
	double* dy;
	for(i=0; i<4; i++){
		if(!tp.specified(i))continue;
		double ts,te;
		if(i%2) {
			tau2=&vtau2; t=&tv; dy=dyv.data(); ntau2=nvtau2;
		}else{
			tau2=&utau2; t=&tu; dy=dyu.data(); ntau2=nutau2;
		}
		ts=(*tau2)[0]; te=(*tau2)[tau2->length()-1];
		size_t ids=i, ide=i+1; ide%=4;
		if(i>=2){ids=ide; ide=i;}
		if(twist[ids].sdim()==0) ECS.set_1st(derivatives[i]->eval(ts,1));
		else ECS.set_1st(twist[ids]);
		if(twist[ide].sdim()==0) ECE.set_1st(derivatives[i]->eval(te,1));
		else ECE.set_1st(twist[ide]);
		MGBPointSeq bp;
		derivatives[i]->eval_line(*tau2,bp);
		double dev=bp(0).len()+bp(ntau2/2).len()+bp(ntau2-1).len();
		dev/=3.;
		dev*=MGTolerance::angle_zero()*.8;
		MGLBRep* lb=new MGLBRep;
		lb->buildSRSmoothedLB_of_1stDeriv(ECS,ECE,*tau2,bp,dy,dev,false);
		derivatives.reset(i,lb);
		//cout<<"driv at perim="<<i<<","<<*(derivatives[i])<<endl;
	}
}

MGLBRep* get_perriderisub(
	const MGLBRep& lb,
	const MGLBRep* deriS,
	const MGLBRep* deriE
){
	double t0=lb.param_s(), t1=lb.param_e();
	MGNDDArray tau(2); tau(0)=t0; tau(1)=t1;
	MGBPointSeq bp(2,3);
	MGBPointSeq bp1(2,3);
	MGLBRepEndC endcS, endcE;
	bp.store_at(0,lb.start_point());
	if(deriS) endcS.set_1st(lb.eval(t0,1));
	bp.store_at(1,lb.end_point());
	if(deriE) endcE.set_1st(lb.eval(t1,1));
	int error;
	return new MGLBRep(endcS,endcE,tau,bp,error);
}

//compute perimeters and derivatives from only 4 corner point data
//of input perimeters and derivatives.
void get_peri2_deri2(
	const MGPvector<MGLBRep>& perimeters,
	const MGPvector<MGLBRep>& derivatives,
	MGPvector<MGLBRep>& perimeters2,
	MGPvector<MGLBRep>& derivatives2
				//array of derivatives[4], must have size 4.
){
	perimeters2.resize(4);
	derivatives2.resize(4);

	size_t m;
///////////perimeter 0 and 2.
	const MGLBRep* deriS=derivatives[3];
	const MGLBRep* deriE=derivatives[1];
	for(m=0; m<=2; m+=2){//for perimeter 0 and 2.
		perimeters2.reset(m,
			get_perriderisub(*(perimeters[m]),deriS,deriE));
		if(derivatives[m]){
			derivatives2.reset(m,
				get_perriderisub(*(derivatives[m]),deriS,deriE));
		}else derivatives2.reset(m);
	}

////////// perimeter 1 and 3.
	deriS=derivatives[0];
	deriE=derivatives[2];
	for(m=1; m<=3; m+=2){//for perimeter 1 and 3.
		perimeters2.reset(m,
			get_perriderisub(*(perimeters[m]),deriS,deriE));
		if(derivatives[m]){
			derivatives2.reset(m,
				get_perriderisub(*(derivatives[m]),deriS,deriE));
		}else derivatives2.reset(m);
	}
}

MGLBRep get1deri_of_peri(
	size_t peri,	//perimeter number
	const MGSBRep& surf
){
	bool alongu=true; if(peri%2) alongu=false;
	const MGKnotVector* t;
	const MGKnotVector* t_other;
	if(alongu){
		t=&(surf.knot_vector_u());
		t_other=&(surf.knot_vector_v());
	}else{
		t=&(surf.knot_vector_v());
		t_other=&(surf.knot_vector_u());
	}
	double t0=t_other->param_s();
	if(peri==1 || peri==2) t0=t_other->param_e();
	size_t len=t->bdim();
	MGBPointSeq dbp(len,3);
	const MGSPointSeq& spoint=surf.surface_bcoef();
	for(size_t j=0; j<len; j++){
		MGLBRep l0(*t_other,MGBPointSeq(!alongu,j,spoint));
		dbp.store_at(j,l0.eval(t0,1));
	}
	return MGLBRep(*t,dbp);
}

MGSBRep* get_1DireSurf(
	bool udire,	//indicates if perimetes[0],[2] or [3],[1] should be used to construct
				//the surface. if udire=true, [0] and [2] are used.
	const MGPvector<MGLBRep>& perimeters,
	const MGPvector<MGLBRep>& derivatives
){
	size_t ncd=3;

	size_t ids,ide,otherS, otherE;
	if(udire){
		ids=0; ide=2;
		otherS=3; otherE=1;
	}else{
		ids=3; ide=1;
		otherS=0; otherE=2;
	}

//Herafter, u and v directions are temporal. When udire==true, u is v, and v is u.
//Construct u-direction data point and knot vector.
	double u0=perimeters[otherS]->param_s(),u1=perimeters[otherS]->param_e();
	MGENDCOND ecu0=MGENDC_1D, ecu1=MGENDC_1D;
	size_t lenu=2;
	if(!derivatives[ids]) ecu0=MGENDC_NO; else lenu++;
	if(!derivatives[ide]) ecu1=MGENDC_NO; else lenu++;
	size_t orderu=4; if(orderu>lenu) orderu=lenu;
	MGNDDArray utau(lenu);
	size_t i=0;
	utau(i++)=u0;
	if(ecu0==MGENDC_1D) utau(i++)=u0;
	if(ecu1==MGENDC_1D) utau(i++)=u1;
	utau(i++)=u1;
	MGKnotVector tu(utau,orderu);

//Construct v-direction data point and knot vector.
	const MGKnotVector& tv=perimeters[ids]->knot_vector();
	MGNDDArray vtau(3,tv.bdim()-2,tv);
	if(!derivatives[otherS]) vtau.add_data((vtau(0)+vtau(1))*.5);
	if(!derivatives[otherE]){
		size_t nvtau2=vtau.length();
		vtau.add_data((vtau(nvtau2-1)+vtau(nvtau2-2))*.5);		
	}
	//cout<<utau<<tu<<endl;cout<<vtau<<tv<<endl;///////////////////

	double* q=new double[lenu*(2*orderu-1)];
	double* wk=new double[lenu*2];

	size_t lenum1=lenu-1, lenum2=lenu-2;
	MGBPointSeq bp(lenu,ncd);
	size_t nvtau=vtau.length();
	std::vector<MGLBRep*> lines;
	int error=2;
	for(size_t j=0; j<nvtau; j++){
		double vtauj=vtau(j);
		bp.store_at(0,perimeters[ids]->eval(vtauj));
		if(ecu0==MGENDC_1D){
			bp.store_at(1,derivatives[ids]->eval(vtauj));
			MGVector der=derivatives[ids]->eval(vtauj);
			//cout<<j<<",tau="<<vtauj<<","<<der<<", len="<<der.len()<<endl;
		}
		if(ecu1==MGENDC_1D){
			bp.store_at(lenum2,derivatives[ide]->eval(vtauj));
			MGVector der=derivatives[ide]->eval(vtauj);
			//cout<<j<<",tau="<<vtauj<<","<<der<<", len="<<der.len()<<endl;
		}
		bp.store_at(lenum1,perimeters[ide]->eval(vtauj));
		MGLBRep* lb=new MGLBRep(lenu,orderu,ncd);
		MGBPointSeq& obp=lb->line_bcoef();//cout<<bp<<endl;////////
		lb->knot_vector()=tu;
		for(size_t k=0; k<ncd; k++){
			blg4sp2_(orderu,&error,ecu0,ecu1,utau.data(),bp.data(0,k)
				,lenu,lenu,1,tu.data(),1,wk,wk+lenu,q,obp.data(0,k));
		}
		lines.push_back(lb);
	}

	MGPvector<MGLBRep> derivatives2(2);
	size_t id[2]={otherS, otherE};
	for(size_t m=0; m<=1; m++){
	const MGLBRep* deri=derivatives[id[m]];
	if(deri){
		bp.store_at(0,deri->eval(u0));
		if(ecu0==MGENDC_1D) bp.store_at(1,deri->eval(u0,1));
		if(ecu1==MGENDC_1D) bp.store_at(lenum2,deri->eval(u1,1));
		bp.store_at(lenum1,deri->eval(u1));
		MGLBRep* deri2=new MGLBRep(lenu,orderu,ncd);
		MGBPointSeq& obp=deri2->line_bcoef();
		deri2->knot_vector()=tu;
		for(size_t k=0; k<ncd; k++){
			blg4sp2_(orderu,&error,ecu0,ecu1,utau.data(),bp.data(0,k)
				,lenu,lenu,1,tu.data(),1,wk,wk+lenu,q,obp.data(0,k));
		}
		derivatives2.reset(m,deri2);//cout<<"get_1DireSurf::"<<(*deri2)<<endl;
	}
	}

	MGSBRep* sb=new MGSBRep(vtau,lines,derivatives2[0],derivatives2[1]);
	if(udire) sb->exchange_uv();

	delete[] q; delete[] wk;
	size_t nlines=lines.size();
	for(i=0; i<nlines; i++) delete lines[i];
	return sb;
}

//Construct Surface B-rep from lines and derivatives.
//Interpolation will be done only along one parameter direction,
//along v.
//tau provides data point sequence along v direction of surface (u,v) parameter
//configuration. deriS and deriE are used to provide the 1st derivative
//B-representation along the perimeter 0 and 2, may be null
//if 1st derivative B-rep is not provided. If derivative
//B-rep is provided, deriS and deriE must have the same knot configuration
//as the one of lines which makes u kont configuration of this surface (u,v)
//parameter. tau[i] is the parameter for lines[i].
MGSBRep::MGSBRep(
	const MGNDDArray& tau,
	const std::vector<MGLBRep*>& lines,
	const MGLBRep* deriS,
	const MGLBRep* deriE
):MGSurface(){
	size_t i,j,k;
	const size_t ncd=3;

	m_uknot=lines[0]->knot_vector();
	MGENDCOND ec0=MGENDC_1D, ec1=MGENDC_1D;
	if(!deriS) ec0=MGENDC_NO;
	if(!deriE) ec1=MGENDC_NO;

// Compute b-rep dimension along v in lenv.
	//v=min and max condition.
	size_t lenu=m_uknot.bdim(), lenv1,lenv;
	lenv1=lenv=lines.size();
	size_t ivs=0;
	if(ec0==MGENDC_1D){lenv+=1; ivs=1;}
	if(ec1==MGENDC_1D) lenv+=1;
	size_t orderv=4; if(orderv>lenv) orderv=lenv;
	MGSPointSeq surface_bcoef(lenv,lenu,ncd);//Temporal spoint seq.

// Prepare data point ordinate.
	// 1. Copy original data.
	size_t js;
	for(j=0; j<lenv1; j++){
		js=ivs+j;
		const MGBPointSeq& bcj=lines[j]->line_bcoef();
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				surface_bcoef(js,i,k)=bcj(i,k);
			}
		}
	}

	// 2. Copy derivative data along perimeter from endc.
	//v=min condition.
	if(ec0==MGENDC_1D){
		const MGBPointSeq& bp1S=deriS->line_bcoef();
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				surface_bcoef(0,i,k)=bp1S(i,k);
			}
		}
	}

	//v=max condition.
	size_t lenvm1=lenv-1, lenvm2=lenv-2;
	if(ec1==MGENDC_1D){
		const MGBPointSeq& bp1E=deriE->line_bcoef();
		for(k=0; k<ncd; k++){
			for(i=0; i<lenu; i++){
				surface_bcoef(lenvm1,i,k)=bp1E(i,k);
			}
		}
	}

	//Exchange positional data for blg4sp2_.
	//Perimeter 0.
	double save;
	if(ec0==MGENDC_1D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=surface_bcoef(1,i,k);
				surface_bcoef(1,i,k)=surface_bcoef(0,i,k);
				surface_bcoef(0,i,k)=save;
			}
	}
	//Perimeter 2.
	if(ec1==MGENDC_1D){
		for(i=0; i<lenu; i++)
			for(k=0; k<ncd; k++){
				save=surface_bcoef(lenvm2,i,k);
				surface_bcoef(lenvm2,i,k)=surface_bcoef(lenvm1,i,k);
				surface_bcoef(lenvm1,i,k)=save;
			}
	}

	size_t lenuv=lenu*lenv;
	double* q=new double[lenv*(2*orderv-1)];
	double* wk=new double[lenv*2];
	MGNDDArray vtau(lenv);
	for(i=0; i<ivs; i++) vtau(i)=tau(0);
	for(j=0; j<lenv1; j++) vtau(i++)=tau(j);
	for(;i<lenv; i++) vtau(i)=tau(lenv1-1);
	vtau.set_length(lenv);
	m_vknot=MGKnotVector(vtau,orderv);

	int error=2;
	m_surface_bcoef.resize(lenu,lenv,ncd);
	for(k=0; k<ncd; k++){
		blg4sp2_(orderv,&error,ec0,ec1,vtau.data(),surface_bcoef.data(0,0,k)
			,lenv,lenv,lenu,m_vknot.data(),lenu,wk,wk+lenv,q,m_surface_bcoef.data()+lenuv*k);
		if(error!=1) break;
	}
	if(error==1) error=0;
	delete[] q; delete[] wk;
}

void bool_sum(
	const MGCurve*	edge_crvl[4],	//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGSBRepTP& tp,			//接続面(パラメータ範囲は境界線と同じ)
	int& error,	//エラーコード
	MGSBRep& surf
){
	double sinmax[4], taumax[4];
	//cout<<"****** non-vector version//input data ********"<<endl;
	//for(size_t ii2=0; ii2<4; ii2++) cout<<"curve "<<ii2<<"="<<*(edge_crvl[ii2])<<endl;
	//cout<<tp;
	MGSBRepTP tp2(tp);
	MGPvector<MGLBRep> perimeters;
	if(rebuild_curve(edge_crvl,tp2,perimeters))
		return;
	//cout<<endl<<"****** after rebuild ********"<<endl;
	//cout<<perimeters<<endl;	cout<<tp2<<endl;
	MGPosition uv(2);//double cur_gap;////////////

	MGKnotVector& tu=perimeters[0]->knot_vector();
	MGKnotVector& tv=perimeters[1]->knot_vector();
	MGSBRep ruled0(true,perimeters);//cout<<ruled0;
	MGSBRep ruled1(false,perimeters);//cout<<ruled1;
	MGSPointSeq sp(ruled0.surface_bcoef()+ruled1.surface_bcoef());
	sp*=.5;
	MGSBRep ruled01(sp,tu, tv);//cout<<ruled01;
	MGPvector<MGLBRep> derivatives(4);
	get_all_derivatives(perimeters,ruled01,tp2,derivatives);
	//cout<<endl<<" Derivatives="<<derivatives;
	/*
	bool twist[4]={false,false,false,false};
	for(size_t i=0; i<4; i++){
		size_t im1=(i+3)%4;
		if(tp.specified(im1) && tp.specified(i)) twist[i]=true;
	}
	cout<<endl<<"Twist::";
	for(i=0; i<4; i++){
		size_t im1=(i+3)%4;
		if(twist[i]){
			cout<<"at "<<i<<"=";
			double t1=derivatives[im1]->param_e(),t2=derivatives[i]->param_s();
			if(im1==3) t1=derivatives[im1]->param_s();
			if(i>=2) t2=derivatives[i]->param_e();
			cout<<derivatives[im1]->eval(t1,1)<<"//"<<derivatives[i]->eval(t2,1)<<endl;
		}
	}
	*/

	//Save the old knot vector.
	MGKnotVector tu2(tu), tv2(tv);

	double error_max;
	//Construct the boolian sum surface.
	MGSBRep* g1=get_1DireSurf(true,perimeters,derivatives);
//cur_gap = g1->eval_gap(*edge_crvl[0],0, uv);/////////////
//cur_gap = g1->eval_gap(*edge_crvl[2],2, uv);/////////////
//cout<<endl<<" g1="<<(*g1);
//cout<<get1deri_of_peri(3,*g1)<<get1deri_of_peri(1,*g1)<<endl;
//	bool eval[4];
//	eval[0]=true;eval[1]=false;eval[2]=true;eval[3]=false;
//	error_max=tp2.get_perimeters_max_sin(*g1,taumax,sinmax,eval);//////////////

	MGSBRep* g2=get_1DireSurf(false,perimeters,derivatives);
//cur_gap = g2->eval_gap(*edge_crvl[1],1, uv);/////////////
//cur_gap = g2->eval_gap(*edge_crvl[3],3, uv);/////////////
//cout<<endl<<" g2="<<(*g2);
	//cout<<get1deri_of_peri(0,*g2)<<get1deri_of_peri(2,*g2)<<endl;
	//eval[0]=false;eval[1]=true;eval[2]=false;eval[3]=true;
	//error_max=tp2.get_perimeters_max_sin(*g2,taumax,sinmax,eval);//////////////

	MGPvector<MGLBRep> derivatives2(4);
	MGPvector<MGLBRep> perimeters2;
	get_peri2_deri2(perimeters,derivatives,perimeters2,derivatives2);
//cout<<perimeters2;
//cout<<derivatives2;
	std::vector<MGLBRep*> lines(2);
	MGNDDArray vtau0(2); vtau0(0)=tv.param_s(); vtau0(1)=tv.param_e();
	lines[0]=perimeters2[0]; lines[1]=perimeters2[2];
	MGSBRep* g12=new MGSBRep(vtau0,lines,derivatives2[0],derivatives2[2]);
//cout<<endl<<" g12="<<(*g12);

	MGBSumSurf g(g1,g2,g12);
	error_max=tp2.get_perimeters_max_sin(g,taumax,sinmax);//////////////
	error_max*=1.005;
//cur_gap = g.eval_gap(edge_crvl, uv);/////////////
	double azero=MGTolerance::angle_zero();
	if(azero>error_max)
		error_max=azero;

	MGNDDArray utau(3,tu2.bdim()-2,tu2);
	MGNDDArray vtau(3,tv2.bdim()-2,tv2);
	MGSPointSeq spoint(utau.length(),vtau.length(),3);
	g.eval_spoint(utau,vtau,spoint);
	MGSBRepEndC endc(utau,vtau,g);
	MGSBRep* surf2=new MGSBRep(endc,utau,vtau,spoint,tu2,tv2,error);
	if(error){
		error = -7; return ;
	}

//cur_gap = surf2->eval_gap(edge_crvl, uv);/////////////

	double line0=MGTolerance::line_zero();
	double line02=line0*2.;
	MGSBRep* surf3=0;
	//Remove knot by the minimum line zero that satisfy the angle continuity.
	for(size_t q=0; q<4; q++){
		delete surf3;
		surf3=new MGSBRep(*surf2);
		line02*=.5;
		MGTolerance::set_line_zero(line02);
		surf3->remove_knot();
//cur_gap = surf3->eval_gap(edge_crvl, uv);/////////////
		double max2=tp2.get_perimeters_max_sin(*surf3,taumax,sinmax);
		if(max2<=error_max) break;
	}
	delete surf2;
	surf2=surf3;
	MGTolerance::set_line_zero(line0);
	surf=*surf2;
	delete surf2;
	//cout<<surf;
}
