/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Pvector.h"
#include "mg/Tolerance.h"
#include "mg/SPointSeq.h"
#include "mg/LBRepEndC.h"
#include "mg/SBRepEndC.h"
#include "mg/SBRepVecTP.h"
#include "mg/SBRep.h"
#include "mg/BSumSurf.h"
using namespace std;

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
//using namespace std;

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

MGSBRep* get_1DireSurf(
	bool udire,	//indicates if perimetes[0],[2] or [3],[1] should be used to construct
				//the surface. if udire=true, [0] and [2] are used.
	const MGPvector<MGLBRep>& perimeters,
	const MGPvector<MGLBRep>& derivatives
);

//compute perimeters and derivatives from only 4 corner point data
//of input perimeters and derivatives.
void get_peri2_deri2(
	const MGPvector<MGLBRep>& perimeters,
	const MGPvector<MGLBRep>& derivatives,
	MGPvector<MGLBRep>& perimeters2,
	MGPvector<MGLBRep>& derivatives2
				//array of derivatives[4], must have size 4.
);

void get_all_derivatives(
 	const MGPvector<MGLBRep>& perimeters,//perimeters.
	const MGSBRep& surf,	//Original surface.
	const MGSBRepVecTP& vectp,	//接続面(パラメータ範囲は境界線と同じ)
	MGPvector<MGLBRep>& derivatives
				//array of derivatives[4], must have size 4.
){
	size_t sd=surf.sdim();
	double azero=MGTolerance::angle_zero()*.05;
	const MGKnotVector& tu=surf.knot_vector_u();
	const MGKnotVector& tv=surf.knot_vector_v();

	size_t m;
	MGVector N;
///////////perimeter 0 and 2.
	size_t nu=surf.bdim_u(); size_t num1=nu-1;
	MGNDDArray utau2; utau2.update_from_knot(tu);
	MGBPointSeq deris_u(nu,sd);
	double v[2]={tv.param_s(),tv.param_e()};

	MGNDDArray utau(nu); utau(0)=utau2[0];
	for(m=0; m<2; m++){//for perimeter 0 and 2.
		size_t perim=m*2;
		if(!vectp.specified(perim)) continue;

		double vm=v[m];
		size_t i2=0;
		MGVector deri=perimeters[3]->eval(vm,1);//cout<<deri;
		deris_u.store_at(i2++,deri);
		for(size_t i=1; i<num1; i++){
			double utaui=utau2[i];
			if(vectp.eval(perim,utaui,N)){//get Normal of tangent plane.
				deri=surf.eval(utaui,vm,0,1);//cout<<deri;
				double vlen=deri.len();
				deri=MGUnit_vector(N*deri*N)*vlen;//cout<<deri<<endl;
				deris_u.store_at(i2,deri);
				utau(i2++)=utaui;
			}
		}
		deri=perimeters[1]->eval(vm,1);//cout<<deri;
		deris_u.store_at(i2,deri);
		utau(i2++)=utau2[num1];
		utau.set_length(i2);
		deris_u.set_length(i2);//cout<<utau<<deris_u<<endl;
		derivatives.reset(perim,new MGLBRep(utau,deris_u));
		//cout<<"driv at perim="<<perim<<","<<*(derivatives[perim])<<endl;
	}

////////// perimeter 1 and 3.
	size_t nv=surf.bdim_v(); size_t nvm1=nv-1;
	MGNDDArray vtau2; vtau2.update_from_knot(tv);
	MGBPointSeq deris_v(nv,sd);
	double u[2]={tu.param_s(),tu.param_e()};
	size_t peri[2]={3,1};

	MGNDDArray vtau(nv); vtau(0)=vtau2[0];
	for(m=0; m<2; m++){//for perimeter 3 and 1.
		size_t perim=peri[m];
		if(!vectp.specified(perim))continue;

		double um=u[m];
		size_t j2=0;
		MGVector deri=perimeters[0]->eval(um,1);//cout<<deri;
		deris_v.store_at(j2++,deri);
		for(size_t j=1; j<nvm1; j++){
			double vtauj=vtau2[j];
			if(vectp.eval(perim,vtauj,N)){//Normal of tangent plane at surf(um,vtau[j]).
				deri=surf.eval(um,vtauj,1,0);//cout<<deri;				
				double vlen=deri.len();
				deri=MGUnit_vector(N*deri*N)*vlen;//cout<<deri<<endl;
				deris_v.store_at(j2,deri);
				vtau(j2++)=vtauj;
			}
		}
		deri=perimeters[2]->eval(um,1);//cout<<deri;
		deris_v.store_at(j2,deri);
		vtau(j2++)=vtau2[nvm1];
		vtau.set_length(j2);
		deris_v.set_length(j2);
		derivatives.reset(perim,new MGLBRep(vtau,deris_v));
		//cout<<"driv at perim="<<perim<<","<<*(derivatives[perim])<<endl;
	}

	size_t i,l;
//*****Twist adjustment. Here dudv and dvdu are only averaged, may be able to improve.
//Adjust the twist vector at the corners.
	MGVector twist[4];
	for(i=0; i<4; i++){
		size_t im1=(i+3)%4;
		if(vectp.specified(im1) && vectp.specified(i)){
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
	MGNDDArray utau3(3,tu.bdim()-2,tu);
	MGNDDArray vtau3(3,tv.bdim()-2,tv);
	size_t ntau2,nutau2=utau3.length(), nvtau2=vtau3.length();
	MGBPointSeq dyu(nutau2,1), dyv(nvtau2,1);
	for(l=0; l<nutau2; l++) dyu(l,0)=1.;
	dyu(0,0)=.01;dyu(nutau2-1,0)=.01;
	for(l=0; l<nvtau2; l++) dyv(l,0)=1.;
	dyv(0,0)=.01;dyv(nvtau2-1,0)=.01;
	MGNDDArray* tau2;
	const MGKnotVector* t;
	double* dy;
	for(i=0; i<4; i++){
		if(!vectp.specified(i))continue;
		double ts,te;
		if(i%2) {
			tau2=&vtau3; t=&tv; dy=dyv.data(); ntau2=nvtau2;
		}else{
			tau2=&utau3; t=&tu; dy=dyu.data(); ntau2=nutau2;
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
		dev*=MGTolerance::angle_zero()*.7;
		MGLBRep* lb=new MGLBRep;
		lb->buildSRSmoothedLB_of_1stDeriv(ECS,ECE,*tau2,bp,dy,dev,false);
		derivatives.reset(i,lb);
		//cout<<"driv at perim="<<i<<","<<*(derivatives[i])<<endl;
	}
}

void bool_sum(
	const MGCurve*	edge_crvl[4],	//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	MGSBRepVecTP& vectp,			//接続面(パラメータ範囲は境界線と同じ)
	int& error,	//エラーコード
	MGSBRep& surf
){
	double sinmax[4], taumax[4];
	//cout<<"****** VECTOR version//input data ********"<<endl;
	//for(size_t ii2=0; ii2<4; ii2++) cout<<"curve "<<ii2<<"="<<*(edge_crvl[ii2])<<endl;
	//cout<<vectp;
	MGSBRepTP tp2;
	MGPvector<MGLBRep> perimeters;
	if(rebuild_curve(edge_crvl,tp2,perimeters)) return;
	//cout<<endl<<"****** after rebuild ********"<<endl;
	//cout<<perimeters<<endl;

	MGKnotVector& tu=perimeters[0]->knot_vector();
	MGKnotVector& tv=perimeters[1]->knot_vector();
	vectp.change_range(true,tu.param_s(),tu.param_e());
	vectp.change_range(false,tv.param_s(),tv.param_e());
	//cout<<vectp<<endl;

	MGSBRep ruled0(true,perimeters);//cout<<ruled0;
	MGSBRep ruled1(false,perimeters);//cout<<ruled1;
	MGSPointSeq sp(ruled0.surface_bcoef()+ruled1.surface_bcoef());
	sp*=.5;
	MGSBRep ruled01(sp,tu, tv);//cout<<ruled01;
	MGPvector<MGLBRep> derivatives(4);
	get_all_derivatives(perimeters,ruled01,vectp,derivatives);
	//cout<<endl<<" Derivatives="<<derivatives;

	////////////////******
	/*bool twist[4]={false,false,false,false};
	for(size_t i=0; i<4; i++){
		size_t im1=(i+3)%4;
		if(vectp.specified(im1) && vectp.specified(i)) twist[i]=true;
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
	}*/
	//********************////////////////////

	//Save the old knot vector.
	MGKnotVector tu2(tu), tv2(tv);

	double error_max;
	//Construct the boolian sum surface.
	MGSBRep* g1=get_1DireSurf(true,perimeters,derivatives);
/*	cout<<endl<<" g1="<<(*g1);cout<<get1deri_of_peri(3,*g1)<<get1deri_of_peri(1,*g1)<<endl;
	bool eval[4];
	eval[0]=true;eval[1]=false;eval[2]=true;eval[3]=false;
	error_max=vectp.get_perimeters_max_sin(*g1,taumax,sinmax,eval);//////////////
*/
	MGSBRep* g2=get_1DireSurf(false,perimeters,derivatives);
	//cout<<endl<<" g2="<<(*g2);
	//cout<<get1deri_of_peri(0,*g2)<<get1deri_of_peri(2,*g2)<<endl;
	//eval[0]=false;eval[1]=true;eval[2]=false;eval[3]=true;
	//error_max=vectp.get_perimeters_max_sin(*g2,taumax,sinmax,eval);//////////////

	MGPvector<MGLBRep> derivatives2(4);
	MGPvector<MGLBRep> perimeters2;
	get_peri2_deri2(perimeters,derivatives,perimeters2,derivatives2);
	//cout<<derivatives2;
	std::vector<MGLBRep*> lines(2);
	MGNDDArray vtau0(2); vtau0(0)=tv.param_s(); vtau0(1)=tv.param_e();
	lines[0]=perimeters2[0]; lines[1]=perimeters2[2];
	MGSBRep* g12=new MGSBRep(vtau0,lines,derivatives2[0],derivatives2[2]);
	//cout<<endl<<" g12="<<(*g12);

	MGBSumSurf g(g1,g2,g12);
	error_max=vectp.get_perimeters_max_sin(g,taumax,sinmax);//////////////
	error_max*=1.005;
	double azero=MGTolerance::angle_zero();
	if(azero>error_max) error_max=azero;

	MGNDDArray utau(3,tu2.bdim()-2,tu2);
	MGNDDArray vtau(3,tv2.bdim()-2,tv2);
	MGSPointSeq spoint(utau.length(),vtau.length(),3);
	g.eval_spoint(utau,vtau,spoint);
	MGSBRepEndC endc(utau,vtau,g);
	MGSBRep* surf2=new MGSBRep(endc,utau,vtau,spoint,tu2,tv2,error);
	if(error){ error = -7; return ;}
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
		double max2=vectp.get_perimeters_max_sin(*surf3,taumax,sinmax);
		if(max2<=error_max) break;
	}
	delete surf2;
	surf2=surf3;

	MGTolerance::set_line_zero(line0);
	surf=*surf2;
	delete surf2;
}
