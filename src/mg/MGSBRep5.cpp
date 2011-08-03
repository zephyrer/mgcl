/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Vector.h"
#include "mg/Unit_vector.h"
#include "mg/Position.h"
#include "mg/Transf.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/Tolerance.h"

extern "C"{
#include "cskernel/Blgi2d.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
using namespace std;

// MGSBRep5.cpp
//
// Implements Surface B-Representation class MGSBRep.

//Compute ratio at s of the span (s0, s1).
//Ratio computation of linear case.
void get_ratio2(
	double s0,double s1,double s,
	double& r0, double& r1//r0:ratio at s0 side, r1:ration at s1 side will be output.
){
	r0=(s1-s)/(s1-s0);
	r1=1.-r0;
}

double cfunc(double x){
		double x42=x*x;	x42*=4.;
		return x42-x*x42;
}
//Ratio computation of cubic case.
void get_ratio3(
	double s0,double s1,double s,
	double& r0, double& r1//r0:ratio at s0 side, r1:ration at s1 side will be output.
){
	double x=(s1-s)/(s1-s0);
	if(x<=.5){
		r0=cfunc(x); r1=1.-r0;
	}else{
		r1=cfunc(1.-x); r0=1.-r1;
	}	
}

void get_angle(
	const MGLBRep& tie,	//tie curve of rail0 and 1.
	bool& is_straight_tie0,
	bool& is_straight_tie1,
		//boolian that indicates if tie's tangent vector at the start(is_straight_tie0) or
		//end(is_straight_tie1) point is parallel to the vector from the start point
		//to the end point of tie will be output.
	double& theta0,
	double& theta1
		//angles between the surface's normal and the tie's principal normal
		//at the tie's start(theta0) or end(theta1) point will be output.
){
	double s0=tie.param_s(), s1=tie.param_e();
	MGUnit_vector E1=tie.eval(s1)-tie.eval(s0);

//at the start of the tie.
	MGUnit_vector T0,N0,B0;
	double crvtr0,torsn0;
	tie.Frenet_frame(s0,T0,N0,B0,crvtr0,torsn0);
//	cout<<T0<<N0<<B0<<",k="<<crvtr0<<",t="<<torsn0<<endl;
	MGUnit_vector E20=T0*E1;//is the normal of the virtual surface.
	is_straight_tie0=false;
	if(MGRZero(crvtr0)){
		if(E1.parallel(T0)) is_straight_tie0=true;
		else B0=T0*E1;
	}
	theta0=E20.angle(B0);

//at the end of the tie.
	MGUnit_vector T1,N1,B1;
	double crvtr1,torsn1;
	tie.Frenet_frame(s1,T1,N1,B1,crvtr1,torsn1);
//	cout<<T1<<N1<<B1<<",k="<<crvtr1<<",t="<<torsn1<<endl;
	MGUnit_vector E21=E1*T1;//is the normal of the virtual surface.
	is_straight_tie1=false;
	if(MGRZero(crvtr1)){
		if(E1.parallel(T1)) is_straight_tie0=true;
		else B1=T1*E1;
	}
	theta1=E21.angle(B1);
}

void get_tie(
	double t,			//param of rail0 and 1(the same).
	const MGLBRep& rail0,
	const MGLBRep& rail1,
	const MGLBRep& tie0,	//initial tie(tie at end point).
	const MGPosition& Q0,	//is tie0.start_point().
	const MGPosition& Q1,	//is tie0.end_point().
	bool parallel0, bool parallel1,
		//indicate if tie0's tangent vector at the start(parallel0) or end(parallel1) point
		//is parallel to the vector from Q0 to Q1.
	double angle_s, double angle_e,
		//are angles at the tie0's start(angle_s) or end(angle_e) point between the surface's normal
		//and the tie0's principal normal.
	MGLBRep& tie
){
	MGPosition C=rail0.eval(t), D=rail1.eval(t);
	tie=tie0*MGTransf(Q0,Q1,C,D);
	if(!parallel0 || !parallel1){
		double theta0, theta1;
		bool para_ui0, para_ui1;
		get_angle(tie,para_ui0,para_ui1,theta0,theta1);
		//cout<<endl<<"Befor theta="<<theta0<<","<<theta1<<endl;
		double dif0, dif1;
		if(parallel0 || para_ui0) dif0=0.; else dif0=theta0-angle_s;
		if(parallel1 || para_ui1) dif1=0.; else dif1=theta1-angle_e;
		double dif=(dif0+dif1)*.5;//cout<<" dif0="<<dif0<<",dif1="<<dif1<<",dif="<<dif<<endl;
		MGTransf tr; tr.set_rotate_3D(D-C,dif,C);
		tie*=tr;
		get_angle(tie,para_ui0,para_ui1,theta0,theta1);
		//cout<<"After theta="<<theta0<<","<<theta1<<endl;
	}
}

//Given two curves(peri0 and 3), get a curve(peri1) that is peri3's slide along peri0
//from the peri0's start point to the end point.
//peri0's start point and peri3's start point must coincide.
void get_peri031(
	MGPvector<MGLBRep>& perimeters2
){
	const MGLBRep& peri0=*(perimeters2[0]);
	const MGLBRep& peri3=*(perimeters2[3]);
	perimeters2.assign(1,new MGLBRep);
	MGLBRep& peri1=*(perimeters2[1]);
	peri1=peri3;
	MGVector V0=peri0.eval(peri0.param_s(),1), V1=peri0.eval(peri0.param_e(),1);
	MGMatrix M; M.set_rotate(V0,V1);
	peri1-=peri3.start_point();
	peri1*=M;
	peri1+=peri0.end_point();
}
//Given two curves(peri0 and 1), get a curve(peri3) that is peri1's slide along peri0
//from the peri0's end point to the start point.
//peri0's end point and peri1's start point must coincide.
void get_peri013(
	MGPvector<MGLBRep>& perimeters2
){
	const MGLBRep& peri0=*(perimeters2[0]);
	const MGLBRep& peri1=*(perimeters2[1]);
	perimeters2.assign(3,new MGLBRep);
	MGLBRep& peri3=*(perimeters2[3]);
	peri3=peri1;
	MGVector V0=peri0.eval(peri0.param_s(),1), V1=peri0.eval(peri0.param_e(),1);
	MGMatrix M; M.set_rotate(V1,V0);
	peri3-=peri1.start_point();
	peri3*=M;
	peri3+=peri0.start_point();
}
//Given two curves(peri1 and 2), get a curve(peri0) that is peri2's slide along peri1
//from the peri1's end point to the start point.
//peri1's end point and peri2's end point must coincide.
void get_peri120(
	MGPvector<MGLBRep>& perimeters2
){
	const MGLBRep& peri1=*(perimeters2[1]);
	const MGLBRep& peri2=*(perimeters2[2]);
	perimeters2.assign(0,new MGLBRep);
	MGLBRep& peri0=*(perimeters2[0]);
	peri0=peri2;
	MGVector V0=peri1.eval(peri1.param_s(),1), V1=peri1.eval(peri1.param_e(),1);
	MGMatrix M; M.set_rotate(V1,V0);
	peri0-=peri2.end_point();
	peri0*=M;
	peri0+=peri1.start_point();
}
//Given two curves(peri2 and 3), get a curve(peri0) that is peri2's slide along peri3
//from the peri3's end point to the start point.
//peri3's end point and peri2's start point must coincide.
void get_peri230(
	MGPvector<MGLBRep>& perimeters2
){
	const MGLBRep& peri3=*(perimeters2[3]);
	const MGLBRep& peri2=*(perimeters2[2]);
	perimeters2.assign(0,new MGLBRep);
	MGLBRep& peri0=*(perimeters2[0]);
	peri0=peri2;
	MGVector V0=peri3.eval(peri3.param_s(),1), V1=peri3.eval(peri3.param_e(),1);
	MGMatrix M; M.set_rotate(V1,V0);
	peri0-=peri2.start_point();
	peri0*=M;
	peri0+=peri3.start_point();
}

//Generalized ruled surface construction.
//Build a surface by ruled surface method. That is, costruct a surface by sliding
//the blending curve(a tie curve of the rails) of perimeters 3 and 1 along
//perimeter 0 and 2(0 and 2 make the rail).
//Or by sliding the blending curve of perimeter 0 and 2
//along perimeter 3 and 1(3 and 1 make the rail).
MGSBRep::MGSBRep(
	bool along_u,	//indicates which perimeters make a rail.
			//if true, perimeter 0 and 2, if false, perimeter 3 and 1 make rail.
	const MGPvector<MGLBRep>& perimeters
			//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
			//perimeters must be the same knot configuration along u(perimeter 0 and 2)
			//and along v(perimeter 3 and1).
){
	const MGLBRep& peri0=*(perimeters[0]);
	const MGLBRep& peri1=*(perimeters[1]);
	const MGLBRep& peri2=*(perimeters[2]);
	const MGLBRep& peri3=*(perimeters[3]);
	size_t sd=peri0.sdim();
	size_t nu=peri0.bdim(), nv=peri1.bdim();
	size_t num1=nu-1, nvm1=nv-1;
	size_t len=nu; if(len<nv) len=nv;
	size_t one=1;

	double u0=peri0.param_s(), u1=peri0.param_e();
	double uspan=u1-u0;
	double v0=peri1.param_s(), v1=peri1.param_e();
	double vspan=v1-v0;
	const MGKnotVector& uknot=peri0.knot_vector();
	const MGKnotVector& vknot=peri1.knot_vector();
	size_t orderu=uknot.order(), orderv=vknot.order();
	size_t order=orderu; if(order<orderv) order=orderv;
	double* q=new double[len*(2*order-1)];
	double* work=new double[len];

	MGPosition P[4]=//Four corner point.
		{peri3.start_point(), peri1.start_point(),peri1.end_point(), peri3.end_point()};		

	double alpha0, alpha1, alpha2, alpha3;
	bool parallel0, parallel1, parallel2,parallel3;
	m_uknot=uknot;m_vknot=vknot;
	m_surface_bcoef.resize(nu,nv,sd);
	double r0,r1;

	if(along_u){
		MGSPointSeq sp(nu,nv,sd);//Temporal spoint.
		size_t usize=sp.capacity_u();
		MGNDDArray utau; utau.update_from_knot(uknot);
		get_angle(peri3,parallel0,parallel3,alpha0,alpha3);
		get_angle(peri1,parallel1,parallel2,alpha1,alpha2);

		sp.store_BP_along_v_at(0,peri3.line_bcoef());
		for(size_t i=1; i<num1; i++){
			double ui=utau[i];
			MGLBRep tiei0;
			get_tie(ui,peri0,peri2,peri3,P[0],P[3],parallel0,parallel3,alpha0,alpha3,tiei0);
			MGLBRep tiei1;
			get_tie(ui,peri0,peri2,peri1,P[1],P[2],parallel1,parallel2,alpha1,alpha2,tiei1);
			get_ratio2(u0,u1,ui,r0,r1);//***********
			sp.store_BP_along_v_at(i,tiei0.line_bcoef()*r0+tiei1.line_bcoef()*r1);
		}
		sp.store_BP_along_v_at(num1,peri1.line_bcoef());
		int error=2;
		for(size_t k=0; k<sd; k++){
			for(size_t j=0; j<nv; j++)
				blgi2d_(&error,utau.data(),sp.data(0,j,k),m_uknot.data(),
				orderu,nu,one,usize,one,work,q,&m_surface_bcoef(0,j,k));
			if(error!=1) break;
		}
	}else{
		MGSPointSeq sp(nv,nu,sd);//Temporal spoint.
		size_t vsize=sp.capacity_u();
		MGNDDArray vtau; vtau.update_from_knot(vknot);
		get_angle(peri0,parallel0,parallel1,alpha0,alpha1);
		get_angle(peri2,parallel3,parallel2,alpha3,alpha2);

		sp.store_BP_along_v_at(0,peri0.line_bcoef());
		for(size_t j=1; j<nvm1; j++){
			double vj=vtau[j];
			MGLBRep tiej0;
			get_tie(vj,peri3,peri1,peri0,P[0],P[1],parallel0,parallel1,alpha0,alpha1,tiej0);
			MGLBRep tiej1;
			get_tie(vj,peri3,peri1,peri2,P[3],P[2],parallel3,parallel2,alpha3,alpha2,tiej1);
			get_ratio2(v0,v1,vj,r0,r1);//********
			sp.store_BP_along_v_at(j,tiej0.line_bcoef()*r0+tiej1.line_bcoef()*r1);
		}
		sp.store_BP_along_v_at(nvm1,peri2.line_bcoef());
		int error=2;
		for(size_t k=0; k<sd; k++){
			blgi2d_(&error,vtau.data(),sp.data(0,0,k),m_vknot.data(),
				orderv,nv,nu,vsize,nu,work,q,&m_surface_bcoef(0,0,k));
			if(error!=1) break;
		}
	}
	delete[]  q; delete[]  work;
}

//Construct 4 perimeters, given at least two of the four.
//Input perimeters may have different knot configuration. In this case they will be updated
//so as to have the same configuration.
//Function's return value indicates which perimeter(s) was missing:
// 10: all of the 4 were input(and knot configurations were updated to have the same).
//  0: only perimeter 0 was missing.
//  1: only perimeter 1 was missing.
//  2: only perimeter 2 was missing.
//  3: only perimeter 3 was missing.
//  4: perimeter 2 and 3 were missing.
//  5: perimeter 1 and 3 were missing.
//  6: perimeter 1 and 2 were missing.
//  7: perimeter 0 and 3 were missing.
//  8: perimeter 0 and 2 were missing.
//  9: perimeter 0 and 1 were missing.
// -2: less than 2 perimeters were provided.
int construct_perimeters(
	const MGPvector<MGCurve>& peris,
			//境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順). Let i be the perimeter number,
			//and the data is missing, perimeter[i] must be null. If perimeter 3 data is missing,
			//perimeters.size() may be 3. If perimeter 2 and 3 data are missing, perimeters.size() may
			//be 2.
			//When perimeters were not the same knot configuration along u(perimeter 0 and 2)
			//or along v(perimeter 3 and1), they will be rebuild to have the same configuration.
	MGPvector<MGLBRep>& perimeters2
			//new perimeters will be output.
){
	size_t n=peris.size();if(n<2) return -2;

	std::vector<const MGCurve*> perimeters(4);
	size_t i;
	for(i=0; i<n; i++) perimeters[i]=peris[i];
	for(i=n; i<4; i++) perimeters[i]=0;

//1. Classify which perimeters are input.
	int type;
	if(!perimeters[0]){//0 missing
		if(!perimeters[3]){//3 missing
			if(!perimeters[1] || !perimeters[2]) return -2;
			type=7;
		}else{//3 provided
			if(perimeters[1]){//1 provided
				if(perimeters[2]) type=0;
				else              type=8;
			}else{
				if(!perimeters[2]) return -2;
				type=9;
			}
		}
	}else{//0 provided
		if(perimeters[1]){//1 provided
			if(perimeters[2]){//2 provided
				if(perimeters[3]) type=10;
				else type=3;
			}else{
				if(perimeters[3]) type=2;
				else type=4;
			}
		}else{//1 missing
			if(perimeters[2]){//2 provided
				if(perimeters[3]) type=1;
				else type=5;
			}else{
				if(perimeters[3]) type=6;
				else return -2;
			}
		}
	}

//2. knot configuration adjustment.
	perimeters2.clear(); perimeters2.resize(4);
	std::vector<const MGCurve*> peri2(2);
	if(type==10 || type==1 || type==3 || type==5){//when 0 and 2 provided.
		peri2[0]=perimeters[0];
		peri2[1]=perimeters[2];
		MGPvector<MGLBRep> lb2=rebuild_knot(peri2);
		perimeters2.assign(0,lb2.release(0));
		perimeters2.assign(2,lb2.release(1));
	}
	if(type==10 || type==0 || type==2 || type==8){//when 1 and 3 provided.
		peri2[0]=perimeters[1];
		peri2[1]=perimeters[3];
		MGPvector<MGLBRep> lb2=rebuild_knot(peri2);
		perimeters2.assign(1,lb2.release(0));
		perimeters2.assign(3,lb2.release(1));
	}
	if(type==10) return 10;

//3. process when parallel two perimeters are provided.
	if(type==5){
		perimeters2.assign(1,new MGLBRep(
			MGStraight(perimeters2[2]->end_point(),perimeters2[0]->end_point())));
		perimeters2.assign(3,new MGLBRep(
			MGStraight(perimeters2[2]->start_point(),perimeters2[0]->start_point())));
		return 5;
	}
	if(type==8){
		perimeters2.assign(0,new MGLBRep(
			MGStraight(perimeters2[1]->start_point(),perimeters2[3]->start_point())));
		perimeters2.assign(2,new MGLBRep(
			MGStraight(perimeters2[1]->end_point(),perimeters2[3]->end_point())));
		return 8;
	}

//4. approximate the perimeter that has not the opposite one.
	for(i=0; i<4; i++){
		if(!perimeters2[i]){
			if(perimeters[i]) perimeters2.assign(i,new MGLBRep(*(perimeters[i])));
		}
	}

//5. Adjust the common corner points.
	MGLBRep& l0=*(perimeters2[0]);
	MGLBRep& l1=*(perimeters2[1]);
	MGLBRep& l2=*(perimeters2[2]);
	MGLBRep& l3=*(perimeters2[3]);

	MGPosition P0, P1, P2, P3;
	if(type==1 || type==2 || type==6) P0=(l0.start_point()+l3.start_point())*.5;
	if(type==2 || type==3 || type==4) P1=(l0.end_point()+l1.start_point())*.5;
	if(type==0 || type==3 || type==7) P2=(l1.end_point()+l2.end_point())*.5;
	if(type==0 || type==1 || type==9) P3=(l3.end_point()+l2.start_point())*.5;

	double fixp[2];
	if(!P0.is_null()){//P0
		fixp[0]=l0.param_e(); l0.move(2,l0.param_s(),P0,fixp);
		fixp[0]=l3.param_e(); l3.move(2,l3.param_s(),P0,fixp);
	}
	if(!P1.is_null()){//P1
		fixp[0]=l0.param_s(); l0.move(2,l0.param_e(),P1,fixp);
		fixp[0]=l1.param_e(); l1.move(2,l1.param_s(),P1,fixp);
	}
	if(!P2.is_null()){//P2
		fixp[0]=l1.param_s(); l1.move(2,l1.param_e(),P2,fixp);
		fixp[0]=l2.param_s(); l2.move(2,l2.param_e(),P2,fixp);
	}
	if(!P3.is_null()){//P3
		fixp[0]=l2.param_e(); l2.move(2,l2.param_s(),P3,fixp);
		fixp[0]=l3.param_s(); l3.move(2,l3.param_e(),P3,fixp);
	}

//6. process when non-parallel two perimeters are provided.
	int type2=type;
	if(type==4){
		get_peri013(perimeters2);type2=2;
	}else if(type==6){
		get_peri031(perimeters2);type2=2;
	}else if(type==7){
		get_peri120(perimeters2);type2=3;
	}else if(type==9){
		get_peri230(perimeters2);type2=1;
	}
	//cout<<perimeters2;

//7. Now 3 of the 4 perimeters are constructed, construct the missing one.
	bool para0,para1;
	double theta0, theta1;
	perimeters2.assign(type2,new MGLBRep);
	MGLBRep& peri=*(perimeters2[type2]);//peri to be constructed.

	MGLBRep& peri_opo=*(perimeters2[(type2+2)%4]);//opposide perimeter of peri.
	size_t is;
	if(type2<=1) is=(type2+3)%4; else is=(type2+1)%4;
	MGLBRep& peri_s=*(perimeters2[is]);//starting side perimeter of peri.
	MGLBRep& peri_e=*(perimeters2[(is+2)%4]);//ending side perimeter of peri.

	get_angle(peri_opo,para0,para1,theta0,theta1);
	MGPosition Q0=peri_opo.start_point(), Q1=peri_opo.end_point();
	double t;
	if(type2==0 || type2==3) t=peri_s.param_s(); else t=peri_s.param_e();
	get_tie(t,peri_s,peri_e,peri_opo,Q0,Q1,para0,para1,theta0,theta1,peri);
	return type;
}
