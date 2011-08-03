/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Vector.h"
#include "mg/Unit_vector.h"
#include "mg/Position.h"
#include "mg/Position_list.h"
#include "mg/LBRep.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/Bluprt.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGSBRep4.cpp
//
// Implements Surface B-Representation class MGSBRep.

//Compute continuity with brep2.
int MGSBRep::continuity(	// Reuturn value is the continuity.
	const MGSBRep& brep2,	// Input second SBRep
	int& which1,	// Outputs which perimeter(which1) of this is
	int& which2,	// connected to which(which2) of brep2.
					// These are valid only when continuity>=0.
	int& opposite,	// Outputs if parameter direction of which2
					// is the same as which1 along common edge.
					// If opposite is true, the direction is opposite.
	double& ratio	// Ratio of 1st derivatives of the two surfaces will
					// be returned.
			// ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
	) const
// Function's return value is:
// -1: G(-1) continuity, i.e. two surfaces are discontinuous.
//  0: G0 continuity, i.e. two surfaces are connected,
//     but tangents are discontinuous
//  1: G1 continuity, i.e. two surfaces are connected,
//     and tangents are also continuous.
//  2: G2 continuity, i.e. two surfaces are connected,
//     and tangents and curvatures are also continuous.
{
	int cont;

	// 1.Test if perimeter 0 or 2 of this is continuous to perimeter 0 or 2 of
	//   brep2.
	// 1.1 Test if knot-vectors are equal.
	if(knot_vector_u()==brep2.knot_vector_u()){
		opposite=0;
		cont=continuity(brep2,1,1,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}
	if(knot_vector_u()==-brep2.knot_vector_u()){
		opposite=1;
		cont=continuity(brep2,1,1,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}

	// 2.Test if perimeter 0 or 2 of this is continuous to perimeter 1 or 3 of
	//   brep2.
	// 2.1 Test if knot-vectors are equal.
	if(knot_vector_u()==brep2.knot_vector_v()){
		opposite=0;
		cont=continuity(brep2,1,0,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}
	if(knot_vector_u()==-brep2.knot_vector_v()){
		opposite=1;
		cont=continuity(brep2,1,0,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}

	// 3.Test if perimeter 1 or 3 of this is continuous to perimeter 0 or 2 of
	//   brep2.
	// 3.1 Test if knot-vectors are equal.
	if(knot_vector_v()==brep2.knot_vector_u()){
		opposite=0;
		cont=continuity(brep2,0,1,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}
	if(knot_vector_v()==-brep2.knot_vector_u()){
		opposite=1;
		cont=continuity(brep2,0,1,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}

	// 4.Test if perimeter 1 or 3 of this is continuous to perimeter 1 or 3 of
	//   brep2.
	// 4.1 Test if knot-vectors are equal.
	if(knot_vector_v()==brep2.knot_vector_v()){
		opposite=0;
		cont=continuity(brep2,0,0,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}
	if(knot_vector_v()==-brep2.knot_vector_v()){
		opposite=1;
		cont=continuity(brep2,0,0,opposite,which1,which2,ratio);
		if(cont>=0) return cont;
	}

	return -1;
}

//Compute continuity with brep2.
int MGSBRep::continuity(	// Reuturn value is the continuity.
	const MGSBRep& brep2,	// Input second SBRep
	int is_u1,		// Input if u-direction of this.
	int is_u2,		// Input if u-direction of brep2.
	int opposite,	// Input if parameter direction of which2 is equal or not.
	int& which1,	// Outputs which perimeter(which1) of this is
	int& which2,	// connected to which(which2) of brep2.
					// These are valid only when continuity>=0.
	double& ratio	// Ratio of 1st derivatives of the two surfaces will
					// be returned.
			// ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
	) const
	// Function's return value is:
	// -1: G(-1) continuity, i.e. two surfaces are discontinuous.
	//  0: G0 continuity, i.e. two surfaces are connected,
	//     but tangents are discontinuous
	//  1: G1 continuity, i.e. two surfaces are connected,
	//     and tangents are also continuous.
	//  2: G2 continuity, i.e. two surfaces are connected,
	//     and tangents and curvatures are also continuous.
{
	size_t i,j,k,i2; int incrmnt;
	double ratio2; int which;
	int cont, contold;
	size_t dim1=sdim(), dim2=brep2.sdim();
	size_t ns1,nt1, ns2,nt2;
	MGLBRep p1a, p1b, p2a, p2b;

	// Test if perimeter of this is continuous to perimeter brep2.
	const MGKnotVector *s1,*s2,*t1,*t2;
	if(is_u1){
		which1=0; 
		s1=&knot_vector_u(); t1=&knot_vector_v();
		p1a=perimeter(0); p1b=perimeter(2);
	}
	else{
		which1=1;
		s1=&knot_vector_v(); t1=&knot_vector_u();
		p1a=perimeter(1); p1b=perimeter(3);
	}
	if(is_u2){
		which2=0;
		s2=&(brep2.knot_vector_u()); t2=&(brep2.knot_vector_v());
		p2a=brep2.perimeter(0); p2b=brep2.perimeter(2);
	}
	else{
		which2=1;
		s2=&(brep2.knot_vector_v()); t2=&(brep2.knot_vector_u());
		p2a=brep2.perimeter(1); p2b=brep2.perimeter(3);
	}
	ns1=(*s1).bdim(); ns2=(*s2).bdim();
	nt1=(*t1).bdim(); nt2=(*t2).bdim();

	cont=0;
	// 1. Test if positional data of two perimeters are equal.
	if(opposite){ p2a.negate(); p2b.negate();}
	if     (p1a.line_bcoef()==p2a.line_bcoef()){
		cont=1;}
	else if(p1a.line_bcoef()==p2b.line_bcoef()){
		which2+=2; cont=1;}
	else if(p1b.line_bcoef()==p2a.line_bcoef()){
		which1+=2; cont=1;}
	else if(p1b.line_bcoef()==p2b.line_bcoef()){
		which1+=2; which2+=2; cont=1;}
	if(cont==0) return -1;

	// There exists a possibility of continuity 1.
	// 2. Test if derivatives along v direction are equal.
	i2=0; incrmnt=1; if(opposite) {i2=ns2-1; incrmnt=-1;}
	MGBPointSeq b1(nt1,dim1), b2(nt2,dim2);
	contold=0;
	for(i=0; i<ns1; i++){
		if(is_u1){
			for(j=0; j<nt1; j++)
				for(k=0; k<dim1; k++) b1(j,k)=coef(i,j,k);
		}
		else{
			for(j=0; j<nt1; j++)
				for(k=0; k<dim1; k++) b1(j,k)=coef(j,i,k);
		}
		if(is_u2){
			for(j=0; j<nt2; j++)
				for(k=0; k<dim2; k++) b2(j,k)=brep2.coef(i2,j,k);
		}
		else{
			for(j=0; j<nt2; j++)
				for(k=0; k<dim2; k++) b2(j,k)=brep2.coef(j,i2,k);
		}
		i2=i2+incrmnt;

		MGLBRep lb1(*t1,b1), lb2(*t2, b2);
		cont=lb1.continuity(lb2,which,ratio2);
		if(cont<=0) return 0;			          //Continuity is C0.
		if(contold==0) {contold=cont; ratio=ratio2;}
		else{
			if(!MGREqual2(ratio2,ratio)) return 0; //Continuity is C0.
			else if(contold>cont) contold=cont;
		}
	}
	return contold;
}

//Test if the surface is planar or not.
//Returned is 0(false) if this is not planar, 1(true) if this is planar.
int MGSBRep::planar(
	MGPlane& plane,		//Plane that might be closest to this.
						//Plane is always output even if not planar.
	double& deviation	//maximum deviation of this from the output plane.
	) const{
	double u0=param_s_u(), u1=param_e_u();
	double v0=param_s_v(), v1=param_e_v();
	MGUnit_vector N=(normal(u0,v0)+normal(u0,v1)+normal(u1,v0)+normal(u1,v1));
	MGPosition P=(eval(u0,v0)+eval(u0,v1)+eval(u1,v0)+eval(u1,v1)
					+eval((u0+u1)*.5,(v0+v1)*.5))/5.;
	plane=MGPlane(N,P);

	size_t m=bdim_u(), n=bdim_v();
	deviation=0.;
	double d=P%N;
	double nx=N.ref(0), ny=N.ref(1), nz=N.ref(2);
	for(size_t i=0; i<m; i++){
	for(size_t j=0; j<n; j++){
		double x=d - (nx*coef(i,j,0)+ny*coef(i,j,1)+nz*coef(i,j,2));
		if(x<0.) x=-x;
		if(deviation<x) deviation=x;
	}
	}
	return MGAZero(deviation);

}

//Test if part of the surface is planar or not within the tolerance tol.
//The part of the surface is input by the surface parameter range uvbox.
//Returned is 0(false) if this is not planar, 1(true) if planar.
int MGSBRep::planar(
	const MGBox& uvbox,//This surface parameter range.
	double tol,	//maximum deviation allowed to regard the sub surface as a plane.
	int* divideU//Direction to subdivide will be output, if this was not planar.
				//=1: u direction, =0: v direction.
	) const{
	MGBox uvb=param_range()&uvbox;
	MGPosition P; 
	MGUnit_vector N;
	int direction;
	if(!flat(uvb,tol,direction,P,N)){
		if(divideU) *divideU=direction;
		return 0;
	};

	const MGInterval& urng=uvbox[0];
	double u0=urng[0].value(), u1=urng[1].value();
	const MGInterval& vrng=uvbox[1];
	double v0=vrng[0].value(), v1=vrng[1].value();
	double um=(u0+u1)*0.5, vm=(v0+v1)*0.5;
	int ncd=sdim();
	unsigned ku=order_u(), kv=order_v();
	const MGKnotVector& tv=knot_vector_v();
	int j0=tv.locate(v0)-kv+1, j1=tv.locate(v1);
		//Coef necessary for v direction part is from j0 to j1
		//(B-rep dimension is j1-j0+1).
		//Knot vector necessary is from j0 to j1+kv.
	int nunew, nvnew=j1-j0+1;

	tol*=1.1;
	size_t nu, nv, sizeu, sizev;
	surface_bcoef().length(nu, nv); surface_bcoef().capacity(sizeu, sizev);
	unsigned kmax=ku; if(kmax<kv) kmax=kv;
	double* work=new double[kmax*kmax];

	MGSPointSeq surf1(nvnew, nu, ncd); size_t irc=nvnew*nu;
	int nmax=nu; if(nmax<nvnew) nmax=nvnew;
	int nvncd=nvnew; if(nvncd<ncd) nvncd=ncd;
	MGBPointSeq temp(nmax,nvncd);
	MGKnotVector t(kmax, nmax);
	double* tpointer=t.data(); double* temppointer=temp.data();

	for(int k=0; k<ncd; k++){
		double Pk=P[k];
		bluprt_(ku, nu, knot_data_u(), coef_data(0,j0,k),
			sizeu,nvnew,u0,u1,nmax,work,&nunew,tpointer,temppointer,1);
		for(int i=0; i<nunew; i++) for(int j=0; j<nvnew; j++) surf1(j,i,k)=temp(i,j)-Pk;
	}
	//surf1.set_length(nvnew,nunew);cout<<surf1<<endl;//////////

	const double* tvnew=knot_data_v()+j0;
	int nvnew2;
	double x;
	int i,j;
	for(i=0; i<nunew; i++){
		bluprt_(kv,nvnew,tvnew,surf1.data(0,i,0),
			irc,ncd,v0,v1,nmax,work,&nvnew2,tpointer,temppointer,1);
		//temp.set_length(nvnew2);cout<<temp<<endl;///////////
		for(j=0; j<nvnew2; j++){
			x=0.;
			for(int k=0; k<ncd; k++)
				x-=temp(j,k)*N[k];
			if(x<0.) x=-x;
			if(x>tol) break;
		}
		if(j<nvnew2) break;
	}

	int retcode=1;
	if(i<nunew || j<nvnew2){
		if(divideU){
//		*divideU=direction;
		MGVector dfdu=eval(um,vm,1,0), dfdv=eval(um,vm,0,1);
		double ulen=dfdu.len()*(u1-u0), vlen=dfdv.len()*(v1-v0);
		if(ulen*5.<vlen) *divideU=0;
		else if(ulen>vlen*5.) *divideU=1;
		else{

		//Compute maximum deviation from the plane along u direction.
		double udevi0, udevi1, vdevi0, vdevi1;
		size_t jm=nvnew2/2;
		x=0.;
		for(int k=0; k<ncd; k++)
			x-=surf1(jm,0,k)*N[k];
		udevi0=udevi1=x;
		for(i=1; i<nunew; i++){
			x=0.;
			for(int k=0; k<ncd; k++)
				x-=surf1(jm,i,k)*N[k];
			if(x<udevi0)
				udevi0=x;
			if(x>udevi1)
				udevi1=x;
		}
		size_t im=nunew/2;
		bluprt_(kv,nvnew,tvnew,surf1.data(0,im,0),
			irc,ncd,v0,v1,nmax,work,&nvnew2,tpointer,temppointer,1);
		x=0.;
		for(int k=0; k<ncd; k++)
			x-=temp(0,k)*N[k];
		vdevi0=vdevi1=x;
		for(j=1; j<nvnew2; j++){
			x=0.;
			for(int k=0; k<ncd; k++)
				x-=temp(j,k)*N[k];
			if(x<vdevi0)
				vdevi0=x;
			if(x>vdevi1)
				vdevi1=x;
		}
		if((udevi1-udevi0)>=(vdevi1-vdevi0))
			*divideU=1;
		else
			*divideU=0;

		}

		}
		retcode=0;
	}
	delete[] work;

	return retcode;

}
