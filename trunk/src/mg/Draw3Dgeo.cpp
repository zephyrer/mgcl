/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Point.h"
#include "mg/Ellipse.h"
#include "mg/Straight.h"
#include "mg/Straight.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/TrimmedCurve.h"
#include "mg/CompositeCurve.h"
#include "mg/SurfCurve.h"
#include "mg/Plane.h"
#include "mg/DrawFunc.h"

extern "C" {
#include "cskernel/bler.h"
#include "cskernel/bpval2.h"
#include "cskernel/bk2fli.h"
#include "cskernel/Blcbpn.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Implements the drawWire functions of all the classes.

//Draw 3D point(vertex) in world coordinates.
//The object is converted to point(s) and is drawn.
//This is valid only for topology objects or MGPoint.
void MGGeometry::draw3DVertex(
) const{;}

//Draw 3D point(vertex) in world coordinates.
//The object is converted to point(s) and is drawn.
//This is valid only for topology objects or MGPoint.
void MGPoint::draw3DVertex(
) const{(*(MGDrawFunc::pfunc()))(ref(0), ref(1), ref(2));}

void MGPoint::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	(*(MGDrawFunc::pfunc()))(ref(0), ref(1), ref(2));
}

void MGEllipse::drawSE(
	double span_length,	//Line segment span length.
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	t0=gp_to_radian(t0); t1=gp_to_radian(t1);
	if(t1<t0){double save=t1; t1=t0; t0=save;}
	double r=major_len();
	double whole_span=t1-t0;
	size_t m=size_t(r*whole_span/span_length);m++;
	double t=t0, dt=whole_span/double(m);
	(*(MGDrawFunc::bfunc()))();

	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
	MGVector V;
	for(size_t i=0; i<m; i++){
		V=eval_in_radian2(t);
		draw_vertex(V[0], V[1], V[2]);
		t+=dt;
	}
	V=eval_in_radian2(t1);
	draw_vertex(V[0], V[1], V[2]);

	(*(MGDrawFunc::efunc()))();
}

//////////////////////////////////////////////

void MGStraight::drawSE(
	double span_length,	//Line segment span length.
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	(*(MGDrawFunc::bfunc()))();
	MGVector V=eval(t0);
	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
	draw_vertex(V[0], V[1], V[2]);
	V=eval(t1);
	draw_vertex(V[0], V[1], V[2]);
	(*(MGDrawFunc::efunc()))();
}

#define INFINITE_LINE_LENGTH 50.
void MGStraight::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
) const{
	double t0=param_s(),t1=param_e();
	if(infinite_above()) t1=INFINITE_LINE_LENGTH;
	if(infinite_below()) t0=-INFINITE_LINE_LENGTH;

	double len=direction_len()*INFINITE_LINE_LENGTH/20.;
	MGUnit_vector X(direction()), Y, Z;
	X.orthonormal(X,Y,Z);
	MGVector Y2(Y*len);
	MGVector P0=eval(t0), P1=eval(t1), V;

	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
	(*(MGDrawFunc::bfunc()))();
	if(infinite_below()){
		V=P0+Y2;
		draw_vertex(V[0], V[1], V[2]);
	}
	draw_vertex(P0[0], P0[1], P0[2]);
	draw_vertex(P1[0], P1[1], P1[2]);
	if(infinite_above()){
		V=P1-Y2-X*len;
		draw_vertex(V[0], V[1], V[2]);
	}
	(*(MGDrawFunc::efunc()))();
}

//////////////////////////////////////////////

void MGLBRep::drawSE(
	double span_length,	//Line segment span length.
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	if(t1<t0){double save=t1; t1=t0; t0=save;}
	size_t k=order(), n=bdim();
	(*(MGDrawFunc::bfunc()))();
	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
	if(k==2){
		const MGKnotVector& t=knot_vector();
		if(t0<t[1]) t0=t[1];
		if(t1>t[n]) t1=t[n];
		MGVector P=eval(t0);
		draw_vertex(P[0],P[1],P[2]);
		size_t i0=t.locate(t0), i1=t.locate(t1);
		const MGBPointSeq& bp=line_bcoef();
		for(size_t i=i0; i<i1; i++)
			draw_vertex(bp(i,0),bp(i,1),bp(i,2));
		P=eval(t1);
		draw_vertex(P[0],P[1],P[2]);
	}else{
		drawgl(span_length, t0, t1);
	}
	(*(MGDrawFunc::efunc()))();
}

//////////////////////////////////////////////
void MGRLBRep::drawSE(
	double span_length,	//Line segment span length.
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	(*(MGDrawFunc::bfunc()))();
	drawgl(span_length, t0, t1);
	(*(MGDrawFunc::efunc()))();
}

//////////////////////////////////////////////
void MGTrimmedCurve::drawSE(
	double span_length,	//Line segment span length.
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	double ts=param_s(), te=param_e();
	if(t0<ts) t0=ts;
	if(t1>te) t1=te;
	m_curve->drawSE(span_length,t0,t1);
}

//////////////////////////////////////////////
void MGCompositeCurve::drawSE(
	double span_length,	//Line segment span length.
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	if(m_composite.size()==0) return;

	double s0,s1;
	if(t0>t1){ s0=t1; s1=t0;}
	else{ s0=t0; s1=t1;}
	size_t i0=find(s0), i1=find(s1);
	if(i0==i1){
		m_composite[i0]->drawSE(span_length,s0,s1);
	}else{
		m_composite[i0]->drawSE(span_length,s0,m_composite[i0]->param_e());
		for(size_t i=i0+1; i<i1; i++)
			m_composite[i]->drawSE(
			span_length,m_composite[i]->param_s(),m_composite[i]->param_e());
		m_composite[i1]->drawSE(
			span_length,m_composite[i1]->param_s(), s1);
	}
}

#define NDIVIDE 3
#define NMAX 150
void MGCurve::drawSE(
	double span_length,	//Line segment span length.
	double t0,			//Start parameter value of the curve.
	double t1			//End parameter value of the curve.
						//Draw will be performed from t0 to t1.
)const{
	size_t i;
	double s0,s1,s[NDIVIDE+1];
	if(t0>t1){ s0=t1; s1=t0;}
	else{ s0=t0; s1=t1;}
	s[0]=range(s0); s[NDIVIDE]=range(s1);
	double ds=(s[NDIVIDE]-s[0])/NDIVIDE;
	for(i=1; i<NDIVIDE; i++) s[i]=s[i-1]+ds;

	(*(MGDrawFunc::bfunc()))();
	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
	MGVector P=eval(s0);
	draw_vertex(P[0], P[1], P[2]);

	MGVector dfdt=eval(s[0]); double vlen1=dfdt.len();
	for(size_t m=0; m<NDIVIDE; m++){
	
	s0=s[m]; s1=s[m+1];
	double s2=(s0+s1)*.5;
	double vlen0=vlen1;
	dfdt=eval(s2,1); double vlen2=dfdt.len();
	dfdt=eval(s1,1); vlen1=dfdt.len();
	double vlen;
	if(vlen0<=vlen1){
		if(vlen1<=vlen2) vlen=vlen0+vlen2;
		else if(vlen0>vlen2) vlen=vlen1+vlen2;
		else vlen=vlen1+vlen0;
	}else{
		if(vlen0<=vlen2) vlen=vlen2+vlen1;
		else if(vlen1>vlen2) vlen=vlen0+vlen2;
		else vlen=vlen0+vlen1;
	}
	vlen*=.5;
	double span=(s1-s0);
	double df=vlen*span;
	size_t n=int(df/span_length); n+=2;
	if(n>NMAX) n=NMAX;
	double dt=span/double(n);
	double t=s0;
	for(size_t i=1; i<n; i++){
		t+=dt;
		P=eval(t);
		draw_vertex(P[0], P[1], P[2]);
	}
	P=eval(s1);
	draw_vertex(P[0], P[1], P[2]);

	}
	(*(MGDrawFunc::efunc()))();
}

//Draw this by converting straight line segments.
void MGLBRep::drawgl(
	double dl,		//approximate line length of the straight line segments.
	double tstart, double tend	//start and end parameter value of this.
)const{
	const size_t c0=0, c1=1;
	const size_t k=order(), n=bdim(), irc=line_bcoef().capacity();
	const double* t=knot_data();
	const double* rcoef=coef_data();

    double x[3]={0.,0.,0.}, dx[3];
	//std::cout<<(*this)<<std::endl;

// ***** INITIAL SET. 
	size_t sd=sdim(); if(sd>3) sd=3;
	double t0=t[k-1], t1=t[n];
	if(tstart>=tend){
		tstart=t0; tend=t1;
	}else{
		if(tstart<t0) tstart=t0;
		if(tend>t1) tend=t1;
	}
	if(t0>=t1)
		return;

// ******* CONVERT TO PP AND DRAW EACH LINE *****
	size_t np1=n+1, is1,is2;
	is1=bk2fli_(np1, t, tstart);
	is2=bk2fli_(np1, t, tend);	while(t[is2-1]==tend) is2--;

	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
// ***** DRAW LINE FROM RW(I) TO RW(I+1) 
//   *** MOVE TO THE FIRST POSITION 
	size_t i;
	for(i=0;i<sd;i++) x[i]=bler_(k,n,t,rcoef+i*irc,tstart,c0);
	draw_vertex(x[0],x[1],x[2]);

//   *** DRAW LINE ONE KNOT SPAN BY ONE
	double *wk1=new double[3*k*k+3*k];
	double* wk1p3k=wk1+3*k;
	for(size_t j=is1; j <= is2; ++j){
		double ts, te;
		double tnow=ts=t[j-1]; te=t[j]; if(ts>=te) continue;
		if(j==is1) ts=tnow=tstart;
		if(j==is2) te=tend;
//     ...CONVERT THE ONE SPAN TO PP-REP
		size_t jmk=j-k; int lpp;
		double brk[2];
		blcbpn_(k,k,t+jmk,rcoef+jmk,irc,sd,c1,wk1p3k,brk,wk1,&lpp);
		double tm=(ts+te)*0.5;
		double vlen=0.;
		for(i=0;i<sd;i++){
			dx[i]=bpval2_(brk,wk1+i*k,k,tm,c1,c1);
			vlen+=dx[i]*dx[i];
		}
	    vlen = sqrt(vlen);//length of 1st deriv at the middle point of the span.
		double span=te-ts;
	    size_t mj = (int)(vlen*span/dl) + 1;
		double dt=span/(double)mj;

		//Draw middle points between the knots BY INCREMENTING DT.
		for(size_t l=1; l<mj; l++){
			tnow+=dt;
			for(i=0;i<sd;i++) x[i]=bpval2_(brk,wk1+i*k,k,tnow,c1,c0);
			draw_vertex(x[0],x[1],x[2]);
		}
		
		//Draw to the end point of the knot.
		for(i=0;i<sd;i++) x[i]=bpval2_(brk,wk1+i*k,k,te,c1,c0);
		draw_vertex(x[0],x[1],x[2]);

	}
	delete[] wk1;
}

//Draw this by converting straight line segments.
void MGRLBRep::drawgl(
	double dl,		//approximate line length of the straight line segments.
	double tstart, double tend	//start and end parameter value of this.
)const{
	const size_t c0=0, c1=1;
	const size_t k=order(), n=bdim(), irc=line_bcoef().capacity();
	const double* t=knot_data();
	const double* rcoef=coef_data();
    double x[3]={0.,0.,0.}, dx[3];

// ***** INITIAL SET.
	double t0=t[k-1], t1=t[n];
	if(tstart>=tend){
		tstart=t0; tend=t1;
	}else{
		if(tstart<t0) tstart=t0;
		if(tend>t1) tend=t1;
	}

// ******* CONVERT TO PP AND DRAW EACH LINE *****
	size_t is1, is2;
	size_t np1=n+1;
	is1=bk2fli_(np1, t, tstart);
	is2=bk2fli_(np1, t, tend);	while(t[is2-1]==tend) is2--;

    size_t i;
	size_t ncd=sdim();
	size_t sd=ncd; if(sd>3) sd=3;
	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
// ***** DRAW LINE FROM RW(I) TO RW(I+1)
//   *** MOVE TO THE FIRST POSITION
	double w=bler_(k,n,t,rcoef+ncd*irc,tstart,c0);
	for(i=0;i<sd;i++) x[i]=bler_(k,n,t,rcoef+i*irc,tstart,c0)/w;
	draw_vertex(x[0],x[1],x[2]);

//   *** DRAW LINE ONE KNOT SPAN BY ONE
	size_t ncdp1=ncd+1;
	double*  pcoef=new double[ncdp1*(k*k+k)];
	double* wpcoef=pcoef+k*ncd;//pcoef for weight of rational b-coef.
	double* scratch=wpcoef+k;//work area.
	for(size_t j=is1; j<=is2; ++j){
		double ts,te;
		double tnow=ts=t[j-1]; te=t[j]; if(ts>=te) continue;
		if(j==is1) ts=tnow=tstart;
		if(j==is2) te=tend;
//     ...CONVERT THE ONE SPAN TO PP-REP
		size_t jmk=j-k; int lpp;
		double brk[2];
	    blcbpn_(k,k,&t[jmk],rcoef+jmk,irc,ncdp1,c1,scratch,brk,pcoef,&lpp);
		double tm=(ts+te)*0.5;
		double vlen=0.;
		for(i=0;i<sd;i++){
			x[i]=bpval2_(brk,pcoef+i*k,k,tm,c1,c0);
			dx[i]=bpval2_(brk,pcoef+i*k,k,tm,c1,c1);
		}
		w=bpval2_(brk,wpcoef,k,tm,c1,c0);
		double dw=bpval2_(brk,wpcoef,k,tm,c1,c1);
		vlen=0.;
		for(i=0;i<sd;i++){
			double vln=(dx[i]-dw*x[i]/w)/w; vlen+=vln*vln;
		}
	    vlen = sqrt(vlen);//length of 1st deriv at the middle point of the span.
		double span=te-ts;
	    size_t mj = size_t(vlen*span/dl) + 1;
		if(mj<2) mj=2;
		double dt=span/double(mj);

		//Draw middle points between the knots BY INCREMENTING DT
		for(size_t l=1; l<mj; l++){
			tnow+=dt;
			w=bpval2_(brk,wpcoef,k,tnow,c1,c0);
			for(i=0;i<sd;i++) x[i]=bpval2_(brk,pcoef+i*k,k,tnow,c1,c0)/w;
			draw_vertex(x[0],x[1],x[2]);
		}

		//Draw to the end point of the knot.
		w=bpval2_(brk,wpcoef,k,te,c1,c0);
		for(i=0;i<sd;i++) x[i]=bpval2_(brk,pcoef+i*k,k,te,c1,c0)/w;
		draw_vertex(x[0],x[1],x[2]);
	}
    delete[] pcoef;
}

void MGPlane::get_uv_display_vector(
	MGVector& u,
	MGVector& v
)const{
	const double L = 10.;
	u = u_deriv(); u *= L / u.len();
	v = v_deriv(); v *= L / v.len();
}

void MGPlane::drawWire(
	double span_length,	//Line segment span length.
	int line_density	//line density to draw a surface in wire mode.
)const{
	struct DrawCross{
		DrawCross(const MGVector& u, const MGVector& v) :
			m_u(u), m_v(v){}

		void operator()(
			const MGVector& start, 
			const MGVector& end,
			bool is_u
			) const{
			MGVector v1, v2;
			if(is_u){
				v1 = start - m_u + m_v;
				v2 = start - m_u - m_v;
			}else{
				v1 = start - m_u - m_v;
				v2 = start + m_u - m_v;
			}
			const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
			const MGDrawFunc::BeginLineStripF begin_line_strip=*(MGDrawFunc::bfunc());
			const MGDrawFunc::EndF end_line_strip=*(MGDrawFunc::efunc());
			begin_line_strip();
				draw_vertex(v1[0], v1[1], v1[2]);
				draw_vertex(start[0], start[1], start[2]);
				draw_vertex(v2[0], v2[1], v2[2]);
			end_line_strip();

			begin_line_strip();
				draw_vertex(start[0], start[1], start[2]);
				draw_vertex(end[0], end[1], end[2]);
			end_line_strip();

			MGVector tmp(end - start);
			v1 += tmp;
			v2 += tmp;

			begin_line_strip();
				draw_vertex(v1[0], v1[1], v1[2]);
				draw_vertex(end[0], end[1], end[2]);
				draw_vertex(v2[0], v2[1], v2[2]);
			end_line_strip();
		}

		const MGVector& m_u;
		const MGVector& m_v;
	};

	MGVector u,v;
	get_uv_display_vector(u,v);

	const MGDrawFunc::VertexF draw_vertex=*(MGDrawFunc::vfunc());
	const MGPosition& cen = center();
	(*(MGDrawFunc::bfunc()))();
		MGVector V(cen + u + v);  // (L, L).
		draw_vertex(V[0], V[1], V[2]);
		V -= 2 * u;
		draw_vertex(V[0], V[1], V[2]);  // (-L, L)
		V -= 2 * v;
		draw_vertex(V[0], V[1], V[2]);  // (-L, -L)
		V += 2 * u;
		draw_vertex(V[0], V[1], V[2]);  // (L, -L)
		V += 2 * v;
		draw_vertex(V[0], V[1], V[2]);  // (L, L)
	(*(MGDrawFunc::efunc()))();

	const double ratio = 4.11;
	u /= ratio;
	v /= ratio;
	DrawCross da(u, v);
	da(cen-3*u, cen+3*u, true);
	da(cen-3*v, cen+3*v, false);
}

MGDrawFunc::BeginLineStripF MGDrawFunc::m_bfunc=0;
MGDrawFunc::VertexF MGDrawFunc::m_vfunc=0;
MGDrawFunc::EndF MGDrawFunc::m_efunc=0;
MGDrawFunc::DrawPoint MGDrawFunc::m_pfunc=0;
