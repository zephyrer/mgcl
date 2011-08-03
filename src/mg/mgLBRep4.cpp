/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/LBRep.h"

extern "C" {
#include "cskernel/Bldrw1.h"
#include "cskernel/Bldrwc.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// MGLBRep.cpp
//
// Implement MGLBRep class.

//Draw this line's 1st and 2nd coordinates in 2D space
//using drawing function moveto( , ) and lineto( , ).
//wind[] is the window of the screen to draw the line in.
//Clipping will be performed about the wind[].
//(wind[0], wind[1]) is the center coordinates of the window.
//wind[2] is width and wind[3] is hight of the window. When wind[2]<=0,
//no clipping is performed. Even when wind[2]<=0, wind[3] is necessary 
//to input to specify the resolution of the line. In this case,
//wind[0] and wind[1] are not referended.
//ynum is the resolution of the line, is the number of
//straight line segments for the curve length of wind[3](height of window).
//***draw_2D does not perform box including judment, always performs clipping
//operation and draws the line. Users must do obvious box inclusion test
//if maximum drawing performance is necessary.
void MGLBRep::draw_2D(
	void (*moveto)(int, int), void (*lineto)(int, int),
	const double wind[4], //window box to draw the line in.
	size_t ynum) const{//Resolution of the line.
	unsigned n=bdim(), k=order();
	const double* rcoef[2]={m_line_bcoef.data(0,0), m_line_bcoef.data(0,1)};
	const double* t=m_knot_vector.data();
	double* work=new double[4*k*k+3*k+2*n];

	int nrw;
	if(wind[2]>0.){
		double x[2], y[2];
		double xw=wind[2]*0.5, yw=wind[3]*.5;
		x[0]=wind[0]-xw+0.6; x[1]=wind[0]+xw-0.6;
		y[0]=wind[1]-yw+0.6; y[1]=wind[1]+yw-0.6;
		bldrw1_(x,y,0,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
	}
	
	bldrwc_(1,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,work,0,k,n,t,rcoef,work+n);
	delete[] work;
}
void MGLBRep::draw_2D(void (*moveto)(float, float), void (*lineto)(float, float),
	const double wind[4],	//window box to draw the line in.
	size_t ynum) const{		//Resolution of the line.
	unsigned n=bdim(), k=order();
	const double* rcoef[2]={m_line_bcoef.data(0,0), m_line_bcoef.data(0,1)};
	const double* t=m_knot_vector.data();
	double* work=new double[4*k*k+3*k+2*n];

	int nrw;
	if(wind[2]>0.){
		double x[2], y[2];
		double xw=wind[2]*0.5, yw=wind[3]*.5;
		x[0]=wind[0]-xw; x[1]=wind[0]+xw;
		y[0]=wind[1]-yw; y[1]=wind[1]+yw;
		bldrw1_(x,y,0,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
	}
	
	bldrwc_(2,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,work,0,k,n,t,rcoef,work+n);
	delete[] work;
}
void MGLBRep::draw_2D(void (*moveto)(double, double), void (*lineto)(double, double),
	const double wind[4],	//window box to draw the line in.
	size_t ynum) const{		//Resolution of the line.
	unsigned n=bdim(), k=order();
	const double* rcoef[2]={m_line_bcoef.data(0,0), m_line_bcoef.data(0,1)};
	const double* t=m_knot_vector.data();
	double* work=new double[4*k*k+3*k+2*n];

	int nrw;
	if(wind[2]>0.){
		double x[2], y[2];
		double xw=wind[2]*0.5, yw=wind[3]*.5;
		x[0]=wind[0]-xw; x[1]=wind[0]+xw;
		y[0]=wind[1]-yw; y[1]=wind[1]+yw;
		bldrw1_(x,y,0,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
	}
	
	bldrwc_(3,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,work,0,k,n,t,rcoef,work+n);
	delete[] work;
}

//Draw this line's coordinate'th coordinate in 2D space as
//(t, LBRep(coordinate)) when t_is_x is true, 
//or as ( LBRep(coordinate),t)  when t_is_x is false,  
//using drawing function moveto(int, int) and lineto(int,int).
//The other behaviours are the same as draw_2D.
void MGLBRep::draw_1D(void (*moveto)(int, int), void (*lineto)(int, int),
	size_t coordinate,		//id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			//=true:t is x coordinate, and false:t is y.
	const double wind[4],	//window box to draw the line in.
	size_t ynum) const{		//Resolution of the line.
	assert(coordinate<sdim());

	unsigned n=bdim(), k=order();
	const double* rcoef[]={m_line_bcoef.data(0,coordinate)};
	const double* t=m_knot_vector.data();
	double* work=new double[4*k*k+3*k+2*n];

	int nrw;
	int klin=2; if(t_is_x) klin=1;
	if(wind[2]>0.){
		double x[2], y[2];
		double xw=wind[2]*0.5, yw=wind[3]*.5;
		x[0]=wind[0]-xw+0.6; x[1]=wind[0]+xw-0.6;
		y[0]=wind[1]-yw+0.6; y[1]=wind[1]+yw-0.6;
		if(t_is_x) bldrw1_(y,x,klin,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
		else bldrw1_(x,y,klin,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
	}
	bldrwc_(0,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,work,klin,k,n,t,rcoef,work+n);
	delete[] work;
}
void MGLBRep::draw_1D(void (*moveto)(float, float), void (*lineto)(float, float),
	size_t coordinate,		//id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			//=true:t is x coordinate, and false:t is y.
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const{		//Resolution of the line.
	assert(coordinate<sdim());

	unsigned n=bdim(), k=order();
	const double* rcoef[]={m_line_bcoef.data(0,coordinate)};
	const double* t=m_knot_vector.data();
	double* work=new double[4*k*k+3*k+2*n];

	int nrw;
	int klin=2; if(t_is_x) klin=1;
	if(wind[2]>0.){
		double x[2], y[2];
		double xw=wind[2]*0.5, yw=wind[3]*.5;
		x[0]=wind[0]-xw; x[1]=wind[0]+xw;
		y[0]=wind[1]-yw; y[1]=wind[1]+yw;
		if(t_is_x) bldrw1_(y,x,klin,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
		else bldrw1_(x,y,klin,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
	}
	bldrwc_(1,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,work,klin,k,n,t,rcoef,work+n);
	delete[] work;
}
void MGLBRep::draw_1D(void (*moveto)(double, double), void (*lineto)(double, double),
	size_t coordinate,		//id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			//=true:t is x coordinate, and false:t is y.
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const{		//Resolution of the line.
	assert(coordinate<sdim());

	unsigned n=bdim(), k=order();
	const double* rcoef[]={m_line_bcoef.data(0,coordinate)};
	const double* t=m_knot_vector.data();
	double* work=new double[4*k*k+3*k+2*n];

	int nrw;
	int klin=2; if(t_is_x) klin=1;
	if(wind[2]>0.){
		double x[2], y[2];
		double xw=wind[2]*0.5, yw=wind[3]*.5;
		x[0]=wind[0]-xw; x[1]=wind[0]+xw;
		y[0]=wind[1]-yw; y[1]=wind[1]+yw;
		if(t_is_x) bldrw1_(y,x,klin,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
		else bldrw1_(x,y,klin,k,n,t,rcoef,work+2*n,work+n,&nrw,work);
	}
	bldrwc_(2,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,work,klin,k,n,t,rcoef,work+n);
	delete[] work;
}
