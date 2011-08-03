/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/CParam_list.h"
#include "mg/RLBRep.h"
#include "mg/Tolerance.h"

extern "C" {
#include "cskernel/Bldrwcr.h"
}

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
// Implement MGRLBRep class.

typedef int (*func)(...);
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
void MGRLBRep::draw_all2D(
	int kfunc,		//Kind of function move and line,
		//1:move(int,int), 2:move(float, float), otherwise:move(double,double).
	int (*moveto)(...), int (*lineto)(...),
	const double wind[4], //window box to draw the line in.
	size_t ynum)const//Resolution of the line.
{
	unsigned n=bdim(), k=order(), wid=sdim();
	const double* rcoef[3]={
		rcoef[0]=m_line.m_line_bcoef.data(0,0),
		rcoef[1]=m_line.m_line_bcoef.data(0,1),
		rcoef[2]=m_line.m_line_bcoef.data(0,wid)
	};
	const double* knotp=m_line.m_knot_vector.data();
	double* work=new double[k*k+3*k];

	std::vector<double> pvector;
	int nrw;
	if(wind[2]>0.){//If clipping is necessary.
		double xwin[2], ywin[2];
		double xw=wind[2], yw=wind[3];
		double error=MGTolerance::rc_zero();
		double xerror=xw*error, yerror=yw*error;
		xw*=0.5; yw*=0.5;
		xwin[0]=wind[0]-xw+xerror; xwin[1]=wind[0]+xw-xerror;
		ywin[0]=wind[1]-yw+yerror; ywin[1]=wind[1]+yw-yerror;
		if(kfunc==1){ xwin[0]+=0.6; xwin[1]-=0.6; ywin[0]+=0.6; ywin[1]-=0.6;}
			//xwin[] , ywin[] are window coordinates.
		MGVector P=start_point();
		if(xwin[0]<=P[0] && P[0]<=xwin[1] && ywin[0]<=P[1] && P[1]<=ywin[1])
			pvector.push_back(param_s());
		MGCParam_list plist=isect_1D(xwin[0],0);
		MGCParam_list::Citerator i, ie;
		ie=plist.end();
		for(i=plist.begin(); i!=ie; i++) pvector.push_back(*i);
		plist=isect_1D(xwin[1],0); ie=plist.end();
		for(i=plist.begin(); i!=ie; i++) pvector.push_back(*i);
		plist=isect_1D(ywin[0],1); ie=plist.end();
		for(i=plist.begin(); i!=ie; i++) pvector.push_back(*i);
		plist=isect_1D(ywin[1],1); ie=plist.end();
		for(i=plist.begin(); i!=ie; i++) pvector.push_back(*i);
		P=end_point();
		if(xwin[0]<=P[0] && P[0]<=xwin[1] && ywin[0]<=P[1] && P[1]<=ywin[1])
			pvector.push_back(param_e());
		//*** sort the parameter value array.
		std::vector<double>::iterator vi=pvector.begin(), ve=pvector.end();
		std::sort(vi,ve);
		nrw=pvector.size();
	}else{
		pvector=std::vector<double>(2);
		pvector[0]=param_s(); pvector[1]=param_e();
		nrw=2;
	}
	double* isparam=&(pvector[0]);
	bldrwcr_(kfunc,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,isparam,0,
			k,n,knotp,rcoef,work);
	delete[] work;
}

void MGRLBRep::draw_2D(
	void (*moveto)(int, int), void (*lineto)(int, int),
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const		//Resolution of the line.
{	draw_all2D(1,(func)moveto,(func)lineto,wind,ynum);}

void MGRLBRep::draw_2D(
	void (*moveto)(float, float), void (*lineto)(float, float),
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const		//Resolution of the line.
{	draw_all2D(2,(func)moveto,(func)lineto,wind,ynum);}

void MGRLBRep::draw_2D(
	void (*moveto)(double, double), void (*lineto)(double, double),
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const		//Resolution of the line.
{	draw_all2D(3,(func)moveto,(func)lineto,wind,ynum);}

////////////////////////////////////////////////////////////
void MGRLBRep::draw_all1D(
	int coordinate,	//indicates coordinate kind to draw.
	bool t_is_x,	//=true: t is x coordinate, and false:t is y.
	int kfunc,		//Kind of function move and line,
		//1:move(int,int), 2:move(float, float), otherwise:move(double,double).
	int (*moveto)(...), int (*lineto)(...),
	const double wind[4], //window box to draw the line in.
	size_t ynum)const//Resolution of the line.
{
	assert(size_t(coordinate)<sdim());
	unsigned n=bdim(), k=order(), wid=sdim();
	const double* rcoef[2]={
		rcoef[0]=m_line.m_line_bcoef.data(0,coordinate),
		rcoef[1]=m_line.m_line_bcoef.data(0,wid)
	};
	const double* knotp=m_line.m_knot_vector.data();
	double* work=new double[k*k+3*k];

	double ts,te,x,y;
	std::vector<double> pvector;
	int nrw;
	if(wind[2]>0.){//If clipping is necessary.
		double xwin[2], ywin[2];
		double xw=wind[2], yw=wind[3];
		double error=MGTolerance::rc_zero();
		double xerror=xw*error, yerror=yw*error;
		xw*=0.5; yw*=0.5;
		xwin[0]=wind[0]-xw+xerror; xwin[1]=wind[0]+xw-xerror;
		ywin[0]=wind[1]-yw+yerror; ywin[1]=wind[1]+yw-yerror;
		if(kfunc==1){ xwin[0]+=0.6; xwin[1]-=0.6; ywin[0]+=0.6; ywin[1]-=0.6;}
			//xwin[] , ywin[] are window coordinates.

		MGVector P=start_point();
		ts=param_s(); y=P[coordinate]; if(t_is_x) x=ts; else{x=y; y=ts;}
		if(xwin[0]<=x && x<=xwin[1] && ywin[0]<=y && y<=ywin[1])
			pvector.push_back(ts);

		MGCParam_list plist; MGCParam_list::Citerator i, ie;

		if(t_is_x) plist=isect_1D(ywin[0],coordinate);
		else plist=isect_1D(xwin[0],coordinate);
		ie=plist.end();
		for(i=plist.begin(); i!=ie; i++) pvector.push_back(*i);

		if(t_is_x) plist=isect_1D(ywin[1],coordinate);
		else plist=isect_1D(xwin[1],coordinate);
		ie=plist.end();
		for(i=plist.begin(); i!=ie; i++) pvector.push_back(*i);

		te=param_e();
		if(t_is_x){
			if(ts<xwin[0] && xwin[0]<te) pvector.push_back(xwin[0]);
			if(ts<xwin[1] && xwin[1]<te) pvector.push_back(xwin[1]);
		}else{
			if(ts<ywin[0] && ywin[0]<te) pvector.push_back(xwin[0]);
			if(ts<ywin[1] && ywin[1]<te) pvector.push_back(xwin[1]);
		}

		P=end_point();
		y=P[coordinate]; if(t_is_x) x=te; else{x=y; y=te;}
		if(xwin[0]<=x && x<=xwin[1] && ywin[0]<=y && y<=ywin[1])
			pvector.push_back(ts);
		//*** sort the parameter value array.
		std::vector<double>::iterator vi=pvector.begin(), ve=pvector.end();
		std::sort(vi,ve);
		nrw=pvector.size();
	}else{
		pvector=std::vector<double>(2);
		pvector[0]=param_s(); pvector[1]=param_e();
		nrw=2;
	}
	double* isparam=&(pvector[0]);
	int klin; if(t_is_x) klin=1; else klin=2;
	bldrwcr_(kfunc,(S_fp)moveto,(S_fp)lineto,ynum,wind,nrw,isparam,klin,
			k,n,knotp,rcoef,work);
	delete[] work;
}

//Draw this line's coordinate'th coordinate in 2D space as
//(t, LBRep(coordinate)) when t_is_x is true, 
//or as ( LBRep(coordinate),t)  when t_is_x is false,  
//using drawing function moveto(int, int) and lineto(int,int).
//The other behaviours are the same as draw_2D.
void MGRLBRep::draw_1D(
	void (*moveto)(int, int), void (*lineto)(int, int),
	size_t coordinate,		//id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			//=true:t is x coordinate, and false:t is y.
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const		//Resolution of the line.
{	draw_all1D(coordinate,t_is_x,1,(func)moveto,(func)lineto,wind,ynum);}

void MGRLBRep::draw_1D(
	void (*moveto)(float, float), void (*lineto)(float, float),
	size_t coordinate,		//id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			//=true:t is x coordinate, and false:t is y.
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const		//Resolution of the line.
{	draw_all1D(coordinate,t_is_x,2,(func)moveto,(func)lineto,wind,ynum);}

void MGRLBRep::draw_1D(
	void (*moveto)(double, double), void (*lineto)(double, double),
	size_t coordinate,		//id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			//=true:t is x coordinate, and false:t is y.
	const double wind[4],	//window box to draw the line in.
	size_t ynum)const		//Resolution of the line.
{	draw_all1D(coordinate,t_is_x,3,(func)moveto,(func)lineto,wind,ynum);}
