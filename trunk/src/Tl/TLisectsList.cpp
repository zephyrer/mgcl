/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position_list.h"
#include "mg/CParam_list.h"
#include "mg/TrimmedCurve.h"
#include "mg/Tolerance.h"
#include "topo/Edge.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "Tl/TLparameter.h"
#include "Tl/TLisects.h"
#include "Tl/TLisectsList.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//mgTLisectsList is a proprietry class for Face tessellation.
//mgTLisectsList holds all the necessary mgTLisects along u or v 
//in the parameter space of a face.
//This class is to eliminate the same value's isect_1D computation.
//Once isect_1D is invoked for some u or v value, the intersections are
//hold in this class.

//////////// Member function ///////////

//Search mgTLisects whose t()==tin. If not found compute intersections
//and add one member of mgTLisects whose t() is tin.
//Function's return value is the iterator of the mgTLisects in mgTLisectsList.
mgTLisectsList::issListItr mgTLisectsList::get(
	mgTLparameter& param,	//tessellation parameter.
	double tin,		//face's parameter u or v according to kcod.
	size_t kcod	//Coordinate kind, =0: tin is u value, =1: tin is v value.
){
	const MGFace& face=param.get_face();
	std::vector<bool>* edgPerim=param.get_edgPerim();

	//m_isectsList includes at least two members(since initialized in mgTLparameter()).
	issListItr is,i=m_isectsList.begin(),ie=m_isectsList.end();
	if(tin<=middle()){
		for(; i!=ie; i++){
			double x=(*i).t();
			if(x==tin) return i;
			else if(x>tin) break;
		}
	}else{
		is=i;
		do{
			i=ie--;
			double x=(*ie).t();
			if(x==tin) return ie;
			else if(x<tin) break;
		}while(ie!=is);
	}

	const MGSurface& sf=param.get_surface();
	//Now tin's intersections not found, and i indicates the place to insert
	//(the position just after the place to insert), generate mgTLisects for tin.
	i=m_isectsList.insert(i,mgTLisects(tin));
	mgTLisects& isects=*i;
	size_t n=face.number_of_loops();
	for(size_t m=0; m<n; m++){
		const MGLoop& lpi=*(face.loop(m));
		mgTLisect1D(param, lpi, edgPerim[m], tin, kcod, isects);
	}
	std::sort(isects.begin(), isects.end());
//	cout<<isects<<endl;
	return i;
}

//Compute intersection points of 1D sub curves of the original loop.
//Parameter values of intersection points(MGLEPoint's) will be returned.
//This is for tessellation and intersection with perimeter boudary edges will
//be excluded.
void mgTLisect1D(
	mgTLparameter& param,
	const MGLoop& lp,	//Target Loop.
	const std::vector<bool>& edgPerimj,//lp's edgPerim[j].
		//if lp is an inner boundary, this is not used.
		//Used only when lp is outer boundary or perimeter boundary.
	double f,		//Coordinate value
	size_t kcod,	//Coordinate kind of the data f, u(kcod=0) or v(=1).
	mgTLisects& iss	//mgTLisect vector will be output.
					//Must be initialized before use.
){
	const MGBox& lpbx=lp.box();
	if(lpbx[kcod]<<f) return;

	const double* prng=param.suf_param_range();
	bool inner=lp.is_inner_boundary();
	double err;
	if(kcod) err=param.isect_verror();
	else     err=param.isect_uerror();
	double errsave=MGTolerance::set_wc_zero(err);
	double tol=err*3.;
	double error2=MGTolerance::rc_zero_sqr();
	MGPosition_list Plist;
		//This is used to avoid the answers of multiple points at a connection.
	MGPosition_list::iterator itr;
	MGComplex::const_pcellItr ei=lp.pcell_begin(), ee=lp.pcell_end();
	size_t id=1; if(kcod) id=0;
	size_t id2=id*2;
	size_t id2p1=id2+1;
	for(size_t i=0; ei!=ee; ei++, i++){
		if(!inner) if(edgPerimj[i]) continue;
		const MGEdge& edg=*(edge_from_iterator(ei));
		const MGInterval& frng=(edg.box())[kcod];
		if((f+tol)<frng[0] || frng[1]<f-tol) continue;
		MGTrimmedCurve tcrv=edg.trimmed_curve();
		MGCParam_list tlist=tcrv.isect_1D(f, kcod);
		MGCParam_list::Citerator ti=tlist.begin(), te=tlist.end();
		for(;ti!=te; ti++){
			MGPosition Q=tcrv.eval(*ti);
			if(!Plist.in(MGBox(Q,tol),itr)){
				Plist.push_back(Q);
				MGVector tangen=tcrv.eval(*ti,1); double len2=tangen%tangen;
				double dt=tangen[kcod]; double dt2=dt*dt;
				int increase=1;
				if(dt2<=error2*len2) increase=0;//Flag of unknown of increase or not.
				else if(dt<0.) increase=-1;		//Flag of decreasing.
				double t=Q[id];
				if(t<prng[id2]) t=prng[id2];
				if(t>prng[id2p1]) t=prng[id2p1];
				iss.m_isects.push_back(mgTLisect(t, MGLEPoint(ei,*ti), increase));
			}
		}
	}
	MGTolerance::set_wc_zero(errsave);
	return;
}

ostream& operator<< (ostream& out, const mgTLisectsList& issl){
	out<<"TLisectsList::m_middle="<<issl.m_middle<<endl;
	mgTLisectsList::CissListItr i=issl.m_isectsList.begin(),
		ie=issl.m_isectsList.end();
	for(; i!=ie; i++) out<<(*i);
	return out;
}
