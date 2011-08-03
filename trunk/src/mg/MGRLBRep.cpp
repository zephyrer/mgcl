/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Matrix.h"
#include "mg/Transf.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/KnotArray.h"
#include "mg/RLBRep.h"
#include "mg/CCisect.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Implements MGRLBRep Class
//
// Defines Rational Line B-Representation.
// This NURBS is of homogeneous form, i.e., B-Coefficients have
// weight included values. 
// When usual NURBS form is (xi, yi, zi, wi) ,
// MGRLBRep form is (xi*wi, yi*wi, zi*wi, wi) for i=0,..., n-1.

//Construct Line NURBS, providing all the member data.
MGRLBRep::MGRLBRep(
	const MGKnotVector& t,		//Knot Vector.
	const MGBPointSeq& bcoef,	//Line B-Coef, each of coefficients
		//includes weight multiplied when homogeneous=1,
		//and not includes when homogeneous =0.
int homogeneous
):MGCurve(), m_line(t,bcoef){
	assert(t.bdim()==bcoef.length());

	if(!homogeneous){
		double weight;
		size_t dim=bcoef.sdim()-1;
		//Multiply each weight to all control polygons.
		for(size_t i=0;i<t.bdim(); i++){
			weight=bcoef(i,dim);
			for(size_t j=0; j<dim; j++)
				m_line.coef(i,j)=m_line.coef(i,j)*weight;
		}
	}
}

//Construct Line NURBS, providing all the member data.
MGRLBRep::MGRLBRep(
	const MGKnotVector& t,	//Knot Vector.
	const MGBPointSeq& bcoef,//Line B-Coef, each of coefficients does not include weights data.
	const std::vector<double>& weights
):MGCurve(){
	assert(t.bdim()==bcoef.length());

	m_line.m_knot_vector=t;
	size_t bdim=t.bdim(), dim=bcoef.sdim();
	m_line.m_line_bcoef.resize(bdim,dim+1);
	double weight;

	//Multiply each weight to all control polygons.
	for(size_t i=0;i<bdim; i++){
		weight=weights[i];
		for(size_t j=0; j<dim; j++)
			m_line.coef(i,j)=bcoef(i,j)*weight;
		m_line.coef(i,dim)=weight;
	}
}

// Construct ellipse NURBS.
MGRLBRep::MGRLBRep(
	const MGEllipse& ellipse	//Original ellipse.
):MGCurve(){
	size_t sd=ellipse.sdim();
	if(ellipse.ellipse_type()==MGELLIPSE_EMPTY) return;
	else if(ellipse.ellipse_type()==MGELLIPSE_SEGMENT){
		const MGVector& major=ellipse.major_axis(), &minor=ellipse.minor_axis();
		double a=major.len(), b=minor.len();
		double angs=ellipse.param_s(), ange=ellipse.param_e();
		angs=ellipse.gp_to_radian(angs); ange=ellipse.gp_to_radian(ange);
		double plen=ange-angs;	//Total parameter length in radian.

		MGRLBRep rlb;
		if(plen<=0.5*mgPAI)		// Case that parameter length<=0.5*PAI.
			rlb=MGRLBRep(a,b,angs,ange);
		else if(plen<=mgPAI){	
			//				Case that parameter 0.5*PAI<=length<=PAI.
			double mid1=angs+plen*0.5;
			rlb=MGRLBRep
				(MGRLBRep(a,b,angs,mid1),0,2,MGRLBRep(a,b,mid1,ange));
		}else if(plen<=1.5*mgPAI){
			//				Case that parameter PAI<=length<=1.5*PAI.
			double mid1=angs+plen/3.;
			double mid2=angs+plen*2./3.;
			MGRLBRep rlb1
				(MGRLBRep(a,b,angs,mid1),0,2,MGRLBRep(a,b,mid1,mid2));
			rlb=MGRLBRep(rlb1,0,2,MGRLBRep(a,b,mid2,ange));
		}else{
			//				Case that parameter 1.5*PAI<=length<2*PAI.
			double mid1=angs+plen*0.25;
			double mid2=angs+plen*0.5;
			double mid3=angs+plen*0.75;
			MGRLBRep rlb1
				(MGRLBRep(a,b,angs,mid1),0,2,MGRLBRep(a,b,mid1,mid2));
			MGRLBRep rlb2
				(MGRLBRep(a,b,mid2,mid3),0,2,MGRLBRep(a,b,mid3,ange));
			rlb=MGRLBRep(rlb1,0,2,rlb2);
		}
		MGMatrix mat; mat.set_xy_vector(major,minor);
		if(ellipse.sdim()==2) mat=MGMatrix(2,mat);
		MGTransf tr(mat,ellipse.center());
		(*this)=rlb*tr;
	}
	else if(ellipse.ellipse_type()==MGELLIPSE_CLOSED){
	//For the case of whole length ellipse(closed ellipse).
		double onebysqrt2=sqrt(2.)/2.;
		MGVector center=ellipse.center();
		const MGVector& major=ellipse.major_axis();
		const MGVector& minor=ellipse.minor_axis();
		MGVector majorPminor(major+minor); 
		MGVector majorMminor(major-minor);
		MGBPointSeq cp(9,sd+1);
		MGPosition P0=center+major;
		cp.store_at(0,P0);cp(0,sd)=1.;
		cp.store_at(1,(center+majorPminor)*onebysqrt2);
						cp(1,sd)=onebysqrt2;
		cp.store_at(2,center+minor);cp(2,sd)=1.;
		cp.store_at(3,(center-majorMminor)*onebysqrt2);
						cp(3,sd)=onebysqrt2;
		cp.store_at(4,center-major);cp(4,sd)=1.;
		cp.store_at(5,(center-majorPminor)*onebysqrt2);
						cp(5,sd)=onebysqrt2;
		cp.store_at(6,center-minor);cp(6,sd)=1.;
		cp.store_at(7,(center+majorMminor)*onebysqrt2);
						cp(7,sd)=onebysqrt2;
		cp.store_at(8,P0);cp(8,sd)=1.;

		MGKnotVector kv(3,9);
		size_t i;
		for(i=0; i<3; i++){kv(i)=0.; kv(9+i)=1.;}
		for(i=0; i<2; i++){kv(i+3)=.25; kv(i+5)=.5; kv(i+7)=.75;}

		m_line=MGLBRep(kv,cp);
	}

	//Adjust parameter length;
	double tlen=eval(0.,1).len();
	m_line.knot_vector()*=tlen;
}

// Construct a conic section NURBS.
MGRLBRep::MGRLBRep(
	const MGPosition& P0, const MGVector& T0,
							//Start point and its tangent
	const MGPosition& P,	//Mid point of the conic section
	const MGPosition& P2, const MGVector& T2
							//End point and its tangent
):MGCurve(){
	MGPosition P1; double w1;
	if(MGRLBRep_ellipse_weight(P0,T0,P,P2,T2,P1,w1)){
		if(w1<=-1.) return;			//Error.

		size_t sd=P0.sdim(), sd1;
		if((sd1=T0.sdim())>sd) sd=sd1; if((sd1=P.sdim())>sd) sd=sd1;
		if((sd1=P2.sdim())>sd) sd=sd1; if((sd1=T2.sdim())>sd) sd=sd1;
		//sd is maximum space dimension.

		MGKnotVector t(3,3,0.,1.);	//Of order 3, and B-Rep dimension 3.
		MGBPointSeq cp(3,sd+1);
		cp.store_at(0,P0);cp(0,sd)=1.;
		cp.store_at(1,P1);cp(1,sd)=w1;
		cp.store_at(2,P2);cp(2,sd)=1.;
		m_line=MGLBRep(t,cp);
		//This conic may contain infinite control point and 
		//w1 may negative.

		double ca;
		if(w1!=0.){
			MGVector P1w=P1/w1;
			MGVector V10(P0-P1w), V12(P2-P1w);
			ca=V10.cangle(V12);
		}
		if(w1<=0. || ca>=.5) split_conic(1);
		if(w1<0. && ca<=0.){
			split_conic(3); split_conic(1);
		}
	}
}

//Construct 2D ellipse RLBRep, whose center is origin,
//starting point is of angle1, and end point is of angle2.
//The ellipse is expressed as below using parameter t.
// x(t)=a*cos(t),  y(t)=b*sin(t),   angle1<=t<=angle2
MGRLBRep::MGRLBRep(
	double a, double b,
	double angle1, double angle2
):MGCurve(){
	double midp=(angle1+angle2)*.5;
	MGPosition P(a*cos(midp), b*sin(midp)); //Mid point of the ellipse.
	double c1=cos(angle1), s1=sin(angle1);  //Starting point data.
	double c2=cos(angle2), s2=sin(angle2);  //End point data.
	MGPosition P0(a*c1,b*s1), P2(a*c2,b*s2), P1;
	MGVector T0(-a*s1, b*c1), T2(-a*s2, b*c2);

	double w1; MGRLBRep_ellipse_weight(P0, T0, P, P2, T2, P1, w1);

	MGKnotVector kt(3,3,0.,1.); MGBPointSeq cp(3,3);
	cp.store_at(0,P0); cp(0,2)=1.;
	cp.store_at(1,P1); cp(1,2)=w1;
	cp.store_at(2,P2); cp(2,2)=1.;
	
	m_line=MGLBRep(kt,cp);
}

//Function to compute control point P1 and weight w1 of rational form of
//an ellipse segment. Pi and Ti are point and tangent of start and end
//for i=0,2. P is mid point of the ellipse.
//Function's output is if obtained(!=0:true) or not(=0:false).
//When obtained, =1:as finite control point, =2:as infinite.
//***Method***
//See "The NURBS Book" of W.Tiller and L.Piegl publised by Springer.
int MGRLBRep_ellipse_weight
	(const MGPosition& P0, const MGVector& T0,
	 const MGPosition& P,
	 const MGPosition& P2, const MGVector& T2,
	 MGPosition& P1, double& w1)
{
	MGStraight
		S0(MGSTRAIGHT_UNLIMIT,T0,P0), S2(MGSTRAIGHT_UNLIMIT,T2,P2);
//S0 is straight line whose start point is P0 and tangent is T0,
//and so on.
	MGCCisect is;
	MGPSRELATION rel=S0.relation(S2,is); //Compute isect P1 of S0 and S2.

	int ret_code=0;
	//ret_code=0 is the return code
	//when T0, T2, P0, and P2 are not in one plane.
	if(rel==MGPSREL_PARALLEL){
		//When T0 and T2 are parallel.
		MGStraight S(MGSTRAIGHT_UNLIMIT,T0,P);
					//Staright that passes P and whose direction is T0.
		MGStraight S02(P2,P0); S02.relation(S,is);
		double r=S02.param_e(), r1=is.param1();
		double tt=sqrt(r1/(r-r1)); double u=tt/(1.+tt);
		double b=2.*u*(1.-u); b=(1.-u)/b;
		P1=(P-S02.eval(r1))*b;
		w1=0.;
		ret_code=2;
	}else if(rel==MGPSREL_ISECT){
		//Normal case(T0 and T1 are not parallel).
		P1=S0.eval(is.param1()); MGVector P1P=P1-P;
		MGStraight S(MGSTRAIGHT_UNLIMIT,P1P,P);
					//Staright that passes P1 and P.
		MGStraight S02(P2,P0); S02.relation(S,is);
		double r=S02.param_e(), r1=is.param1();
		double tt=sqrt(r1/(r-r1)); double u=tt/(1.+tt);
		double oneu=1.-u;
		w1=oneu/u*((P-P0)%P1P)+u/oneu*((P-P2)%P1P);
		w1/=(2.*(P1P%P1P));
		P1*=w1;
		ret_code=1;
	}
	return ret_code;
}

//**** 3.Conversion Constructor.****

//Approximate an original NURBS by a new knot configuration.
//The new knot config must be inside the range of the original NURBS
//parameter. However new knots may be coarse or fine.
MGRLBRep::MGRLBRep(
	const MGRLBRep& old_brep,//Original NURBS.
	const MGKnotVector& t,	//knot vector
	int &error)			//Error flag.
:MGCurve(old_brep), m_line(old_brep.m_line, t,error){
	update_mark();
}

//Gets new NURBS by adding knots to an original NURBS.
MGRLBRep::MGRLBRep(
	const MGRLBRep& old_brep,	//Original NURBS.
	const MGKnotArray& knots)	//Knots to add.
:MGCurve(old_brep), m_line(old_brep.m_line,knots){
	update_mark();
}

// Gets new NURBS by computing a part of the original. New one is exactly
// the same as the original except that it is partial.
//If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1
//will be guaranteed. Here, n=bdim() and k=order().
//Both t1 and t2 must be inside te range of old_brep.
MGRLBRep::MGRLBRep(
	double t1, double t2, //New parameter range. t1 must be less than t2.
	const MGRLBRep& old_brep,	//Original NURBS.
	int multiple) //Indicates if start and end knot multiplicities
						//are necessary. =0:unnecessary, !=0:necessary.
:MGCurve(old_brep), m_line(t1,t2,old_brep.m_line,multiple){
	update_mark();
}

// Convert from Non ratoinal to Rational form.
MGRLBRep::MGRLBRep(
	const MGLBRep& brep,	//Original LBRep. This can be ordinary LBRep, or 
		//homogeneous form of MGRLBRep. When homogeneous form,
		//the last space dimension elements are weights.
	int homogeneous)		//true(non zero): homogeneous form,
							//false(zero):ordinary SBRep.
:MGCurve(brep), m_line(brep){
	update_mark();
	if(!homogeneous){
		size_t n=brep.bdim(), dimm1=brep.sdim();
		MGBPointSeq cp(n,dimm1+1);
		size_t i,j;
		for(i=0; i<n; i++){
			for(j=0; j<dimm1; j++) cp(i,j)=brep.coef(i,j);
			cp(i,dimm1)=1.;
		}
		m_line=MGLBRep(brep.knot_vector(), cp);
	}
}

//Gets new NURBS by connecting two NURBS to one.
MGRLBRep::MGRLBRep(
	const MGRLBRep& brep1,	//NURBS 1.
	int continuity,			//continuity.
	int which,				//which point of brep1 to which of brep2.
							//meaingfull when continuity>=0.
				// =0: start of this and start of brep1.
				// =1: start of this and end of brep1.
				// =2: end of this and start of brep1.
				// =3: end of this and end of brep1.
				// continuity and which can be obtained using continuity.
	const MGRLBRep& brep2)	//NURBS 2.
:MGCurve(brep1){
	update_mark();
//Change connecting points' weight to 1. for both brep.
	double t1,t2;
	switch(which){
	case 0:	//Start of brep1 and start of brep2.
		t1=brep1.param_s(); t2=brep2.param_s(); break;
	case 1: //Start of brep1 and end of brep2.
		t1=brep1.param_s(); t2=brep2.param_e(); break;
	case 3: //End of brep1 and end of brep2.
		t1=brep1.param_e(); t2=brep2.param_e(); break;
	default: //Otherwise(2) End of brep1 and start of brep2.
		t1=brep1.param_e(); t2=brep2.param_s(); break;
	}

	size_t idw1=brep1.sdim(), idw2=brep2.sdim();
	double weight1=brep1.m_line.eval(t1)(idw1);
	double weight2=brep2.m_line.eval(t2)(idw2);

//Compute temporary LBReps that have the same weights at the
//connecting points.
	MGLBRep line1(brep1.m_line), line2(brep2.m_line);
	line1.line_bcoef()  /=weight1;
	line2.line_bcoef()  /=weight2;
//Connect two.
	m_line=MGLBRep(line1,continuity,which, line2);
}

//Construct a Line NURBS by changing space dimension and ordering of
//coordinates.
MGRLBRep::MGRLBRep(
	size_t dim,			// New space dimension.
	const MGRLBRep& lbrep,	// Original Line B-rep.
	size_t start1, 		// Destination order of new line.
	size_t start2) 		// Source order of original line.
:MGCurve(lbrep){
	update_mark();
	size_t dim0=lbrep.sdim();
	MGBPointSeq cp1(dim0,lbrep.line_bcoef());//Exclude weights.
	MGBPointSeq cp2(dim,cp1,start1,start2);  //Change order of coordinates.
	MGBPointSeq cp3(dim+1,cp2);			     //Get area for weights.
	for(size_t i=0; i<cp3.length(); i++)
		cp3(i,dim)=lbrep.line_bcoef()(i,dim0);//Set weights.
	m_line=MGLBRep(lbrep.knot_vector(),cp3);
}

//Member Function

//Return minimum box that includes the whole line.
MGBox* MGRLBRep::compute_box() const
{	return m_line.line_bcoef().non_homogeneous().compute_box();}
	
// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
//Return minimum box that includes the partial line.
MGBox MGRLBRep::box_limitted(
	const MGInterval& l
	) const
{
	MGInterval range(param_s(), param_e());
	if(range<<l) return box_unlimit();
	MGInterval prange=l&range;
	double t1=prange.low_point(), t2=prange.high_point();
	if(MGREqual_base(t1,t2,knot_vector().param_span())){
		return MGBox(eval(t1),eval(t2));
	}
	MGLBRep lb(t1,t2,m_line);
	return lb.line_bcoef().non_homogeneous().box();
}

//Return minimum box that includes the whole line.
MGBox MGRLBRep::box_unlimit() const{
	return m_line.line_bcoef().non_homogeneous().box();
}

//Changing this object's space dimension.
MGRLBRep& MGRLBRep::change_dimension(
	size_t dim,		// new space dimension
	size_t start1, 		// Destination order of new object.
	size_t start2) 		// Source order of this object.
{
	size_t dim0=sdim();
	MGBPointSeq cp1(dim0,line_bcoef());		//Exclude weights.
	MGBPointSeq cp2(dim,cp1,start1,start2); //Change order of coordinates.
	MGBPointSeq cp3(dim+1,cp2);			    //Get area for weights.
	const MGBPointSeq& bp=line_bcoef();
	for(size_t i=0; i<cp3.length(); i++) cp3(i,dim)=bp(i,dim0);//Set weights.
	m_line=MGLBRep(knot_vector(),cp3);
	update_mark();
	return *this;
}

//Connect brep2 to this brep to make one B-Representation.
//This parameter range will not be changed, instead brep2's range
//will be so changed that brep2 has the same 1st derivative magnitude
//as the original this brep's at the connecting point(start or end point of
//this).
//continuity and which can be obtained using the fucntion continuity().
void MGRLBRep::connect(
	int continuity,	//continuity. must be>=0.
	int which,		//which point of this to which of brep2.
				// =0: start of this and start of brep2.
				// =1: start of this and end of brep2.
				// =2: end of this and start of brep2.
				// =3: end of this and end of brep2.
	const MGRLBRep& brep2)	//B-Rep 2.
{
	size_t sd1=sdim(), sd2=brep2.sdim();
	if(sd1<sd2){
		sd1=sd2;
		change_dimension(sd1);
	}
	MGRLBRep rlb2(sd1,brep2);

//Compute temporary LBReps that have the same weights at the
//connecting points.
	double t1,t2;
	switch(which){
		case 0: t1=param_s(); t2=brep2.param_s();break;
		case 1: t1=param_s(); t2=brep2.param_e();break;
		case 2: t1=param_e(); t2=brep2.param_s();break;
		case 3: t1=param_e(); t2=brep2.param_e();break;
	}

	MGLBRep& line2=rlb2.m_line;
	double weight1=m_line.eval(t1)(sd1), weight2=line2.eval(t2)(sd1);
	double w=weight1/weight2;
	line2.line_bcoef()  *=w;

	m_line.connect(continuity,which,line2);
}

//Compute continuity with brep2.
int MGRLBRep::continuity(
	const MGRLBRep& brep2,
	int& which,	//Indicates which point of this is connected
				// to which of brep2, is meaingfull when continuity>=0.
				// =0: start of this to start of brep2.
				// =1: start of this to end of brep2.
				// =2: end of this to start of brep2.
				// =3: end of this to end of brep2.
	double& ratio // Ratio of 1st derivatives of the two line will
				// be returned when G0 continuity.
				// ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
	) const
// Function's return value is:
// -1: G(-1) continuity, i.e. two lines are discontinuous.
//  0: G0 continuity, i.e. two lines are connected,
//     but tangents are discontinuous
{
	double start1=param_s(), end1=param_e();
	double start2=brep2.param_s(), end2=brep2.param_e();

//Determine which point of brep1 is continuous to which
//point of brep2.
	MGPosition p1[2]={eval(start1), eval(end1)};
	MGPosition p2[2]={brep2.eval(start2), brep2.eval(end2)}; 

	which=2; ratio=1.;	//These values are default.

	double t1,t2;		//Used for the parameters at connecting points.
	int cn=0;
// Check of C-0 continuity.
	if(p1[0]==p2[0]) {
		which=0;		//which=0 means start of brep1 connects
						// to start of brep2.
		t1=start1; t2=start2;
	}
	else if(p1[0]==p2[1]) {
		which=1;		//which=1 means start of brep1 connects
						// to end of brep2.
		t1=start1; t2=end2;
	}
	else if(p1[1]==p2[0]) {
		which=2;		//which=2 means end of brep1 connects
						// to start of brep2.
		t1=end1; t2=start2;
	}
	else if(p1[1]==p2[1]) {
		which=3;		//which=3 means end of brep1 connects
						// to end of brep2.
		t1=end1; t2=end2;
	}
	else return -1;

//Compute temporary LBReps that have the same weights at the
//connecting points.
	size_t idw1=sdim(), idw2=brep2.sdim();
	MGLBRep line1(m_line), line2(brep2.m_line);
	double weight1=line1.eval(t1)(idw1);
	double weight2=line2.eval(t2)(idw2);
	line1.line_bcoef()  /=weight1;
	line2.line_bcoef()  /=weight2;
	int which2; double ratio2;
	int cn2=line1.continuity(line2,which2,ratio2);
	if(cn2>cn) cn=cn2;

// Compute ratio;
	double dlen1=eval(t1, 1).len();
	double dlen2=brep2.eval(t2, 1).len();
	if(!MGMZero(dlen1)) ratio=dlen2/dlen1;
	return cn;
}

// Evaluate right continuous n'th derivative data.
// nderiv=0 means positional data evaluation.
MGVector MGRLBRep::eval(
	double tau,		// Parameter value.
	size_t nderiv,	// Order of Derivative.
	int left		//Left continuous(left=true)
					//or right continuous(left=false).
)const{
	size_t m;
	size_t dim=sdim();
	MGVector result(dim);
	if(nderiv==0){
		MGVector data=m_line.eval(tau,0,left);
		double weight=data[dim];
		for(m=0; m<dim; m++) result(m)=data[m]/weight;
	}else{
		double* deriv=new double[dim*(nderiv+1)];
		eval_all(tau,nderiv,deriv,left);
		result=MGVector(dim,deriv+dim*nderiv);
		delete[] deriv;
	}
	return result;
}

//Compute position, 1st and 2nd derivatives.
// パラメータ値を与えて位置、一次微分値、二次微分値をもとめる。
void MGRLBRep::eval_all(
	double tau,			//Input parameter value(パラメータ値)
	MGPosition& P,		//Position(位置)
	MGVector& V1,		//1st derivative(1次微分値)
	MGVector& V2		//2nd derivative(2次微分値)
)const{
	size_t sd=sdim();
	double* data=new double[sd*3];
	eval_all(tau,2,data);
	P=MGPosition(sd,data);
	V1=MGVector(sd,data+sd);
	V2=MGVector(sd,data+2*sd);
	delete[] data;
}

//Evaluate all of i'th derivative data for 0<=i<=nderiv.
//Output will be put on deriv[j+i*sdim()]
//for 0<=i<=nderiv and 0<=j<sdim(), i.e. 
//deriv[j+i*sdim()] is i-th derivative data for 0<=j<sdim(). 
void MGRLBRep::eval_all(
	double tau,		// Parameter value to evaluate.
	size_t nderiv,	// Order of Derivative.
	double* deriv,	// Output area of size (nderiv+1)*sdim().
	int left		//Left continuous(left=true)
					//or right continuous(left=false).
)const{
	size_t dim=sdim(); size_t dimp1=dim+1;
	size_t ndp1=nderiv+1;

	//Prepare necessary local area bcp_pre, bcp, and data.
	double bc1[4], bc2[4]; double *bcp_pre, *bcp, *bcp_save, *bcp_temp;
	if(ndp1<=4){ bcp_pre=bc1; bcp=bc2;}
	else{ bcp_save=bcp_pre=new double[ndp1*2]; bcp=bcp_pre+ndp1;}
	//Actually bcp and bcp_pre are array of bcp[ndp1], bcp_pre[ndp1],
	//used for binominal coefficients area.
	double data_area[16]; double *data;
	size_t data_len=dimp1*ndp1;
	if(data_len<=16) data=data_area;
	else data=new double[data_len];
	//Actually data is two dimensional array of data[ndp][dimp1],
	//i.e., data[m+dimp1*i]=data[i][m];

	m_line.eval_all(tau,nderiv,data,left);
	// data[m+dimp1*i]=data[i][m] is the data of i-th derivative
	// of m-th space dimension element. For m=dim, wieght.

	double v, weight=data[dim];
	size_t m, i, j;
	for(m=0; m<dim; m++){
		deriv[m]=data[m]/weight;
		for(j=1; j<=nderiv;j++){
			bcp[0]=1.;
			v=data[m+dimp1*j];
			for(i=1; i<j; i++){
				bcp[i]=bcp_pre[i-1]+bcp_pre[i];
				v=v-bcp[i]*data[dim+dimp1*i]*deriv[m+dim*(j-i)];
			}
			v=v-data[dim+dimp1*j]*deriv[m];
			deriv[m+dim*j]=v/weight;
			bcp[j]=1.;
			bcp_temp=bcp_pre; bcp_pre=bcp; bcp=bcp_temp;
		}
	}
	if(data_len>16) delete[] data;
	if(ndp1>4) delete[] bcp_save;
}

// 自身に指定したパラメータ範囲のlimitをつける。
//Get the sub interval line of the original line.
MGRLBRep& MGRLBRep::limit(const MGInterval& itvl){
	update_mark();
	m_line.limit(itvl);
	return *this;
}

//Return non_homogeneous B-Coefficients with weights of
//the rational B-Spline. This MGBPointSeq includes weights.
MGBPointSeq MGRLBRep::non_homogeneous_bcoef()const{
	size_t i,j, sd=sdim(), n=bdim();
	MGBPointSeq cp(n,sd+1);
	for(i=0; i<n; i++){
		double weight=coef(i,sd);
		for(j=0; j<sd; j++) cp(i,j)=coef(i,j)/weight;
		cp(i,sd)=weight;
	}
	return cp;
}

//Test if this is actually non_rational, i.e. , all of the weights are
//same values.
int MGRLBRep::non_rational()const{
	size_t sd=sdim();
	double weight=coef(0,sd);
	for(size_t i=1; i<bdim(); i++)
		if(!MGREqual2(weight,coef(i,sd))) return 0;
	return 1;
}

//Compute part of this curve from parameter t1 to t2.
//Returned is the pointer to newed object, and so should be deleted
//by calling program, or memory leaked.
MGRLBRep* MGRLBRep::part(double t1, double t2, int multiple) const{
	MGRLBRep*rlb=new MGRLBRep(*this);
	rlb->limit(MGInterval(t1,t2));
	return rlb;
}

//Check if the line B-rep is planar.
//Funtion's return value is;
// 0: Not planar, nor a point, nor straight line.
// 1: NURBS is a point.		2: NURBS is a straight line.
// 3: NURBS is planar.
int MGRLBRep::planar(
	MGPlane& plane			//When Brep is not straight line nor a point,
							// plane is returned.
	//Even when not planar(return value is 0), plane nearest is returned.
	, MGStraight& line		//When Brep is a line, line is returned.
	, MGPosition& point		//When Brep is a point, point is returned.
)const{
	return m_line.line_bcoef().non_homogeneous().planar(plane,line,point);
}

//Split conic RLBRep at the i-th control point.
//This RLBRep must be conic section.
MGRLBRep& MGRLBRep::split_conic(size_t i){
	size_t sd=sdim();
	MGVector P0=coef(i-1), P1=coef(i), P2=coef(i+1);
	double w1=P1.ref(sd); double onepw1=1.+w1;
	double weight=sqrt(onepw1*0.5);
	double weight2=weight*2.;
	MGVector Q=(P0+P1)/weight2, R=(P1+P2)/weight2;
	Q(sd)=R(sd)=weight;
	MGVector S=(P0+2.*P1+P2)/onepw1/2.; S(sd)=1.;
	double tmid=(m_line.knot(i+1)+m_line.knot(i+2))*0.5;
	m_line.knot_vector().add_data(MGKnot(tmid,2),3);
	m_line.line_bcoef().store_at(i,R);
	m_line.line_bcoef().insert_at(i,S);
	m_line.line_bcoef().insert_at(i,Q);
	return *this;
}

//Operator overload

//Assignment.
//When the leaf object of this and obj2 are not equal, this assignment
//does nothing.
MGRLBRep& MGRLBRep::operator=(const MGRLBRep& gel2){
	if(this==&gel2)
		return *this;

	MGCurve::operator=(gel2);
	m_line=gel2.m_line;
	return *this;
}
MGRLBRep& MGRLBRep::operator=(const MGGel& gel2){
	const MGRLBRep* gel2_is_this=dynamic_cast<const MGRLBRep*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

// 曲線の平行移動を行いオブジェクトを生成する。
//Translation of the curve.
MGRLBRep MGRLBRep::operator+(const MGVector& v)const{
	MGBPointSeq bc(m_line.line_bcoef());
	bc.homogeneous_transform(v);
	return MGRLBRep(m_line.knot_vector(),bc);
}

//Translation of the curve.
MGRLBRep operator+(const MGVector& v, const MGRLBRep& lb){
	return lb+v;
}

// 与ベクトルだけ曲線を平行移動して自身とする。
//Translation of the curve.
MGRLBRep& MGRLBRep::operator+= (const MGVector& v){
	m_line.line_bcoef().homogeneous_transform(v);
	if(m_box) (*m_box)+=v;
	return *this;
}

// 曲線の逆方向に平行移動を行いオブジェクトを生成する。
//Translation of the curve.
MGRLBRep MGRLBRep::operator- (const MGVector& v) const{
	MGBPointSeq bc(m_line.line_bcoef());
	bc.homogeneous_transform(-v);
	return MGRLBRep(m_line.knot_vector(),bc);
}

// 与ベクトルだけ曲線をマイナス方向に平行移動して自身とする。
//Translation of the curve.
MGRLBRep& MGRLBRep::operator-= (const MGVector& v){
	m_line.line_bcoef().homogeneous_transform(-v);
	if(m_box) (*m_box)-=v;
	return *this;
}

// 与えられたスケールをかけオブジェクトを生成する。
//generate line by scaling.
MGRLBRep MGRLBRep::operator* (double s) const{
	MGBPointSeq bc(m_line.line_bcoef());
	bc.homogeneous_transform(s);
	return MGRLBRep(m_line.knot_vector(),bc);
}

// 与えられたスケールをかけオブジェクトを生成する。
//generate line by scaling.
MGRLBRep operator* (double scale, const MGRLBRep& lb){
	return lb*scale;
}

// 自身の曲線に与えられたスケールをかける。
//Scale the curve.
MGRLBRep& MGRLBRep::operator*= (double s){
	m_line.line_bcoef().homogeneous_transform(s);
	update_mark();
	return *this;
}

// 与えられた変換で曲線の変換を行いオブジェクトを生成する。
//Matrix transformation of the curve.
MGRLBRep MGRLBRep::operator* (const MGMatrix& m) const{
	MGBPointSeq bc(m_line.line_bcoef());
	bc.homogeneous_transform(m);
	return MGRLBRep(m_line.knot_vector(),bc);
}

// 与えられた変換で曲線の変換を行い自身の曲線とする。
//Matrix transformation of the curve.
MGRLBRep& MGRLBRep::operator*= (const MGMatrix& m){
	m_line.line_bcoef().homogeneous_transform(m);
	update_mark();
	return *this;
}

// 与えられた変換で曲線のトランスフォームを行いオブジェクトを生成する。
//General transformation of the curve.
MGRLBRep MGRLBRep::operator* (const MGTransf& t) const{
	MGBPointSeq bc(m_line.line_bcoef());
	bc.homogeneous_transform(t);
	return MGRLBRep(m_line.knot_vector(),bc);
}

// 与えられた変換で曲線のトランスフォームを行い自身とする。
//General transformation of the curve.
MGRLBRep& MGRLBRep::operator*= (const MGTransf& t){
	m_line.line_bcoef().homogeneous_transform(t);
	update_mark();
	return *this;
}

// 論理演算子の多重定義
// 自身とCurveが等しいかどうか比較し判定する。
// 与曲線と自身が等しいかの比較判定を行う。
//Two curves comparison.
bool MGRLBRep::operator==(const MGLBRep& gel2)const{
	if(sdim()!=gel2.sdim())
		return 0;//Check of space dimension.
	if(order()!=gel2.order())
		return 0;	//Check of order.

	size_t bdm=bdim();
	if(bdm!=gel2.bdim())
		return 0;		//Check of B-Rep dimension.
	if(bdm<=0)
		return 1;
	if(!non_rational())
		return 0;		//Check of rationality.
	if(knot_vector() != gel2.knot_vector())
		return 0;//Check of knot vector.

	//Finally, check of control polygon.
	return m_line.line_bcoef().non_homogeneous()==gel2.line_bcoef() ;
}
bool MGRLBRep::operator==(const MGRLBRep& rlb2)const{
	size_t bd1=bdim(), bd2=rlb2.bdim(); if(bd1!=bd2) return 0;
	if(bd1<=0) return 1;
	size_t sd1=sdim(), sd2=rlb2.sdim(); if(sd1!=sd2) return 0;
	if(order()!=rlb2.order()) return 0;	//Check of order.
	if(knot_vector() != rlb2.knot_vector()) return 0;

	double ratio=rlb2.coef(0,sd2)/coef(0,sd1);
	//Check if weights are equal.
	for(size_t i=1; i<bd1; i++)
		if(!MGREqual2(ratio,rlb2.coef(i,sd2)/coef(i,sd1))) return 0;

	return m_line.line_bcoef().non_homogeneous()
			==rlb2.m_line.line_bcoef().non_homogeneous() ;
}
bool MGRLBRep::operator<(const MGRLBRep& gel2)const{
	return m_line<gel2.m_line;
}
bool MGRLBRep::operator==(const MGGel& gel2)const{
	const MGRLBRep* gel2_is_this=dynamic_cast<const MGRLBRep*>(&gel2);
	if(gel2_is_this)
		return operator==(*gel2_is_this);
	return false;
}
bool MGRLBRep::operator<(const MGGel& gel2)const{
	const MGRLBRep* gel2_is_this=dynamic_cast<const MGRLBRep*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}
