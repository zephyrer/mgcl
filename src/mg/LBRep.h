/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGLBRep_HH_
#define _MGLBRep_HH_

#include "mg/MGCL.h"
#include "mg/Vector.h"
#include "mg/KnotVector.h"
#include "mg/BPointSeq.h"
#include "mg/Curve.h"

// MGLBRep.h
// Forward Declaration
class  MGInterval;
class  MGNDDArray;
class  MGPosition;
class  MGOscuCircle;
class  MGKnotArray;
class  MGCParam_list;
class  MGPosition_list;
class  MGLBRepEndC;
class  MGRLBRep;
class  MGPPRep;
class  MGSBRepTP;
class  MGIfstream;
class  MGOfstream;

// Defines Line B-Representation.

/** @addtogroup GEO
 *  @{
 */

///MGLBRep is a class for B-SPline representation.
///For a general NURBS(non uniform rational B-Spline) is MGRLBRep.
///LBRep abbrebiates Line B-Representation. MGLBRep consists of a knot vector(MGKnotVector)
///and a control polygon(MGBPointSeq) whose B-representaiton dimension are the same.
class MGCLASS MGLBRep: public MGCurve{

public:

///translation by a vector.
MGDECL friend MGLBRep operator+ (const MGVector& v, const MGLBRep& lb);

///scaling by a scalar.
MGDECL friend MGLBRep operator* (double scale, const MGLBRep&);

///曲線列のノットベクトルを再構築する
///入力された複数曲線列を指定オーダーで再構築する。トレランスはline_zero()を使用している。
///オーダーが指定されていないとき曲線列のうちで最も大きいオーダーを使用する。このとき、
///Ellipse, Straightのオーダーは4として考える。
///パラメータ範囲は1次微分値の大きさが１になるようにしたときの長さの平均を使用している。
///戻り値は再構築後の曲線列が返却される。エラーのときヌルが返却される。
MGDECL friend MGPvector<MGLBRep> rebuild_knot(
	const std::vector<const MGCurve*>& brepl,///<入力曲線列
	size_t order = 0,	///<指定オーダー
	MGLBRep**	tp=0	///<接続面
);

///Same as above, except that the input is MGPvector<>.
MGDECL friend MGPvector<MGLBRep> rebuild_knot(
	const MGPvector<MGCurve>& brepl,	///<入力曲線列
	size_t order = 0,	///<指定オーダー
	MGLBRep**	tp=0	///<接続面
);				

///複数カーブの共通で削除できるノットを削除する。
///ただし、入力カーブは同じノットベクトルを持つものとする。
MGDECL friend void remove_knot_curves(
	MGPvector<MGLBRep>& brepList,	///<曲線列
	MGLBRep**		tp=0,	///<接続面	input and output.
		///<if tp[i] for crvl[i] was not null, converted new tp will be output.
	double tp_length=0.
);

///////<< Constructor >>////////

///Default(dummy) constructor.
MGLBRep():MGCurve(),m_line_bcoef(),m_knot_vector(){;}

///Dummy constructor that specifies area length.
MGLBRep(size_t bdim,	///< b-rep dimension of the lbrep.
		size_t order,	///<order of the lbrep.
		size_t sdim		///<space dimension of the lbrep
):MGCurve(),m_line_bcoef(bdim,sdim),m_knot_vector(order, bdim){;}

///Construct Line B-Representation, providing all the member data.
///***** This is the fundamental constructor.*****
MGLBRep(
	const MGKnotVector& t,	///<Knot Vector.
	const MGBPointSeq& bcoef///<Line B-Coef.
);

// **** 1. Interpolation Constructor ****

/// Construct Line B-rep by intepolation from Point data only.
///If circular is true, start and end points are connected smoothly.
///If circular is true, order will be always 4, and input order is neglected.
MGLBRep(
	const MGBPointSeq& points,	///<Point seq data
	int& error,				///<Error flag.
	unsigned order=4,		///<Order
	int circular=0		///<Circular flag
);

/// Construct Line B-rep of a specified order, given data point abscissa and
///the ordinates.
MGLBRep(
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& points,	///<Point seq data
	size_t order=4,				///<order
	double ratio=-1.			///<Maximum of data point ratio of pre and after spans.
			///< Let d(i)=tau[i]-tau[i-1], then if d(i)/d(i-1)>ratio or
			///< d(i-1)/d(i)>ratio, either tau[i] or tau[i-1] will be removed.
			///< This is done to prevent control polygon computation error.
			///< When ratio<0. no data point removal will be done.
);

/// Construct Line B-rep of any order number by interpolation 
///from Point data only with knot vector.
MGLBRep(
	const MGNDDArray& tau,		///<Data point abscissa.
	const MGBPointSeq& points,	///<Point seq data(data point ordinate).
	const MGKnotVector& t,		///<knot vector.
	int &error				///<Error flag.
);

///Construct Line B-rep of order 4 by interpolation from Point data
///and end condition. (tau(i), value(i,.)) for 0<=i<=n(the length of value).
///For the start and end point, tau does not have multiplicity. However,
///if tau has multiplicity at inner point, this means 1st derivative data is
///provided for the associated value, i.e.:
///If tau has multiplicity 2 as tau(i)=tau(i+1), value(i,.) is 1st derivative
///at tau(i) and value(i+1,.) is positional data at tau(i)(=tau(i+1)).
///If tau has multiplicity 3 as tau(i)=tau(i+1)=tau(i+2),
///value(i,.) is 1st derivative at tau(i)- ,
///value(i+1,.) is positional data at tau(i)(=tau(i+1)), 
///value(i+2,.) is 1st derivative at tau(i)+.
///Maximum multiplicity allowed is 3.
MGLBRep(
	const MGLBRepEndC& begin,	///<Begin end condition
	const MGLBRepEndC& end,		///<End end conditoion
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& value,	///<Data point ordinate
	int &error				///<Error flag.
);

///Construct Line B-rep of any order by interpolation from Point data
///and end condition. (tau(i), value(i,.)) for 0<=i<=n(the length of value).
///For the start and end point, tau does not have multiplicity. However,
///if tau has multiplicity at inner point, this means 1st derivative data is
///provided for the associated value, i.e.:
///If tau has multiplicity 2 as tau(i)=tau(i+1), value(i,.) is 1st derivative
///at tau(i) and value(i+1,.) is positional data at tau(i)(=tau(i+1)).
///If tau has multiplicity 3 as tau(i)=tau(i+1)=tau(i+2),
///value(i,.) is 1st derivative at tau(i)- ,
///value(i+1,.) is positional data at tau(i)(=tau(i+1)), 
///value(i+2,.) is 1st derivative at tau(i)+.
///Maximum multiplicity allowed is 3.
MGLBRep(
	size_t order,				///<Order of the MGLBRep
	const MGLBRepEndC& begin,	///<Begin end condition
	const MGLBRepEndC& end,		///<End end conditoion
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& value,	///<Data point ordinate
	int &error				///<Error flag.
);

///Construct Line B-rep of any order by interpolation from Point data
///with end condition and the knot vector for the B-rep to construct.
/// (tau(i), value(i,.)) for 0<=i<=n(the length of value).
///tau(i) and knot vector t must satisfy Shoenberg's variation diminishing
///constraint.
///For the start and end point, tau does not have multiplicity. However,
///if tau has multiplicity at inner point, this means 1st derivative data is
///provided for the associated value, i.e.:
///If tau has multiplicity 2 as tau(i)=tau(i+1), value(i,.) is 1st derivative
///at tau(i) and value(i+1,.) is positional data at tau(i)(=tau(i+1)).
///If tau has multiplicity 3 as tau(i)=tau(i+1)=tau(i+2),
///value(i,.) is 1st derivative at tau(i)- ,
///value(i+1,.) is positional data at tau(i)(=tau(i+1)), 
///value(i+2,.) is 1st derivative at tau(i)+.
///Maximum multiplicity allowed is 3.
MGLBRep(
	const MGLBRepEndC& begin,	///<Begin end condition
	const MGLBRepEndC& end,		///<End end conditoion
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& value,	///<Data point ordinate
	const MGKnotVector& t,		///<knot vector.
	int &error					///<Error flag.
);

/// Construct Line B-rep of order 4 from point and point-kind followed by
///osculating circle data.
///point_kind[i] is point kind of the point points(i,.):
/// =0:G2 point, =1:G0 point, =2:G1 point.
///If two consecutive points are 1 or 2,
///the span is a straight line. point_kind 2 is a start of G2 curve.
///If two straight line span meet at points(i), osculating circle can be generated
///at this point by providing circle data at this point.
MGLBRep(
	const MGLBRepEndC& begin,	///<Begin end condition
	const MGLBRepEndC& end,		///<End end conditoion
	const MGBPointSeq& points,	///<Point seq data
	const int* point_kind,		///<Point kind of above point.
	const MGOscuCircle& circle,	///<Provides osculating circle data.
	int &error					///<Error flag.
);

// **** 2. Approximation Constructor ****

///Construct Curve B-Rep.
///This is an approximation, and the tolerance is MGTolerance::line_zero().
explicit MGLBRep(
	const MGCurve& crv,	///<Original Curve.
	size_t order=0///<Order. When order=0 is input, and crv was a MGLBRep,
		///<the original order will be used. Otherwise(order=0 and crv was not an MGLBRep)
		///<order will be set to 4.
);

///Approximate an original B-Rep by a new knot configuration.
///The new knot config must be inside the range of the original B-Rep
///parameter. However new knots may be coarse or fine.
MGLBRep(
	const MGLBRep& old_brep,///<Original B-Rep.
	const MGKnotVector& t,	///<knot vector
	int &error			///<Error flag.
);

///Gets new B-Rep by a new knots.
///The parameter range of t must be inside the one of old_curve.
///The constructed MGLBRep's knot vector is t, and the MGLBRep is the approximation
///of old_curve's parameter range from t.param_s() to t.param_e();
MGLBRep(
	const MGCurve& old_curve,	///<Original curve.
	const MGKnotVector& t	///<knot vector
);

/// Construct 3D B-Rep by mixing two 2D B-Rep.
///The two 2D B-Rep's directions and start and end points must be the same.
///Second 2D B-Rep can be girth representaion, ie,
///Let brep1 is f(t)=(f1(t),f2(t)) and brep2 is g(s)=(g1(s),g2(s)), where
///f1,f2 are two coordinates, g1 is parameter t of f(t), 
///g2(s) is the missing coordinate of f(t). Given parameter s,
/// ( f1(g1(s)), f2(g1(s)), g2(s)) is a 3D space point.
MGLBRep(
	unsigned coordinate1,	///<Missing oordinate kind of the brep1 
							///< 0:x, 1:y, 2:z.
	const MGLBRep& brep1,	///<Original 2D B-Rep1. Coordinates are
		///<(y,z), (z,x), (x,y) according to coordinate1.
	unsigned coordinate2,	///<Missing coordinate kind of the brep2.
							///< 0:x, 1:y, 2:z, and 3:girth rep.
	const MGLBRep& brep2	///<Original 2D B-Rep2. Coordinates are
		///<(y,z), (z,x), (x,y) and (t, g2) according to coordinate2.
		///<t is parameter of brep1 and g2 is x, y, or z according to coordinate1.
);

/// Construct Line B-rep of any order number by least square approximation
///from Point data with approximation weights and knot vector of B-Rep.
///weight[i] is for points(i,.). For detail information of the approximation
///method, see "A Practical Guide to Splines by Carl de Boor" Springer-Verlag.
MGLBRep(
	const MGNDDArray& tau,		///<Data point abscissa
	const MGBPointSeq& points,	///<Point seq data
	const double* weight,		///<Weights for each points 
	const MGKnotVector& t		///<knot vector
);

// **** 3.Conversion Constructor.****

///Convert PP-Rep to B-rep.
MGLBRep(const MGPPRep& pprep);

///This constructor constructs B-Rep, converting from PP-Rep.
///Knot Vector is input. Each knot of the knot vector is break point of
///pprep. The continuities at all the break points must be C(k-2) where
///k is the order of pprep.
MGLBRep(
	const MGPPRep& pprep,	///<PP-rep
	const MGKnotVector t	///<Knot Vector
);

///Gets new B-Rep by adding knots to an original B-Rep.
MGLBRep(
	const MGLBRep& old_brep,	///<Original B-Rep.
	const MGKnotArray& knots	///<Knots to add.
);

///Construct LBRep by connecting brep1 and brep2 to make one B-Representation.
///brep1's parameter range will not be changed, instead brep2's range
///will be so modified that brep2 has the same 1st derivative magnitude
///as the original brep1's at the connecting point
///(start or end point of brep1).
///continuity and which can be obtained using the fucntion continuity().
MGLBRep(
	const MGLBRep& brep1,	///<B-Rep 1.
	int continuity,			///<continuity. must be>=0.
	int which,				///<which point of brep1 to which of brep2.
							///<meaingfull when continuity>=0.
				///< =0: start of brep1 and start of brep2.
				///< =1: start of brep1 and end of brep2.
				///< =2: end of brep1 and start of brep2.
				///< =3: end of brep1 and end of brep2.
	const MGLBRep& brep2	///<B-Rep 2.
);

///Gets new B-Rep by computing a part of the original. New one is exactly
///the same as the original except that it is partial.
///If multiple==true(!=0), knot(i)=t1 and knot(n+i)=t2 for i=0,..., k-1
///will be guaranteed. Here, n=bdim() and k=order().
///Both t1 and t2 must be inside te range of old_brep.
MGLBRep(
	double t1, double t2, ///<New parameter range. t1 must be less than t2.
	const MGLBRep& old_brep,	///<Original B-Rep.
	int multiple=0///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

///Construct a Line B-Rep by changing space dimension and ordering of
///coordinates.
MGLBRep(
	size_t dim,				///< New space dimension.
	const MGLBRep& lbrep,	///< Original Line B-rep.
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of original line.
);

//Copy constructor.
//MGLBRep(const MGLBRep& lb2); ///We can use default copy constructor.

//Destructor
//	~MGLBRep();	///We can use default destructor.

////////Operator overload////////

//Assignment.
//MGLBRep& operator= (const MGLBRep& lb2);//We can use default assignment.

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGLBRep& operator=(const MGGel& gel2);

///Assignment.
MGLBRep& operator=(const MGLBRep& el2);

///Transformation object construction
MGLBRep operator+ (const MGVector& ) const;
MGLBRep operator- (const MGVector& ) const;
MGLBRep operator* (double) const;
MGLBRep operator* (const MGMatrix& ) const;
MGLBRep operator* (const MGTransf& ) const;

///Object transformation.
MGLBRep& operator+=(const MGVector& v);
MGLBRep& operator-=(const MGVector& v);
MGLBRep& operator*=(double scale);
MGLBRep& operator*=(const MGMatrix& mat);
MGLBRep& operator*=(const MGTransf& tr);

///comparison
bool operator==(const MGLBRep& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGLBRep& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator==(const MGRLBRep& gel2)const;

/////////Member Function////////

///Returns B-Rep Dimension.
size_t bdim() const{return m_line_bcoef.length();}
	
/// 入力のパラメータ範囲の曲線部分を囲むボックスを返す。
///Return minimum box that includes the partial line.
MGBox box_limitted(const MGInterval& l) const;

///Return minimum box that includes the whole line.
const MGBox& box_unlimit() const;

///Build line B-Rep by Schoenberg and Reinsch smoothing function, supposing the end
///conditions are free end conditions, given
///data points (tau,y), weights dy at data points, and a deviation.
///If dy[i] gets larger, deviation at tau(i) gets larger.
///n can be any number greater than or equal to 2.
///***End conditions are free end condition.***
void buildSRSmoothedLB_of_FreeEnd(
	const MGNDDArray& tau,	///<Data point abscissa
	const MGBPointSeq& y,	///<Data point ordinates.
	const double* dy,///<dy[i] is the weights  at tau[i] for i=0,..., tau.length()-1.
	double deviation,///<if dev_is_sum is true,
		///<deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		///<if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		///<dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		///<dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		///<for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum=false
);

///Build line B-Rep by Schoenberg and Reinsch smoothing function, given
///1st derivatives on the start and end points, data points (tau,y),
///weights dy at data points, and a mean deviation deviation.
///If dy[i] gets larger, deviation at tau(i) gets larger.
///n can be any number greater than or equal to 2.
void buildSRSmoothedLB_of_1stDeriv(
	const MGLBRepEndC& begin,///<Begin end condition
	const MGLBRepEndC& end,	///<End end conditoion.
		///<begin.cond() and end.cond() must be MGENDC_1D or MGENDC_12D.
	const MGNDDArray& tau,	///<Data point abscissa
	const MGBPointSeq& y,	///<Data point ordinates.
	const double* dy,///< dy[i] is the weight  at tau[i] for i=0,..., tau.length()-1.
	double deviation,///<if dev_is_sum is true,
		///<deviation is the upper bound of Sum(((points(i)-pout(i))/dp[i])**2.
		///<if dev_is_sum is false, deviation is max_deviation of each point at tau[i],i.e.,
		///<dev_is_sum=true: deviation>=Sum(((points(i)-pout(i))/dp[i])**2),
		///<dev_is_sum=false:deviation>=Max((points(i)-pout(i))**2),
		///<for i=0,...,n-1. Here pout(i) is the this->eval(tau(i)).
	bool dev_is_sum=false
);

///Changing this object's space dimension.
MGLBRep& change_dimension(
	size_t sdim,		///< new space dimension
	size_t start1=0,	///< Destination order of new object.
	size_t start2=0 	///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
);

///Change order of the B-Rep. When new order is greater than the original,
///new B-rep is guaranteed to be the same line as the original. However,
///if new order is less than the original one, new line is not the same
///in general.
MGLBRep& change_order(
	unsigned order	///<New order number. 
);

///Change order of the B-Rep by approximation.
MGLBRep& change_order_by_approximation(
	unsigned ordr	///<New order number. 
);

///Access to (i,j)th element of coef
///( left-hand side version)
double& coef(size_t i, size_t j){return m_line_bcoef(i,j);}

///Access to (i,j)th element of coef
///(right hand side version)
double coef(size_t i, size_t j) const{return m_line_bcoef(i,j);}

///Extract (i,j)element for 0<=j<sdim().
MGVector coef(size_t i) const{return m_line_bcoef(i);}

///Returns a pointer to the line b-coef data.
const double* coef_data(size_t i=0, size_t j=0)const{return m_line_bcoef.data(i,j);}
double* coef_data(size_t i=0, size_t j=0){return m_line_bcoef.data(i,j);}

///Connect brep2 to this brep to make one B-Representation.
///This parameter range will not be changed, instead brep2's range
///will be so changed that brep2 has the same 1st derivative magnitude
///as the original this brep's at the connecting point(start or end point of
///this).
///continuity and which can be obtained using the fucntion continuity().
void connect(
	int continuity,	///<continuity. must be>=0.
	int which,	///<which point of this to which of brep2.
				///< =0: start of this and start of brep2.
				///< =1: start of this and end of brep2.
				///< =2: end of this and start of brep2.
				///< =3: end of this and end of brep2.
	const MGLBRep& brep2	///<B-Rep 2.
);

///Compute continuity with brep2.
/// Function's return value is:
/// -1: G(-1) continuity, i.e. two lines are discontinuous.
///  0: G0 continuity, i.e. two lines are connected,
///     but tangents are discontinuous
///  1: G1 continuity, i.e. two lines are connected,
///     and tangents are also continuous.
///  2: G2 continuity, i.e. two lines are connected,
///     and tangents and curvatures are also continuous.
int continuity(
	const MGLBRep& brep2,
	int& which,	///<Indicates which point of this is connected
				///< to which of brep2, is meaingfull when continuity>=0.
				///< =0: start of this to start of brep2.
				///< =1: start of this to end of brep2.
				///< =2: end of this to start of brep2.
				///< =3: end of this to end of brep2.
	double& ratio///< Ratio of 1st derivatives of the two line will be returned.
			///< ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
) const;

///Exchange ordering of the coordinates.
///Exchange coordinates (j1) and (j2).
MGLBRep& coordinate_exchange(size_t j1, size_t j2);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGLBRep* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const{return clone();};

///copy as a newed curve. The new curve will be MGLBRep.
///Returned object must be deleted.
MGLBRep* copy_as_LBRep() const;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGLBRep* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

///Construct new curve object by copying to newed area,
///and limitting the parameter range to prange.
///Returned is newed object and must be deleted.
MGCurve* copy_limitted(const MGInterval& prange) const;

///Compute curvilinear integral of the 1st two coordinates.
///This integral can be used to compute area sorounded by the curve.
///(線積分）を求める。
double curvilinear_integral(double t1, double t2) const;
#ifdef __sgi
	double curvilinear_integral()const
	{return curvilinear_integral(param_s(), param_e());}
#endif

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	MGPvector<MGCurve>& crv_list,	///<divided curves will be appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at.
						///<When multiplicity<=0, order()-1 is assumed.
						///<When multiplicity>=order(), order() is assumed.
) const;

///Display control polygons using mgGDL::MGDrawPointSeq()
void display_control_polygon()const;

///Draw this line's 1st and 2nd coordinates in 2D space
///using drawing function moveto( , ) and lineto( , ).
///wind[] is the window of the screen to draw the line in.
///Clipping will be performed about the wind[].
///(wind[0], wind[1]) is the center coordinates of the window.
///wind[2] is width and wind[3] is hight of the window. When wind[2]<=0,
///no clipping is performed. Even when wind[2]<=0, wind[3] is necessary 
///to input to specify the resolution of the line. In this case,
///wind[0] and wind[1] are not referended.
///ynum is the resolution of the line, is the number of
///straight line segments for the curve length of wind[3](height of window).
///***draw_2D does not perform box including judment, always performs clipping
///operation and draws the line. Users must do obvious box inclusion test
///if maximum drawing performance is necessary.
void draw_2D(void (*moveto)(int, int), void (*lineto)(int, int),
	const double wind[4],	///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

void draw_2D(void (*moveto)(float, float), void (*lineto)(float, float),
	const double wind[4],	///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

void draw_2D(void (*moveto)(double, double), void (*lineto)(double, double),
	const double wind[4],///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

///Draw this line's coordinate'th coordinate in 2D space as
///(t, LBRep(coordinate)) when t_is_x is true, 
///or as ( LBRep(coordinate),t)  when t_is_x is false,  
///using drawing function moveto(int, int) and lineto(int,int).
///Here t is the parameter of the LBRep.
///The other behaviours are the same as draw_2D.
void draw_1D(void (*moveto)(int, int), void (*lineto)(int, int),
	size_t coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

void draw_1D(void (*moveto)(float, float), void (*lineto)(float, float),
	size_t coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

void draw_1D(void (*moveto)(double, double), void (*lineto)(double, double),
	size_t coordinate,		///<id of coordinate, that is =0:x, =1:y, and so on.
	bool t_is_x,			///<=true:t is x coordinate, and false:t is y.
	const double wind[4],	///<window box to draw the line in.
	size_t ynum		///<Resolution of the line.
)const;

void drawSE(
	double span_length,	///<Line segment span length.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve.
						///<Draw will be performed from t0 to t1.
)const;

/// Evaluate right continuous n'th derivative data.
/// nderiv=0 means positional data evaluation.
MGVector eval(
	double,			///< Parameter value.
	size_t nderiv=0,///< Order of Derivative.
	int leftcon=0	///< Left continuous(leftcon=true)
					///< or right continuous(leftcon=false).
)const;

///Compute position, 1st and 2nd derivatives.
/// パラメータ値を与えて位置、一次微分値、二次微分値をもとめる。
void eval_all(
	double,		///< Input parameter value(パラメータ値)
	MGPosition&,///< Position(位置)
	MGVector&,	///< 1st derivative(1次微分値)
	MGVector&	///< 2nd derivative(2次微分値)
)const;

///Evaluate all of i'th derivative data for 0<=i<=nderiv.
///Output will be put on deriv[j+i*sdim()]
///for 0<=i<=nderiv and 0<=j<sdim(), i.e. 
///deriv[j+i*sdim()] is i-th derivative data for 0<=j<sdim(). 
void eval_all(
	double tau,		///< Parameter value to evaluate.
	size_t nderiv,	///< Order of Derivative.
	double* deriv,	///< Output area of size (nderiv+1)*sdim().
	int leftcon=0	///<Left continuous(leftcon=true) or right continuous(leftcon=false).
)const;

///Evaluate line data at data point seq.(BLELIN)
void eval_line(
	MGENDCOND begin,		///<Begin end condition
	MGENDCOND end,			///<End end conditoion 
	const MGNDDArray& tau,	///<Data points.
	MGBPointSeq& value		///<Values evaluated.
)const;

///Evaluate line data at data point seq.(BLELIN)
void eval_line(
	const MGNDDArray& tau,	///<Data points.
	MGBPointSeq& value		///<Values evaluated.
)const{eval_line(MGENDC_NO,MGENDC_NO,tau,value);};

///Extrapolate the curve by the chord length.
///The extrapolation is C2 continuous if the order >=4.
MGLBRep& extend(
	int start,		///<Flag of start or end poit of the line.
					///<If start is true extend on the start point.
	double length,	///<chord length to extend. 
	double dk     ///<Coefficient of how curvature should vary at the connecting point.
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point.
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
);

///Extrapolate the curve by the parameter value.
MGLBRep& extend_with_parameter(
	double tau,	///<The parameter value at the end of extended point.
				///<When tau<param_s(), extension will be done at the starting point.
				///<When tau>param_e(), extension will be done at the end point.
	double dk	///<Coefficient of how curvature should vary at the connecting point.
				///<See extend();
);

///Extracts control points.
///Fucntion's return value is 
///true if control points was obtained, false if not.
bool get_control_points(
	MGBPointSeq& cpoints	///<Control points will be output.
)const;

/// Return This object's typeID
long identify_type() const;

///Test if this cure is co-planar with the 2nd curve curve2.
///MGPlane expression will be out to plane if this is co-planar.
///Function's return value is true if co-planar.
bool is_coplanar(const MGCurve& curve2, MGPlane& plane)const;

///Test if this cure is planar or not.
///MGPlane expression will be out to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

/// Spline と Curve の交点を求める。
///Intersection point of spline and curve.
MGCCisect_list isect(const MGCurve& curve2) const;
MGCCisect_list isect(const MGStraight& curve2)const;
MGCCisect_list isect(const MGRLBRep& curve2)const;
MGCCisect_list isect(const MGEllipse& curve2)const;
MGCCisect_list isect(const MGLBRep& curve2)const;
MGCCisect_list isect(const MGSurfCurve& curve2)const;
MGCCisect_list isect(const MGBSumCurve& curve2)const;

///Intersection of MGLBRep and Surface
MGCSisect_list isect(const MGSurface& surf) const;
MGCSisect_list isect(const MGPlane& surf) const;
MGCSisect_list isect(const MGSphere& surf)const;
MGCSisect_list isect(const MGCylinder& surf)const;
MGCSisect_list isect(const MGSBRep& surf)const;
MGCSisect_list isect(const MGRSBRep& surf)const;
MGCSisect_list isect(const MGBSumSurf& surf)const;

///Access to i-th element of knot.
///( left-hand side version)
double& knot(size_t i){return m_knot_vector(i);}

///Access to i-th element of knot.
///(right hand side version)
double knot(size_t i) const{return m_knot_vector(i);}			

///Returns a pointer to
/// the knot vector data.
const double* knot_data() const{return m_knot_vector.data();}
double* knot_data(){return m_knot_vector.data();}

///Returns the knot vector.
///(RHS version)
const MGKnotVector& knot_vector() const{return m_knot_vector;}

///Returns the knot vector.
///(LHS version)
MGKnotVector& knot_vector(){return m_knot_vector;}

/// 自身に指定したパラメータ範囲のlimitをつける。
///Get the sub interval line of the original line.
MGLBRep& limit(const MGInterval& );

///Returns the B-coef's(RHS version).
const MGBPointSeq& line_bcoef() const{return m_line_bcoef;}

///Returns the B-coef's(LHS version).
MGBPointSeq& line_bcoef(){
	update_mark();
	return m_line_bcoef;}

///Modify the original line by moving move_point to to_point. fix_point can be
///applied according to move_kind.
/// move_kind=1: Start and end point of the line are fixed. The line is modified
///              linearly so that move_point_param point is the maximum move.
///          =2: The point fix_point[0] is fixed and the other end of
///				the move_point_param side is moved. In this case, maximum move
///              is the end point of the line.
///          =3: fix_point[0]<move_point_param<fix_point[1], and two point
///              fix_point[.] are fixed. The line is modified
///              linearly so that move_point_param point is the maximum move.
///   otherwise: Two fix point fix_point[.] are computed so that the modify
///	            range is the minimum. Other move is same as move_kind=3.
/// Restriction: For the case move_kind=3, actual fix point is wider than
///  specified range. The range is the smallest one possible including
///  fix_point[].
MGLBRep& move(
	int move_kind,			///<Indicates how to move line.
	double move_point_param,///<indicate object point to move by the parameter value.
	const MGPosition& to_point,	///<destination point of the abve source point.
	const double fix_point[2]
);

///Change direction of the line.
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///Returns the order.
unsigned order() const{return m_knot_vector.order();}

/// Return ending parameter value.
double param_e() const;

/// Return starting parameter value.
double param_s() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;
	
///Compute part of this curve from parameter t1 to t2.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGCurve* part(
	double t1, double t2,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

/// 与ポイントから曲線へ下ろした垂線の足の，曲線のパラメータ値を
/// すべて求める。
///Return all the foots of the  straight lines that is perpendicular
///to the line.
MGCParam_list perps(
	const MGPosition& point	///< 与ポイント
)const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const;
MGPosition_list perps(const MGStraight& crv2)const;
MGPosition_list perps(const MGRLBRep& crv2)const;
MGPosition_list perps(const MGEllipse& crv2)const;
MGPosition_list perps(const MGLBRep& crv2)const;
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const;

///Check if the line B-rep is planar.
///Funtion's return value is;
/// 0: Not planar, nor a point, nor straight line.
/// 1: B-Rep is a point.		2: B-Rep is a straight line.
/// 3: B-Rep is planar.
int planar(
	MGPlane& plane	///<When Brep is not straight line nor a point, plane is returned.
	///<Even when not planar(return value is 0), plane nearest is returned.
	, MGStraight& line	///<When Brep is a line, line is returned.
	, MGPosition& point	///<When Brep is a point, point is returned.
)const;

///Change the B-Rep by decreasing B-Rep dimension by ndec. This is
///an approximation of the origimal B-Rep. Return value is error flag.
int reduce(
	int ndec	///<Number of B-rep dimension to decrease 
);

///Change an original B-Rep to new one with subdivided knot configuration.
///Knots t must be subdivided knots.
MGLBRep& refine(				///<BLUNK
		const MGKnotVector& t	///<knot vector
);

///ノット削除関数
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///removal knot. line_zero tolerance is used.
void remove_knot(){remove_knot(0,0);};

///Remove knot if removed line has the difference less than line_zero();
///The difference is checked only for the space id of coef(.,j+k)
///of j=0, ..., snum-1. When snum=0, snum is set as sdim();
void remove_knot(size_t j, size_t snum);

///ノット削除関数(1つのノット)
///戻り値は、削除したノットの数
///When snum!=0, tolerance of totalTol is checked only for coef(.,sid+j),
///where j=0, ..., snum-1. When snum=0, snum is set as sdim();
int remove_knot_one(
	double line0,	///<Tolerance allowed for the knot removal.
					///<When line0=<0., the knot will be removed unconditionally.
	size_t	nKnot,	///<削除しようとするノットの番号
	double& totalTol,///<誤差合計
	size_t& num_knot,///<Remained knot number at knot(id) after removed.
	size_t sid=0,	///<Space dimension start id of this LBRep's B-coef.
	size_t snum=0	///<Num of space dimension for the totalTol tolerance check.
);

///Returns the space dimension.
size_t sdim() const{return m_line_bcoef.sdim();}

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist				///<distance to end edge.
)const;

/// 曲線のタイプをを返す。
///Return the curve type.
MGCURVE_TYPE type() const{return MGCURVE_SPLINE;}

/// ｌｉｍｉｔをはずす。
MGCurve& unlimit(){return *this;};

///Unlimit parameter range of the curve to the end point direction
///(終点方向にlimitをはずす)
MGCurve& unlimit_end(){return *this;};

///Unlimit parameter range of the curve to the start point direction
///(始点方向にlimitをはずす)
MGCurve& unlimit_start(){return *this;};

///Returns the knot vector.
///MGKnotVector knot_vector_real() const{return knot_vector();}

///IGES output function. PD126.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Debug Function
std::ostream& out(std::ostream&) const;

std::string whoami()const{return "LBRep";};

protected:

///Compute intersection point of 1D sub B-Rep of original B-rep.(BLIPP)
///Parameter values of intersection point will be returned.
///isect_1D covers this LBRep's C0 ontinuity.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	size_t coordinate=0	///< Coordinate kind of the data f(from 0).
						///< Id of m_line_bcoef.
)const;	

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> oneD(
	const double g[4]	///<Plane expression(a,b,c,d) where ax+by+cz=d.
)const;

///メンバデータを読み出す関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
void ReadMembers(MGIfstream& buf);

///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
void WriteMembers(MGOfstream& buf)const;

private:

///Construct MGLBRep that interpolates data points (tau(i), rlb.eval(tau(i)))
/// for 0<=i<tau.length(). tau is output of data point used to construct
/// MGLBRep.
MGLBRep(const MGRLBRep& rlb, const MGNDDArray& tau);

///Function for BLUMIX constructor. Given 3D point F,
///compute correct point of F that is closest to F.
MGPosition closest_mix(
	unsigned coordinate1,	///<Missing coordinate kind of this line.
	unsigned coordinate2,	///<Missing coordinate kind of Point P.
	double tau,				///<Parameter value of P of other line,
							///<used only when coordinate2=3.
	const MGPosition& P,	///<Point of 2nd line.
	const MGPosition& F		///<3D point, used to coose closest point to F
							///<when more than one point are found.
)const;

///Compute whole box of the curve. Retured is a pointer of a newed MGBox.
MGBox* compute_box() const;

///Draw this by converting straight line segments.
void drawgl(
	double dl,		///<approximate line length of the straight line segments.
	double tstart, double tend	///<start and end parameter value of this.
)const;

///Provide divide number of curve span for function intersect.
size_t intersect_dnum()const;

///compute isects by splitting this curve to sub-curves that do not have
///C0 points in it.
MGCCisect_list isect_by_split_to_C1(const MGCurve& crv2)const;

///isect for this LBRep that does not include C0 continuity points in it.
///For general intersection computation, use isect.
MGCCisect_list C1isect(const MGCurve& crv2) const;
MGCCisect_list C1isect(const MGStraight& crv2) const;

///perps for this LBRep that does not include C0 continuity points in it.
///lb2 may include them. For general perps computation, use perps.
MGPosition_list C1perps(const MGCurve& crv2)const;

///isect_order2 is a private function for isect, computes intersections of
///LBRep(this) of order 2, i.e. polyline B-Rep, with crv or surface.
MGCCisect_list isect_order2(const MGCurve& crv2) const;
MGCSisect_list isect_order2(const MGSurface& crv2) const;

///isect with SurfCurve whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGCCisect_list isect_with_noCompoSC(const MGSurfCurve& curve2)const;

///Compute intersections with MGLBRep curve2 that does not have C0 continuity in it.
MGCCisect_list isect_withC1LB(const MGLBRep& crv2)const;

///compute isects by splitting this curve to sub-curves that do not have
///C0 points in it.
MGCSisect_list isect_by_split_to_C1(const MGSurface& surf)const;

///perps_by_split_to_C1 is a private function for perps, computes perpendicular
///points of LBRep(this) with crv.
MGPosition_list perps_by_split_to_C1(const MGCurve& crv2)const;

///Perpendicular points with C1 conitnuity LBRep lbC1.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv's parameter value.
MGPosition_list perps_withC1LB(
   const MGLBRep& lbC1
)const;

///Perpendicular points with SurfCurve
///whose m_curve is not a MGTrimmedCurve of MGCompositeCurve.
MGPosition_list perps_with_noCompoSC(const MGSurfCurve& curve2)const;

///perps_order2 is a private function for isect, computes intersection of
///LBRep(this) of order 2, i.e. polyline B-Rep, with crv.
MGPosition_list perps_order2(const MGCurve& crv)const;

///isect for this LBRep that does not include C0 continuity points in it.
///lb2 may include them. For general perpendicular computation, use isect.
MGPosition_list perps_1span(const MGLBRep& lb2)const;

///Member Data
private:
	MGKnotVector m_knot_vector;	///< Knot Vector.
	MGBPointSeq  m_line_bcoef;	///< Line B-Coef.

friend class MGCurve;
friend class MGStraight;
friend class MGRLBRep;

};

/** @} */ // end of GEO group
#endif
