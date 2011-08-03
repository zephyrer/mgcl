/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSBRep_HH_
#define _MGSBRep_HH_

#include <assert.h>
#include "mg/MGCL.h"
#include "mg/Vector.h"
#include "mg/SPointSeq.h"
#include "mg/KnotVector.h"
#include "mg/Surface.h"
#include "mg/CSisect_list.h"
#include "mg/SSisect_list.h"

// MGSBRep.h
//
// Forward Declaration
class  MGBox;
class  MGNDDArray;
class  MGPosition;
class  MGBPointSeq;
class  MGKnotArray;
class  MGMatrix;
class  MGTransf;
class  MGLBRep;
class  MGSBRepEndC;
class  MGSBRepTP;
class  MGSBRepVecTP;
class  MGStraight;
class  MGEllipse;
class  MGPlane;
class  MGIfstream;
class  MGOfstream;

/** @addtogroup GEO
 *  @{
 */

///Defines Surface B-Representation, that is , B-Spline surface.
///MGSBRep is a tensor product surface representation of MGLBRep.
///It has:
///(1) Surface control polygon MGSPointSeq surface_bcoef.
///(2) U direciton knot vector MGKnotVector uknot.
///(3) V direction knot vector MGKnotVector vknot.
class MGCLASS MGSBRep: public MGSurface {

public:

/// 与えられたスケーリングで曲面の変換を行いオブジェクトを生成する。
///Scaling.
MGDECL friend MGSBRep operator* (double scale, const MGSBRep& sb);

///build_SBRep_1direction constructs MGSBRep brep, by interpolating only u-direction.
///Input (utaui, tu) will make brep's v direction parameter, and (vtaui, tv) will make
/// udirection parameter. That is, u and v will be exchanged for output brep.
friend void build_SBRep_1direction(	///<Derivative Inf and Knot Vector
	MGSBRepEndC& endc,			///<end condition, only perimeter 3 and 1 are used.
	const MGNDDArray& utaui,	///<Data point of u-direction	of value
	const MGNDDArray& vtaui,	///<Data point of v-direction	of value
	const MGSPointSeq& value,	///<Data point ordinate
	const MGKnotVector& tu,	///<knot vector of u-direction
	const MGKnotVector& tv,	///<knot vector of v-direction
	MGSBRep& brep,			///<constructed b-rep will be output
	int &error			///<Error flag.
);

//////////// Constructor ////////////

///Default constructor(dummy surface brep).
MGSBRep();

/// Construct Surface B-rep by providing raw data of Surface B-Rep.
/// ***** This is the fundamental constructor.*****
MGSBRep(
	const MGSPointSeq& vertex,	///<Control Vertex of Surface B-Rep
	const MGKnotVector& tu,		///<knot vector of u-direction
	const MGKnotVector& tv	///<knot vector of v-direction
);

//**** 1. Interpolation Constructor ****

/// Construct Surface B-rep by intepolation from Point data only.
MGSBRep(
	const MGSPointSeq& points,	///<Point seq data
	int& error,				///<Error flag.
	unsigned orderu=4,		///< Order of u-direction
	unsigned orderv=4		///< Order of v-direction
);

/// Construct Surface B-rep of any order number by interpolation 
///from data point, point data , and knot vector.
MGSBRep(
	const MGNDDArray& utau,		///<Data point abscissa of u-direction
	const MGNDDArray& vtau,		///<Data point abscissa of v-direction
	const MGSPointSeq& points,	///<Point seq data
	const MGKnotVector& tu,		///<knot vector of u-direction
	const MGKnotVector& tv,		///<knot vector of v-direction
	int &error				///<Error flag.
);

/// Construct Surface B-rep of order 4 by interpolation from Point data
/// and tangent plane end condition.
/// Inner point must be homogeneous, should not include first derivative inf.
MGSBRep(				///<Tangent Plane
	const MGSBRepTP& tp,		///<end condition	of Tangent Plane
	const MGNDDArray& utau,		///<Data point of u-direction
	const MGNDDArray& vtau,		///<Data point of v-direction
	const MGVector   uvec[4],	///<Tangent vector of u-direction for 
			///< v=min and v=max boundary line.
			///<uvec[0], uvec[1]: start and end of v=min line
			///<uvec[2], uvec[3]: end and start of v=max line.
	const MGVector   vvec[4],	///<Tangent vector of v-direction for 
			///< u=min and u=max boundary line.
			///<vvec[0]:start of u=min line, vvec[1]:start of u=max line
			///<vvec[2]: end  of u=max line, vvec[3]: end  of u=min.
	///< It is allowed that uvec[i] or vvec[j] is null vector(i.e. sdim==0).
	///<When this is the case, it means no inf provided
	///<about the tangent vector.
	const MGSPointSeq& points,	///<Data point ordinate
	int &error				///<Error flag.
);

/// Construct Surface B-rep of order 4 by interpolation from Point data
///and end condition.
///Inner point may include first derivative inf if the corresponding data
///points are multiple.
///For perimeters, utaui and vtaui do not have multiplicity. However,
///if utaui or vtaui has multiplicity at inner point, this means 1st derivative
///data is provided for the associated value, i.e.:
///If utaui has multiplicity 2 as utaui(i)=utaui(i+1), value(i,.,.) is
///1st derivative along u direction at utaui(i) and value(i+1,.,.) is positional
///data at utaui(i)(=utaui(i+1)).
///If utaui has multiplicity 3 as utaui(i)=utaui(i+1)=utaui(i+2),
///value(i,.,.) is 1st derivative along u direction at utaui(i)- ,
///value(i+1,.,.) is positional data at utaui(i)(=utaui(i+1)), 
///value(i+2,.,.) is 1st derivative along u direction at utaui(i)+.
///Maximum multiplicity allowed is 3.
///About vtaui, the same.
MGSBRep(				///<Derivative Inf.
	MGSBRepEndC& endc,		///<end condition
	const MGNDDArray& utaui,///<Data point of u-direction of value
	const MGNDDArray& vtaui,///<Data point of v-direction	of value
	const MGSPointSeq& value,///<Data point ordinate
	int &error			///<Error flag.
);

/// Construct Surface B-rep of specified order by interpolation from Point data
///and end condition.
///Inner points must not include first derivative inf and the data
///points must not be multiple if the order is not 4.
MGSBRep(				///<Derivative Inf.
	MGSBRepEndC& endc,			///<end condition
	size_t orderu, size_t orderv,///<order along u or v direction,
								///<must be greater than  or equal to 4.
	const MGNDDArray& utaui,	///<Data point of u-direction of value
	const MGNDDArray& vtaui,	///<Data point of v-direction	of value
	const MGSPointSeq& value,	///<Data point ordinate
	int &error			///<Error flag.
);

/// Construct Surface B-rep by interpolation from Point data
/// and end condition with knot vector.
///This is a general version of:
///		MGSBRep(				///Derivative Inf.
///		MGSBRepEndC& endc,			///end condition
///		const MGNDDArray& utaui,	///Data point of u-direction of value
///		const MGNDDArray& vtaui,	///Data point of v-direction	of value
///		const MGSPointSeq& value,	///Data point ordinate
///		int &error);			///Error flag.
///The difference is to provide or not knot vector tu and tv.
///(utaui(i),vtaui(j),value(i,j)) is one pair of data,
///      for 0<=i<=utaui.lenght(), 0<=j<=vtaui.lenghth().
///Knot vector tu and tv must satisfy Schoenberg-Whitney consdition, i.e.
///t(i)<tau(i)<t(i+k) for all i, where k=order.
///If not, error=2 will be returned.
///(t=tu, tau=utaui+multiplicity included at both ends according to endc,
///for udirection. About v, the same)
MGSBRep(			///<Derivative Inf and Knot Vector
	MGSBRepEndC& endc,			///<end condition
	const MGNDDArray& utaui,	///<Data point of u-direction	of value
	const MGNDDArray& vtaui,	///<Data point of v-direction	of value
	const MGSPointSeq& value,	///<Data point ordinate
	const MGKnotVector& tu,	///<knot vector of u-direction
	const MGKnotVector& tv,	///<knot vector of v-direction
	int &error			///<Error flag.
);

///Construct Surface B-rep by sweep B-Rep and sweep length.
///The sweep surface is defined as:
///lbrep(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSBRep(
	const MGLBRep& lbrep,		///<Sweep crv.
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist			///<distance to end edge.
);

/// Construct Surface B-rep of specified order by interpolation from Point data.
///This is an easy-to-version.
MGSBRep(				///<Derivative Inf.
	const MGNDDArray& utau,	///<Data point of u-direction of value
	const MGNDDArray& vtau,	///<Data point of v-direction	of value
	const MGSPointSeq& value,	///<Data point ordinate
	size_t orderu=4, size_t orderv=4///<order along u or v direction.
);

///Construct Surface B-rep by sweep Straight and sweep length.
///The sweep surface is defined as:
///st(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSBRep(
	const MGStraight& st,		///<Sweep crv.
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist			///<distance to end edge.
);

///Construct Surface B-rep from lines and derivatives.
///Interpolation will be done only along one parameter direction,
///along v.
///tau provides data point sequence along v direction of surface (u,v) parameter
///configuration. deriS and deriE are used to provide the 1st derivative
///B-representation along the perimeter 0 and 2, may be null
///if 1st derivative B-rep is not provided. If derivative
///B-rep is provided, deriS and deriE must have the same knot configuration
///as the one of lines which makes u kont configuration of this surface (u,v)
///parameter. tau[i] is the parameter for lines[i].
MGSBRep(
	const MGNDDArray& tau,
	const std::vector<MGLBRep*>& lines,
	const MGLBRep* deriS,
	const MGLBRep* deriE
);

//**** 2. Approximation Constructor ****

///Approximate an original B-Rep by a new knot configuration.
///The new knot config must be inside the range of the original B-Rep
///parameter. However new knots may be coarse or fine.
MGSBRep(
	const MGSBRep& old_brep,///<Original B-Rep.
	const MGKnotVector& ut,	///<knot vector of u-direction
	const MGKnotVector& vt,	///<knot vector of v-direction
	int &error			///<Error flag.
);

///Obtain Surface B-Rep by Least square approximation.
MGSBRep(				
	const MGNDDArray& utau,	///<Data point abscissa of u-direction
	const MGNDDArray& vtau,	///<Data point abscissa of v-direction
	const MGSPointSeq& points,	///<Point seq data
	const MGKnotVector& tu,	///<knot vector of u-direction
	const MGKnotVector& tv,	///<knot vector of v-direction
	const MGSPointSeq& weight	///<Weight for each point(space dimension 1).
	///<weight(i,j,0) is the weight for points(i,j,.).
	///<When weight(i,j,0) is large compared with other weights, the deviation
	///<from the surface at (utau(i),vtau(j)) is small.
);

///Obtain Surface B-Rep by Shoenberg and Reinch's smoothing function
///If dp(i,j,0) becomes bigger, deviation at (utau(i),vtau(j)) becomes bigger.
MGSBRep(				
	const MGNDDArray& utau,	///<Data point abscissa of u-direction
	const MGNDDArray& vtau,	///<Data point abscissa of v-direction
	const MGSPointSeq& points,	///<Point seq data
	const MGSPointSeq& dp,	///<Error estimate for each point(space dimension 1).
	double deviation		///<Mean deviation of each point.
);

///4本の境界線、ブレンド関数、接続面を与えて面を生成する。
///境界線はvmin,umax,vmax,uminの順で、vmin,vmaxの向きをuminからumaxの方向に
///umin,umaxの向きをvminからvmaxの方向になっているものとする。境界線のノットベクトル
///をあわせるときの誤差はline_zero()を使用している。ブレンド曲線はパラメータ範囲0,1
///で値域も0,1である。接続面(MGSBRepTP)のパラメータ範囲は各境界線と同じとする。
///		エラーコード：
///		0:			正常終了
///		-1:			境界線がC1連続でなかった
///		-2:			接続面のパラメータ範囲が境界線と違った
///		-3:			境界線のノットベクトルをあわせられなかった
///		-4:			接続面のパラメータ範囲を変更できなかった
///		-5:			点列が求まらなかった
///		-6:			端末条件が求まらなかった
///		-7:			面が生成できなかった
MGSBRep(
	const MGCurve*	edge_crvl[4],///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGCurve&	blendCrvU,	///<空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve&	blendCrvV,	///<空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGSBRepTP&tp,			///<接続面(パラメータ範囲は境界線と同じ)
	int& error		///<エラーコード
);

///Easy to use version of the above constructor.
///When blendCCrvU,V were null, straight line from 0. to 1. will be used.
MGSBRep(
	MGPvector<MGLBRep>& edges,	///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	int&			error	,	///<エラーコードが出力される（0:正常終了）
	const MGCurve*	blendCrvU=0,///<空間次元1のu方向のブレンド曲線(パラメータ、値域ともに0,1)
	const MGCurve*	blendCrvV=0 ///<空間次元1のv方向のブレンド曲線(パラメータ、値域ともに0,1)
);

MGSBRep(
	const MGCurve*	edge_crvl[4],///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGSBRepTP&tp,			///<接続面(パラメータ範囲は境界線と同じ)
	int& error,					///<エラーコード
	const double* alpha=0		///<derivative magnitude coefficients for
		///<perimeter i in alpha[i].
		///<=1. is normal, if<1., curvature will be made large, and
		///<if>1. curvature will be made small.
		///<If alpha=null, all of the coefs are assumed to be 1.
);

///Generalized ruled surface construction.
///Build a surface by ruled surface method. That is, costruct a surface by sliding
///the blending curve(a tie curve of the rails) of perimeters 3 and 1 along
///perimeter 0 and 2(0 and 2 make the rail).
///Or by sliding the blending curve of perimeter 0 and 2
///along perimeter 3 and 1(3 and 1 make the rail).
MGSBRep(
	bool along_u,	///<indicates which perimeters make a rail.
					///<if true, perimeter 0 and 2, if false, perimeter 3 and 1 make rail.
	const MGPvector<MGLBRep>& perimeters
			///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
			///<perimeters must be the same knot configuration along u(perimeter 0 and 2)
			///<and along v(perimeter 3 and1).
);

//**** 3.Conversion Constructor.****

///Gets new B-Rep by adding knots to an original B-Rep.
MGSBRep(
	const MGSBRep& old_brep,	///<Original B-Rep.
	const MGKnotArray& uknots,	///<Knots to add for u-direction
	const MGKnotArray& vknots	///<Knots to add for v-direction.
);

///Gets new B-Rep by connecting two B-Rep to one.
///The parameter (which1,continuity,which2,opposite) can be obtained by
///public function continuity.
MGSBRep(
	const MGSBRep& brep1,	///<B-Rep 1.
	int which1,				///<which perimeter of brep1.
	int continuity,			///<continuity. Must >=0.
	const MGSBRep& brep2,	///<B-Rep 2.
	int which2,				///<which perimeter of brep2.
	int opposite			///< Input if parameter direction of which2
					///< is the same as which1 along common edge.
					///< If opposite is true, the direction is opposite.
);

///Obtain partial Surface B-Rep restricted by sub interval of u and v parameter
///range.
///If multiple==true(!=0), knot_u(i)=t1 and knot_u(n+i)=t2 for i=0,..., k-1
///will be guaranteed. Here, n=bdim_u(), k=order_u(),
///t1=uvrange(0).low_point(), and t2=uvrange(0).high_point().
///About knot_v(j), the same.
/// Both u-range and v-range must be inside the range of old_brep.
/// Gets new B-Rep by computing a part of the original. New one is exactly
/// the same as the original except that it is partial.
MGSBRep(
	const MGBox& uvrange,		///<u and v parameter range.
	const MGSBRep& old_brep,	///<Original B-Rep.
	int multiple=0///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

/// Construct a Surface B-Rep by changing space dimension and ordering
///of coordinates.
MGSBRep(
	size_t dim,				///< New space dimension.
	const MGSBRep& sbrep,	///< Original Surface B-rep.
	size_t start1=0, 		///< Destination order of new Surface.
	size_t start2=0); 		///< Source order of original Surface.

/// Construct a Surface B-Rep (order = 2 ) from Plane and Parameter ranges.
MGSBRep(
	const MGPlane& plane,	///< Original Plane
	const MGBox& prange	///< parameter range of new Surface.
);

//	MGSBRep(const MGSBRep&);  ///Copy constructor.
//  We can use default copy constructor.

///リブ曲線列から面を作成する
///作成する面のノットベクトルはリブ曲線の向きをu,リブ列方向をvとする
///This constructor only generates MGSBRep even if curves are MGRLBRep of the same
///knot configuration. To avoid this, use createSurfaceFromRibs() that generates
///MGRSBRep when curves are MGRLBRep of the same knot configuration.
///Let v0=start parameter value, v1=terminate parameter value along v, then
///v=v0 const parameter line is curves[0], and v=v1 const parameter line is
///curves[n-1], where n=curves.size(). n must be greater or equal to 2.
///When n==2, the surface is a ruled surface(that is, this->order_u() is 2).
///When direction_adjustment=true, curves[i] directions are adjusted to line up
///to the same direciton as the curves[i-1]'s.
///When direction_adjustment=false, no curves of curves[i] are negated and
///the directions of curves[i] are the direcitons of ribs.
MGSBRep(const std::vector<const MGCurve*>& curves,bool direction_adjustment=true);
MGSBRep(const MGPvector<MGCurve>& curves,bool direction_adjustment=true);

///Destructor
///	~MGSBRep();	///We can use default destructor.

//////////// Operator overload ////////////

///Assignment.
///When the leaf object of this and srf2 are not equal, this assignment
///does nothing.
MGSBRep& operator=(const MGGel& gel2);
MGSBRep& operator=(const MGSBRep& gel2);

/// 曲面の平行移動を行いオブジェクトを生成する。
///Translation.
MGSBRep operator+ ( const MGVector & ) const;

/// 曲面の逆方向に平行移動を行いオブジェクトを生成する。
///Translation.
MGSBRep operator- ( const MGVector & ) const;

/// 与えられたスケーリングで曲面の変換を行いオブジェクトを生成する。
///Scaling.
MGSBRep operator* (double) const;

/// 与えられた変換で曲面の変換を行いオブジェクトを生成する。
///Matrix transformation.
MGSBRep operator* ( const MGMatrix& ) const;

/// 与えられた変換で曲面のトランスフォームを行いオブジェクトを生成する。
///General transformation.
MGSBRep operator*( const MGTransf & ) const;

///Object transformation.
MGSBRep& operator+=(const MGVector& v);
MGSBRep& operator-=(const MGVector& v);
MGSBRep& operator*=(double scale);
MGSBRep& operator*=(const MGMatrix& mat);
MGSBRep& operator*=(const MGTransf& tr);

///Comparison of two objects.
bool operator==(const MGSBRep& gel2)const;
bool operator<(const MGSBRep& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator<(const MGGel& gel2)const;
bool operator!=(const MGGel& gel2)const{return !(gel2==(*this));};
bool operator!=(const MGSBRep& gel2)const{return !(gel2==(*this));};
bool operator==(const MGRSBRep& gel2)const;

//////////// Member Function ////////////

/// 自身の曲面の全体の面積を返却する。
///Compute total surface area.
///	double area () const;

/// 与えられたパラメータ範囲の曲面の面積を返す。
///Compute surface area limitted by parameter range box.
///	double area(const MGBox& box) const;

///Returns B-Rep Dimension of u.
size_t bdim_u() const{return m_surface_bcoef.length_u();}

///Returns B-Rep Dimension of v.
size_t bdim_v() const{return m_surface_bcoef.length_v();}

///Compute minimum box that includes the linitted parameter range surface.
MGBox box_limitted(const MGBox& uvrange) const;

///Changing this object's space dimension.
MGSBRep& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
MGSBRep& change_range(
	int is_u,	///<if true, (t1,t2) are u-value. if not, v.
	double t1,	///<Parameter value for the start of original. 
	double t2	///<Parameter value for the end of original. 
);

///Access to (i,j)th element of coef(left-hand side version).
double& coef(size_t i, size_t j, size_t k){
	assert(i<bdim_u() && j<bdim_v() &&k<sdim());
	update_mark();
	return m_surface_bcoef(i,j,k);
}		

///Access to (i,j)th element of coef(right-hand side version).
double coef(size_t i, size_t j, size_t k) const{return m_surface_bcoef(i,j,k);}

///Extract (i,j,k) elements for 0<=k<sdim() as a vector.
MGVector coef(size_t i, size_t j) const{return m_surface_bcoef(i,j);}

///Returns a pointer to the surface b-coef data.
const double* coef_data(size_t i=0, size_t j=0, size_t k=0) const
{return m_surface_bcoef.data(i,j,k);}

///Compute continuity with brep2.
/// Function's return value is:
/// -1: G(-1) continuity, i.e. two surfaces are discontinuous.
///  0: G0 continuity, i.e. two surfaces are connected,
///     but tangents are discontinuous
///  1: G1 continuity, i.e. two surfaces are connected,
///     and tangents are also continuous.
///  2: G2 continuity, i.e. two surfaces are connected,
///     and tangents and curvatures are also continuous.
/// Reuturn value is the continuity.
int continuity(
	const MGSBRep& brep2,	///< Input second SBRep
	int& which1,	///< Outputs which perimeter(which1) of this is
	int& which2,	///< connected to which(which2) of brep2.
					///< These are valid only when continuity>=0.
	int& opposite,	///< Outputs if parameter direction of which2
					///< is the same as which1 along common edge.
					///< If opposite is true, the direction is opposite.
	double& ratio	///< Ratio of 1st derivatives of the two surfaces will
					///< be returned.
			///< ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
) const;

///Construct new surface object by copying to newed area.
///User must delete this copied object by "delete".
MGSBRep* clone() const;

///Construct new surface object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGSBRep* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0 		///< Source order of this line.
)const;

///Display control polygons using mgGDL::MGDrawPointSeq(sp)
void display_control_polygon()const;

///uまたはv方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot(
    MGPvector<MGSurface>& srfl     ///<分割した曲面リスト
) const;

///Evaluate right continuous ndu'th and ndv'th derivative data.
///Function's return value is (d(ndu+ndv)f(u,v))/(du**ndu*dv**ndv).
/// ndu=0 and ndv=0 means positional data evaluation.
MGVector eval(
	double u, double v,	///< Parameter value of the surface.
	size_t ndu=0,	///< Order of Derivative along u.
	size_t ndv=0	///< Order of Derivative along v.
) const;

///Evaluate surface data.
MGVector eval(
	const MGPosition& uv	///< Parameter value of the surface.
	, size_t ndu=0			///< Order of derivative along u.
	, size_t ndv=0			///< Order of derivative along v.
) const;

///Evaluate right continuous surface data.
///Evaluate all positional data and 1st and 2nd derivatives.
void eval_all(
	double u, double v,		///< Parameter value of the surface.
	MGPosition& f,			///< Positional data.
	MGVector&   fu,			///< df(u,v)/du
	MGVector&   fv,			///< df/dv
	MGVector&   fuv,		///< d**2f/(du*dv)
	MGVector&   fuu,		///< d**2f/(du**2)
	MGVector&   fvv			///< d**2f/(dv**2)
)const;

///Evaluate all of derivative data (d(i+j)f(u,v))/(du**i*dv**j),
///for 0<=i<=ndu and 0<=j<=ndv.
void eval_all(
	double u, double v,	///< Parameter value of the surface.
	size_t ndu,		///< Order of Derivative along u.
	size_t ndv,		///< Order of Derivative along v.
	double* deriv	///<Output. (d(i+j)f(u,v))/(du**i*dv**j) in
					///<deriv[r+j*dim+i*(ndv+1)*dim] for 0<=r<dim=sdim().
					///<for 0<=i<=ndu and 0<=j<=ndv.
///<deriv is an array of deriv[ndu+1][ndv+1][r],
///<(d(i+j)f(u,v))/(du**i*dv**j) is returned in deriv[i][j][r].
)const;

///Exchange parameter u and v.
MGSurface& exchange_uv();

///Modify the original Surface by extrapolating the specified perimeter.
///The extrapolation is C2 continuous if the order >=4.
///The extrapolation is done so that extrapolating length is "length"
///at the position of the parameter value "param" of the perimeter.
MGSBRep& extend(
	int perimeter,	///<perimeter number of the Surface.
					///< =0:v=min, =1:u=max, =2:v=max, =3:u=min.
	double param,	///< parameter value of above perimeter.
	double length,	///<chord length to extend at the parameter param of the perimeter.
	double dk=0.  ///<Coefficient of how curvature should vary at
///<    extrapolation start point. When dk=0, curvature keeps same, i.e.
///<    dK/dS=0. When dk=1, curvature becomes zero at length extrapolated point,
///<    i.e. dK/dS=-K/length at extrapolation start point.
///<    (S=parameter of arc length, K=Curvature at start point)
///<    That is, when dk reaches to 1 from 0, curve changes to flat.
);

/// Return This object's typeID
long identify_type() const;

///Test if input parameter value is inside parameter range of the surface.
bool in_range(double u, double v) const;
bool in_range(const MGPosition& uv) const;

/// Surface と Curve の交点を求める。
///Compute curve and surface intersection point(s)
MGCSisect_list isect(const MGCurve& curve)const;
MGCSisect_list isect(const MGStraight& sl)const{return isectSl(sl);};
MGCSisect_list isect(const MGRLBRep& curve)const;
MGCSisect_list isect(const MGEllipse& curve)const;
MGCSisect_list isect(const MGLBRep& curve)const;
MGCSisect_list isect(const MGSurfCurve& curve)const;
MGCSisect_list isect(const MGBSumCurve& curve)const;

///Surface と Surface の交線を求める。
///Surface and Surface intersection.
///Compute intersectio line(s) of two surface.
///Restriction:Currently if two surface do not have intersection on
///any of 4 perimeters, this function does not compute surface to surface
///intersection.
MGSSisect_list isect(const MGSurface& srf2)const;
MGSSisect_list isect(const MGPlane& srf2)const;
MGSSisect_list isect(const MGSphere& srf2)const;
MGSSisect_list isect(const MGCylinder& srf2)const;
MGSSisect_list isect(const MGSBRep& srf2)const;
MGSSisect_list isect(const MGRSBRep& srf2)const;
MGSSisect_list isect(const MGBSumSurf& srf2)const;

///Access to i-th element of u knot
/////(right-hand side version)
double& knot_u(size_t i){return m_uknot(i);}

///Access to i-th element of u knot
///( left-hand side version)
double knot_u(size_t i) const{return m_uknot(i);}

///Access to i-th element of v knot
///( left-hand side version)
double& knot_v(size_t i){return m_vknot(i);}			

///Access to i-th element of v knot
///(right hand side version)
double knot_v(size_t i) const{return m_vknot(i);}

///Returns a pointer to the u knot vector data.
const double* knot_data_u() const{return m_uknot.data();}

///Returns a pointer to the v knot vector data.
const double* knot_data_v() const{return m_vknot.data();}

///Returns the u knot vector.
const MGKnotVector& knot_vector_u() const{return m_uknot;}
MGKnotVector& knot_vector_u(){return m_uknot;};

///Returns the v knot vector.
const MGKnotVector& knot_vector_v() const{return m_vknot;}
MGKnotVector& knot_vector_v(){return m_vknot;};

///Compare two parameter values. If uv1 is less than uv2, return true.
///Comparison is done after prjected to i-th perimeter of the surface.
bool less_than(
	size_t i,	///perimeter number.
	const MGPosition& uv1,
	const MGPosition& uv2) const;

/// 自身に指定したパラメータ範囲のlimitをつける。
///Update this by limitting the parameter range.
///uvrange is parameter value range of (umin, vmin) to (umax, vmax).
MGSBRep& limit(const MGBox& uvrange);

///Modify the original Surface by moving move_point to to_point.
///fix_point can be applied according to move_kind.
///move_kind_u and _v mean as below for u and v direction:
/// move_kind=1: Start and end perimeter of the surface are fixed. The surface
///              is modified linearly so that move_point_param line is
///              the maximum move.
///          =2: The parameter line fix_point[0] is fixed and the other end of
///				the move_point_param side is moved. In this case,
///				maximum move is the end perimeter of the surface.
///          =3: fix_point[0]<move_point_param<fix_point[1], and two parameter
///              line fix_point[.] are fixed. The surface is modified
///              linearly so that move_point_param point is the maximum move.
///   otherwise: Two fix parameter line fix_point[.] are computed so that
///	            the modify range is the minimum. Other move is same as
///				move_kind=3.
/// Restriction: For the case move_kind=3, actual fix parameter line is wider
/// than specified range. The range is the smallest one possible including
/// fix_point[].
MGSBRep& move(
	int move_kind_u,	///<Indicates how to move Surface for u direction.
	int move_kind_v,	///<Indicates how to move Surface for v direction.
	const MGPosition& move_point_param, ///<indicate object point to move
						///< by the (u,v) parameter value.
	const MGPosition& to_point,	///<destination point(x,y,z coordinates) of
						///< the abve source point.
	const MGPosition fix_point[2]///<(u,v) parameter value pair.
);

///Change direction of the surface.
void negate(			
	int is_u	///< Negate along u-direction if is_u is ture, else along v-direction.
);

///Obtain parameter value if this surface is negated by "negate()".
/// Negate along u-direction if is_u is ture,
/// else along v-direction.
MGPosition negate_param(const MGPosition& uv, int is_u=1)const;

///Returns the B-Rep order(u-direction).
unsigned order_u() const{return m_uknot.order();}

///Returns the B-Rep order(v-direction).
unsigned order_v() const{return m_vknot.order();}

/// Return ending parameter value.
MGPosition param_e() const;
double param_e_u() const;
double param_e_v() const;

/// Compute parameter curve.
///Returned is newed area pointer, and must be freed by delete.
MGCurve* parameter_curve(
	int is_u			///<Indicates x is u-value if is_u is true.
	, double x			///<Parameter value.
						///<The value is u or v according to is_u.
) const;

/// Compute parameter line.
MGLBRep parameter_line(
	int is_u			///<Indicates x is u-value if is_u is true.
	, double x			///<Parameter value.
						///<The value is u or v according to is_u.
	, unsigned nderiv=0	///<Order of derivative.
)const;

/// パラメータ範囲を返す。
///Return parameter range.
MGBox param_range() const;

/// Return starting parameter value.
MGPosition param_s() const;
double param_s_u() const;
double param_s_v() const;

///Compute part of the surface limitted by the parameter range uvbx.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
///Retured is newed object, must be deleted.
MGSBRep* part(
	const MGBox& uvbx,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
)const;

/// Compute perimeter line B-Rep.
/// i is perimeter number:
/// =0: v=min line, =1: u=max line, =2: v=max line, =3: u=min line
///perimeter_curve returs newed object, must be deleted.
MGLBRep perimeter(size_t i) const;
MGCurve* perimeter_curve(size_t i) const;

///Return how many perimeters this surface has.
size_t perimeter_num() const{return 4;};

///Test if the surface is planar or not.
///Returned is 0(false) if this is not planar, 1(true) if this is planar.
int planar(
	MGPlane& plane,		///<Plane that might be closest to this.
						///<Plane is always output even if not planar.
	double& deviation	///<maximum deviation of this from the output plane.
) const;

///Test if part of the surface is planar or not within the tolerance tol.
///The part of the surface is input by the surface parameter range uvbox.
///Returned is 0(false) if this is not planar, 1(true) if planar.
int planar(
	const MGBox& uvbox,///<This surface parameter range.
	double tol,	///<maximum deviation allowed to regard the sub surface as a plane.
	int* divideU=0///<Direction to subdivide will be output, if this was not planar,
				///<=1: u direction, =0: v direction.
	) const;

///Change the B-Rep by decreasing B-Rep dimension by ndec. This is
///an approximation of the origimal B-Rep. Return value is error flag.
int reduce(
	int is_u,		///<if true, reduce b-rep dimension of u-direction.
	int ndec);		///<Number of B-rep dimension to decrease 

///Change an original B-Rep to new one with subdivided knot configuration.
///Knots t must be subdivided knots.
MGSBRep& refine(
	const MGKnotVector& uknot,	///< new knot of u-direction
	const MGKnotVector& vknot	///< new knot of v-direction
);

///ノット削除関数
///トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
///removal knot. line_zero tolerance is used.
void remove_knot();

///uノットを一つ削除する
///関数の戻り値は削除したノットの数
int remove_knot_one(
	double line0,		///<Tolerance allowed for the knot removal. 
						///<When line0<=0., removal will be done uncoditionally.
	size_t	id,	///<削除しようとするノットの番号
	double& tol,///<削除後の誤差が出力される
	bool u_knot=true	///<削除対象が（u,v)のいずれのknot vectorかを入力する
				///<=trueのとき、u-knot_vectorを削除
);

///Returns the space dimension.
size_t sdim() const{return m_surface_bcoef.sdim();}

///Shrink this surface to the part limitted by the parameter range of uvbx.
///New parameter range uvbx2 is so determined that uvbx2 is the smallest
///box tha includes uvbx, and all of the u or v values of uvbx2 is one of 
///the values of u or v knots of the surface knotvector.
///uvbx(0) is the parameter (us,ue) and uvbx(1) is (vs,ve).
///That is u range is from us to ue , and so on.
void shrink_to_knot(
	const MGBox& uvbx,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
);

///Returns the B-coef's.
const MGSPointSeq& surface_bcoef() const{return m_surface_bcoef;}

///Returns the B-coef's.
MGSPointSeq& surface_bcoef(){update_mark(); return m_surface_bcoef;}

///Compute surface integral of the 1st two coordinates.
///(面積分）を求める。
///This integral can be used to compute volume sorounded by the surface.
///double surface_integral(const MGBox&) const;

/// 曲面のタイプをを返す。
///Return the surface type.
MGSURFACE_TYPE type() const{return MGSURFACE_SPLINE;}

/// ｌｉｍｉｔをはずす。
///Unlimit the parameter range. Return the same.
MGSurface& unlimit(){return *this;};

///IGES output function. PD128(NURBS Surface).
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

///Debug Function
std::ostream& out(std::ostream&) const;

std::string whoami()const{return "SBRep";};

protected:

///Intersection of Surface and a straight line.
MGCSisect_list isectSl(
	const MGStraight& sl,
	const MGBox& uvbox=mgNULL_BOX ///<indicates if this surface is restrictied to the parameter
					///<range of uvbox. If uvbox.is_null(), no restriction.
) const;

///メンバデータを読み込む関数
/// 戻り値boolは正常に読み出しが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void ReadMembers(MGIfstream& buf);
	
///メンバデータを書き込む関数
/// 戻り値boolは正常に書き込みが出来ればtrue、失敗すればfalseになる
/// ここでは処理対象となるデータメンバが無いので何も処理をしない。
void WriteMembers(MGOfstream& buf) const;

private:

///Member Data
	MGKnotVector m_uknot;			///< Knot Vector of u-direction
	MGKnotVector m_vknot;			///< Knot Vector of v-direction
	MGSPointSeq  m_surface_bcoef;	///< Surface B-Coef.

///Return minimum box that includes whole of the surface.
///Returned is a newed object pointer.
MGBox* compute_box() const;

///Compute m_uknot and m_vknot from point sequence.
void compute_knot(
	const MGSPointSeq& points,	///<Point seq data
	unsigned orderu,			///< Order of u-direction
	unsigned orderv,			///< Order of v-direction
	MGNDDArray& utau,			///< u data point will be output.
	MGNDDArray& vtau			///< v data point will be output.
);

///Compute continuity with brep2.
/// Function's return value is:
/// -1: G(-1) continuity, i.e. two surfaces are discontinuous.
///  0: G0 continuity, i.e. two surfaces are connected,
///     but tangents are discontinuous
///  1: G1 continuity, i.e. two surfaces are connected,
///     and tangents are also continuous.
///  2: G2 continuity, i.e. two surfaces are connected,
///     and tangents and curvatures are also continuous.
/// Reuturn value is the continuity.
int continuity(
	const MGSBRep& brep2,	///< Input second SBRep
	int is_u1,		///< Input if u-direction of this.
	int is_u2,		///< Input if u-direction of brep2.
	int opposite,	///< Input if parameter direction of which2 is equal or not.
	int& which1,	///< Outputs which perimeter(which1) of this is
	int& which2,	///< connected to which(which2) of brep2.
					///< These are valid only when continuity>=0.
	double& ratio	///< Ratio of 1st derivatives of the two surfaces will
					///< be returned.
			///< ratio= d2/d1, where d1=1st deriv of this and d2=of brep2
)const;

///The following two function will be used in perps or isect
///to decide how many division of the surface along u or v direction
///should be applied before using perp_guess or isect_guess.
size_t intersect_dnum_u() const;
size_t intersect_dnum_v() const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
MGCurve* isect_incr_pline(
	const MGPosition& uv,	///<last intersection point.
	int kdt,				///<Input if u=const v-parameter line or not.
							///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,				///<next u value will be output
	double& v,				///<next v value will be output
	size_t incr=0			///<Incremental valuse of B-coef's id.
) const;

///"isect1_incr_pline" is a dedicated function of isect_start_incr, will get
/// shortest parameter line necessary to compute intersection.
void isect_incr_pline2(
	const MGPosition& uv,///<last intersection point.
	int kdt,			///<Input if u=const v-parameter line or not.
						///< true:u=const, false:v=const.
	double du, double dv,///<Incremental parameter length.
	double& u,			///<next u value will be output
	double& v,			///<next v value will be output
	size_t incr,		///<Incremental valuse of B-coef's id.
	MGLBRep& pline		///<parameter line will be output.
) const;

///"isect_sub_interval" is a dedicated function of isect_incr_pline.
///Computes subinterval id and length of parameter line necessary for
///intersection computation.
size_t isect_sub_interval(
	size_t kdt,	///<indicates if u=const(kdt=1) or v=const parameter line is used
				///<or not to get the intersection
	double u, double v,	///< u and v parameter value.
	double du, double dv,///< incremental parameter length so far.
	size_t& index,			///<index of B-Rep coef that is the start of
	///<sub line b-rep(parameter line) for the computation of intersection
	size_t incr			///<Incremental valuse of B-coef's id.
)const;

///Return order of intersection line order of MGLBRep.
///The default is 4.
size_t isect_order() const;

///uノットを削除する
int remove_knot_u();

///uノットを一つ削除する
///関数の戻り値は削除したノットの数
int remove_knot_u_one(
	double line0,		///<Tolerance allowed for the knot removal. 
						///<When line0<=0., removal will be done uncoditionally.
	size_t	id,			///<削除しようとするノットの番号
	double& totalTol,	///<最初は０を入力する。あとはremove_knot_u_oneが管理する。
	///<削除を繰り返しているときは誤差が加算される。これ以上削除できない時０クリアされる。
	size_t& num_knot	///<Remained knot number at knot(id) after removed.
);

///Obtain 1D surface rep. of this surf which can be used for
///isect(const MGPlane& pl). This surf1D is used in isect for
///the argument of isect_startPlane, which will use surf1D to compute isect(pl).
///surf1D=0.(intersection with x=0. plane) is the intersection lines.
MGSBRep* surf1D(const MGPlane& pl)const;

///u方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot_u(
    MGPvector<MGSBRep>& srfl      ///<分割した曲面リスト
) const;

///v方向に折れ(マルチノット)があるとき面を分割する
///戻り値は、分割数を返却する
int divide_multi_knot_v(
    MGPvector<MGSurface>& srfl      ///<分割した曲面リスト
) const;

friend class MGSurface;
friend class MGRSBRep;

};

///Construct 4 perimeters, given at least two of the four.
///Input perimeters may have different knot configuration. In this case they will be updated
///so as to have the same configuration.
///Function's return value indicates which perimeter(s) was missing:
/// 10: all of the 4 were input(and knot configurations were updated to have the same).
///  0: only perimeter 0 was missing.
///  1: only perimeter 1 was missing.
///  2: only perimeter 2 was missing.
///  3: only perimeter 3 was missing.
///  4: perimeter 2 and 3 were missing.
///  5: perimeter 1 and 3 were missing.
///  6: perimeter 1 and 2 were missing.
///  7: perimeter 0 and 3 were missing.
///  8: perimeter 0 and 2 were missing.
///  9: perimeter 0 and 1 were missing.
/// -2: less than 2 perimeters were provided.
int construct_perimeters(
	const MGPvector<MGCurve>& peris,
			///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順). Let i be the perimeter number,
			///<and the data is missing, perimeter[i] must be null. If perimeter 3 data is missing,
			///<perimeters.size() may be 3. If perimeter 2 and 3 data are missing, perimeters.size() may
			///<be 2.
			///<When perimeters were not the same knot configuration along u(perimeter 0 and 2)
			///<or along v(perimeter 3 and1), they will be rebuild to have the same configuration.
	MGPvector<MGLBRep>& perimeters2
			///<new perimeters will be output.
);

///construct boolean sum surface from the four perimeters and their tangent planes.
void bool_sum(
	const MGCurve*	edge_crvl[4],	///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	const MGSBRepTP& tp,			///<接続面(パラメータ範囲は境界線と同じ)
	int& error,	///<エラーコード
	MGSBRep& surf
);

void bool_sum(
	const MGCurve*	edge_crvl[4],	///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
	MGSBRepVecTP& vectp,			///<接続面(パラメータ範囲は境界線と同じ)
	int& error,	///<エラーコード
	MGSBRep& surf
);

/** @} */ // end of GEO group
#endif
