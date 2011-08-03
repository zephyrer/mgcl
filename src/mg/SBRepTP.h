/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGSBRepTP_HH_
#define _MGSBRepTP_HH_

#include <assert.h>
#include <memory>
#include "mg/LBRep.h"

// MGSBRepTP.h
//

class MGSurface;
class MGLBRep;
class MGOfstream;
class MGIfstream;

/** @addtogroup GEORelated
 *  @{
 */

///Defines Tangent Plane Line B-Representation Class.
///Tangent plane is a line b-representation of (unit)normal vector of
///tangent plane along surface perimeter.
class MGCLASS MGSBRepTP{

public:

///String stream Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGSBRepTP& );

//////////// Constructor ////////////

///Default Constructor, will be set as no TPs' are specified.
MGSBRepTP();

///Copy Constructor.
MGSBRepTP(const MGSBRepTP&);

///Compute TP of four boundaries of a Surface B-Rep.
MGSBRepTP(const MGSurface& srf);
	  
//////////// Destructor ////////////

~MGSBRepTP();

//////////// Operator overload. ////////////

///Assignment.
MGSBRepTP& operator=(const MGSBRepTP&);

//////////// Member Function ////////////

/// Compute the maximum (absolute) cos value of between vector deris[i](t) 
/// and vector this->TP(i)(t) for i=0,1,2,3, where t is a common
/// parameter of the data point obtained from deris[i]'s knot vector.
///Function's return value is the max out of cosmax[.].
double get_perimeters_max_cos(
	const MGPvector<MGLBRep>& deris,     ///< the size must be 4.
	double taumax[4], ///< parameter on which the maximum value attains will be stored.
	double cosmax[4]  ///< the maximum value will be stored.
)const;

/// Compute the maximum (absolute) sin value of between vector srf.normal(uv(t))
/// and vector this->TP(i)(t) for i=0,1,2,3, where perim[i] is 
/// the same as srf.perimeter_curve(i), and t is a common parameter
/// of deris[i] and TP(i).
///Function's return value is the max out of sinmax[.].
double get_perimeters_max_sin(
	const MGSurface& srf,    /// surface which must corresponds to this object.
	double         taumax[4],/// parameters on which the maximum value attains will be stored.
	double         sinmax[4],/// the maximum value will be stored.
	bool*          eval=0	///indicates perimeters to evalate if eval!=null,
			///When eval[i] is true, perimeter i is evaluated for 0<=i<=3.
)const;

///Compute maximun abs(cons(theta)), where theta=angle of TP(i) and  corresponding 
///edge_crvl[i]'s start and end points' tangent vector.
double max_cos(
	const MGCurve*	perimeter[4]///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
)const;

///Compute maximun abs(cons(theta)), where theta=angle of TP(i) and  corresponding 
///edge_crvl[i]'s start and end points' tangent vector.
double max_cos(
	const MGPvector<MGLBRep>& perimeters///<境界線リスト(vmin,umax,vmax,uminの順,辺番号0,1,2,3の順)
)const;

///Return if i-th perimeter's TP specified(true) or not.
///i=0, 2 are v=min and max u-parameter line.
///i=1, 3 are u=max and min v-parameter line.
bool specified(size_t i) const{assert(i<4); return m_TP[i]!=0;};

///Set i-th perimeter's TP as a null, as an unspecified one.
void set_TP_null(size_t i);

///Set i-th perimeter's TP(copy version).
void set_TP(size_t i, const MGLBRep& tp);

///Set i-th perimeter's TP(auto_ptr version).
void set_TP(size_t i, std::auto_ptr<MGLBRep>& tp);

///Return i-th perimeter's TP.
const MGLBRep& TP(size_t i) const{assert(i<4);return *(m_TP[i]);}
MGLBRep& TP(size_t i) {assert(i<4);return *(m_TP[i]);}

MGLBRep** TP(){return m_TP;};

private:
//////////// Member Data ////////////
	MGLBRep* m_TP[4];	///<Tangent Plane will be stored.
				///<Tangent plane m_TP[i] is a line b-representation of
				///<(unit)normal vector along the i-th perimeter.
				///<Parameter range of the TP is the same as u or v parameter
				///<range of the corresponding surface representation.
	///< m_TP[0]: v=min boundary line, m_TP[1]: u=max boundary line
	///< m_TP[2]: v=max boundary line, m_TP[3]: u=min boundary line

};

///Construct a tangent plane LBRep of a curve wcrv. wcrv is a world coordinate
///curve that lies on a surface srf.
///TP_from_world_curve is for the users who do not know the parameter
///representation of wcrv of the surface srf. If you know the parameter
///representation of wcrv, you should use TP_from_parameter_curve.
///TP_from_parameter_curve's computation time is much smaller
///than TP_from_world_curve's.
///Let f(t) be wcrv whose parameter is t. Then the parameter of this function's
///output is the same as t. Let g(t) be the output, then g(t)=normal at srf(u,v).
///where srf(u,v)=f(t).
MGDECL MGLBRep TP_from_world_curve(
	const MGSurface& srf,	///<接続面
	const MGCurve& wcrv,		///<接続曲線
	size_t order,			///<作成するカーブのオーダ
	int& error);	///<エラーコード  0:OK, !=0 error code when constructing LBRep.

///Construct a tangent plane LBRep of a curve pcrv. pcrv is a parameter coordinate
///curve of a surface srf. wcrv is a target curve that lies on the surface srf, is
///world coordinate representation of pcrv. Output MGLBRep's parameter is 
///Let f(t) be wcrv whose parameter is t. Then the parameter of this function's
///output is the same as t. Let g(t) be the output, then g(t)=normal at srf(u,v).
///where srf(u,v)=f(t).
///pcrvは接続しようとしている曲面の(u,v)表現による２次元の曲線で、
///曲線が２曲線に分かれているような場合、接続して一つの曲線として、
///次のwcrvとおなじ範囲を表現する曲線に変換して入力します。
///wcrvは４辺面の１辺となる曲線で(x,y,z)の世界座標表現です。
MGDECL MGLBRep TP_from_parameter_curve(
	const MGSurface& srf,	///<接続面
	const MGCurve& pcrv,	///<接続曲線(接続面上のパラメータカーブ)
	const MGCurve& wcrv,	///<接続曲線(世界座標上のカーブ)
	size_t order,			///<作成するカーブのオーダ
	int& error);	///<エラーコード  0:OK, !=0 error code when constructing LBRep.

/** @} */ // end of GEORelated group

#endif
