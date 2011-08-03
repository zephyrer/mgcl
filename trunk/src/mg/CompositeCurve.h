/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGCompositeCurve_HH_
#define _MGCompositeCurve_HH_

#include "mg/Interval.h"
#include "mg/Position_list.h"
#include "mg/Straight.h"
#include "mg/RLBRep.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/SurfCurve.h"
#include "mg/BSumCurve.h"

#if defined(MGCL_DLL)
#include "mg/DequeProxy.h"
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS MGDequeProxy<MGCurve*>;
#pragma warning( pop )
#else
#include <deque>
#endif

class MGPlane;
class MGIgesDirectoryEntry;
class MGSurfCurve;
class MGStraight;
class MGCompositeCurve;

/** @addtogroup GEO
 *  @{
 */

///MGCompositeCurve is a composite of other leaf curves.
///Assumedly they are connected as C0 continuity. However, MGCompositeCurve
///does not check their continuity, but only put first or last as the user says
/// (in connect_to_end or in connect_to_start).
///Parameter ranges of the member curves are always continuous. param_s() of the
///1st curve to param_e() of the last form MGCompositeCurve's paramter range.
///number_of_curves() indicates the number of leaf curves(in other words,
///element curves to construct the MGCompositeCurve).
class MGCLASS MGCompositeCurve:public MGCurve{

public:
#if defined(MGCL_DLL)
	typedef MGDequeProxy<MGCurve*> container_type;
#else
	typedef std::deque<MGCurve*> container_type;
#endif

typedef container_type::iterator iterator;
typedef container_type::const_iterator const_iterator;
typedef container_type::reverse_iterator reverse_iterator;
typedef container_type::const_reverse_iterator const_reverse_iterator;

MGDECL friend MGCompositeCurve operator+ (const MGVector& v, const MGCompositeCurve& lb);
MGDECL friend MGCompositeCurve operator* (double scale, const MGCompositeCurve&);

//////////////Constructor///////////////

///Void constructor(�������Ȃ��ŃI�u�W�F�N�g���쐬����B)
MGCompositeCurve(){;};				///Constructor of no curves.

///Constructor of one curve.
///crv is a newed object pointer and MGCompositeCurve takes the ownership.
explicit MGCompositeCurve(MGCurve* crv);

///Copy constructor.
MGCompositeCurve(const MGCompositeCurve&);

//////////// Destructor ////////////
~MGCompositeCurve();

//////////// Operator overload(���Z�q���d��`) ////////////

///Assignment.
///When the leaf object of this and crv2 are not equal, this assignment
///does nothing.
MGCompositeCurve& operator=(const MGGel& gel2);
MGCompositeCurve& operator=(const MGCompositeCurve& original);

///Transformation object construction
MGCompositeCurve operator+ (const MGVector& v) const;
MGCompositeCurve operator- (const MGVector& v) const;
MGCompositeCurve operator* (double scale) const;
MGCompositeCurve operator* (const MGMatrix& mat) const;
MGCompositeCurve operator* (const MGTransf& tr) const;

///Object transformation.
MGCompositeCurve& operator+=(const MGVector& v);
MGCompositeCurve& operator-=(const MGVector& v);
MGCompositeCurve& operator*=(double scale);
MGCompositeCurve& operator*=(const MGMatrix& mat);
MGCompositeCurve& operator*=(const MGTransf& tr);

///Comparison of two curves.
bool is_same_curve(const MGCurve& curve2)const;
bool operator==(const MGCompositeCurve& gel2)const;
bool operator==(const MGGel& gel2)const;
bool operator==(const MGTrimmedCurve& gel2)const;
bool operator<(const MGCompositeCurve& gel2)const;
bool operator<(const MGGel& gel2)const;

//////////// Member Function ////////////

///Returns B-Rep Dimension.
///bdim of MGCompositeCurve is the sum of each member curves.
size_t bdim() const;

const_iterator begin()const{return m_composite.begin();};
iterator begin(){return m_composite.begin();};

///Return minimum box that includes the curve of parameter interval.
/// ���͂̃p�����[�^�͈͂̋Ȑ��������͂ރ{�b�N�X��Ԃ��B
MGBox box_limitted(
	const MGInterval& ///< Parameter Range of the curve.
) const;

///Changing this object's space dimension.
MGCompositeCurve& change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new object.
	size_t start2=0 		///< Source order of this object.
);

///Change parameter range, be able to change the direction by providing
///t1 greater than t2.
void change_range(
	double t0,	///<Parameter value for the start of original. 
	double t1	///<Parameter value for the end of original. 
);

///Compute the closest point parameter value of this curve from a point.
double closest(const MGPosition& point) const;

///�֐����F
///		common
///�ړI�F
///		�^����ꂽ�Ȑ��Ǝ��g�̌�_�������͋��ʕ��������邩�ǂ������ׂ�B
///�����F
///		const MGCurve&			crv2,		(I/ )	�^������Ȑ�
///		std::vector<double>&	vec_param	( /O)	���ʕ����̃p�����[�^�͈�
///		MGCCisect_list&			isect		( /O)	��_
///				 4n�̔z��ŁAt(4*i+0),t(4*i+1)�����g�̃p�����[�^�͈�(t(4*i+0) < t(4*i+1))�A
///							 t(4*i+2),t(4*i+3)���^�Ȑ��̃p�����[�^�͈�(f(t(4*i+0))=f(t(4*i+2))
///�߂�l�F
///		3:��_�����ʕ��������܂���
///		2:��_�݂̂����܂���
///		1:���ʕ����݂̂����܂���
///		0:��_�����ʕ������Ȃ�����
///		-1:���ʃG�b�W�̎����v�Z�G���[
///		-2:���ʃG�b�W���S�ȏ㋁�܂���(�̂��Ă��Ȃ��ƌ��Ȃ�)
///�ǋL�F
///	�Ȑ������ʂ��ǂ����̌덷�ɂ�line_zero()�A���p�����[�^�͈͂̎����v�Z��
///	�덷�ɂ́A�p�����[�^�͈�*rc_zero()���g�p����
int common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan,
	MGCCisect_list& isect
)const;

///�֐����F
///		common
///�ړI�F
///		�^����ꂽ�Ȑ��Ǝ��g�̋��ʕ��������邩�ǂ������ׂ�B
///�����F
///		const MGCurve&			crv2,		(I/ )	�^������Ȑ�
///		std::vector<double>&	vec_param	( /O)	���ʕ����̃p�����[�^�͈�
///				 4n�̔z��ŁAt(4*i+0),t(4*i+1)�����g�̃p�����[�^�͈�(t(4*i+0) < t(4*i+1))�A
///							 t(4*i+2),t(4*i+3)���^�Ȑ��̃p�����[�^�͈�(f(t(4*i+0))=f(t(4*i+2))
///�߂�l�F
///		���ʕ����̐�:	���ʕ��������܂���
///		0:				���ʕ������Ȃ�����
///		-1:				���ʃG�b�W�̎����v�Z�G���[
///		-2:				���ʃG�b�W���S�ȏ㋁�܂���(�̂��Ă��Ȃ��ƌ��Ȃ�)
///�ǋL�F
///	�Ȑ������ʂ��ǂ����̌덷�ɂ�line_zero()���A�p�����[�^�͈͂̎����v�Z�̌덷�ɂ́A
///  �p�����[�^�͈�*rc_zero()���g�p����
int common(
	const MGCurve& crv2,
	std::vector<double>& vecComSpan
)const;

///Connect the input curve to the end(to_end) or start(to_start) of this curve.
///(1) End(start) point of this curve is assumedly the same as the start(end) point
///of add_curve. However, connect_to_xxx does not check the continuity except when
///connecting two curves are both MGLBRep.
///(2) add_curve must be a newed object pointer and connect_to_end takes the ownership.
///(3) add_curve can be a MGCompositeCurve. In this case each curves in add_curve
///become members of this MGCompositeCurve.
///(4) When add_curve is a SurfCurve, or a TrimmedCurve, it is changed to non SurfCurve
///or non TrimmedCurve. Thus the original surface or curve can be deleted or modified.
///(5) When MGCompositeCurve was not an empty curve, the original part of
///the MGCompositeCurve's parameter range will not be changed by this connect_to_end.
///Instead add_curve's parameter range will be so modified that the magnitude of
///1st derivatives of the two curves are the same at the connecting point and
///the modified add_curve's parameter range is continuous to the original.
///(6) connect_to_end(start) will change add_curve's direction if necessary.
///Function's return value is the new parameter range of add_curve after added.
///connect() connects to either start or end of this depending the distance
///of the start or end points of this curve to the start or end points of add_curve.
///connect_to_end() connects either start or end points of add_curve to the end to this curve.
///connect_to_start() connects either start or end points of add_curve to the start to this curve.
MGInterval connect(MGCurve* add_curve);
MGInterval connect_to_end(MGCurve* add_curve);
MGInterval connect_to_start(MGCurve* add_curve);

///Exchange ordering of the coordinates.
///Exchange coordinates (i) and (j).
MGCompositeCurve& coordinate_exchange(size_t i, size_t j);

///Construct new curve object by copying to newed area.
///User must delete this copied object by "delete".
MGCompositeCurve* clone() const;

///copy as a newed curve. The new curve will be MGLBRep or MGRLBRep.
///When original curve was a MGRLBRep, the new curve will be a MGRLBRep.
///Otherwise,  the new curve will be a MGLBRep.
///Returned object must be deleted.
MGCurve* copy_as_nurbs() const;

///Construct new curve object by changing
///the original object's space dimension.
///User must delete this copied object by "delete".
MGCompositeCurve* copy_change_dimension(
	size_t sdim,			///< new space dimension
	size_t start1=0, 		///< Destination order of new line.
	size_t start2=0		///< Source order of this line.
)const;

///Get the i-th curve in this CompositeCurve.
const MGCurve& curve(size_t i)const{return *(m_composite[i]);};
MGCurve& curve(size_t i){return *(m_composite[i]);};

///Compute curvilinear integral of the 1st two coordinates from parameter t0
///to t1.
///This integral can be used to compute area sorrounded by the curve.
///The sum of all the member curve's curvilinear_integral.
double curvilinear_integral(double t1, double t2) const;

///Compute curvilinear integral of the 1st two coordinates.
///(All the parameter range of the curve.)
///This integral can be used to compute area sorrounded by the curve.
///The sum of all the member curve's curvilinear_integral.
double curvilinear_integral() const;

void display_break_points()const;
void display_control_polygon()const;
void display_curvatures(
	double	scale,	///<scaling of the graph.
	int		density,///<densitiy of the graph.
	bool	use_radius///<true:radius display, false:curvature display.
)const;

///Divide this curve at the designated knot multiplicity point.
///Function's return value is the number of the curves after divided.
int divide_multi(
	MGPvector<MGCurve>& crv_list,	///<divided curves will be appended.
	int multiplicity=-1	///<designates the multiplicity of the knot to divide at.
						///<When multiplicity<=0, order()-1 is assumed.
						///<When multiplicity>=order(), order() is assumed.
) const;

void drawSE(
	double span_length,	///<Line segment span length.
	double t0,			///<Start parameter value of the curve.
	double t1			///<End parameter value of the curve.
						///<Draw will be performed from t0 to t1.
)const;

const_iterator end()const{return m_composite.end();};
iterator end(){return m_composite.end();};

/// Evaluate n'th derivative data. n=0 means positional data evaluation.
MGVector eval(
	double t,			///< Parameter value.
	size_t nderiv=0,	///< Order of Derivative.
	int left=0			///<Left continuous(left=true)
						///<or right continuous(left=false).
)const;

///Extrapolate this curve by an (approximate) chord length.
///The extrapolation is C2 continuous.
void extend(
	double length,	///<approximate chord length to extend. 
	bool start=false///<Flag of which point to extend, start or end point of the line.
					///<If start is true extend on the start point.
);

///Find which curve parameter range input t belongs to.
///Function's return value m is the id of m_composite which satisfy
///m_composite[m]<=t<m_composite[m+1];
///m_composite.size() should be >=1, or find will abort.
size_t find(double t, bool right_continuous=true) const;

/// Return This object's typeID
long identify_type() const;

///Test if this cure is planar or not.
///MGPlane expression will be output to plane if this is planar.
///Function's return value is true if planar.
bool is_planar(MGPlane& plane)const;

///Intersection of CompositeCurve and another curve crv2.
MGCCisect_list isect(const MGCurve& crv2) const;
MGCCisect_list isect(const MGStraight& curve2)const;
MGCCisect_list isect(const MGRLBRep& curve2)const;
MGCCisect_list isect(const MGEllipse& curve2)const;
MGCCisect_list isect(const MGLBRep& curve2)const;
MGCCisect_list isect(const MGSurfCurve& curve2)const;
MGCCisect_list isect(const MGBSumCurve& curve2)const;
MGCCisect_list isect(const MGTrimmedCurve& curve2)const;
MGCCisect_list isect(const MGCompositeCurve& curve2)const;

///Intersection with a Surface
MGCSisect_list isect(const MGSurface& surf) const;
MGCSisect_list isect(const MGPlane& surf) const;
MGCSisect_list isect(const MGSphere& surf)const;
MGCSisect_list isect(const MGCylinder& surf)const;
MGCSisect_list isect(const MGSBRep& surf)const;
MGCSisect_list isect(const MGRSBRep& surf)const;
MGCSisect_list isect(const MGBSumSurf& surf)const;

///Access to i-th element of knot.
double knot(size_t i) const;

///Returns the knot vector of the curve.
///This should not be used.
const MGKnotVector& knot_vector() const;

///Cmpute curve length of the interval.
///If t1 is greater than t2, return negative value.
/// �p�����[�^�������ŗ^����ꂽ�Ƃ��͐��l�A�~���̂Ƃ��͕��l��Ԃ��B
///The sum of all the member curve's length.
double length(double t1, double t2) const;
double length()const{return MGCurve::length();}

///Update this by limiting the parameter range of the curve.
/// ���g�Ɏw�肵���p�����[�^�͈͂̂���������������B
MGCompositeCurve& limit(const MGInterval& rng);

///Negate the curve direction(�Ȑ��̕����𔽓]����)
void negate();

///Obtain parameter value if this curve is negated by "negate()".
double negate_param(double t)const;

///Get the number of included curves in this CompositeCurve.
size_t number_of_curves()const{return m_composite.size();};

//****Offset of CompositeCurve is CompositeCurve of the offsets of each
//member curves********

///���I�t�Z�b�g�֐�
///�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
///�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
///�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
///costant offset curve. if the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
MGPvector<MGCurve> offset(
	double ofs_value,							///<�I�t�Z�b�g��
	const MGVector& norm_vector = mgNULL_VEC	///<�@���x�N�g��
) const;

///�σI�t�Z�b�g�֐�
///�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
///�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
///�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
///�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����X�g���ԋp�����B
///valuable offset curve. if the norm_vector is given, the positive offset direction decide
///to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
MGPvector<MGCurve> offset(
	const MGLBRep& ofs_value_lb,					///<��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector = mgNULL_VEC	///<�@���x�N�g��
)const;

///C2�A���Ȑ��̈��I�t�Z�b�g�֐�
///�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
///�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
///�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
///costant offset curve of C2 continuous curve. if the norm_vector is given, the positive offset direction
///decide to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
MGLBRep offset_c2(
	double ofs_value,								///<�I�t�Z�b�g��
	const MGVector& norm_vector = mgNULL_VEC	///<�@���x�N�g��
) const;

///C2�A���Ȑ��̉σI�t�Z�b�g�֐�
///�I�t�Z�b�g�ʂ͋�Ԏ���1�̐�B�\���ŗ^������B
///�I�t�Z�b�g�����́A�@���������猩�ē��͋Ȑ��̐i�s���������𐳂Ƃ���B
///�@���x�N�g�����k���̏ꍇ�A�n�_�ɂ����ċȗ����S�����𐳂Ƃ���B�������A�ȗ����S�֋ȗ����a�ȏ�̃I�t�Z�b�g
///�͍s��Ȃ��B�g�������X��line_zero()���g�p���Ă���B�߂�l�́A�I�t�Z�b�g�Ȑ����ԋp�����B
///valuable offset curveof C2 continuous curve. if the norm_vector is given, the positive offset direction
///decide to left hand side from ahead, or the direction to center of curvature at start parameter.
///the offset value is less than radius of curvature. line_zero() is used.
MGLBRep offset_c2(
	const MGLBRep& ofs_value_lb,					///<��Ԏ����P�̐�B�\���Ŏ������I�t�Z�b�g��
	const MGVector& norm_vector = mgNULL_VEC	///<�@���x�N�g��
) const;

///Test if given point is on the curve or not. If yes, return parameter
///value of the curve. Even if not, return nearest point's parameter.
/// �w��_�����g��ɂ��邩�𒲂ׂ�B�Ȑ���ɂ���΁C���̃p�����[�^�[�l���C
/// �Ȃ��Ă��ŋߖT�_�̃p�����[�^�l��Ԃ��B
/// Function's return value is >0 if the point is on the curve,
/// and 0 if the point is not on the curve.
bool on(
	const MGPosition&,	///<Point(�w��_)
	double&				///<Parameter of the curve(�p�����[�^)
)const;

///Returns the order.
///Returns the maximum order among the curves.
unsigned order() const;

/// Return ending parameter value.
double param_e() const;

///Normalize parameter value t to the nearest knot if their distance is
///within tolerance.
double param_normalize(double t) const;

/// Return starting parameter value.
double param_s() const;
	
///Compute part of this curve from parameter t0 to t1.
///Returned is the pointer to newed object, and so should be deleted
///by calling program, or memory leaked.
MGCurve* part(
	double t0, double t1,
	int multiple=0	///<Indicates if start and end knot multiplicities
					///<are necessary. =0:unnecessary, !=0:necessary.
) const;

///Return perpendicular point from a point P,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
int perp_guess(
	double t0, double t1,	///<parameter range of this.
							///<(t0>=t1) indicates no range specified.
	const MGPosition& P,	///<Point(�w��_)
	double tg,				///<Guess parameter values of this curve.
	double& t				///<Output parameter
) const;

///Return perpendicular points of two curves,
///given guess starting paramter values.
///Function's return value is:
///   perp_guess=true if perpendicular points obtained,
///   perp_guess=false if perpendicular points not obtained,
int perp_guess(
	double s0, double s1,		///<parameter range of this.
					///<When s0>=s1, no limit for this parameter range.
	const MGCurve& curve2,		///<2nd curve.
	double t0, double t1,		///<parameter range of curve2.
					///<When t0>=t1, no limit for curve2 parameter range.
	double sg, double tg,		///<Guess parameter values of the two curves
	///<sg: this curve's parameter, tg:curve2's parameter.
	MGPosition& st				///<perpendicular points' parameter values
								///<will be output.
	///<st(0): this curve's parameter, st(1):curve2's parameter.
) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps(const MGCurve& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGStraight& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGRLBRep& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGEllipse& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGLBRep& crv2)const{return perps_1by1(crv2);};
MGPosition_list perps(const MGSurfCurve& crv2)const;
MGPosition_list perps(const MGBSumCurve& crv2)const{return perps_1by1(crv2);};

///Compute all foot points of the perpendicular line from point to
///the curve.
/// �^�|�C���g����Ȑ��։��낵�������̑��́C�Ȑ��̃p�����[�^�l��
/// ���ׂċ��߂�B
MGCParam_list perps(
	const MGPosition& P		///Point(�w��_)
) const;

///Obtain the projected curve of this onto the surface srf.
///The direction of the projection is along the vector vec if the vec is not
///NULL, and normal to the surface if the vec is NULL.
///Output of 'project' is two kind of curves:
///one is general world coordinate curves('vec_crv'), and the other is (u,v) curves of
///the parameter space of the surfaces(vec_crv_uv).
///vec_crv_uv.size() is equal to vec_crv.size(). Let the size be n, then
/// (vec_crv_uv[i], vec_crv[i]) is one pair for 0<=i<n.
int project_onto_surface(
	const MGFSurface& srf,
	MGPvector<MGCurve>& vec_crv_uv,	
		///<Projected curve(surface parameter (u,v) representation) will be appended.
	MGPvector<MGCurve>& vec_crv,
		///<Projected curve(world coordinate(x,y,z) representation) will be appended.
	const MGVector& vec
)const;

///Release the pointer of the last curve.
///Returned will be the released MGCurve pointer.
MGCurve* release_back();

///Release the pointer of the 1st curve.
///Returned will be the released MGCurve pointer.
MGCurve* release_front();

///�m�b�g�폜�֐�(B�\���Ȑ��̂�)
///�g�������X��line_zero���g�p����B���̃m�b�g���ׂ������̂قǍ폜���₷��
///Remove redundant knot, and reduce the b-rep dimension.
///The tolerance used is MGTolerance::line_zero().
void remove_knot();

///Return the space dimension. It is the maximum space dimension of the
///member curves.
size_t sdim() const;

///Return sweep surface from crv
///Returned is a newed MGSurface, must be deleted.
///The sweep surface is defined as:
///This curve(say c(t)) is the rail and the straight line segments from
///C(t)+start_dist*uvec to C(t)+end_dist*uvec are the generatrix.
MGSurface* sweep(
	const MGUnit_vector& uvec,	///<Sweep Direction.
	double start_dist,			///<distance to start edge.
	double end_dist				///<distance to end edge.
) const;	

///Return curve type(�Ȑ��̃^�C�v��Ԃ�)
MGCURVE_TYPE type() const{return MGCURVE_COMPOSITE;};

///Unlimit parameter range of the curve(limit���͂���)
MGCurve& unlimit();

///Unlimit parameter range of the curve to the end point direction
///(�I�_������limit���͂���)
MGCurve& unlimit_end();

///Unlimit parameter range of the curve to the start point direction
///(�n�_������limit���͂���)
MGCurve& unlimit_start();

///Output function.
std::ostream& out(std::ostream&) const;

///Output to IGES stream file.
int out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch=0
)const;

protected:

///Compute intersection point of 1D sub curve of original curve.
///Parameter values of intersection point will be returned.
MGCParam_list intersect_1D(						
	double f,			///< Coordinate value
	size_t coordinate=0	///< Coordinate kind of the data f(from 0).
) const;	

///Obtain so transformed 1D curve expression of this curve that
///f(t)={sum(xi(t)*g[i]) for i=0(x), 1(y), 2(z)}-g[3], where f(t) is the output
///of oneD and xi(t) is i-th coordinate expression of this curve.
///This is used to compute intersections with a plane g[4].
std::auto_ptr<MGCurve> oneD(
	const double g[4]			///<Plane expression(a,b,c,d) where ax+by+cz=d.
) const;

///�����o�f�[�^��ǂݏo���֐�
/// �߂�lbool�͐���ɓǂݏo�����o�����true�A���s�����false�ɂȂ�
/// �����ł͏����ΏۂƂȂ�f�[�^�����o�������̂ŉ������������Ȃ��B
void ReadMembers(MGIfstream& buf);

///�����o�f�[�^���������ފ֐�
/// �߂�lbool�͐���ɏ������݂��o�����true�A���s�����false�ɂȂ�
/// �����ł͏����ΏۂƂȂ�f�[�^�����o�������̂ŉ������������Ȃ��B
void WriteMembers(MGOfstream& buf) const;

private:

	container_type m_composite;	///<m_composite consists of newed MGCurve pointer.

///Return minimum box that includes whole of the curve.
///�Ȑ��������͂ރ{�b�N�X��Ԃ��B
///Or'ed box of all the member curves.
MGBox* compute_box() const;

///Provide divide number of curve span for function intersect.
///****This is declared here to prohibit the use****.
size_t intersect_dnum() const;

///Intersection of CompositeCurve and another curve crv2.
MGCCisect_list isect_1by1(const MGCurve& crv2) const;

///Intersection with a surface. Obtains the isects by computing isects
///of each elemet curve.
MGCSisect_list isect_1by1(const MGSurface& surf) const;

///Compute all the perpendicular points of this curve and the second one.
///That is, if f(s) and g(t) are the points of the two curves f and g,
///then obtains points where the following conditions are satisfied:
///  fs*(f-g)=0.    gt*(g-f)=0.
///Here fs and gt are 1st derivatives at s and t of f and g.
///MGPosition P in the MGPosition_list contains this and crv's parameter
///as:     P(0)=this curve's parameter, P(1)=crv2's parameter value.
MGPosition_list perps_1by1(
	const MGCurve& crv2		///<The second curve
)const;

friend class MGHHisect;

};

/** @} */ // end of GEO group
#endif
