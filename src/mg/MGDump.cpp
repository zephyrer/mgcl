/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Unit_vector.h"
#include "mg/Position_list.h"
#include "mg/KnotVector.h"
#include "mg/KnotArray.h"
#include "mg/BPointSeq.h"
#include "mg/Transf.h"
#include "mg/Geometry.h"
#include "mg/Point.h"
#include "mg/Curve.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/LBRepEndC.h"
#include "mg/SPointSeq.h"
#include "mg/Surface.h"
#include "mg/Plane.h"
#include "mg/SBRep.h"
#include "mg/OscuCircle.h"
#include "mg/OscuCircleData.h"
#include "mg/RLBRep.h"
#include "mg/RSBRep.h"
#include "mg/PPRep.h"

#include "mg/Ifstream.h"
#include "mg/Ofstream.h"
#include "mg/TrimmedCurve.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// ***MGEReal***
//Dumped size calculation
size_t MGEReal::dump_size()const{
	return sizeof m_value;
}
// Dump Function
int MGEReal::dump(MGOfstream& buf)const{
	int len = -1;
	try{
		buf<<m_value;
		len = dump_size();
	}
	catch(...){
		throw;
	}

	return len;
}

// Restore & Create New Object
int MGEReal::restore(MGIfstream& buf){
	int len = -1;
	try{
		double value;
		buf>>value;
		m_value = value;
		len = dump_size();
	}
	catch (...){
		throw;
	}

	return len;
}

// ***MGInterval***
// Dumped Size calculation
size_t MGInterval::dump_size() const
{
	size_t len = 0;
	len = m_high.dump_size() + m_low.dump_size();
	return len;
}

//Dump Function
int MGInterval::dump(MGOfstream& buf)const{
	int len = -1;
	try
	{
		m_high.dump(buf);
		m_low.dump(buf);
		len = dump_size();
	}
	catch(...)
	{
		throw;
	}
	return len;
}

//Restore & Create Object
int MGInterval::restore(MGIfstream& buf){
	int len = -1;
	try
	{
		m_high.restore(buf);
		m_low.restore(buf);
		dump_size();
	}
	catch (...){
		throw;
	}
	return len;
}

// ***MGBox***
//Calculate dump size
size_t MGBox::dump_size() const{
	size_t len = 0;
	len = (sizeof sdim());
	for (size_t i=0; i<sdim(); i++)	{
		len += (m_range[i].dump_size());
	}
	return len;
}

// Dump Functions
int MGBox::dump(MGOfstream & buf) const{
	int len = -1;
	try	{
		size_t dim = sdim();
		buf << dim;
		for(size_t i=0; i<sdim(); i++)		{
			m_range[i].dump(buf);
		}
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

int MGBox::restore(MGIfstream& buf){
	int len = -1;
	try{
		int dim;
		buf >> dim;
		resize(dim);
		for (int i = 0; i < dim ; i++){
			m_range[i].restore(buf);
		}
		len = dump_size();
	}
	catch (...){
		throw;
	}
	return len;
}

// ***NDDArray***
//Calculate dump size
size_t MGNDDArray::dump_size() const{
	size_t len = 0;
	len = sizeof m_length;
	len += (sizeof ref(0))*length(); // ref(i)は全てdoubleなので。
	return len;
}
//Dump Function
// length() == size() となるようにサイズを変更してダンプする。
int MGNDDArray::dump(MGOfstream & buf) const{
	int len = -1;
	try	{
		buf << m_length;
		for(size_t i=0; i<m_length; i++)
			buf<<m_element[i];
		len = dump_size();
	}
	catch (...)	{
		throw;
	}
	return len;
}

//Restore Function
int MGNDDArray::restore(MGIfstream& buf){
	int len = -1;
	try	{
		buf >> len;
		resize(len);
		for(int i=0; i<len; i++)
			buf>>m_element[i];
		len = dump_size();
	}
	catch(...){
		throw;
	}
	return len;
}

// ***MGVector***
//Calculate dump size
size_t MGVector::dump_size() const{//Dumped Length;
// !!! In stead of to USE STL Vector, BUT NOT to USE Iterator
	size_t len = 0;
	len = (sizeof sdim()) + (sizeof ref(0)) * sdim();
	return len;
}
//Dump Function
// Dumped Vector Dimension and each element
int MGVector::dump(MGOfstream& buf) const{
	int len = -1;
	try{
		int dim = sdim();
		buf<<dim;
		for(int i=0; i<dim; i++)
			buf << m_element[i];
		len = dump_size();
	}
	catch(...){
		throw;
	}
	return len;
}

//Restore Function
int MGVector::restore(MGIfstream& buf){  // 0:OK !0:NG -99:file not open
	int len = -1;
	try{
		int dim;
		buf>>dim;
		resize(dim);
		for(int i=0; i<dim; i++)
			buf>>m_element[i];
		len = dump_size();
	}
	catch (...)	{
		throw;
	}
	return len;
}

// ***MGMatrix***
//Calculate dump size
size_t MGMatrix::dump_size() const{
// !!! In stead of to USE STL Vector, BUT NOT to USE Iterator
	size_t len = 0;
	len = (sizeof sdim()) + (sizeof ref(0,0)) * sdim() * sdim();
	return len;
}

//Dump Function
int MGMatrix::dump(MGOfstream& buf) const{
	int len = 0;
	try{
		size_t dim = sdim();
		size_t mtx_len = dim*dim;
		buf << dim;
		for(size_t i=0; i<mtx_len; i++)
			buf<<m_matrix[i];
		len = dump_size();
	}
	catch (...){
		throw;
	}
	return len;
}

//Restore Function
int MGMatrix::restore(MGIfstream& buf){
	int len = -1;
	try{
		size_t sdim;
		buf >> sdim;
		resize(sdim);
		size_t mtx_len=sdim*sdim;
		for(size_t i=0; i<mtx_len; i++)
			buf>>m_matrix[i];
		len = dump_size();
	}
	catch (...)	{
		throw;
	}
	return len;
}

// ***MGKnot***
//Calculate dump size
size_t MGKnot::dump_size() const{//Dumped Length;

	size_t len;
	len = (sizeof m_value) + (sizeof m_multiplicity);
	return len;
}

//Dump Function
int MGKnot::dump(MGOfstream& buf) const{ // 0:OK !0:NG -99: file not open

	int len = -1;
	try
	{
		buf << m_value << m_multiplicity;
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGKnot::restore(MGIfstream& buf){
	int len = -1;
	try
	{
		buf >> m_value >> m_multiplicity;
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

// ***MGKnotArray***
//Calculate dump size
size_t MGKnotArray::dump_size() const{//Dumped Length;

	size_t len;
	len = (sizeof size()) + (front().dump_size()) * size();
	return len;
}

//Dump Function
int MGKnotArray::dump(MGOfstream& buf) const{
	int len = -1;
	try{
		std::vector<MGKnot>::const_iterator Citr = begin();
		size_t dim = size();
		buf << dim;
		for(size_t i=0; i<dim; i++,Citr++){
			(*Citr).dump(buf);
		}
		len = dump_size();
	}
	catch (...){
		throw;
	}
	return len;
}

//Restore Function
int MGKnotArray::restore(MGIfstream& buf){
	int len = -1;
	try{
		size_t size;
		buf>>size;
		for(size_t i = 0; i < size; i++){
			MGKnot knot;
			knot.restore(buf);
			push_back(knot);
		}
		len=dump_size();
	}
	catch(...){
		throw;
	}
	return len;
}

// ***MGKnotVector*** (class MGKnotVector: public MGNDDArray)
//Calculate dump size
size_t MGKnotVector::dump_size() const{
	size_t len;
	len = MGNDDArray::dump_size();
	len += sizeof m_order;
	return len;
}

//Dump Function
int MGKnotVector::dump(MGOfstream& buf) const{
	int len = -1;
	try
	{
		MGNDDArray::dump(buf);
		buf << m_order;
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGKnotVector::restore(MGIfstream& buf){
	int len = -1;
	try
	{
		MGNDDArray::restore(buf);
		buf >> m_order; if(m_order) m_current=m_order-1;
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

// ***MGBPointSeq***
//Calculate dump size
size_t MGBPointSeq::dump_size() const{
	size_t len = (sizeof m_sdim) + (sizeof m_length)
			+ (sizeof m_bpoint[0]) * m_sdim * m_length;
	return len;
}

//Dump Function
int MGBPointSeq::dump(MGOfstream& buf) const{
	int len = -1;
	try
	{
		buf << m_sdim << m_length;
		for (size_t j = 0; j < m_sdim; j++)
			for (size_t i = 0; i < m_length; i++)
			{
				double val = ref(i,j);
				buf << val;
			}
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGBPointSeq::restore(MGIfstream& buf){
	int len = -1;
	try
	{
		size_t sdim, len;
		buf >> sdim >> len;
		resize(len, sdim);
		for (size_t j = 0; j < sdim; j++)
			for (size_t i = 0; i <len; i++)
			{
				double val = 0.;
				buf >> val;
				(*this)(i,j) = val;
			}
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

// ***MGTransf***
//Calculate dump size
size_t MGTransf::dump_size() const{
	size_t len;
	len = m_affine.dump_size() + m_translation.dump_size();
	return len;
}

//Dump Function
int MGTransf::dump(MGOfstream& buf) const{
	int len = -1;
	try
	{
		m_affine.dump(buf);
		m_translation.dump(buf);
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGTransf::restore(MGIfstream& buf){
	int len = -1;
	try
	{
		m_affine.restore(buf);
		m_translation.restore(buf);
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}
// ***MGSPointSeq***
// !! In stead of To USE STL Vector, BUT NOT to USE iterator
//Calculate dump size
size_t MGSPointSeq::dump_size() const{
	size_t len;
	len = (sizeof m_lengthu) + (sizeof m_lengthv) + (sizeof m_sdim);
	len += (sizeof m_spoint[0]) * (m_lengthu * m_lengthv * m_sdim);
	return len;
}

//Dump Function
int MGSPointSeq::dump(MGOfstream& buf) const{
	int len = -1;
	try{
		buf << m_lengthu << m_lengthv << m_sdim;
		for(size_t k = 0; k < m_sdim; k++)
			for(size_t j = 0; j < m_lengthv; j++)
				for(size_t i = 0; i < m_lengthu; i++){
					double val = ref(i,j,k);
					buf << val;
				}

		len = dump_size();
	}
	catch (...)	{
		throw;
	}
	return len;
}

//Restore Function
int MGSPointSeq::restore(MGIfstream& buf){
	int len = -1;
	try{
		size_t lenu, lenv, dim;
		buf >> lenu >> lenv >> dim;
		resize(lenu, lenv, dim);
		size_t vec_len=lenu*lenv*dim;
		for (size_t i=0; i<vec_len; i++) buf>>m_spoint[i];
		len = dump_size();
	}
	catch (...)	{
		throw;
	}
	return len;
}
// ***MGLBRepEndC***
//Calculate dump size
size_t MGLBRepEndC::dump_size() const{
	size_t len;
	len = (sizeof m_cond) + m_1deriv.dump_size() + m_2deriv.dump_size();
	return len;
}

//Dump Function
int MGLBRepEndC::dump(MGOfstream& buf) const{
	int len = -1;
	try
	{
		buf << m_cond;
		m_1deriv.dump(buf);
		m_2deriv.dump(buf);
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGLBRepEndC::restore(MGIfstream& buf){
	int len = 0;
	try
	{
		buf >> (int&)m_cond;
		m_1deriv.restore(buf);
		m_2deriv.restore(buf);
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

/*// ***MGSBRepEndC***
//Calculate dump size
size_t MGSBRepEndC::dump_size() const
{
	size_t len = (sizeof m_sdim) + (sizeof m_nu) + (sizeof m_nv);
	for (int i = 0; i < 4; i++)
	//MGVector, MGPointSeqは可変長なので
	{
		len += (sizeof m_cond[i])
			+ m_1deriv[i].dump_size() + m_2deriv[i].dump_size()
			+ m_11d[i].dump_size() + m_12d[i].dump_size()
			+ m_21d[i].dump_size() + m_22d[i].dump_size();
	}
	return len;
}

//Dump Function
int MGSBRepEndC::dump(MGOfstream& buf) const
{
	int len = -1;
	try
	{
		//大きさの変わらないものを先に入れる
		buf << m_sdim << m_nu << m_nv;
		for (int i = 0; i < 4; i++)
		{
			buf << m_cond[i];
			m_1deriv[i].dump(buf);
			m_2deriv[i].dump(buf);
			m_11d[i].dump(buf);
			m_12d[i].dump(buf);
			m_21d[i].dump(buf);
			m_22d[i].dump(buf);
		}
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGSBRepEndC::restore(MGIfstream& buf)
{
	int len = -1;
	try
	{
		//大きさの変わらないものを先に入れる
		buf >> m_sdim >> m_nu >> m_nv;
		for (int i = 0; i < 4; i++)
		{
			buf >> (int&)m_cond[i];
			m_1deriv[i].restore(buf);
			m_2deriv[i].restore(buf);
			m_11d[i].restore(buf);
			m_12d[i].restore(buf);
			m_21d[i].restore(buf);
			m_22d[i].restore(buf);
		}
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

// ***MGSBRepTP***
//Calculate dump size
size_t MGSBRepTP::dump_size() const
{
	size_t len = 0;
	for (int i = 0; i < 4; i++)
	{
		const MGLBRep& tp = *(m_TP[i]);
		len += (tp.knot_vector().dump_size() + tp.line_bcoef().dump_size()
			+ 8);
			// ヘッダ部の８をプラスする。
	}

	return len;
}

//Dump Function
int MGSBRepTP::dump(MGOfstream& buf) const
{
	int len = -1;
	try
	{
		for(int i = 0; i < 4; i++) *(m_TP[i]).Write(buf);
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGSBRepTP::restore(MGIfstream& buf)
{
	int len = -1;
	try
	{
		for(int i = 0; i < 4; i++)
		{
			*(m_TP[i].Read(buf);
		}
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}
*/

// ***MGTolerance***
//Calculate dump size
size_t MGTolerance::dump_size() const{
	size_t len =
//			(sizeof m_mach_zero) + (sizeof m_wc_zero)
			(sizeof mach_zero()) + (sizeof wc_zero())
//			+ (sizeof m_wc_zero_sqr) + (sizeof m_rc_zero)
			+ (sizeof wc_zero_sqr()) + (sizeof rc_zero())
//			+ (sizeof m_rc_zero_sqr) + (sizeof m_angle_zero)
			+ (sizeof rc_zero_sqr()) + (sizeof angle_zero())
//			+ (sizeof m_line_zero) + (sizeof m_max_knot_ratio);
			+ (sizeof line_zero()) + (sizeof max_knot_ratio());
	return len;
}

//Dump Function
int MGTolerance::dump(MGOfstream& buf) const{ // 0:OK !0:NG -99: file not open

	int len = -1;
	try
	{
//		buf << m_mach_zero << m_wc_zero << m_wc_zero_sqr 
		buf << mach_zero() << wc_zero() << wc_zero_sqr() 
//			<< m_rc_zero << m_rc_zero_sqr << m_angle_zero
			<< rc_zero() << rc_zero_sqr() << angle_zero()
//			<< m_line_zero << m_max_knot_ratio;
			<< line_zero() << max_knot_ratio();
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGTolerance::restore(MGIfstream& buf){
	int len = -1;
	try
	{
		double m_mach_zero,m_wc_zero,m_wc_zero_sqr,
			m_rc_zero,m_rc_zero_sqr,m_angle_zero,
			m_line_zero,m_max_knot_ratio;
		buf >> m_mach_zero >> m_wc_zero >> m_wc_zero_sqr 
			>> m_rc_zero >> m_rc_zero_sqr >> m_angle_zero
			>> m_line_zero >> m_max_knot_ratio;
		len = dump_size();

		set_mach_zero(m_mach_zero);
		set_wc_zero(m_wc_zero);
		set_rc_zero(m_rc_zero);
		set_angle_zero(m_angle_zero);
		set_line_zero(m_line_zero);
		set_max_knot_ratio(m_max_knot_ratio);
	}
	catch (...)
	{
		throw;
	}
	return len;
}

// ***MGPosition_list
//Calculate dump size
size_t MGPosition_list::dump_size() const{
	const_iterator i;
	size_t len = sizeof entries();
	for (i = begin(); i != end(); i++)
	{
		len += (*i).dump_size();
	}
	return len;
}

//Dump Function
int MGPosition_list::dump(MGOfstream& buf) const{
	int len = -1;
	try
	{
		size_t len = entries();
		buf << len;
		const_iterator i;
		for(i = begin(); i != end(); i++)
		{
			(*i).dump(buf);
		}
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGPosition_list::restore(MGIfstream& buf){
	int len = -1;
	try{
		size_t len;
		buf >> len;
		for(size_t i = 0; i < len; i++){
			MGPosition Q;
			Q.restore(buf);
			push_back(Q);
		}
		len = dump_size();
	}
	catch (...){
		throw;
	}
	return len;
}

// ***MGOscuCircleData***
//Calculate dump size
size_t MGOscuCircleData::dump_size() const{
	size_t len = (sizeof m_index) + (sizeof m_radius);
	return len;
}

//Dump Function
int MGOscuCircleData::dump(MGOfstream& buf) const{
	int len = -1;
	try
	{
		buf << m_index << m_radius;
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

//Restore Function
int MGOscuCircleData::restore(MGIfstream& buf){
	int len = 0;
	try
	{
		buf >> m_index >> m_radius;
		len = dump_size();
	}
	catch (...)
	{
		throw;
	}
	return len;
}

// ***MGOscuCircle***
//Calculate dump size
size_t MGOscuCircle::dump_size() const{
	std::vector<MGOscuCircleData>::const_iterator i = m_circle.begin();
	size_t len = sizeof m_n;
	for(size_t iloop = 0; iloop<m_n; iloop++, i++){
		len += (*i).dump_size();
	}
	return len;
}

//Dump Function
int MGOscuCircle::dump(MGOfstream& buf) const{
	int len = -1;
	try{
		std::vector<MGOscuCircleData>::const_iterator i = m_circle.begin();
		buf<<m_n;
		for(size_t iloop = 0; iloop < m_n; iloop++, i++){
			(*i).dump(buf);
		}
		len = dump_size();
	}
	catch (...)	{
		throw;
	}
	return len;
}

//Restore Function
int MGOscuCircle::restore(MGIfstream& buf){
	int len = -1;
	try{
		buf >> m_n;
		for(size_t i = 0; i < m_n; i++){
			MGOscuCircleData Q;
			Q.restore(buf);
			m_circle.push_back(Q);
		}
		len = dump_size();
	}
	catch (...)	{
		throw;
	}
	return len;
}

// MGPPRep
size_t MGPPRep::dump_size() const{
	size_t len = (sizeof m_order) + (sizeof m_nbreak) + (sizeof m_sdim)
		+ m_break_point.dump_size();
	len+=(m_order*m_nbreak*m_sdim)*(sizeof m_coef[0]);
	return len;
}

int MGPPRep::dump(MGOfstream& buf) const{
	int len = -1;
	try{
		buf<<m_order<<m_nbreak<<m_sdim;
		m_break_point.dump(buf);
		size_t dim = m_order * m_nbreak * m_sdim;
		for(size_t i=0; i<dim; i++) buf<<m_coef[i];
		len = dump_size();
	}
	catch(...){
		throw;
	}
	return len;
}

int MGPPRep::restore(MGIfstream& buf){
	int len = -1;
	try{
		size_t order, nbreak, sdim;
		buf>>order>>nbreak>>sdim;
		resize(order, nbreak, sdim);
		m_break_point.restore(buf);
		size_t dim = m_order * m_nbreak * m_sdim;
		for(size_t i = 0; i < dim ; i++) buf>>m_coef[i];
		len = dump_size();
	}
	catch(...){
		throw;
	}
	return len;
}

// ***MGPosition***
//Calculate dump size
size_t MGPosition::dump_size() const{//Dumped Length;
	size_t len = 0;
	len = m_element.dump_size();
	return len;
}
int MGPosition::dump(MGOfstream& buf) const{
	int len = -1;
	try{
		m_element.dump(buf);
		len = dump_size();
	}
	catch(...){
		throw;
	}
	return len;
}

//Restore Function
int MGPosition::restore(MGIfstream& buf){  // 0:OK !0:NG -99:file not open
	int len = -1;
	try	{
		m_element.restore(buf);
		len = dump_size();
	}
	catch (...){
		throw;
	}
	return len;
}
