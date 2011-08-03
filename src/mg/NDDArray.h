/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGNDDArray_HH_
#define _MGNDDArray_HH_
/** @addtogroup BASE
 *  @{
 */
#include "mg/MGCL.h"

// MGNDDArray.h
//

// Forward Declaration
class MGBPointSeq;
class MGKnotVector;
class MGKnot;
class MGIfstream;
class MGOfstream;

///Defines non-decreasing double data array.
///Used for data point abscissa, or knot vector, etc.
///MGNDDArray has size and effective data length.
class MGCLASS MGNDDArray{

public:

///Friend Function
MGDECL friend MGNDDArray operator* (double scale, const MGNDDArray& nd);

///String stream Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGNDDArray&);

//////////// Constructor ////////////

///void constructor.
MGNDDArray():m_capacity(0),m_length(0),m_current(0),m_element(0){;}; 

///Constructor MGNDDArray of size n and lenght=n.
explicit MGNDDArray(
	size_t n,		///<size of this.
	const double* data=0///<data array of length n if data!=NULL.
);

/// Construct MGNDDArray so that initial data is init and
/// is incremented by increment.
MGNDDArray(size_t n, double init, double increment=1.0);

/// From Data Point ordinate, obtain data point seq abscissa.
explicit MGNDDArray(const MGBPointSeq&);

/// From Data Point ordinate with End condition, obtain data point
/// seq abscissa.
MGNDDArray(MGENDCOND begin, MGENDCOND end, const MGBPointSeq&);

/// From knot vector with End conditions, obtain data point
/// seq abscissa. This data point can be input of MGSBRep constructor
///of the following form:
///MGSBRep::MGSBRep(
///MGSBRepEndC& endc,///end condition
///const MGNDDArray& utaui,	///Data point of u-direction
///const MGNDDArray& vtaui,	///Data point of v-direction
///const MGSPointSeq& value,	///Data point ordinate
///const MGKnotVector& tu,	///knot vector of u-direction, of order 4
///const MGKnotVector& tv,	///knot vector of v-direction, of order 4
///int &error)				///Error flag.
///That is, generate utaui or vtaui from tu or tv each, taking endc 
///into account.
///The order is assumed to be 4.
MGNDDArray(MGENDCOND begin, MGENDCOND end, const MGKnotVector& t);

/// From data point, obtain data point of updated number(nnew).
///If original data point has multiplicities and nnew>=length(),
///original data point parameters and the multiplicities are preserved.
///Same as change_number().
MGNDDArray(const MGNDDArray&, size_t nnew);

/// From data point, obtain data point of updated value range.
///Update so that (*this)(0)=ts, and (*this)(lenght()-1)=te.
///Must be ts<te.
MGNDDArray(const MGNDDArray&, double ts, double te);

/// Construct by extracting sub interval of array2.
MGNDDArray(
	size_t start_id,			///<Start id of array2(from 0).
	size_t num,					///<new array length.
	const MGNDDArray& array2	///<Original NDDArray.
);

///Construct by mixing two arrays.
///Mixing is so done as that no too close points are included.
///Although data point multiplicities of array1 are preserved,
///multiplicities of array2 are not.
///DATA POINT MULTIPLICITY IS ALLOWED ONLY IN array1.
MGNDDArray(
	size_t id1,			///<Start id of array1(from 0).
	size_t num1,		///<new array length to use of array1.
	const MGNDDArray& array1,///<Original NDDArray1.
	size_t id2,			///<Start id of array2(from 0).
	size_t num2,		///<new array length to use of array2..
	const MGNDDArray& array2///<Original NDDArray2.
);

///Copy constructor.
MGNDDArray(const MGNDDArray& nd);

//////////// Destructor ////////////
virtual ~MGNDDArray(){if(m_element) delete[] m_element;};

//////////// Operator Overload ////////////
	
///Access to i-th element.
double operator[] (size_t i) const{return m_element[i];};
double operator() (size_t i) const{return m_element[i];};

///Access to i-th element.
double& operator[] (size_t i){return m_element[i];};
double& operator() (size_t i){return m_element[i];};

///Assignment
MGNDDArray& operator= (const MGNDDArray& vec2);

///Addition and subtraction of real number.
///All of the elements will be added or subtracted.
MGNDDArray operator+(double) const;
MGNDDArray& operator+=(double);
MGNDDArray operator-(double) const;
MGNDDArray& operator-=(double);

/// 単項マイナス。
///Unary minus. Reverse the ordering of elements by changing all of the
///signs.
MGNDDArray operator-() const;

///Scaling
MGNDDArray operator*(double scale) const;
MGNDDArray& operator*= (double scale);

///Copmarison operator.
bool operator== (const MGNDDArray& t2) const;
bool operator!= (const MGNDDArray& t2) const{return !operator==(t2);}

//////////// Member Function ////////////

///Add data of multiplicity 1 into data points. mult_max is the maximum
///multiplicity allowed for NDDArray.
///Return value is number of data actually added.
int add_data(double value, size_t mult_max=1);

///Add data with multiplicity into data points. mult_max is the maximum
///multiplicity allowed.	Return value is number of data actually added.
int add_data(const MGKnot& knot, size_t mult_max);

///Copy data points from array by removing the mltiple knot.
void copy_removing_multi(
	size_t start_id,			///<Start id of array2(from 0).
	size_t num,					///<new array length from start_id
	const MGNDDArray& array2	///<Original NDDArray.
);

///Return a pointer to raw data of MGNDDArray.
const double* data(size_t i=0) const{ return &m_element[i]; };

///Return a pointer to raw data of MGNDDArray.
double* data(size_t i=0) { return &m_element[i]; };

///Delete one data at index.
///Return value is new total number of data in the array, generally
///is original_length()-1.
int del_data(size_t index);

///Test if this is null.
bool is_null()const{return m_length==0;};

/// Return the length of MGNDDArray.
size_t length() const{ return m_length; };

///Set the length of effective data.
void set_length(size_t length);

///Set this as a null NDDArray.
void set_null();

///Finds index where tau is located in MGNDDArray as an index of knot:
/// tau < (*this)(0) : index=-1.
/// (*this)(0) <= tau< (*this)(n-1) : 0<= index <n-1, such that
///                         (*this)(index) <= tau < (*this)(index+1)
/// (*this)(n-1) <= tau : index=n-1.
///Here n=lenght().
virtual int locate(double tau) const;

///Locate where data of multiplicity of multi is after start.
///index is the starting point index of this found first after start.
///index>=start if index>=0.
///Function's return value locate_multi is actual multiplicity at the
///index, i.e. locate_multi>=multi.
///If position of the multiplicity is not found to the end,
///index=(lenght()-1) (index of the last element) and locate_multi=0
///will be returned.
///multi must be >=1.
size_t locate_multi(size_t start, size_t multi, size_t& index) const;

///Update so that (*this)(0)=ts, and (*this)(lenght()-1)=te.
///Must be ts<te.
virtual void change_range(double ts, double te);

///Update array length.
///Updated array is so generated that the original proportions of
///neighbors hold as much as possible. 
///If original data point has multiplicities and nnew>=length(),
///original data point parameters and the multiplicities are preserved.
MGNDDArray& change_number(size_t nnew);

///Reference to i-th element.
double ref(size_t i) const{ return m_element[i];};

///Remove too near data points. Removal will be done with the ordinates.
void remove_too_near(
	MGBPointSeq& ordinates,	///<ordinate,
							///<ordinates.length() must be equal to this->length().
	bool allow_multi=false,	///<indicates if multiple data point is allowed or not,
							///<when allow_multi=false, multiple data points will be removed,
							///<when allow_multi=true, will not be removed.
	double ratio=6.	///<maximum ratio allowed for neighboring span,
			///<let ti=(*this)[i], then
			///<if (t(i+1)-ti)/(ti-t(i-1))>ratio or (t(i+1)-ti)/(ti-t(i-1))<1/ratio,
			///<a data point will be removed(along with the ordinates).
);

///Change the size. 
///start is to indicate from which location of new area to start storing.
///Although size can be less than original length, some of end data will be
///lost in this case. When 'start'>0, first 'start' data will be garbage.
void reshape(size_t size, size_t start=0);

///Resize the array. Result will contain garbages.
///length() and size() will have nsize.
void resize(size_t nsize);

/// Return the size of MGNDDArray.
size_t capacity() const { return m_capacity; };

///Obtain data point from KnotVector t and replace own with it.
MGNDDArray& update_from_knot(const MGKnotVector& t);

///Dump Functions.
///Calculate dump size
size_t dump_size() const;

///Dump Function
int dump(MGOfstream& ) const;

///Restore Function
int restore(MGIfstream& );
	
//////////// Member data /////////
protected:
	mutable int m_current;///<current interval of the array is held.

private:
	size_t m_capacity;	///<Size of m_element.
	size_t m_length;	///<Length of effective data in m_element.
	double* m_element;

};

/** @} */ // end of BASE group
#endif
