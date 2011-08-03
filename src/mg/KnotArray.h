/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGKnotArray_HH_
#define _MGKnotArray_HH_
/** @addtogroup BASE
 *  @{
 */

#include <vector>
#include "mg/Knot.h"

#if defined(MGCL_DLL)
#pragma warning( push )
#pragma warning( disable : 4231 )
MGTEMPLATE template class MGCLASS std::vector<MGKnot>;
#pragma warning( pop )
#endif

// MGKnotArray.h
//

//Forward Declaration
class MGIfstream;
class MGOfstream;

/// Defines Array of Knots.
///This is not knot vector for B-Rep.
///See MGKnotVector.
class MGCLASS MGKnotArray{

///Member Data

public:

std::vector<MGKnot> m_ktArray;

///String stream Function.
MGDECL friend std::ostream& operator<< (std::ostream&, const MGKnotArray&);

////////////Constructor////////////

/// Dummy constructor.
MGKnotArray(){;};

///From Knot.
explicit MGKnotArray(const MGKnot&);

///From knot and the multiplicity.
MGKnotArray(double knot, int mult);

//	MGKnotArray(const MGKnotArray&);	//Copy Constructor.

///Destructor
	~MGKnotArray(){;};

////////////Operator overload.////////////

///[ ] reference.
const MGKnot& operator[](size_t pos) const
{return m_ktArray.operator[](pos);}
MGKnot& operator[](size_t pos)
{return m_ktArray.operator[](pos);}

////////////Member Function////////////

///Add to the end of list.
MGKnotArray& add(const MGKnot&);
MGKnotArray& add(double knot, int mult);

///Get the last element in the vector.
const MGKnot& back() const{return m_ktArray.back();};
MGKnot& back(){return m_ktArray.back();};

///Obtain the iterator of the first element.
std::vector<MGKnot>::iterator begin(){return m_ktArray.begin();};
std::vector<MGKnot>::const_iterator begin() const{return m_ktArray.begin();};

///Clear all the elements in m_CCilist.
void clear(){m_ktArray.clear();};

///Obtain the iterator of the first element.
std::vector<MGKnot>::iterator end(){return m_ktArray.end();};
std::vector<MGKnot>::const_iterator end() const{return m_ktArray.end();};

///Get the first element in the vector.
const MGKnot& front() const{return m_ktArray.front();};
MGKnot& front(){return m_ktArray.front();};

void push_back(const MGKnot& knot){m_ktArray.push_back(knot);};

/// Return the number of items that are in the list.
size_t length() const {return m_ktArray.size();};
size_t size() const{return m_ktArray.size();};

///Operator overload.
///	MGKnotArray& operator =(MGKnotArray&);///Assignment operator overload.

///Dump Functions.
	///Calculate dump size
	size_t dump_size() const;

	///Dump Function
	int dump(MGOfstream& ) const;

	///Restore Function
	int restore(MGIfstream& );

};

/** @} */ // end of BASE group
#endif
