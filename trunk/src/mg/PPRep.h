/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGPPRep_HH_
#define _MGPPRep_HH_
/** @addtogroup GEORelated
 *  @{
 */

#include "mg/NDDArray.h"
#include "mg/Vector.h"

// MGPPRep.h
//

// Forward Declaration
class MGLBRep;
class MGIfstream;
class MGOfstream;

///Defines PP-Represetation.
///The indices are:
///MGPPRep(i,j,k) where i is order, j is break point id, and k is space dimensio id.
///When tau(j)<=t<tau(j+1)(tau(j) is break_point(j)), using coef(i,j,k), the PP-representation
///f(t) is: f(t)=sum of i(coef(i,j,k)/factorial(i)*(power(i) of (t-tau(j)))) for i=0,...,order-1.
///In other words, the coefficinet of (t-tau(j))**i is coef(i,j,k)/factorial(i) for k=0,...,sdim()-1,
///and i=0, ... ,order-1.
///Here, ** i indicates power i.
class MGCLASS MGPPRep {
 
public:

///String stream output Function
MGDECL friend std::ostream& operator<< (std::ostream&, const MGPPRep& );

//////////// Constructor ////////////

///Default constructor.
MGPPRep();

///Constructor of dummy PP-Rep of specified size.
MGPPRep(unsigned order, size_t nbreak, size_t sdim);

///Constructor of dummy PP-Rep(no data except data points) of specified size.
MGPPRep(unsigned order, size_t sdim, const MGNDDArray& tau);

///Constructor to convert from Line B-Representation.
MGPPRep(const MGLBRep& lbrep);

///Constructor to change order of original PP-Representation.
///New order may be greater or less than the original one. However,
///if new one is less than the original, PP-Rep constructed may not
///be able to hold the same shape.
MGPPRep(unsigned order, const MGPPRep& pprep);

///Copy constructor.
MGPPRep(const MGPPRep& rhs);

////////////Destructor/////////

~MGPPRep(){if(m_coef) delete[] m_coef;};
									
//////////// Operator Overload ////////////

MGPPRep& operator=(const MGPPRep&); ///Assignment

///Return i-th element of the position.
double operator()(size_t i,size_t j, size_t k) const{return coef(i,j,k);}
double& operator()(size_t i,size_t j, size_t k){return coef(i,j,k);}

///Extract (i,j,k)elements for 0<=k<sdim() as a vector.
MGVector operator()(size_t i, size_t j)const{return coef(i,j);};

//////////// Member Function ////////////

///Returns break point sequence.
const MGNDDArray& break_point()const{return m_break_point;}

///Returns i-th break point value.
double break_point(size_t i)const{return m_break_point(i);}

///Returns a pointer to i-th break_point to access.
double& break_point(size_t i){return m_break_point(i);}

///Returns a pointer to the break_point data.
const double* break_point_data()const{return m_break_point.data();}		

///Returns (i,j,k)-th coef value.
double coef(size_t i, size_t j, size_t k)const;

///Returns a pointer to coef(i,j,k) to access.
double& coef(size_t i, size_t j, size_t k);

///Extract (i,j,k)elements for 0<=k<sdim() as a vector.
MGVector coef(size_t i, size_t j)const;

///Returns a pointer to the PPCoef data.
const double* coef_data(size_t i=0, size_t j=0, size_t k=0)const;

///Evaluate right continuous n'th derivative(BPVAL)
MGVector eval(		
	double t,	///<Parameter value to evaluate.
	size_t n=0	///<Degree of derivative. When n=0, compute positional data.
)const;

///Evaluate i-th span's n'th derivative(BPVAL)
MGVector eval_i(	
	size_t i,		///<span number(from 0) of the PP-Rep.
	double t,		///<Parameter value to evaluate.
	size_t n		///<Dgree of derivative. When n=0, compute positional data.
)const;		

///Returns the number of Break points including the end point.
size_t nbreak()const{return m_nbreak;}

///"normalize" normalizes the PP-Rep, i.e. changes break point data
///and pp-coefficients so as that length of first derivatives
/// from left and right at each break point are the same.
MGPPRep& normalize();

///Returns the order of PPRep.
unsigned order()const{return m_order;}

///Return the data of (i,j,k).
double ref(size_t i, size_t j, size_t k)const;	

///Change size. Change of sdim not allowed.
///Stored data so far will be guarateed to hold in the same id of coef(i,j,k).
void reshape(size_t nbreak); 

///Resize to (order, nbreak, dim).
///Resutl will contain only garbages.
void resize(
	size_t order,	///<Order number.
	size_t nbreak,	///<number of break points, includes last points.
	size_t dim	///<Space dimension.
);

///Returns the space dimension.
size_t sdim()const{return m_sdim;}

///Store the vector v at coef(i,j).
void store_at(size_t i, size_t j, const MGVector& v);

///Dump Functions
size_t dump_size()const;

///Dump Function
int dump(MGOfstream& )const;

///Restore Function
int restore(MGIfstream& );

//////////// Member Data ////////////

private:

	size_t m_order;		///< Order of PP-Rep.
	size_t m_nbreak;	///< number of Break points(Interval number=m_nbreak-1).
						///< m_nbreak==m_break_point.length()
	size_t m_sdim;		///< Space Dimension.
	MGNDDArray m_break_point;///<Break point sequence.
	double* m_coef;		///<PP coef.

};

/** @} */ // end of GEORelated group
#endif
