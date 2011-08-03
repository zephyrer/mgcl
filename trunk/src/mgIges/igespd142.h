/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD142_H__)
#define __MGIGESPD142_H__

#include "mgiges/IgesIfstream.h"
#include "mgiges/IgesOfstream.h"
#include "mgiges/IgesPD.h"
class MGLoop;

///MGIgesPD142 is the class for Iges parameter data type 142(Curve on parameteric space).
class MGIgesPD142: public MGIgesPD{
public:
	/// Constructors.

	/// Constructs an object of class MGIgesPD142.
	MGIgesPD142(MGIgesDirectoryEntry* DEpointer=0);

	/// Constructs an object of class MGIgesPD142.
	MGIgesPD142(
		const MGLoop& loop,	///loop to make PD142. This is a loop of the face.
		int surface_DE,	///the base surface. The surface must be output to IGES file first.
		MGIgesOfstream& igesfile///Iges file to output.
	);

	///Destructor;
	~MGIgesPD142(){;};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);
	
	///Obtain both the parametric space curve of the surface and the model space curve.
	void trim_face(
		const MGIgesIfstream& igesifstrm,
		std::auto_ptr<MGFace>& face,///<Face to be trimmed by this boundary MGIgesPD142.
		bool outer=true	///<True if this be the outer boundary.
	)const;

	///Write out this PD as MGIgesParamLine's(into plines).
	///Except for string data, one integer or double data is output
	///into one MGIgesParamLine, not striding over more than one line.
	///Only when string data is output(to Holleris string), the data
	///may stride over more than one lines.
	///plines[i] for 0<=i<plines.size() are valid.
	void write_out_into_string(
		const MGIgesGSec& gsec,	///<Input gsec to input delimeter_param and delimeter_record;
		MGPvector<std::string>& plines ///<output plines.
	)const;

public:
//Member data.

	short m_created_way;///<Indicates the way the curve on the surface has been created:
			///< =0: unspecified,
			///< =1: projection of a given curve on the surface,
			///< =2: intersection of two surfaces.
			///< =3: isoparametric curve, either a u or v-parameter curve.
	short m_prefered;///<indicates prefered representation:
			///< =0: unspecified,
			///< =1:S(m_param_curve_DE(t)) is prefered,
			///< =2: m_modelcurve_DE is prefered.
			///< =3: m_param_curve_DE and m_model_curve_DE are equally prefered.
	int m_surface_DE;///<Directory entry of the surface on which the curve lies.
	int m_param_curve_DE;///<Directory entry of the parametric space curve of the surface.
	int m_model_curve_DE;///<Directory entry of the curve(in the model space).
};

#endif // __MGIGESPD142_H__
