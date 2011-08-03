/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __MGIGESPD141_H__)
#define __MGIGESPD141_H__

#include <vector>
#include "mgiges/IgesPD.h"

///@cond
//A Private class for MGIges141Edge, express one IGES edge for IGES boundary.
class MGIges141Edge{
public:
	MGIges141Edge(
		int curve_DE=0, 
		int sense=1
	):m_curve_DE(curve_DE), m_sense(sense){;};

	void push_back_pcurve(int pcurve_DE){m_pcurves.push_back(pcurve_DE);};

public:
	int m_curve_DE;	///<DE pointer of the model space curve of MGIgesPD141.
	int m_sense;	///<ordering flag.
		///<=1: the direction of the model space curve does not require reversal.
		///<=2: the direction of the model space curve need to be reversed.
		///<	  the parameter space curve and model space curve orientations disagree.
	std::vector<int> m_pcurves;///<vector of parameter space curve directory entries of
		///<the model space curve m_curve_DE.
};
///@endcond

///MGIgesPD141 is the class for Iges parameter data type 141(BOUNDARY entity).
class MGIgesPD141: public MGIgesPD{
public:
	// Constructors.

	/// Constructs an object of class MGIgesPD141.
	MGIgesPD141(MGIgesDirectoryEntry* DEpointer=0);

	///Destructor;
	~MGIgesPD141(){;};

	///Read in parameter data from string stream data.
	void read_in(
		char pDelimeter,
		std::istringstream& pdstream
	);

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
///Member data.

	short m_type;///<Type of bounded surface representation,
			///<=0: the boundary shall reference only model space trimming curves.
			///<=1: the boundary shall reference model space and associated parameter space
			///<curve collections.
	short m_prefered;///<Prefered representation of the trimming curves.
			///<=0:Unspecified, =1:Model space, =2: Parameter space, =3:of equal preference.
	int m_surface_DE;///<Directory entry of the untrimmed surface(base surface).
	std::vector<MGIges141Edge> m_edges;///<vector of MGIges141Edge.
};

#endif // __MGIGESPD141_H__