/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTrimLoop_HH_
#define _MGTrimLoop_HH_

#include <iosfwd>
#include "topo/LEPoint.h"

class MGLoop;

/** @addtogroup TOPORelated
 *  @{
 */

///MGTrimLoop is a private and utility class to implement trimming of MGFSurface.
///MGTrimLoop expresses a loop and the points of the start and end that are connected to
///a legacy(loops before trimming is done) boundary.
class MGTrimLoop{

public:

	///Stream output of the content.
	friend std::ostream& operator<< (std::ostream&, const MGTrimLoop&);

	MGTrimLoop():m_loop(0){;}
	MGTrimLoop(const MGTrimLoop& linf2);
	MGTrimLoop(
		MGLoop* loop,int star_loop_id,MGLEPoint& start_lep, int end_loop_id,MGLEPoint& end_lep
	);

	~MGTrimLoop();
	MGTrimLoop& operator=(const MGTrimLoop& loop2);
	bool is_null()const{return m_loop==0;};
	MGLoop* loop(){return m_loop;};
	MGLoop* release_loop();

	bool start_is_on_boundary(){return m_start_loopid==4 || m_start_loopid<0;};
	MGLEPoint start_lep()const{return m_start;};

	///Valid only when start_is_on_boundary();
	size_t start_loopid()const;

	bool end_is_on_boundary(){return m_end_loopid==4 || m_end_loopid<0;};
	MGLEPoint end_lep()const{return m_end;};

	///Valid only when start_is_on_boundary();
	size_t end_loopid()const;

	void set_null();

	///Set used flag of used_loops@as true@for the both end loops id.
	void set_used_loop_flag(std::vector<bool>& used_loops)const;

private:
	MGLoop* m_loop;///<Newed object of MGLoop
	int m_start_loopid;
		///<loop id(output of MGFace::in_range_with_on()) which start point of m_loop
		///<is connected to.
	MGLEPoint m_start;///<Only when m_start_loopid<0(on an inner loop), or 
		///<m_start_loopid=4(on the outer boundary loop), m_start is valid.

	int m_end_loopid;
		///<loop id(output of MGFace::in_range_with_on()) which end pointof m_loop
		///<is connectetd to.
	MGLEPoint m_end;///<Only when m_end_loopid<0(on an inner loop), or 
		///<m_end_loopid=4(on the outer boundary loop), m_end is valid.

};

/** @} */ // end of TOPORelated group
#endif
