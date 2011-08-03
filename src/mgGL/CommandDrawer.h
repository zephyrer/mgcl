/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#if !defined( __COMMANDDRAWER_H__)
#define __COMMANDDRAWER_H__
#include "mg/MGCL.h"
#include "mgGL/Color.h"
#include "mg/Position.h"
#include <vector>

/** @addtogroup ALGORITHM
 *  @{
 */

///MGCommandDrawer is a utility class fo MGOpenGLView class.
///Using MGOpenGLView's attach_drawer(), users can draw some pictures in
///the MGOpenGLView as highest priority ones. That is the pictures are drawn
///last after hilighted pictures. MGCommandDrawer is designed for GL6 command
///as a tool to draw some pictures of the commands.
class MGCLASS MGCommandDrawer{
protected:
	const MGPosition* m_cursor;///Mouse cursor position of the interaction.
	const std::vector<MGPosition>* m_points;///Input positions so far for the command.

public:
	MGCommandDrawer(
		const MGPosition* cursor=0,
		const std::vector<MGPosition>* points=0
		):m_cursor(cursor),m_points(points){;};
	virtual ~MGCommandDrawer(){;};
	virtual void draw(){;};///This function must be overided.

	const MGPosition& cursor()const{return *m_cursor;};
	const std::vector<MGPosition>& points()const{return *m_points;};

	///Points drawer service program.
	virtual void draw_points(
		const MGColor& color,
		const std::vector<MGPosition>& ipos
	)const;

	///Polyline drawer service program.
	virtual void draw_lines(
		const MGColor& color,
		const std::vector<MGPosition>& ipos
	)const;

	///Polyline and points  drawer service program.
	virtual void draw_points_lines(
		const MGColor& pcolor,	///<color of points(boundary)
		const MGColor& lcolor,	///<color of line.
		const std::vector<MGPosition>& ipos
	)const;

	void set_points(const std::vector<MGPosition>* points){m_points=points;};
	void set_cursor(const MGPosition* cursor){m_cursor=cursor;};
};

/** @} */ // end of ALGORITHM group

#endif //__COMMANDDRAWER_H__
