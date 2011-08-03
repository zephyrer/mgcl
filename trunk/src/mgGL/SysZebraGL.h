/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
/// SysZebraGL.h: mgSysZebraGL クラスのインターフェイス
#if !defined(AFX_SysZebraGL_H__1EC223A1_C0D3_4AA7_80BD_89392C970098__INCLUDED_)
#define AFX_SysZebraGL_H__1EC223A1_C0D3_4AA7_80BD_89392C970098__INCLUDED_

#include "mgGL/SysGL.h"
#include "Tl/TLInputParam.h"

/** @addtogroup DisplayHandling
 *  @{
 */

///Defines a tool to display zebra map on a surface.
class mgSysZebraGL: public mgSysGL{

public:
	
enum{TEXGEN_WIDTH = 16};
union zebraData{
	char BYTE[4];
	int INTGER;
};

mgSysZebraGL(
	size_t function_code,	///<Function code.
	const MGGel* gel,	///<This must be MGSurface, MGFace, or MGShell(mgAll_2Manifold).
	const mgTLInputParam& tlparam,	///<Tessellation parameter.
	unsigned char color[4],	///<zebra color
	float thickness,	///<zebra thickness.
	bool is_vertical	///<direction of the zera stripes.
);

///Construct new object by copying to newed area.
///User must delete this copied object by "delete".
mgSysZebraGL* clone()const;

///Draw this Sysgl.
///This draw is used to draw the pictures for Undo(, Redo) operations.
///When this draw is invoked, all of the functions of mgGDL can be
///used in that context.
void draw(MGOpenGLView& glv)const;

///Process necessary before transformation in MGOpenGLView::DrawScene.
virtual void pre_transform_process()const;

/// Output virtual function.
///Output to stream file:メンバデータを標準出力に出力する。
virtual std::ostream& out(std::ostream& ostrm) const;

private:
	mgTLInputParam m_tlparam;
	float m_thickness;
	bool m_vertical;///<vertical or horizontal.
	unsigned char m_color[4];///<zebra color
};

/** @} */ // end of DisplayHandling group
#endif // !defined(AFX_SysZebraGL_H__1EC223A1_C0D3_4AA7_80BD_89392C970098__INCLUDED_)
