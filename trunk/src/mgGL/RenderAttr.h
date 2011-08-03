/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGRenderAttr_HH_
#define _MGRenderAttr_HH_

#include "mgGL/ShadeModel.h"
#include "mgGL/TranspMode.h"

class MGMaterial;
class MGTexture;
class MGOfstream;
class MGIfstream;

//
//Define MGRenderAttr Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGRenderAttr defines the attributes of rendering attributes.
///These attrubutes are not used for drawing(line drawing mode).
///Regarding to m_material, m_back_material, and m_texture, MGRenderAttr behaves 
///just like auto_ptr. That is, these three pointers are newed object pointers,
///and when used in copy constructor or assignment, the ownerships are transfered to
///the new MGRenderAttr object.
class MGCLASS MGRenderAttr:public MGGLAttrib{

public:

enum RENDERSIDE{
	UNDEFINED=MGGLAttrib::UNDEFINED,
	DISABLED=MGGLAttrib::DISABLED,
	FRONT=GL_FRONT,		///render front side
	BACK=GL_BACK,		///render back side
	FRONT_AND_BACK=GL_FRONT_AND_BACK	///render both sides.
};

friend GLenum GLrender_side(RENDERSIDE rs){return static_cast<GLenum>(rs);};

MGRenderAttr(
	RENDERSIDE rs=FRONT,
	MGShadeModel sm=MGShadeModel::SMOOTH
):MGGLAttrib(static_cast<int>(rs)),m_shade_model(sm),m_material(0),
m_back_material(0),m_texture(0){;};

///copy constructor.
MGRenderAttr(const MGRenderAttr& attr);

///Destructor.
~MGRenderAttr();

///Assignment
MGRenderAttr& operator=(const MGGel& gel2);
MGRenderAttr& operator=(const MGRenderAttr& gel2);

///comparison
bool operator<(const MGRenderAttr& gel2)const;
bool operator<(const MGGel& gel2)const;

////////////Member Function////////////

///Generate a newed clone object.
MGRenderAttr* clone()const;

///Invoke appropriate OpenGL fucntion to this attribute.
void exec()const;

///ShadeModel
MGShadeModel::MODEL get_shade_model()const{return m_shade_model.get_model();};
void set_shade_model(MGShadeModel::MODEL m=MGShadeModel::SMOOTH){
	m_shade_model.set_model(m);
}

///Material
bool material_defined()const;
const MGMaterial* material()const{return m_material;};

///Set the material. When rs=FRONT_AND_BACK and different material for the back side
///is used, set_back_material must be invoked after invoking set_material.
///Else the same material will be appllied for the both sides.
void set_material(
	RENDERSIDE rs,
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess=0.,
	float transparency=0.
);

///Set the back side material. Invoking set_back_material means two sided material
///and setting different material to the back side.
///Before use of set_back_material, set_material must be invoked first.
///set_back_material will set two sided material.
void set_back_material(
	const float ambient[3],
	const float diffuse[3],
	const float specular[3],
	const float emission[3],
	float shininess=0.,
	float transparency=0.
);

///Set the render side.
void set_render_side(RENDERSIDE side=FRONT){m_flag=static_cast<int>(side);};

///Return the material pointer of the back side when this has the two sided rendering
///(FRONT_AND_BACK) and the different materials are applied to each side.
///NULL will be returned when render_side()!=FRONT_AND_BACK. 
///Even when render_side()==FRONT_AND_BACK, if m_back_material==null(the same material
///is applied to both side), null will be returned.
const MGMaterial* back_material()const;

///Get the render side.
///FRONT:render front side, BACK:render back side, FRONT_AND_BACK: render both sides.
RENDERSIDE render_side()const{return static_cast<RENDERSIDE>(m_flag);};
GLenum GLrender_side()const{return static_cast<GLenum>(render_side());};

///Texture
const MGTexture* texture()const{return m_texture;};
void set_texture(MGTexture* texture){m_texture=texture;};

///draw GLAttribute process.
void drawAttrib(
	bool no_color=false	///<if true, color attribute will be neglected.
)const{;};

///render GLAttribute process.
void render()const{exec();};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_draw_attrib_mask(unsigned int& mask)const{;};

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_draw_attrib_mask(unsigned int& mask)const{;};

///Turn on the appropriate mask bit for this attribute. See glPushAttrib().
void set_render_attrib_mask(unsigned int& mask)const;

///Turn off the appropriate mask bit for this attribute. See glPushAttrib().
void reset_render_attrib_mask(unsigned int& mask)const;

/// Return This object's typeID
long identify_type() const{return MGRENDER_ATTR_TID;};

std::string whoami()const{return "RenderAttr";};

///Read all member data.
void ReadMembers(MGIfstream& buf);
///Write all member data
void WriteMembers(MGOfstream& buf)const;

/// Output function.
std::ostream& out(std::ostream&) const;

private:

	MGShadeModel m_shade_model;///< <<lighting>>
	MGMaterial* m_material;		///<redering side and material definition.<<lighting>>
		///<m_material->render_side() defines the redering side.
		///<When m_material->render_side()==FRONT_AND_BACK and m_back_material is not null,
		///<different materials are applied to each side(m_material for front and m_back_material
		///<for back).
		///<When m_material->render_side()==FRONT_AND_BACK and m_back_material is null,
		///<the same material(m_material) is applied to each side.
		///<When m_material->render_side()==BACK, m_material is used and m_back_material is NOT
		///<used for the back side rendering.
	MGMaterial* m_back_material;///<back material.<<lighting>>
	MGTexture* m_texture;		///<texture<<texture>>

	friend class MGAppearance;

};

/** @} */ // end of GLAttrib group
#endif
