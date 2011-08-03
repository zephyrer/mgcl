/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTextureImage_HH_
#define _MGTextureImage_HH_

#include <iosfwd>
#include <memory>
#include "mg/MGCL.h"
#include "mg/Position.h"
#include "mg/Box.h"
#include "mgGL/image.h"

//Define MGTextureImage Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGTextureImage defines the image data for OpenGL's glTexImage2D.
///That is, the image size will be rounded to 2's power.
class MGCLASS MGTextureImage{

public:

MGTextureImage(){;};

MGTextureImage(
	Gdiplus::Bitmap& bitmap,	///<texture image data.
	double left, double bottom,	///<texture coordinates of left and bottom point of the image.
	double right, double top,	///<texture coordinates of right and top point of the image.
	bool scalable=true			///<indicates if the image can be scalabe or not.
			///<if scalable=false, (right, top) will be changed according to the expanded image area.
);

MGTextureImage(
	MGImage& bitmap,		///<image data.
	double left, double bottom,///<texture coordinates of left and bottom point of the image.
	double right, double top,///<texture coordinates of right and top point of the image.
	bool scalable=true		///<indicates if the image can be scalabe or not.
			///<if scalable=false, (right, top) will be changed according to the expanded image area.
);

///copy constructor.
MGTextureImage(const MGTextureImage& txtr);

///Destructor.
~MGTextureImage(){;};

///Assignment.
MGTextureImage& operator=(const MGTextureImage& txtr);

///render process.
void bind_texture_object()const;
void apply_texture_object()const;

///Return the box.
const MGBox& box()const{return m_box;};
MGBox& box(){return m_box;};

int nwidth()const{return m_image.width();};
int nheight()const{return m_image.height();};

double left()const{return m_box[0][0].value();};
double bottom()const{return m_box[1][0].value();};
double right()const{return m_box[0][1].value();};
double top()const{return m_box[1][1].value();};
double width()const{return m_box[0].length().value();};
double height()const{return m_box[1].length().value();};
MGPosition center()const{return m_box.mid();};
size_t dlist_name()const{return size_t(this);};
const MGImage::pixel* image()const{return m_image.image();};

///Set texture coordinates of this.
void set_address(
	double left, double bottom,
	double right, double top
);

protected:
	MGImage m_image;

private:
	MGBox m_box;
			///<Texture coordinates from (m_box[0], m_box[1]) to (m_box[2], m_box[3]) in the world coordinates.
			///<(m_box[2]-m_box[0])/m_width is one pixel width.

void round(bool scalable);

};

/** @} */ // end of GLAttrib group
#endif
