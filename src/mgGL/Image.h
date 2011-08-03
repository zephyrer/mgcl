/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGImage_HH_
#define _MGImage_HH_

#include "mg/MGCL.h"

//class MGOfstream;
//class MGIfstream;

/** @addtogroup DisplayHandling
 *  @{
 */

///Define MGImage Class.
///MGImage defines the attributes of Image.
class MGCLASS MGImage{

public:

///Union to treat one pixel as a unit and RGBA elements.
///One pixel consists of 4 bytes, which are (red, green , blue, alpha).
///If uint is used, the image's pixel is treated as a image,
///and if uchar[i] is used each element of the pixel can be handled.
union pixel{
	size_t uint;
	unsigned char uchar[4];
};

MGImage():m_width(0), m_height(0){;};

///Conversion constructor from Gdiplus::Bitmap.
MGImage(Gdiplus::Bitmap& bitmap);

///Extract a part of bitmap.
MGImage(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	///<left bottom address of bitmap.
	int width, int height
);

int width()const{return m_width;};
int height()const{return m_height;};
pixel* image(){return m_image.get();};
const pixel* image()const{return m_image.get();};
pixel& operator()(int i, int j);
const pixel& operator()(int i, int j)const;

///resize the image size to (width,height).
///Scaling the whole image to the size (width, height).
void resize(
	int width,
	int height
);

///resize the image size to (width,height) filling the color to
///the extra part for the size(width,height).
///resize_with_fill_color() does not perform scaling to the image.
///width and height can be less than the original length. In this case,
///image trimming will be done.
void resize_with_fill_color(
	int width,
	int height,
	const pixel& pdata
);

///Add border color(0,0,0,0)=(alfa=0) of pixel size2 for each perimeter.
void resize_and_add_zero_border(int nwidth2, int nheight2);

private:
	int m_width, m_height;		///<number of pixels about the width and the height.
	std::auto_ptr<pixel> m_image;///<array of pixels of size m_width*m_height;

///Extract a part of bitmap into this, from(x,y) to (x+width, y+height).
void extract(
	Gdiplus::Bitmap& bitmap,
	int x,  int y,	///<left bottom address of bitmap.
	int width, int height
);

};

///Compute 2's power of width and height.
void MGImageCompute_2spower(
	int width, int height,
	int& width2, int& height2///<The smallest 2's power of width and height will be output.
);

/** @} */ // end of DisplayHandling group
#endif
