/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGTextureImages_HH_
#define _MGTextureImages_HH_

#include <vector>
#include "mg/Box.h"
#include "mgGL/TextureImage.h"

//Define MGTextureImages Class.

/** @addtogroup GLAttrib
 *  @{
 */

///MGTextureImages defines the attributes of TextureImages.
class MGCLASS MGTextureImages{

public:

MGTextureImages()
:m_total_pixel_w(0), m_total_pixel_h(0),
	m_texture_size(0), m_remainder_w(0), m_remainder_h(0),
	m_repeatS(true), m_repeatT(true){;};

MGTextureImages(
	Gdiplus::Bitmap& bitmap,
	int texture_size,	///<texture size to use.
	double pixel_width, double pixel_height,///<One pixel size, width and height.
	bool repeatS,bool repeatT///<repeat to s or t coordinate.
);

MGTextureImages(
	Gdiplus::Bitmap& bitmap,
	int texture_size,	///<texture size to use.
	double pixel_length,///<One pixel size in world coordinate for width and height,
						///<which are the same.
	bool repeat			///<repeat to s or t coordinate, which are the same.
);

/*
///copy constructor.
MGTextureImages(const MGTextureImages& txtr);

///Destructor.
~MGTextureImages();

///Assignment.
MGTextureImages& operator=(const MGTextureImages& attr);
*/

///Refer to (i,j)th MGTextureImage.
const MGTextureImage& operator()(int i, int j)const;

///Refer to (i,j)th MGTextureImage.
MGTextureImage& operator()(int i, int j);

///compute the center of the image in the world (s,t) coordinage.
MGPosition center()const;

///Given stbox, find which MGTextureImage this stbox belongs to.
///Function's return value is:
///0: stbox is outside all the MGTextureImage's.
///1: stbox is inside from (*this)(i,j) to (*this)(i+nwidth-1, j+nheight-1).
///2: stbox is inside from (*this)(i,j) to (*this)(i+nwidth-1, j+nheight-1) when
///	the output base is subtracted from stbox.
///3: a part of stbox is inside total_box and a part is outside total_box.
int compute_image_range(
	const MGBox& stbox,	///<Box of the target triangle.
	int& iTotal, int& jTotal,///<id of the rectangle of this MGTextureImages.
	int& nwTotal, int& nhTotal,///<number of repitition of this MGTextureImages.
	int& i, int&j,		///<id in the 1st MGTextureImages will be output.
	int& nwidth, int& nheight,///<number of repitition of MGTextureImage in the 1st MGTextureImage.
		///<i,j, nwidth, nheight are valid only when return=1 or 2.
	MGVector& base		///<Base(left, bottom) address of the (i,j) MGTextureImage.
)const;

///delete texture objects made by make_texture_object();
void delete_texture_object();

///render process.
void make_texture_object()const;

double pixel_width()const;
double pixel_height()const;

double total_width()const;
double total_height()const;

bool repeatS()const{return m_repeatS;};
bool repeatT()const{return m_repeatT;};

///Append one MGTextureImage.
void push_back(MGTextureImage& image){m_images.push_back(image);};

private:
	bool m_repeatS, m_repeatT;	///<repeat flag for the texture coordinate (s,t) each.
			///<If true, GL_REPEAT will be applied, if false, GL_CLAMP be applied.
	double m_total_width, m_total_height;
			///<The box the texture coordinates of m_images belong to.
			///<From (0,0) to (m_total_width, m_total_height) is the box.
			///<Left bottom address of the m_images must be (0,0).
	int m_total_pixel_w, m_total_pixel_h;///<pixel number of m_total_width and m_total_height;
	int m_texture_size;///<texture size of all the m_images except the remainders.
			///<If m_texture_size==0, m_images is irregular, else m_images is regular.
			///<When m_images is regular, (*this)(i,j)'s nwidth() and nheight() are m_texture_size
			///<for 0<=i<=m_n_columns-2 and 0<=j<=m_n_rows-2.
			///<And (*this)(m_n_columns-1,.)'s nwidth() is m_remainder_w,
			///<(*this)(.,m_n_rows-1)'s nheight() is m_remainder_h.
	double m_1width, m_1height;///<valid only when m_texture_size>0;
			///<indicates m_texture_size's texture_coordinate width and height.
	int m_remainder_w, m_remainder_h;///<m_remainder_w=m_total_pixel_w%m_texture_size. If m_remainder_w!=0,
			///<m_images[m_n_columns-1+m_n_columns*x]'s number of wide pixels are
			///<m_remainder_w for x=0,..., m_n_rows-1.
			///<For m_remainder_h, the same.
	double m_remainder_wd, m_remainder_hd;
			///<m_remainder_wd=m_remainder_w*pixel_width(), m_remainder_hd=m_remainder_h*pixel_height();
	int m_n_columns, m_n_rows;
	std::vector<MGTextureImage> m_images;///<array of MGTextureImage;
			///<m_images.size()=m_n_columns*m_n_rows;

///Build textures after m_repeatS, m_repeatT, and m_texture_size are set up.
void build_textures(
	Gdiplus::Bitmap& bitmap,
	double pixelw, double pixelh///One pixel size, width and height.
);

	friend class mgTLData;
};

void MGTextureImage_get_ijwh(
  double width, double height,
  double s0,double s1,
  double t0,double t1,
  int&i, int&j,
  int& nwidth, int&nheight
 );

/** @} */ // end of GLAttrib group
#endif
