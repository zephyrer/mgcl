/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include <bitset>
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/Color.h"
#include "mgGL/RenderAttr.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// Color data.
// This is completely compatible to GDI+ color enumuration.
const static size_t mgColorValue[MGColor::endID]={
    0xFFFFFFFF,//WhiteV        
    0xFF000000,//BlackV        
    0xFFFF0000,//RedV          
    0xFF008000,//GreenV        
    0xFF0000FF,//BlueV         
    0xFFFFFF00,//YellowV       
    0xFFFF00FF,//MagentaV      
    0xFF00FFFF,//CyanV         
    0xFFFFFFFF,//WhiteV        
    0xFFF0F8FF,//AliceBlueV    
    0xFFFAEBD7,//AntiqueWhiteV 
    0xFF00FFFF,//AquaV         
    0xFF7FFFD4,//AquamarineV   
    0xFFF0FFFF,//AzureV        
    0xFFF5F5DC,//BeigeV        
    0xFFFFE4C4,//BisqueV       
    0xFFFFEBCD,//BlanchedAlmondV 
    0xFF8A2BE2,//BlueVioletV     
    0xFFA52A2A,//BrownV          
    0xFFDEB887,//BurlyWoodV      
    0xFF5F9EA0,//CadetBlueV      
    0xFF7FFF00,//ChartreuseV     
    0xFFD2691E,//ChocolateV      
    0xFFFF7F50,//CoralV          
    0xFF6495ED,//CornflowerBlueV 
    0xFFFFF8DC,//CornsilkV       
    0xFFDC143C,//CrimsonV        
    0xFF00008B,//DarkBlueV       
    0xFF008B8B,//DarkCyanV       
    0xFFB8860B,//DarkGoldenrodV  
    0xFFA9A9A9,//DarkGrayV             Gray1
    0xFF006400,//DarkGreenV      
    0xFFBDB76B,//DarkKhakiV      
    0xFF8B008B,//DarkMagentaV    
    0xFF556B2F,//DarkOliveGreenV 
    0xFFFF8C00,//DarkOrangeV     
    0xFF9932CC,//DarkOrchidV     
    0xFF8B0000,//DarkRedV        
    0xFFE9967A,//DarkSalmonV     
    0xFF8FBC8B,//DarkSeaGreenV   
    0xFF483D8B,//DarkSlateBlueV  
    0xFF2F4F4F,//DarkSlateGrayV  
    0xFF00CED1,//DarkTurquoiseV  
    0xFF9400D3,//DarkVioletV     
    0xFFFF1493,//DeepPinkV       
    0xFF00BFFF,//DeepSkyBlueV    
    0xFF696969,//DimGrayV        
    0xFF1E90FF,//DodgerBlueV     
    0xFFB22222,//FirebrickV      
    0xFFFFFAF0,//FloralWhiteV    
    0xFF228B22,//ForestGreenV    
    0xFFFF00FF,//FuchsiaV        
    0xFFDCDCDC,//Gainsboro			//Gray2
    0xFFF8F8FF,//GhostWhiteV     
    0xFFFFD700,//GoldV           
    0xFFDAA520,//GoldenrodV      
    0xFF808080,////Gray3GrayV    
    0xFFADFF2F,//GreenYellowV    
    0xFFF0FFF0,//HoneydewV       
    0xFFFF69B4,//HotPinkV        
    0xFFCD5C5C,//IndianRedV      
    0xFF4B0082,//IndigoV         
    0xFFFFFFF0,//IvoryV          
    0xFFF0E68C,//KhakiV          
    0xFFE6E6FA,//LavenderV       
    0xFFFFF0F5,//LavenderBlushV  
    0xFF7CFC00,//LawnGreenV      
    0xFFFFFACD,//LemonChiffonV   
    0xFFADD8E6,//LightBlueV      
    0xFFF08080,//LightCoralV     
    0xFFE0FFFF,//LightCyanV      
    0xFFFAFAD2,//LightGoldenrodYe
    0xFFD3D3D3,//LightGrayV      
    0xFF90EE90,//LightGreenV     
    0xFFFFB6C1,//LightPinkV      
    0xFFFFA07A,//LightSalmonV    
    0xFF20B2AA,//LightSeaGreenV  
    0xFF87CEFA,//LightSkyBlueV   
    0xFF778899,//LightSlateGrayV 
    0xFFB0C4DE,//LightSteelBlueV 
    0xFFFFFFE0,//LightYellowV    
    0xFF00FF00,//LimeV           
    0xFF32CD32,//LimeGreenV      
    0xFFFAF0E6,//LinenV          
    0xFF800000,//MaroonV         
    0xFF66CDAA,//MediumAquamarine
    0xFF0000CD,//MediumBlueV     
    0xFFBA55D3,//MediumOrchidV   
    0xFF9370DB,//MediumPurpleV   
    0xFF3CB371,//MediumSeaGreenV 
    0xFF7B68EE,//MediumSlateBlueV
    0xFF00FA9A,//MediumSpringGree
    0xFF48D1CC,//MediumTurquoiseV
    0xFFC71585,//MediumVioletRedV
    0xFF191970,//MidnightBlueV   
    0xFFF5FFFA,//MintCreamV      
    0xFFFFE4E1,//MistyRoseV      
    0xFFFFE4B5,//MoccasinV       
    0xFFFFDEAD,//NavajoWhiteV    
    0xFF000080,//NavyV           
    0xFFFDF5E6,//OldLaceV        
    0xFF808000,//OliveV          
    0xFF6B8E23,//OliveDrabV      
    0xFFFFA500,//OrangeV         
    0xFFFF4500,//OrangeRedV      
    0xFFDA70D6,//OrchidV         
    0xFFEEE8AA,//PaleGoldenrodV  
    0xFF98FB98,//PaleGreenV      
    0xFFAFEEEE,//PaleTurquoiseV  
    0xFFDB7093,//PaleVioletRedV  
    0xFFFFEFD5,//PapayaWhipV     
    0xFFFFDAB9,//PeachPuffV      
    0xFFCD853F,//PeruV           
    0xFFFFC0CB,//PinkV           
    0xFFDDA0DD,//PlumV           
    0xFFB0E0E6,//PowderBlueV     
    0xFF800080,//PurpleV         
    0xFFBC8F8F,//RosyBrownV      
    0xFF4169E1,//RoyalBlueV      
    0xFF8B4513,//SaddleBrownV    
    0xFFFA8072,//SalmonV         
    0xFFF4A460,//SandyBrownV     
    0xFF2E8B57,//SeaGreenV       
    0xFFFFF5EE,//SeaShellV       
    0xFFA0522D,//SiennaV         
    0xFFC0C0C0,//SilverV         
    0xFF87CEEB,//SkyBlueV        
    0xFF6A5ACD,//SlateBlueV      
    0xFF708090,//SlateGrayV      
    0xFFFFFAFA,//SnowV           
    0xFF00FF7F,//SpringGreenV    
    0xFF4682B4,//SteelBlueV      
    0xFFD2B48C,//TanV            
    0xFF008080,//TealV           
    0xFFD8BFD8,//ThistleV        
    0xFFFF6347,//TomatoV         
    0x00FFFFFF,//TransparentV    
    0xFF40E0D0,//TurquoiseV      
    0xFFEE82EE,//VioletV         
    0xFFF5DEB3,//WheatV          
    0xFFF5F5F5,//WhiteSmokeV		//gray4
    0xFF9ACD32 //YellowGreenV        
};
			   
MGColors::MGColors(){
	MGColors::m_colors=new MGColor*[MGColor::endID+1];
	MGColors::m_colors[0]=new MGColor(mgColorValue[MGColor::White]);
	for(size_t i=1; i<=MGColor::endID; i++)
		MGColors::m_colors[i]=new MGColor(mgColorValue[i]);
}

MGColors::~MGColors(){
	for(int i=0; i<=MGColor::endID; i++)
		delete m_colors[i];
	delete[] m_colors;
}

//get the color id of GDI+ color enumeration ID.
//If not found, 0 will be returned.
int MGColor::get_ColorID(int maxID)const{
	if(maxID<=0 || maxID>=endID)
		maxID=endID-1;

	size_t R=size_t(m_color[0]*255.);R<<=16;
	size_t G=size_t(m_color[1]*255.);G<<=8;
	size_t B=size_t(m_color[2]*255.);
	size_t IColor=0xFF000000 | R | G | B;

	for(int i=1; i<=maxID; i++){
		if(mgColorValue[i] == IColor)
			return i;
	}
	return 0;
}

//////////////////MGColor//////////////////

//Get the color instance reference by the color id.
const MGColor& MGColor::get_instance(MGColor::ColorID id){
	assert(id<MGColor::endID);
	return m_colors.color(id);
}

//Get unsigned integer valur of a color id.
//Function's return value is (A, R, G, B) data of GDI+.
size_t MGColor::get_ARGBinstance(MGColor::ColorID id){
	assert(id<MGColor::endID);
	return mgColorValue[id];
}

MGColor::MGColor(float red, float green, float blue, float alpha)
:MGGLAttrib(ENABLED){
	m_color[0]=red;
	m_color[1]=green;
	m_color[2]=blue;
	m_color[3]=alpha;
}

// Construct a color from ARGB value. In argb, each 8bits of (a,r,g,b) is the value of
//the range 0 to 255 .
MGColor::MGColor(unsigned int argb):MGGLAttrib(ENABLED){
	m_color[0] = (((unsigned char)(argb >> 16))) / 255.f;//R
	m_color[1] = (((unsigned char)(argb >> 8))) / 255.f;//G
	m_color[2] = (((unsigned char)(argb >> 0))) / 255.f;//B
	m_color[3] = (((unsigned char)(argb >> 24))) / 255.f;//A
}

MGColor* MGColor::clone()const{
	return new MGColor(*this);
}
	
//assignment
MGColor& MGColor::operator=(const MGColor& gel2){
	if(this==&gel2)
		return *this;

	MGGLAttrib::operator=(gel2);
	for(size_t i=0; i<4; i++)
		m_color[i]=gel2.m_color[i];
	return *this;
}
MGColor& MGColor::operator=(const MGGel& gel2){
	const MGColor* gel2_is_this=dynamic_cast<const MGColor*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGColor::operator<(const MGColor& gel2)const{
	if(m_flag==gel2.m_flag){
		if(m_color[0]==gel2.m_color[0]){
			if(m_color[1]==gel2.m_color[1]){
				if(m_color[2]==gel2.m_color[2])
					return m_color[3]<gel2.m_color[3];
				else
					return m_color[2]<gel2.m_color[2];
			}else
					return m_color[1]<gel2.m_color[1];
		}else
			return m_color[0]<gel2.m_color[0];
	}else
		return m_flag<gel2.m_flag;
}
bool MGColor::operator<(const MGGel& gel2)const{
	const MGColor* gel2_is_this=dynamic_cast<const MGColor*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//Scaling the color values by the factor scale.
//All of the elements except trancparency element will be multiplied
//by the scale and be clamped between 0. and 1.
MGColor& MGColor::operator*=(float scale){
	for(size_t i=0; i<3; i++){
		float& colori=m_color[i];
		colori*=scale;
		if(colori<0.)
			colori=0.;
		if(colori>1.)
			colori=1.;
	}
	return *this;
}

//Add a color values to RGB data..
//All of the elements except trancparency element will be added by value
//and be clamped between 0. and 1.
MGColor& MGColor::operator+=(float value){
	for(size_t i=0; i<3; i++){
		float& colori=m_color[i];
		colori+=value;
		if(colori<0.)
			colori=0.;
		if(colori>1.)
			colori=1.;
	}
	return *this;
}

// Output function.
std::ostream& MGColor::out(std::ostream& ostrm) const{
	ostrm<<"Color="; MGGLAttrib::out(ostrm);
	if(undefined()) return ostrm;

	ostrm<<"=";
	const float* colr=color();
	ostrm<<"color=["<<colr[0];
	for(size_t i=1; i<4; i++) ostrm<<","<<colr[i];
	ostrm<<"]";
	return ostrm;
}

//draw GLAttribute process.
void MGColor::drawAttrib(
	bool no_color	//if true, color attribute will be neglected.
)const{
	if(no_color)
		return;
	exec();
}

void MGColor::exec()const{
	if(undefined())
		return;
	glColor4fv(color());
}

void MGColor::set_color(const float color[4]){
	for(size_t i=0; i<4; i++) m_color[i]=color[i];
	m_flag=	ENABLED;
}

void MGColor::set_color(float red, float green, float blue, float alpha){
	m_color[0]=red;
	m_color[1]=green;
	m_color[2]=blue;
	m_color[3]=alpha;
	m_flag=	ENABLED;
}

void MGColor::get_color(float color[4])const{
	for(size_t i=0; i<4; i++) color[i]=m_color[i];
}
void MGColor::get_color(float& red, float& green, float& blue, float& alpha)const{
	red=m_color[0];
	green=m_color[1];
	blue=m_color[2];
	alpha=m_color[3];
}

void MGColor::WriteMembers(MGOfstream& buf)const{
	MGGLAttrib::WriteMembers(buf);
	if(undefined()) return;
	for(size_t i=0; i<4; i++) buf<<m_color[i];
}
void MGColor::ReadMembers(MGIfstream& buf){
	MGGLAttrib::ReadMembers(buf);
	if(undefined()) return;
	for(size_t i=0; i<4; i++) buf>>m_color[i];
}
