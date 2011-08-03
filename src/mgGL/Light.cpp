/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Position.h"
#include "mg/Transf.h"
#include "mg/Ofstream.h"
#include "mg/Ifstream.h"
#include "mgGL/DirectionalLight.h"
#include "mgGL/SpotLight.h"
#include "mgGL/LightEnable.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

////////////////////////////////////////////////////////////////////

MGLight::MGLight()
:MGGLAttrib(0), m_intensity(1.), m_ambientIntensity(0.){
	for(size_t i=0; i<3; i++) m_color[i]=1.;
}
MGLight::MGLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3]	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
):MGGLAttrib(0), m_intensity(intensity),
 m_ambientIntensity(ambientIntensity){
	for(size_t i=0; i<3; i++) m_color[i]=color[i];
}

//assignment
MGLight& MGLight::operator=(const MGLight& gel2){
	set_light(gel2);
	return *this;
}
MGLight& MGLight::operator=(const MGGel& gel2){
	const MGLight* gel2_is_this=dynamic_cast<const MGLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

MGLight* MGLight::clone()const{
	return new MGLight(*this);
}

//assignment
MGLight& MGLight::set_light(const MGLight& gel2){
	MGGLAttrib::operator=(gel2);

	m_intensity=gel2.m_intensity;
	m_ambientIntensity=gel2.m_ambientIntensity;
	for(size_t i=0; i<3;i++)
		m_color[i]=gel2.m_color[i];
	return *this;
}

bool MGLight::operator<(const MGLight& gel2)const{
	if(m_intensity==gel2.m_intensity)
		if(m_ambientIntensity==gel2.m_ambientIntensity)
			return m_color[0]<gel2.m_color[0];
		else
			return m_ambientIntensity<gel2.m_ambientIntensity;
	else
		return m_intensity<gel2.m_intensity;
}
bool MGLight::operator<(const MGGel& gel2)const{
	const MGLight* gel2_is_this=dynamic_cast<const MGLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

size_t MGLight::exec()const{
	GLenum light=get_light_num();
	size_t i;
	float color[4];color[3]=1.;
	for(i=0; i<3; i++) color[i]=m_color[i]*m_ambientIntensity;
	glLightfv(light,GL_AMBIENT,color);
	for(i=0; i<3; i++) color[i]=m_color[i]*m_intensity;
	glLightfv(light,GL_DIFFUSE,color);
	glLightfv(light,GL_SPECULAR,color);
	glEnable(light);
	return light;
}

// Output function.
std::ostream& MGLight::out(std::ostream& ostrm) const{
	ostrm<<"Lgiht=";
	ostrm<<",intensity="<<m_intensity<<",ambientIntensity="<<m_ambientIntensity;
	ostrm<<std::endl;
	ostrm<<"color=["<<m_color[0]<<","<<m_color[1]<<","<<m_color[2]<<"]";
	return ostrm;
}

std::ostream& operator<< (std::ostream& ostrm, const MGLight& light){
	return light.out(ostrm);
}

/////////// MGDirectionalLight ///////////
/////////// Constructors & destructor ///////////

MGDirectionalLight::MGDirectionalLight()
:MGLight(){
	m_direction[0]=0.;
	m_direction[1]=0.;
	m_direction[2]=1.;
}
MGDirectionalLight::MGDirectionalLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3],	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGVector& direction
):MGLight(intensity,ambientIntensity,color){
	for(size_t i=0;i<3;i++) m_direction[i]=float(direction[i]);
}

//MGDirectionalLight::~MGDirectionalLight(){;}

MGDirectionalLight* MGDirectionalLight::clone()const{
	return new MGDirectionalLight(*this);
}

//assignment
MGDirectionalLight& MGDirectionalLight::operator=(const MGDirectionalLight& gel2){
	if(this==&gel2)
		return *this;

	MGLight::operator=(gel2);
	for(size_t i=0; i<3;)
		m_direction[i]=gel2.m_direction[i];
	return *this;
}
MGDirectionalLight& MGDirectionalLight::operator=(const MGGel& gel2){
	const MGDirectionalLight* gel2_is_this=dynamic_cast<const MGDirectionalLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGDirectionalLight::operator<(const MGDirectionalLight& gel2)const{
	if(MGLight::operator<(gel2))
		return true;
	return m_direction[0]<gel2.m_direction[0];
}
bool MGDirectionalLight::operator<(const MGGel& gel2)const{
	const MGDirectionalLight* gel2_is_this=dynamic_cast<const MGDirectionalLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//exec GLAttribute process.
size_t MGDirectionalLight::exec()const{
	GLenum light=MGLight::exec();
	float position[4];position[3]=0.;
	for(size_t i=0; i<3;i++) position[i]=-m_direction[i];
	glLightfv(light,GL_POSITION,position);
	glLightfv(light,GL_SPOT_DIRECTION,m_direction);
	glLightf(light,GL_SPOT_CUTOFF,180.);
	return light;
}

//////// Set/Get  ///////////

void MGDirectionalLight::setDirection(const MGVector& direction){
	for(size_t i=0; i<3; i++) m_direction[i]=float(direction[i]);
}

void MGDirectionalLight::getDirection(MGVector& direction)const{
	direction.resize(3);
	for(size_t i=0; i<3; i++) direction(i)=m_direction[i];
}

// Output function.
std::ostream& MGDirectionalLight::out(std::ostream& ostrm) const{
	ostrm<<std::endl<<"DirectionalLight="; MGLight::out(ostrm);
	ostrm<<",m_direction=["<<m_direction[0];
	for(size_t i=1; i<3; i++) ostrm<<","<<m_direction[i];
	ostrm<<"]";
	return ostrm;
}

//Transform the gel by the argument.
void MGDirectionalLight::transform(const MGVector& v)//translation
{
	for(size_t i=0; i<3; i++) m_direction[i]+=float(v[i]);
}
void MGDirectionalLight::transform(double scale)//scaling.
{
	for(size_t i=0; i<3; i++) m_direction[i]+=float(scale);
}
void MGDirectionalLight::transform(const MGMatrix& mat)//matrix transformation.
{
	size_t i;
	MGPosition V(3);for(i=0; i<3; i++) V(i)=m_direction[i];
	V*=mat;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}
void MGDirectionalLight::transform(const MGTransf& tr)//general transformation.
{
	size_t i;
	MGPosition V(4);for(i=0; i<4; i++) V(i)=m_direction[i];
	V*=tr;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}

/////////// MGPointLight ///////////
//////// Constructors  /////////////
MGPointLight::MGPointLight()
:MGLight(), m_radius(1.){
	m_location[0]=0.;
	m_location[1]=0.;
	m_location[2]=1.;
	m_location[3]=1.;

	m_attenuation[0]=1.;
	m_attenuation[1]=0.;
	m_attenuation[2]=0.;
}

MGPointLight::MGPointLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3],	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGPosition& location,//location of the light.
	float radius,
	const float attenuation[3]//[0]=GL_CONSTANT_ATTENUATION,
							//[1]=GL_LINEAR_ATTENUATION
							//[2]=GL_QUADRATIC_ATTENUATION
):MGLight(intensity,ambientIntensity,color),m_radius(radius){
	for(size_t i=0;i<3;i++){
		m_location[i]=float(location[i]);
		m_attenuation[i]=float(attenuation[i]);
	}
	m_location[3]=1.;
}

//MGPointLight::~MGPointLight(){;}

//assignment
MGPointLight& MGPointLight::operator=(const MGPointLight& gel2){
	if(this==&gel2)
		return *this;

	MGLight::operator=(gel2);
	for(size_t i=0; i<4;)
		m_location[i]=gel2.m_location[i];
	m_radius=gel2.m_radius;
	for(size_t i=0; i<3;)
		m_attenuation[i]=gel2.m_attenuation[i];
	return *this;
}
MGPointLight& MGPointLight::operator=(const MGGel& gel2){
	const MGPointLight* gel2_is_this=dynamic_cast<const MGPointLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGPointLight::operator<(const MGPointLight& gel2)const{
	if(MGLight::operator<(gel2))
		return true;
	return m_radius<gel2.m_radius;
}
bool MGPointLight::operator<(const MGGel& gel2)const{
	const MGPointLight* gel2_is_this=dynamic_cast<const MGPointLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

MGPointLight* MGPointLight::clone()const{
	return new MGPointLight(*this);
}

//exec GLAttribute process.
size_t MGPointLight::exec()const{
	GLenum light=MGLight::exec();
	glLightfv(light,GL_POSITION,m_location);
	glLightf(light,GL_SPOT_CUTOFF,180.);
	glLightf(light,GL_CONSTANT_ATTENUATION,m_attenuation[0]);
	glLightf(light,GL_LINEAR_ATTENUATION,m_attenuation[1]);
	glLightf(light,GL_QUADRATIC_ATTENUATION,m_attenuation[2]);
	return light;
}

// Output function.
std::ostream& MGPointLight::out(std::ostream& ostrm) const{
	ostrm<<"PointLight=";
	MGLight::out(ostrm);
	ostrm<<std::endl<<
		"m_location=["<<m_location[0]<<","<<m_location[1]<<","<<m_location[2]<<"]";
	ostrm<<", m_attenuation=["<<
		m_attenuation[0]<<","<<m_attenuation[1]<<","<<m_attenuation[2]<<"]";
	ostrm<<", m_radius="<<m_radius;
	return ostrm;
}

//Transform the gel by the argument.
void MGPointLight::transform(const MGVector& v)//translation
{
	for(size_t i=0; i<4; i++) m_location[i]+=float(v[i]);
}
void MGPointLight::transform(double scale)//scaling.
{
	for(size_t i=0; i<4; i++) m_location[i]+=float(scale);
	m_radius*=float(scale);
}
void MGPointLight::transform(const MGMatrix& mat)//matrix transformation.
{
	size_t i;
	MGPosition P(4);for(i=0; i<4; i++) P(i)=m_location[i];
	P*=mat;
	for(i=0; i<4; i++) m_location[i]=float(P[i]);
	m_radius*=float(mat.scale());
}
void MGPointLight::transform(const MGTransf& tr)//general transformation.
{
	size_t i;
	MGPosition P(4);for(i=0; i<4; i++) P(i)=m_location[i];
	P*=tr;
	for(i=0; i<4; i++) m_location[i]=float(P[i]);
	m_radius*=float(tr.affine().scale());
}

////////// Set/Get ///////////

void MGPointLight::setLocation(const MGPosition& location){
	for(size_t i=0; i<3; i++) m_location[i]=float(location[i]);
}
void MGPointLight::getLocation(MGPosition& location)const{
	location.resize(3);
	for(size_t i=0; i<3; i++) location(i)=m_location[i];
}

/////////// MGSpotLight ///////////

////////// Constructors //////////////

MGSpotLight::MGSpotLight()
:m_cutOffAngle(180.), m_exponent(0.){
	m_direction[0]=0.;
	m_direction[1]=0.;
	m_direction[2]=1.;
}

MGSpotLight::MGSpotLight(
    float intensity,		//applied to GL_DIFFUSE and GL_SPECULAR
    float ambientIntensity,	//applied to GL_AMBIENT
    const float color[3],	//applied to GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR
	const MGPosition& location,
	float radius,
	const float attenuation[3],//[0]=GL_CONSTANT_ATTENUATION,
							//[1]=GL_LINEAR_ATTENUATION
							//[2]=GL_QUADRATIC_ATTENUATION
    const MGVector& direction,//GL_SPOT_DIRECTION
    float exponent,		//GL_SPOT_EXPONENT
    float cutOffAngle	//GL_SPOT_CUTOFF
):MGPointLight(intensity,ambientIntensity,color,location,radius,attenuation)
,m_exponent(exponent),m_cutOffAngle(cutOffAngle){
	for(size_t i=0; i<3; i++) m_direction[i]=float(direction[i]);
}

//MGSpotLight::~MGSpotLight(){;}
	
MGSpotLight* MGSpotLight::clone()const{
	return new MGSpotLight(*this);
}

///////////// Set/Get  ///////////

void MGSpotLight::setDirection(const MGVector& direction){
	for(size_t i=0; i<3; i++) m_direction[i]=float(direction[i]);
}
void MGSpotLight::getDirection(MGVector& direction){
	direction.resize(3);
	for(size_t i=0; i<3; i++) direction(i)=m_direction[i];
}

//assignment
MGSpotLight& MGSpotLight::operator=(const MGSpotLight& gel2){
	if(this==&gel2)
		return *this;

	MGLight::operator=(gel2);
	for(size_t i=0; i<3;)
		m_direction[i]=gel2.m_direction[i];
	m_exponent=gel2.m_exponent;
	m_cutOffAngle=gel2.m_cutOffAngle;
	return *this;
}
MGSpotLight& MGSpotLight::operator=(const MGGel& gel2){
	const MGSpotLight* gel2_is_this=dynamic_cast<const MGSpotLight*>(&gel2);
	if(gel2_is_this)
		operator=(*gel2_is_this);
	return *this;
}

bool MGSpotLight::operator<(const MGSpotLight& gel2)const{
	if(MGLight::operator<(gel2))
		return true;
	return m_exponent<gel2.m_exponent;
}
bool MGSpotLight::operator<(const MGGel& gel2)const{
	const MGSpotLight* gel2_is_this=dynamic_cast<const MGSpotLight*>(&gel2);
	if(gel2_is_this)
		return operator<(*gel2_is_this);
	return false;
}

//exec GLAttribute process.
size_t MGSpotLight::exec()const{
	GLenum light=MGPointLight::exec();
	glLightfv(light,GL_SPOT_DIRECTION,m_direction);
	glLightf(light,GL_SPOT_CUTOFF,m_cutOffAngle);
	glLightf(light,GL_SPOT_EXPONENT,m_exponent);
	return light;
}

// Output function.
std::ostream& MGSpotLight::out(std::ostream& ostrm) const{
	ostrm<<"SpotLight=";
	MGPointLight::out(ostrm);
	ostrm<<std::endl<<
		"m_direction=["<<m_direction[0]<<","<<m_direction[1]<<","<<m_direction[2]<<"]";
	ostrm<<", m_exponent="<<m_exponent<<", m_cutOffAngle="<<m_cutOffAngle;
	return ostrm;
}

//Transform the gel by the argument.
void MGSpotLight::transform(const MGVector& v)//translation
{
	MGPointLight::transform(v);
	for(size_t i=0; i<3; i++) m_direction[i]+=float(v[i]);
}
void MGSpotLight::transform(double scale)//scaling.
{
	MGPointLight::transform(scale);
	for(size_t i=0; i<3; i++) m_direction[i]+=float(scale);
}
void MGSpotLight::transform(const MGMatrix& mat)//matrix transformation.
{
	MGPointLight::transform(mat);
	size_t i;
	MGPosition V(3);for(i=0; i<3; i++) V(i)=m_direction[i];
	V*=mat;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}
void MGSpotLight::transform(const MGTransf& tr)//general transformation.
{
	MGPointLight::transform(tr);
	size_t i;
	MGPosition V(4);for(i=0; i<4; i++) V(i)=m_direction[i];
	V*=tr;
	for(i=0; i<4; i++) m_direction[i]=float(V[i]);
}

// Serialization fucntion.
void MGLight::WriteMembers(MGOfstream& buf)const{
	buf<<m_intensity;
	buf<<m_ambientIntensity;
	for(size_t i=0; i<3; i++) buf<<m_color[i];
}
void MGLight::ReadMembers(MGIfstream& buf){
	buf>>m_intensity;
	buf>>m_ambientIntensity;
	for(size_t i=0; i<3; i++) buf>>m_color[i];
}

void MGDirectionalLight::WriteMembers(MGOfstream& buf)const{
	MGLight::WriteMembers(buf);
	for(size_t i=0; i<3; i++) buf<<m_direction[i];
}
void MGDirectionalLight::ReadMembers(MGIfstream& buf){
	MGLight::ReadMembers(buf);
	for(size_t i=0; i<3; i++) buf>>m_direction[i];
}

void MGPointLight::WriteMembers(MGOfstream& buf)const{
	MGLight::WriteMembers(buf);
	size_t i;
	for(i=0; i<4; i++) buf<<m_location[i];
	buf<<m_radius;
	for(i=0; i<3; i++) buf<<m_attenuation[i];
}
void MGPointLight::ReadMembers(MGIfstream& buf){
	MGLight::ReadMembers(buf);
	size_t i;
	for(i=0; i<4; i++) buf>>m_location[i];
	buf>>m_radius;
	for(i=0; i<3; i++) buf>>m_attenuation[i];
}

void MGSpotLight::WriteMembers(MGOfstream& buf)const{
	MGPointLight::WriteMembers(buf);
	for(size_t i=0; i<3; i++) buf<<m_direction[i];
	buf<<m_exponent;
	buf<<m_cutOffAngle;
}
void MGSpotLight::ReadMembers(MGIfstream& buf){
	MGPointLight::ReadMembers(buf);
	for(size_t i=0; i<3; i++) buf>>m_direction[i];
	buf>>m_exponent;
	buf>>m_cutOffAngle;
}

//Get light number available next.
GLenum MGLight::get_light_num()const{
	return GL_LIGHT0+data();
}
