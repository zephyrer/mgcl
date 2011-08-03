/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Implementaion for class MGIgesOfstream.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mg/Point.h"
#include "mg/Straight.h"
#include "mg/Ellipse.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/CompositeCurve.h"
#include "mg/TrimmedCurve.h"
#include "mg/BSumCurve.h"
#include "mg/Plane.h"
#include "mg/Cylinder.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/SurfCurve.h"
#include "topo/PVertex.h"
#include "topo/BVertex.h"
#include "topo/Loop.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "mg/AttribedGel.h"
#include "mg/Group.h"
#include "mgGL/Appearance.h"
#include "mgGL/LineStipple.h"
#include "mgGL/LineWidth.h"
#include "mgGL/Color.h"
#include "mgIges/IgesOfstream.h"
#include "mgIges/IgesPD100.h"
#include "mgIges/IgesPD102.h"
#include "mgIges/IgesPD104.h"
#include "mgIges/IgesPD190.h"
#include "mgIges/IgesPD110.h"
#include "mgIges/IgesPD116.h"
#include "mgIges/IgesPD122.h"
#include "mgIges/IgesPD123.h"
#include "mgIges/IgesPD124.h"
#include "mgIges/IgesPD126.h"
#include "mgIges/IgesPD128.h"
#include "mgIges/IgesPD142.h"
#include "mgIges/IgesPD144.h"
#include "mgIges/IgesPD196.h"
#include "mgIges/IgesPD314.h"
#include "mgIges/IgesPD402.h"
#include "mgIges/IgesPD502.h"
#include "mgIges/IgesPD504.h"
#include "mgIges/IgesPD508.h"
#include "mgIges/IgesPD510.h"
#include "mgIges/IgesPD514.h"
using namespace std;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//SESwitch of Physically and Logically dependent.
static const MGIgesDEStatusNumber::SESwitch PLD=MGIgesDEStatusNumber::PLDependent;

//MGShellIgesOutDataContainer is a private class to output MGShell IGES data.
//Contains MGBVertex(BVertices_of_ModelEdges) pointers, and
//map of VERTEX id in the vector and the pointer.
//Also contains MGEdge pointer and the map data to generate an edge list entity.
class MGShellIgesOutDataContainer{
public:

	///Constructor.
	MGShellIgesOutDataContainer(
		const MGShell& shell,
		MGIgesOfstream& igesfile
	);

	int edgeListDE()const{return m_edgesDE;};
	MGIgesOfstream& igesfile()const{return m_igesfile;};
	int get_Vid(const MGBVertex* bv)const;//get the VERTEX id from m_bvertexmap.
	int get_Eid(const MGEdge* be)const;//get the edge id from m_edgemap.
	
	//Create de of Iges type 502(Vertices list DE) and output to the igesfile.
	void create_BVeticesDE();

	//Create de of Iges type 510(Edges list DE) and output to the igesfile.
	void create_EdgesDE();

private:
	const MGShell& m_shell;//Target shell;
	MGIgesOfstream& m_igesfile;//Output iges file.

	int m_vetrticesDE;//Vertices list DE(PD502) of m_shell(only one for m_shell).
	std::vector<const MGBVertex*> m_bvertices;//to create id of PD502(vertices list).
		//m_bvertices[0] will be dummy since Iges data is indexed so that the 1st id is 1.
	std::map<const MGBVertex*,int> m_bvertexmap;//to get id of MGBVertex* in m_bvertices.

	int m_edgesDE;//Edges list DE(PD504) of m_shell(only one for m_shell).
	std::map<const MGEdge*,int> m_edgemap;
		//to get id of MGEdge* in MGComplex's bcells of this MGShell.
};

//make DE of type 508(LOOP entity).
//Function's return value is th DE number of the entity.
int make_DE508(
   	const MGLoop* loop,	//loop to output
	MGShellIgesOutDataContainer& shellIgesOdata
){
	MGIgesPD508* pd508=new MGIgesPD508();
	MGIgesOfstream& file=shellIgesOdata.igesfile(); //output IGES file.

	int nedges=loop->number_of_edges();
	int el=shellIgesOdata.edgeListDE();
	for(int i=0; i<nedges; i++){
		const MGEdge& ei=*(loop->edge(i));//i-th parameter edge.
		const MGEdge* bei=ei.binder_edge();//the binder of ei.
		int beiDE=shellIgesOdata.get_Eid(bei);
		short orientation=ei.equal_direction_to_binder() ? 1:0;
		MGIges508Edge* edgei=new MGIges508Edge(0,el,beiDE,orientation);
		edgei->m_isoparameterics.push_back(false);
		
		std::auto_ptr<MGCurve> eiCurve(ei.curve_limitted());
		int eiCurveDE=eiCurve->out_to_IGES(file,PLD);
		edgei->m_pcurves.push_back(eiCurveDE);
		pd508->push_back(edgei);
	}
	return file.create_de(pd508,"LOOP",PLD);
}

//make DE of type 510(FACE entity).
//Function's return value is th DE number of the entity.
int make_DE510(
   	const MGFace* face,	//face to output
	MGShellIgesOutDataContainer& shellIgesOdata
){
	//std::cout<<*face<<std::endl;
	MGIgesPD510* pd510=new MGIgesPD510();
	MGIgesOfstream& file=shellIgesOdata.igesfile(); //output IGES file.

	const MGSurface& srf=*(face->surface());
	pd510->m_surface_DE=srf.out_to_IGES(file,PLD);
	int nloops=face->number_of_loops();
	for(int i=0; i<nloops; i++){
		pd510->m_loops.push_back(make_DE508(face->loop(i),shellIgesOdata));
	}
	return file.create_de(pd510,"FACE",PLD);
}

void MGIgesFstream::set_initial_StartSection(){
	m_StartSection="MGCL Version:";
	m_StartSection.append(MGCL_Version());
	m_StartSection+="/IGES 5.03";
	m_StartSection+="; DG Technologies, Inc. Kiriyama Buildng 2F, 2-1-1 Sarugaku-cho, ";
	m_StartSection+="Chiyoda-ku, Tokyo 101-0064, JAPAN";
}

// Constructors.

// Creates an object of class MGIgesOfstream with a filename.
MGIgesOfstream::MGIgesOfstream(const char* filename){
	open(filename);
}

// Destroys an object of class MGIgesOfstream.
MGIgesOfstream::~MGIgesOfstream(){
	if(is_open())
		close();
}

// Writes all objects of the MGGroup stored in group
// as IGES objects.
MGIgesOfstream& MGIgesOfstream::operator<<(const MGGel& gel){
	 gel.out_to_IGES(*this);
	 return *this;
}

// Writes all objects of the MGGroup stored in group
// as IGES objects.
MGIgesOfstream& MGIgesOfstream::operator<<(const MGGroup& group){
	MGGroup::const_iterator itr = group.begin();
	for(; itr != group.end(); ++itr){
		(*this)<<(**itr);
	}
	return *this;
}

//Open the file. When this is opened already, this will be closed, then will
//be opened.
void MGIgesOfstream::open(const char* filename){
	if(is_open())
		close();
	if(filename)
		m_ofstream.open(filename);
	initialize(filename);
}

//Initialize all the member data to the default value.
void MGIgesOfstream::initialize(const char* filename){
	MGIgesFstream::initialize(filename);
	set_initial_StartSection();

	//1. Write out start section.
	write_out_start_section();

	//2. Write out global section.
	m_nlineGSec=m_GSection.write_out(*this);
}

void MGIgesOfstream::close(){
	write_out_DE_PD_lines();
	write_out_terminate_section();
	m_ofstream.close();
}

//Create de as a newed object, and push back to this MGIgesOfstream.
//The ownership of pd is transfered to the created DE, and the ownership of
//the created DE will be transfered to this MGIgesOfstream.
//Function's return value is the DE id created. If DE pointer is necessary,
//directoryEntry() is the one.
int MGIgesOfstream::create_de(
	MGIgesPD* pd,//newed MGIgesPD object.
	const std::string& EntityLabel,
	int ses,//Subordinate entity switch.
	const MGAttribedGel* gel,//when gel=0 is input, De will be set as dependent.
	int FormNumber	//Form number
){
	int ColorNumber=0;//When negated value, is a pointer to the directory entry.
	int LineWeightNumber=0;
	int LineFontPattern=0;//When negated value, is a pointer to the directory entry.
	MGIgesDEStatusNumber sn;
	if(gel){
		const MGAppearance* appr=gel->appearance();
		if(appr){
			if(appr->no_display())
				sn.set_as_blank();
			MGAppearance::const_iterator ie=appr->end();

			//Color
			MGAppearance::const_iterator icolor=appr->search_by_id(MGCOLOR_TID);
			if(icolor!=ie){//If found.
				const MGColor* clr=static_cast<MGColor*>(*icolor);
				int colorID=clr->get_ColorID(8);
				if(colorID)
					ColorNumber=colorID;
				else{
					ColorNumber=clr->out_to_IGES(*this,PLD);
					ColorNumber*=-1;
				}
			}

			//Line width
			MGAppearance::const_iterator ilwidth=appr->search_by_id(MGLINE_WIDTH_TID);
			if(ilwidth!=ie){//If found.
				const MGLineWidth* width=static_cast<const MGLineWidth*>(*ilwidth);
				double w=width->get_width();
				LineWeightNumber=int(w+.4);
			}

			//Line Stipple
			MGAppearance::const_iterator ilstipple=appr->search_by_id(MGLINE_STIPPLE_TID);
			if(ilstipple!=ie){//If found.
				const MGLineStipple* font=static_cast<const MGLineStipple*>(*ilstipple);
				MGLineStipple::LineFont lf=font->get_font_number();
				if(lf!=MGLineStipple::UndefinedFont)
					LineFontPattern=int(lf);
			}
		}
	}
	sn.set_SubordinateEntitySwitch(MGIgesDEStatusNumber::SESwitch(ses));//set subordinate entity switch.

	int EntityTypeNumber=pd->type_number();
	MGIgesDirectoryEntry* de=new MGIgesDirectoryEntry(
		EntityTypeNumber,EntityLabel,sn,ColorNumber,LineWeightNumber,LineFontPattern,0,FormNumber);
	de->setPD(auto_ptr<MGIgesPD>(pd));
	return push_back_DE(de);
}

//PD116=POINT.
//Function's return value is the directory entry id created.
int MGPosition::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD116* pd116=new MGIgesPD116(*this);
	return igesfile.create_de(pd116,"POSITION",PLD);
}

//PD123=DIRECTION.
//Function's return value is the directory entry id created.
int MGVector::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD123* pd123=new MGIgesPD123(*this);
	return igesfile.create_de(pd123,"VECTOR",PLD);
}

//PD124=Transformation matrix.
//Function's return value is the directory entry id created.
int MGTransf::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	return igesfile.create_de(new MGIgesPD124(*this), "TRANSF",PLD);
}

//PD116=POINT.
//Function's return value is the directory entry id created.
int MGPoint::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD116* pd116=new MGIgesPD116(*this);
	return igesfile.create_de(pd116,"POINT",SubordinateEntitySwitch,this);
}

//PD110=Straight line.
//Function's return value is the directory entry id created.
int MGStraight::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGPosition P1=start_point(), P2=end_point();
	int formnumber=1;
	if(infinite_above()){
		if(infinite_below()){
			P1=root_point();
			P2=P1+MGUnit_vector(direction());
			formnumber=2;
		}else{
			P2=P1+MGUnit_vector(direction());
		}
	}else if(infinite_below()){
		P1=P2;
		P2=P1-MGUnit_vector(direction());
	}else
		formnumber=0;

	MGIgesPD110* pd110=new MGIgesPD110(P1,P2);
	return igesfile.create_de(pd110,"STRAIGHT",	SubordinateEntitySwitch,this,formnumber);
}

//PD100=circular arc.
//Function's return value is the directory entry id created.
int MGEllipse::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGMatrix mat; mat.set_xy_vector(major_axis(),minor_axis());
	MGTransf tran(mat,center());
	int id_tr_de=tran.out_to_IGES(igesfile,PLD);

	double ts=gp_to_radian(param_s());
	double cosvs=cos(ts), sinvs=sin(ts);
	double te=gp_to_radian(param_e());
	double cosve=cos(te), sinve=sin(te);
	MGIgesPD* pd;
	if(circle()){
		double r=radius();
		double center[2]={0.,0.};
		double start[2]={cosvs*r,sinvs*r};
		double terminate[2]={cosve*r,sinve*r};
		pd=new MGIgesPD100(center,start,terminate);
	}else{
		double a=m_m.len(), b=m_n.len();
		double coef[6]={1./(a*a), 0., 1./(b*b), 0., 0., -1.};
		double start[2]={cosvs*a,sinvs*b};
		double terminate[2]={cosve*a,sinve*b};
		pd=new MGIgesPD104(coef,0.,start,terminate);
	}
	int deID=igesfile.create_de(pd,"ELLIPSE",SubordinateEntitySwitch,this);
	igesfile.directoryEntry(deID)->setTransformID(id_tr_de);
	return deID;
}

//IGES output function. PD126=LBRep.
//Function's return value is the directory entry id created.
int MGLBRep::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD126* pd126=new MGIgesPD126(*this);
	return igesfile.create_de(pd126,"LBREP",SubordinateEntitySwitch,this);
}

//IGES output function. PD126=RLBRep.
//Function's return value is the directory entry id created.
int MGRLBRep::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD126* pd126=new MGIgesPD126(*this);
	return igesfile.create_de(pd126,"RLBREP",SubordinateEntitySwitch,this);
}

//PD102=Composite curve.
//Function's return value is the directory entry id created.
int MGCompositeCurve::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	int n=number_of_curves();
	MGIgesPD102* pd102=new MGIgesPD102();
	for(int i=0; i<n; i++){
		const MGCurve& ci=curve(i);
		int cide_id=ci.out_to_IGES(igesfile,PLD);
		if(!cide_id)
			continue;
		pd102->append_curve(cide_id);
	}
	return igesfile.create_de(pd102,"COMPOSIT",SubordinateEntitySwitch,this);
}

//IGES output function.
//Function's return value is the directory entry id created.
int MGTrimmedCurve::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	std::auto_ptr<MGCurve> curve(copy_limitted());
	return curve->out_to_IGES(igesfile,SubordinateEntitySwitch);																																																																			
}

//IGES output function. PD126(approximated as a NURBS line).
//Function's return value is the directory entry id created.
int MGSurfCurve::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGLBRep lb(*this);
	MGIgesPD126* pd126=new MGIgesPD126(lb);
	return igesfile.create_de(pd126,"SRFCURVE",SubordinateEntitySwitch,this);
}

//Output to IGES stream file.
//BSumCurve is approximated as MGLBRep and output as PD126.
//Function's return value is the directory entry id created.
int MGBSumCurve::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGLBRep lb(*this);
	MGIgesPD126* pd126=new MGIgesPD126(lb);
	return igesfile.create_de(pd126,"BSUMCURV",SubordinateEntitySwitch,this);
}

//Output to IGES stream file(PLANE=PD190).
//Function's return value is the directory entry id created.
int MGPlane::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	//make origin DE.
	const MGPosition& location=root_point();
	int locationDE=location.out_to_IGES(igesfile,PLD);

	//make normal DE.
	const MGVector& norm=normal();
	int normalDE=norm.out_to_IGES(igesfile,PLD);

	//make reference DE.
	const MGVector& udir=u_deriv();
	int uderiDE=udir.out_to_IGES(igesfile,PLD);

	//make plane(of pd190) DE.
	MGIgesPD190* pd190=new MGIgesPD190(locationDE,normalDE,uderiDE);
	return igesfile.create_de(pd190,"PLANE",SubordinateEntitySwitch,this,1);//Form number is 1.
}

//PD122=Cylinder.
//Function's return value is the directory entry id created.
int MGCylinder::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	double us=param_s_u(), ve=param_e_v();
	MGPosition Pusve=eval(us,ve);
	std::auto_ptr<MGCurve> vs_curve(perimeter_curve(0));//curve of perimeter 0, i.e. v=v-min.
	int vs_curve_de=vs_curve->out_to_IGES(igesfile,PLD);
	MGIgesPD122* pd122=new MGIgesPD122(vs_curve_de,Pusve.data());//Tabulated Cylinder.
	return igesfile.create_de(pd122, "CYLINDER",SubordinateEntitySwitch,this);
}

//IGES output function. PD128(NURBS Surface).
//Function's return value is the directory entry id created.
int MGSBRep::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD128* pd128=new MGIgesPD128(*this);
	return igesfile.create_de(pd128,"SBREP",SubordinateEntitySwitch,this);
}

//IGES output function. PD128(NURBS Surface).
//Function's return value is the directory entry id created.
int MGRSBRep::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD128* pd128=new MGIgesPD128(*this);
	return igesfile.create_de(pd128,"RSBREP",SubordinateEntitySwitch,this);
}

//Make PD142(Curve on a surface). Function's return value is
//the DE number.
int make_curve_on_surface(
	const MGLoop& loop,	//loop to make curve on a surface.
	int baseSurfaceDE,	//The base surface loop lies on.
	MGIgesOfstream& igesfile//Iges file to output.
){
	MGIgesPD142* pd142=new MGIgesPD142(loop,baseSurfaceDE,igesfile);
	return igesfile.create_de(pd142,"CRVONSRF",PLD);
}

//PD196=Spherical surface(parameterized).
//Function's return value is the directory entry id created.
int MGSphere::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	if(m_ellipseu.circle() && m_ellipsev.circle()
		&& m_ellipseu.is_whole_ellipse() && m_ellipsev.is_whole_ellipse()){

		//make origin DE.
		int centerDE=sphere_center().out_to_IGES(igesfile,PLD);

		//make axis DE.
		int axisDE=B().out_to_IGES(igesfile,PLD);

		//make reference vector DE.
		int refdirDE=M().out_to_IGES(igesfile,PLD);

		//make Spherical surface(of pd196) DE.
		MGIgesPD196* pd196=new MGIgesPD196(centerDE,radius(),axisDE,refdirDE);
		return igesfile.create_de(pd196,"SPHERE",SubordinateEntitySwitch,this);
	}

	return 0;
}

//PD144=MGFace. Output to PD144(Trimmed surface).
//Function's return value is the directory entry id created.
int MGFace::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	const MGFace* f=this;
	std::auto_ptr<MGFace> f2;
	if(hasInactiveLoop() || !hasOuterBoundaryLoop()){
		f2=std::auto_ptr<MGFace>(new MGFace(*this));
		f2->remove_inactive_loops();
		f2->make_outer_boundary();
		f=f2.get();
	}

	//std::cout<<(*f)<<std::endl;
	const MGSurface& srf=*(f->surface());
	int surfaceDE=srf.out_to_IGES(igesfile,PLD);
	const MGLoop& oloop=*(f->loop(size_t(0)));
	int outloopDE=make_curve_on_surface(oloop,surfaceDE,igesfile);

	MGIgesPD144* pd144=new MGIgesPD144(surfaceDE,outloopDE);
	int ninner=f->number_of_inner_boundaries();
	for(int i=0; i<ninner; i++){
		int inneri=make_curve_on_surface(*(f->loop(i+1)),surfaceDE,igesfile);
		pd144->append_inner_boundary(inneri);
	}

	return igesfile.create_de(pd144,"FACE",SubordinateEntitySwitch,this);
}

MGShellIgesOutDataContainer::MGShellIgesOutDataContainer(
	const MGShell& shell,
	MGIgesOfstream& igesfile
):m_shell(shell),m_igesfile(igesfile),m_bvertices(1){
	shell.ensure_BVertices_of_ModelEdges();
	create_BVeticesDE();
	create_EdgesDE();
}

int MGShellIgesOutDataContainer::get_Vid(
	const MGBVertex* bv
)const{
	std::map<const MGBVertex*,int>::const_iterator i=m_bvertexmap.find(bv);
	assert(i!=m_bvertexmap.end());
	if(i==m_bvertexmap.end())
		return 0;
	return i->second;
}

int MGShellIgesOutDataContainer::get_Eid(
	const MGEdge* be
)const{
	std::map<const MGEdge*,int>::const_iterator i=m_edgemap.find(be);
	assert(i!=m_edgemap.end());
	if(i==m_edgemap.end())
		return 0;
	return i->second;
}

//Create de of Iges type 502(Vertices list DE) and output to the igesfile.
void MGShellIgesOutDataContainer::create_BVeticesDE(
){
	MGIgesPD502* pd502=new MGIgesPD502();
	int nedges=m_shell.number_of_bcells();
	MGShell::const_bcellItr j=m_shell.bcell_begin(), jend=m_shell.bcell_end();
	int idV=1;//idV=1, because Iges id starts from 1.
	for(; j!=jend; ++j){//idV=1, because Iges id starts from 1.
		const MGCellNB* cj=*j;
		const MGBVertex* vertj=dynamic_cast<const MGBVertex*>(cj);
		if(!vertj)
			continue;

		pd502->push_back(vertj->position());
		pair<std::map<const MGBVertex*,int>::iterator, bool> insertR=
			m_bvertexmap.insert(make_pair(vertj,idV++));
		assert(insertR.second==true);
	}
	m_vetrticesDE=m_igesfile.create_de(pd502,"BVERTCES",PLD,0,1);
}

void MGShellIgesOutDataContainer::create_EdgesDE(
){
	MGIgesPD504* pd504=new MGIgesPD504();
	int nedges=m_shell.number_of_bcells();
	MGShell::const_bcellItr j=m_shell.bcell_begin(), jend=m_shell.bcell_end();
	int idE=1;//idE=1, because Iges id starts from 1.
	for(; j!=jend; ++j){//idE=1, because Iges id starts from 1.
		const MGCellNB* cj=*j;
		const MGEdge* edgej=dynamic_cast<const MGEdge*>(cj);
		if(!edgej)
			continue;
		
		pair<std::map<const MGEdge*,int>::iterator, bool> insertR=
			m_edgemap.insert(make_pair(edgej,idE++));
		assert(insertR.second==true);

		std::auto_ptr<MGCurve> ecurve(edgej->curve_limitted());
		int curveDE=ecurve->out_to_IGES(m_igesfile,PLD);

		const MGBVertex* sbv=edgej->vertex_start()->binder_vertex();
		int sv=get_Vid(sbv);
		const MGBVertex* ebv=edgej->vertex_end()->binder_vertex();
		int tv=get_Vid(ebv);
		MGIges504Edge e504(curveDE,m_vetrticesDE,sv,m_vetrticesDE,tv);
		pd504->push_back(e504);
	}
	m_edgesDE=m_igesfile.create_de(pd504,"EDGES",PLD,0,1);
}

//PD514=MGShell. Output to PD514(SHELL).
//Function's return value is the directory entry id created.
int MGShell::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGShellIgesOutDataContainer shellIgesOdata(*this,igesfile);
	MGIgesPD514* pd514=new MGIgesPD514();

	int nfaces=number_of_faces();
	for(int i=0; i<nfaces; i++)
		pd514->push_back(make_DE510(face(i),shellIgesOdata));
	return igesfile.create_de(pd514,"SHELL",SubordinateEntitySwitch,this);
}

//IGES output function. PD402=Group.
//Function's return value is the directory entry id created.
int MGGroup::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	MGIgesPD402* pd402=new MGIgesPD402();

	MGGroup::const_iterator itr = begin();
	for(; itr != end(); ++itr){
		const MGGel& geli=**itr;
		int dei=geli.out_to_IGES(igesfile,MGIgesDEStatusNumber::PDependent);
		pd402->append_DE(dei);
	}
	return igesfile.create_de(pd402,"GROUP",SubordinateEntitySwitch,this);
}

//Output to IGES stream file(Color=PD314).
//Function's return value is the directory entry id created.
int MGColor::out_to_IGES(
	MGIgesOfstream& igesfile,
	int SubordinateEntitySwitch
)const{
	//make color(of pd314) DE.
	MGIgesPD314* pd314=new MGIgesPD314(*this);
	return igesfile.create_de(pd314,"COLOR",PLD);
}
