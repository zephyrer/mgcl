/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Tolerance.h"
#include "mg/Position.h"
#include "Tl/TLFan.h"
#include "Tl/TLFans.h"
#include "Tl/TLTriangle.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

//////////// mgTLFans///////////

////////// private class //////////

//エッジクラス
class mgTLEdge{
private:
	size_t m_start, m_end;
	//Edge's start and end ids of mgTLPoints(mgTLparamter.m_points).
public:
	mgTLEdge(){;};
	mgTLEdge(size_t start, size_t end):m_start(start),m_end(end){;};
	size_t start()const {return m_start;};
	size_t end()const {return m_end;};
};

class mgTLEdges{
private:
	std::deque<mgTLEdge> m_edges;
public:
	typedef std::deque<mgTLEdge>::const_iterator eitr;
	mgTLEdges(){;};
	eitr begin()const{return m_edges.begin();};
	eitr end()const{return m_edges.end();};
	bool empty()const{return m_edges.empty();};
	void pop_back(){m_edges.pop_back();};
	void push_back(size_t start, size_t end){
		m_edges.push_back(mgTLEdge(start,end));
	};
	void push_front(size_t start, size_t end){
		m_edges.push_front(mgTLEdge(start,end));
	};
	void push_front(const mgTLEdge& edge){m_edges.push_front(edge);};
};

std::ostream& operator<< (std::ostream& out, const mgTLEdges& edges){
	mgTLEdges::eitr i=edges.begin(), ie=edges.end();
	for(size_t j=0; i!=ie; i++, j++){
		out<<" "<<j<<"("<<(*i).start()<<","<<(*i).end()<<")";
	}
	return out;
}

//ファンを作成する際に使用するステータス
enum mgTRIANG_STATUS{
	UNKNOWN,	//0
	TWOTOUCH,	//1
	NOTOUCH,	//2
	REGULAR,	//3
	RIGHTTOUCH,	//4
	LEFTTOUCH	//5
};

///////////constructor////////////

mgTLFans::mgTLFans(
	double uerror, double verror,//Error along u and v in the parameter space.
	const mgTLPoints& tlPoints,
	mgTLTriangle& polygon)
:m_uerror(uerror),m_verror(verror),m_points(tlPoints),m_polygon(polygon){
	size_t num_vert = polygon.size(), i=0;//頂点の数を数える
	if(num_vert < 3)
		return;	//一応チェックしておく

	//多角形の辺をスタックに積む
	mgTLEdges edgeStack;		//スタックエッジ
	init(edgeStack);
	//cout<<edgeStack<<endl;//////////////
	if(num_vert==3)
		return;

	//スタックがからになるまで、頂点と周辺の頂点リストのベクトルを求める
	size_t nv2=num_vert*2;
	for(; !edgeStack.empty() && i<nv2;i++){
		//スタックからポップする
		mgTLEdges::eitr endEdgeIter = edgeStack.end();
		const mgTLEdge edge = *(--endEdgeIter); edgeStack.pop_back();
		size_t alpha=edge.start(), beta=edge.end();
		if(is_boundary(alpha,beta)){
			if(used(alpha,beta)) continue;
		}

		//3点目の頂点を求める
		size_t status;
		size_t gamma = find3rdV(alpha, beta, status);
		//if(alpha==73 && beta==143 && gamma==71)
		//	cout<<edgeStack<<endl;//////////////
		//cout << "status = " << status << " 3-> "<< gamma << endl << endl;
		//求まったステータスに応じた処理を行う
		if(status==UNKNOWN)continue;						//UNKNOWN==0
		switch(status){
			case TWOTOUCH: continue;						//TWOTOUCH==1
			case NOTOUCH: edgeStack.push_front(edge); break;//NOTOUCH==2
			case REGULAR:									//REGULAR==3
				push1Vaft(alpha, beta, gamma);
				push1Vbefore(beta, gamma, alpha);
				push2V(gamma, alpha, beta);
				set_edge_used(alpha, beta);
				set_edge_used(beta, gamma);
				set_edge_used(gamma, alpha);
				break;
			case RIGHTTOUCH:								//RIGHTTOUCH==4
				push1Vbefore(beta, gamma, alpha);
				push1Vaft(gamma, alpha, beta);
				set_edge_used(alpha, beta);
				set_edge_used(beta, gamma);
				break;
			case LEFTTOUCH:									//LEFTTOUCH==5
				push1Vaft(alpha, beta, gamma);
				push1Vbefore(gamma, alpha, beta);
				set_edge_used(alpha, beta);
				set_edge_used(gamma, alpha);
				break;
			default: break;//求まらなかったとき				
		}
		if(!is_boundary(alpha,gamma) && (status!=RIGHTTOUCH)){
			//cout<<" Stacking("<<alpha<<","<<gamma<<")"<<endl;
			edgeStack.push_back(alpha, gamma);
		}
		if(!is_boundary(gamma,beta) && (status!=LEFTTOUCH)){
			//cout<<" Stacking("<<gamma<<","<<beta<<")"<<endl;
			edgeStack.push_back(gamma, beta);
		}
	}
}

//3点目の頂点(id of m_fans)を求める
size_t mgTLFans::find3rdV(
	size_t		alpha,	//エッジの始点(id of m_fans)
	size_t		beta,	//エッジの終点(id of m_fans)
	size_t&		status	//ステータス
){
	const size_t zcoord=2;
	const MGPosition& sPos=uv(alpha);
	const MGPosition& ePos=uv(beta);
	MGVector cedge(ePos-sPos);

	size_t nv=m_polygon.size();
	//多角形をループして3点目を求める
	bool gamma_is_used, right_touch, left_touch;
	double maxCang=2.0;
	size_t gamma=0;
	for(size_t i=0; i<nv; i++){
		if((alpha==i) || (beta==i))
			continue;//始終点と等しいときパスする

		//外積が負ではいけない
		const MGPosition &p3 = uv(i);//3点目の座標
		if((cedge*MGVector(p3-sPos))[zcoord] <= 0.0)
			continue;

		bool i_is_used, iright_touch=false, ileft_touch=false;
		double cang((sPos-p3).cangle(ePos-p3));
		if(cang < maxCang){//When the angle at i is larger than the before,
			//check if the edges (i,alpha) and (i,beta) do not have
			//any intersections with existing edges.
			i_is_used=used(i);
			if(i_is_used){
				iright_touch=used(i,alpha);
				if(!iright_touch){//共通エッジがあれば交点がない
					if(has_isect(i,alpha))
						continue;
				}
				ileft_touch=used(beta,i);
				if(!ileft_touch){//共通エッジがあれば交点がない
					if(has_isect(beta,i))
						continue;
				}
			}else{
				if(has_isect(i,alpha))
					continue;
				if(has_isect(beta,i))
					continue;
			}
			maxCang = cang;
			gamma_is_used=i_is_used;
			right_touch = iright_touch;
			left_touch = ileft_touch;
			gamma=i;
		}
	}
	if(maxCang > 1.5){//3点目が見つからなかった
		status=UNKNOWN;
		return gamma;
	}

	//ステータスの更新を行う
	if(!gamma_is_used){//未使用の点のときREGULAR
		status=REGULAR;
	}else if(right_touch){
		if(left_touch) status=TWOTOUCH;
		else status = RIGHTTOUCH;
	}else{
		if(left_touch) status=LEFTTOUCH;
		else status=NOTOUCH;
	}
	return gamma;
}

//線分(v1,v2)が既存のedgeと交点があるかどうかを調べる
//Function returns true if (v1,v2) had an isect.
bool mgTLFans::has_isect(
		size_t 	v1,	//線分1の点
		size_t 	v2	//線分1の点
)const{
	if(is_boundary(v1,v2))
 return false;

	double error=MGTolerance::rc_zero();
	const MGPosition& p1=uv(v1), &p2=uv(v2);
	const MGVector dir1(p2-p1);
	double udir1=dir1[0], vdir1=dir1[1];
	size_t nfan=size()-1;
	for(size_t i=0; i<nfan; i++){//Loop over fans(m_fans).
	if(v1==i || v2==i)
		continue;

	const MGPosition& p3=uv(i);//Center of the fani.
	const mgTLFan& fani=*(m_fans[i]);
	size_t nv=fani.size();
	MGVector dir31(p3-p1), dir32(p3-p2);
	double udir31=dir31[0], vdir31=dir31[1];
	double udir32=dir32[0], vdir32=dir32[1];
	double z13=vdir31*udir1-udir31*vdir1;
		//z value of vector product of dir31 and dir1
	for(size_t j=0; j<nv; j++){//Loop over vertices on the fani.
		size_t vj=fani[j];
		if(vj<=i)
			continue;
		if(vj==v1 || vj==v2)
			continue;
//		if(is_boundary(i,vj)) continue;
		if(!is_boundary(i,vj) && !used(i,vj))
			continue;

		const MGPosition& p4=uv(vj);
		MGVector dir41(p4-p1);
		double z14=dir41[1]*udir1-dir41[0]*vdir1;
			//z value of vector product of dir41 and dir1
		double z134=z13*z14;
		if(z134>=0.)
			continue;//This means p3 and p4 are located at the same side about the straight line (p1, p2).
//		if(-z134<=error) continue;///////////

		const MGVector dir2(p4 - p3);
		double udir2=dir2[0], vdir2=dir2[1];
		double z23=vdir31*udir2-udir31*vdir2;
		double z24=vdir32*udir2-udir32*vdir2;
		double z234=z23*z24;
		if(z234>=0.)
			continue;//This means p1 and p2 are located at the same side about the
							//straight line (p3, p4).
//		if(-z234<=error) continue;///////////

		return true;//In this case, (p1,p2) and (p3,p4) has and intersection.
	}

	}
	return false;
}

//Fan生成に必要な変数の準備を行う
//edgeStackにedgeのstackを積む
void mgTLFans::init(
	mgTLEdges& edgeStack
){
	//多角形の辺をスタックに積む、このとき頂点を未使用にしておく
	int i,ip1,nv=m_polygon.size();
	
	//Remove too close points in m_polygons.
	for(i=nv-1; i>=0; i--){
		ip1=i+1;
		if(ip1==m_polygon.size()) ip1=0;
		MGVector dif(uv(i)-uv(ip1));
		if(fabs(dif[0])<m_uerror && fabs(dif[1])<m_verror){
			if(ip1>i) m_polygon.erase(ip1);
			else m_polygon.erase(i);
		}
	}

	//push edges on the stack edgeStack.
	nv=m_polygon.size();
	m_fans.resize(nv);
	int im1=nv-1;
	for(i=0; i<nv; i++){
		ip1=i+1;
		if(ip1==nv) ip1=0;
		m_fans.reset(i,new mgTLFan(ip1,im1));
		edgeStack.push_back(i,ip1);
		im1=i;
	}
}

//Test if the edge(alpha, beta) is boundary or not.
bool mgTLFans::is_boundary(size_t alpha, size_t beta) const{
	if(alpha>beta){
		size_t temp=alpha; alpha=beta; beta=temp;
	}
	size_t n=size();
	if(alpha==0 && beta==n-1)
		return true;
	size_t alpha_n=(alpha+1)%n;
	return alpha_n==beta;
}

//目的：中心点(center)がalphaのmgTLFanに対し､頂点gammaを基準betaの
//後ろに追加する
void mgTLFans::push1Vaft(
	size_t	alpha,	//中心点のインデックス
	size_t	beta,	//基準となる頂点のインデックス
	size_t	gamma	//追加する頂点のインデックス
){
//	cout<<"===push1Vaft("<<alpha<<","<<beta<<","<<gamma<<")"<<endl;
	mgTLFan& fan=*(m_fans[alpha]);
//	cout<<"before:"<<alpha<<"|";fan.print_indices(cout,m_polygon);cout<<endl;
	mgTLFan::IndexItr j=fan.find_aft(beta);
	if(j!=fan.end()){//If not found
		j++;
		if(j!=fan.end()&&(*j)==gamma) return;
	}
	fan.insert(j,gamma);
	fan.set_vertex_used();
//	cout<<"aft:"<<alpha<<"|";fan.print_indices(cout,m_polygon);cout<<endl;
}

//目的：中心点(center)がalphaのmgTLFanに対し､頂点betaを基準gammaの
//前に追加する
void mgTLFans::push1Vbefore(
	size_t	alpha,	//中心点のインデックス
	size_t	beta,	//追加する頂点のインデックス
	size_t	gamma	//基準となる頂点のインデックス
){
//	cout<<"===push1Vbefore("<<alpha<<","<<beta<<","<<gamma<<")"<<endl;
	mgTLFan& fan=*(m_fans[alpha]);
//	cout<<"before:"<<alpha<<"|";fan.print_indices(cout,m_polygon);cout<<endl;
	mgTLFan::IndexItr j=fan.find(gamma);
	if(j!=fan.end()){//If not found
		if(j!=fan.begin()){
			mgTLFan::IndexItr jm1=j;jm1--;
			if((*jm1)==beta) return;
		}
	}else{
		j=fan.begin();
	}
	fan.insert(j,beta);
	fan.set_vertex_used();
//	cout<<"aft:"<<alpha<<"|";fan.print_indices(cout,m_polygon);cout<<endl;
}

//目的：中心点(center)がalphaのfanに頂点(beta,gamma)を新規に作成する
//push2V is invoked only for unused vertices.
void mgTLFans::push2V(
	size_t	alpha,	//γのインデックス
	size_t	beta,	//αのインデックス
	size_t	gamma	//βのインデックス
){
//	cout<<"===push2V("<<alpha<<","<<beta<<","<<gamma<<")"<<endl;
	mgTLFan& fan=*(m_fans[alpha]);
//	cout<<"before:"<<alpha<<"|";fan.print_indices(cout,m_polygon);cout<<endl;
	mgTLFan::IndexItr j=fan.begin();
	if((*j++)!=beta){
		fan.insert(j,beta);
	}
	mgTLFan::IndexItr je=fan.end(); j=je;j--;
	if(*j!=gamma) fan.insert(j,gamma);
	fan.set_vertex_used();
//	cout<<"aft:"<<alpha<<"| ";fan.print_indices(cout,m_polygon);cout<<endl;
}

//Set edge(i,j) as used.
void mgTLFans::set_edge_used(size_t alpha, size_t beta){
	if(alpha>beta){
		size_t temp=alpha; alpha=beta; beta=temp;
	}
	mgTLFan* fan=m_fans[alpha];
	fan->set_edge_used(beta);
}

//check if vertex(alpha) is used or not.
bool mgTLFans::used(size_t alpha) const{
	return m_fans[alpha]->vertex_is_used();
}

//check if edge(alpha, beta) is used or not.
bool mgTLFans::used(size_t alpha, size_t beta) const{
	if(alpha>beta){
		size_t temp=alpha; alpha=beta; beta=temp;
	}
	const mgTLFan* fan=m_fans[alpha];
	return fan->edge_is_used(beta);
}

const MGPosition& mgTLFans::uv(size_t i)const{
	size_t j=m_polygon[i];
	return m_points[j];
}

std::ostream& operator<< (std::ostream& out, const mgTLFans& fans){
	out<<"mgTLFans::num of fans="<<fans.size()<<std::endl;
	mgTLFans::const_iterator i=fans.begin(), ie=fans.end();
	for(size_t j=0; i!=ie; i++,j++){
		out<<j<<"("<<fans.m_polygon[j]<<")"<<"| ";
		(**i).print_indices(out,fans.m_polygon);
		out<<endl;
	}
	return out;
}
