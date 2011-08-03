/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/BPointSeq.h"
#include "mg/Position.h"
#include "mg/LBRep.h"
#include "mg/RLBRep.h"
#include "mg/Straight.h"
#include "mg/SPointSeq.h"
#include "mg/SBRep.h"
#include "mg/RSBRep.h"
#include "mg/Tolerance.h"

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

using namespace std;

// Implementation of knot remove.

//ノット削除関数(B表現曲線のみ)
//トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
//removal knot. line_zero tolerance is used.
void MGCurve::remove_knot(){;}

#define START_DELNUM 4000
//Remove knot if removed line has the difference less than line_zero();
//The difference is checked only for the space id of coef(.,j+k)
//for k=0, ..., nd-1.
void MGLBRep::remove_knot(size_t j, size_t nd){
	int k=order();
	int km1 = k-1;	//次元数
	int i, n=bdim();
	int mid=(km1+n)/2;
	if(mid<=km1) mid=k;
	else if(mid>START_DELNUM){mid=START_DELNUM;}
	double totalTol = 0.0;
	const double line0=MGTolerance::line_zero();

	size_t num_remained;
	for(i=n-1; i>=mid;){
		size_t ndel = remove_knot_one(line0,i,totalTol,num_remained,j,nd);
		i-=int(ndel+num_remained);
	}
//	cout<<(*this)<<endl;
	totalTol = 0.0;
	i=k;
//	int nt=int(bdim())-mid;
	int nt=int(bdim())-mid-1;if(nt<0) nt=0;
	while(int(bdim())-i>nt){
		remove_knot_one(line0,i,totalTol,num_remained,j,nd);
		i+=num_remained;
	}
//	cout<<(*this)<<endl;
}

//ノット削除関数(1つのノット)
//戻り値は、削除したノットの数
//When snum!=0, tolerance of totalTol is checked only for coef(.,sid+j),
//of j=0, ..., snum-1. When snum=0, snum is set as sdim();
int MGLBRep::remove_knot_one(
	double line0,		//Tolerance allowed for the knot removal.
	size_t	id,			//削除しようとするノットの番号
	double& totalTol,	//誤差合計
	size_t& num_knot,	//Remained knot number at knot(id) after removed.
	size_t sid,			//Space dimension start id of this LBRep's B-coef.
	size_t snum			//Num of space dimension for the totalTol tolerance check.
){
	size_t sd=sdim();
	if(!snum)
		snum=sd;
	MGPosition P(snum), Q(snum), R(snum);
	const MGBPointSeq& bcoef=line_bcoef();

	const int k = order();	//オーダー
	const int km1 = k - 1;	//次数
	const int n=bdim();	//B表現次元
	const double tau=knot(id);

	//Get the multiplicity of knot(id).
	//We need the end knot id of the multiplicity in 'id'.
	int m, idold=id;
	int multi=1;
	for(m=id+1; m<n && tau==knot(m); m++){
		multi++; id++;
	}
	for(m=idold-1; m>=k && tau==knot(m); m--)
		multi++;
	num_knot=multi;

	//const double line0=MGTolerance::line_zero();
	MGBPointSeq temp_bp(2 * km1 + 1, sd);

	int	last = id - multi;
	int first = id - km1;
	int ndel;				//削除したノットの数
	if(k == 2){
	//オーダー2のときの処理
		const int off = first - 1;
		bcoef.point(off,sid,snum,P);
		bcoef.point(last+1,sid,snum,Q);
		MGStraight st1(P, Q);
		bcoef.point(first,sid,snum,R);
		double tmpTol = st1.distance(R);
		totalTol += tmpTol;
		if(line0>0. && totalTol>line0)
			ndel=0;
			//誤差合計がトレランス以上の場合ノット削除を行わない
		else ndel=1;
	}else{

	//オーダー3以上のときの処理
	for(ndel = 0; ndel<multi; ndel++){
		const int off = first - 1;
		temp_bp.store_at(0, coef(off));
		temp_bp.store_at(last + 1 - off, coef(last + 1));
		int i = first;	int j = last;
		int ii = 1;		int jj = last - off;
		while((j-i) > ndel){
			double ti=knot(i), tjmndel=knot(j - ndel);
			double alfi = (tau-ti) / (knot(i+k+ndel) - ti);
			double alfj = (tau - tjmndel) / (knot(j+k) - tjmndel);
			temp_bp.store_at(ii, (coef(i) - (1.0 - alfi) * temp_bp(ii-1)) / alfi);
			temp_bp.store_at(jj, (coef(j) - alfj * temp_bp(jj+1)) / (1.0 - alfj));
			i++;	j--;	ii++;	jj--;
		}

		if((j-i) < ndel){
			temp_bp.point(ii - 1,sid,snum,P);
			temp_bp.point(jj + 1,sid,snum,Q);
			totalTol += (P - Q).len();
		}else{
			double ti=knot(i);
			double alfi = (tau - ti) / (knot(i+k+ndel) - ti);
			temp_bp.point(ii+ndel+1,sid,snum,P);
			temp_bp.point(ii-1,sid,snum,Q);
			bcoef.point(i,sid,snum,R);
			totalTol += (R - (alfi*P + (1.0 - alfi)*Q)).len();
		}
		//誤差合計がトレランス以上の場合ノット削除を行わない
		if(line0>0. && totalTol > line0) break;

		i = first;	j = last;
		while((j-i) > ndel){
			m_line_bcoef.store_at(i, temp_bp(i - off));
			m_line_bcoef.store_at(j, temp_bp(j - off));
			i++;	j--;
		}
		first--; last++;
	}

	}

	if(!ndel){totalTol=0.; return 0;}
	int l = id + 1;
	const int num_end_knot = n + km1;	//最終ノットの番号
	for(; l <= num_end_knot; l++) knot(l-ndel) = knot(l);
	int j = (2*id - multi - km1)/2;//最初のコントロールポイント id
	int i = j;
	for(l = 1; l < ndel; l++){
		if((l % 2) == 1)i++; else j--;
	}
	for(l = i+1; l<n; l++, j++) m_line_bcoef.store_at(j, coef(l));

	//曲線を更新する
	size_t newnbd=n - ndel;
	m_knot_vector.set_bdim(newnbd);
	m_line_bcoef.set_length(newnbd);
	num_knot=multi - ndel;
	if(num_knot) totalTol = 0.0;//ノット削除が行われないとき誤差合計をクリアする
	return ndel;
}

//ノット削除関数
void MGRLBRep::remove_knot(){
	double Pmax = 0.0, Wmin = 1.0;
	size_t sd=sdim();
	for(size_t i = 0;i < bdim(); i++){
		MGPosition pos = eval_position(knot(i));
		double P = pos.len();
		double W = coef(i, sd);		//最終次元が重みである
		if(P > Pmax) Pmax = P;
		if(W < Wmin) Wmin = W;
	}
	double save=MGTolerance::line_zero();
	double tol =  save* Wmin /(1 + Pmax);
	MGTolerance::set_line_zero(tol);
	m_line.remove_knot();
	MGTolerance::set_line_zero(save);
}

//ノット削除関数
//トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
void MGSBRep::remove_knot(){
	remove_knot_u();
	exchange_uv();
	remove_knot_u();
	exchange_uv();
}

//uノットを削除する
int MGSBRep::remove_knot_u(){
	int k = order_u(), n = bdim_u();//order and bdim along u.
	int km1=k-1;
	int mid=(km1+n)/2;
	if(mid<=km1)
		mid=k;
	size_t num_remained;

	int nDel = 0;		//削除したノットの数
	double totalTol = 0.0;
	double line0=MGTolerance::line_zero();
	int i;
	for(i=n-1; i>=mid;){
		size_t nDelOne = remove_knot_u_one(line0,i, totalTol, num_remained);
		i-=(nDelOne+num_remained);
		nDel+=nDelOne;
	}
	i=k;
	int nt=bdim_u()-mid-1;
	if(nt<0) nt=0;

	totalTol=0.;
	while(int(bdim_u())-i>nt){
		size_t nDelOne=remove_knot_u_one(line0,i,totalTol,num_remained);
		i+=num_remained;
		nDel+=nDelOne;
	}
	return nDel;
}

//uノットを削除する
//関数の戻り値は削除したノットの数
int MGSBRep::remove_knot_u_one(
	double line0,		//Tolerance allowed for the knot removal. 
						//When line0<=0., removal will be done uncoditionally.
	size_t	id,			//削除しようとするノットの番号
	double& totalTol,	//最初は０を入力する。あとはremove_knot_u_oneが管理する。
	//削除を繰り返しているときは誤差が加算される。これ以上削除できない時０クリアされる。
	size_t& num_knot	//Remained knot number at knot(id) after removed.
){
	const int k = order_u();	//オーダー
	const size_t km1 = k - 1;		//次数
	const int m = bdim_u(), n = bdim_v();
	const size_t sd=sdim();
	const double tau=knot_u(id);

	//Get the multiplicity of knot(id).
	//We need the end knot id of the multiplicity in 'id'.
	int a, multi=1, idold=id;
	for(a=id+1; a<n && tau==knot_u(a); a++){
		multi++; id++;
	}
	for(a=idold-1; a>=k && tau==knot_u(a); a--)
		multi++;
	num_knot=multi;

	int	last = id - multi;
	int first = id - km1;

	int l;
	int ndel=0;			//削除したノットの数
	if(k==2){

//オーダー2のときの処理
	const int off = first - 1;
	double maxTol = 0.0;
	int l;
	for(l=0; l<n; l++){
		MGStraight st1(coef(off, l), coef(last + 1, l));
		double tmpTol = st1.distance(coef(first, l));
		if(line0>0. && tmpTol>line0){ndel= 0; break;}//ノットが削除できなかった
		if(tmpTol>maxTol) maxTol = tmpTol;
	}
	if(l>=n){//When possible to delete knot.
		ndel=1;totalTol += maxTol;
	}

	}else{
//オーダー3以上のときの処理
	for(ndel = 0; ndel < multi; ndel++){

	MGSPointSeq temp(2*km1+1,n,sd);
	const int off = first - 1;
	int i, j, ii, jj;
	double maxTol = 0.0;
	for(l=0; l<n; l++){
		i = first;	j = last;
		ii = 1;		jj = last-off;
		temp.store_at(0, l, coef(off,l));
		temp.store_at(last+1-off,l,coef(last+1,l));
		double alfi, alfj, tmpTol = 0.0;
		while((j-i) > ndel){
			double ti=knot_u(i), tjmndel=knot_u(j-ndel);
			alfi = (tau-ti) / (knot_u(i+k+ndel) - ti);
			alfj = (tau-tjmndel) / (knot_u(j+k) - tjmndel);
			temp.store_at(ii,l,(coef(i,l) - (1.0 - alfi)*temp(ii-1,l)) / alfi);
			temp.store_at(jj, l, (coef(j,l) - alfj * temp(jj+1,l)) / (1.-alfj));
			i++;	j--;	ii++;	jj--;
		}
		if((j-i) < ndel){
			tmpTol = (temp(ii-1,l) - temp(jj+1,l)).len();
			if(line0>0. && tmpTol>line0){ maxTol=tmpTol; break;}
		}else{
			double ti=knot_u(i);
			alfi = (tau-ti) / (knot_u(i+k+ndel)-ti);
			tmpTol=(coef(i,l)-(alfi * temp(ii+ndel+1,l)+(1.-alfi)*temp(ii-1,l))).len();
			if(line0>0. && tmpTol>line0){ maxTol=tmpTol; break;}
		}
		if(tmpTol > maxTol) maxTol=tmpTol;
	}
	totalTol += maxTol;
	if(line0>0. && totalTol>line0) break;
	i = first;	j = last;
	while((j-i) > ndel){
		for(int vk=0; vk<n; vk++){
			m_surface_bcoef.store_at(i,vk,temp(i-off,vk));
			m_surface_bcoef.store_at(j,vk,temp(j-off,vk));
		}
		i++; j--;
	}
	first--;	last++;

	}
	}

	if(!ndel){totalTol=0.; return 0;}
	const int num_end_uknot = m + km1;		//最終uノットの番号
	for(l=id+1; l<=num_end_uknot; l++) knot_u(l-ndel)=knot_u(l);
	int j=(2*id-multi-km1)/2;	//最初のコントロールポイント
	int i=j;
	for(l=1; l<ndel; l++){
		if((l % 2) == 1) i++; else j--;
	}
	for(l=i+1; l<m; l++){
		for(int vk = 0; vk < n; vk++) m_surface_bcoef.store_at(j,vk,coef(l,vk));
		j++;
	}
	//曲面を更新する
	size_t newnbd=m-ndel;
	m_uknot.set_bdim(newnbd);
	m_surface_bcoef.set_length(newnbd, n);
	num_knot=multi - ndel;
	if(num_knot) totalTol = 0.0;//ノット削除が行われないとき誤差合計をクリアする
	return ndel;
}

//uノットを一つ削除する
//関数の戻り値は削除したノットの数
int MGSBRep::remove_knot_one(
	double line0,		//Tolerance allowed for the knot removal. 
						//When line0<=0., removal will be done uncoditionally.
	size_t	id,	//削除しようとするノットの番号
	double& tol,//削除後の誤差が出力される
	bool u_knot	//削除対象が（u,v)のいずれのknot vectorかを入力する
				//=trueのとき、u-knot_vectorを削除
){
	size_t num_knot;
	int delnum;
	if(u_knot) delnum=remove_knot_u_one(line0,id,tol,num_knot);
	else{
		exchange_uv();
		delnum=remove_knot_u_one(line0,id,tol,num_knot);
		exchange_uv();
	}
	return delnum;
}

//ノット削除関数
//トレランスはline_zeroを使用する。元のノットが細かいものほど削除しやすい
void MGRSBRep::remove_knot(){
	size_t i = 0, j = 0;
	double Pmax = 0.0, Wmin = 1.0;
	for(;i < bdim_u(); i++){
		for(j = 0; j < bdim_v(); j++){
			MGPosition pos = eval(knot_u(i), knot_v(j), 0);
			double P = pos.len();
			double W = coef(i, j, sdim());		//最終次元が重みである
			if(P > Pmax)Pmax = P;
			if(W < Wmin)Wmin = W;
		}
	}
	double save=MGTolerance::line_zero();
	double tol =  save* Wmin /(1 + Pmax);
	MGTolerance::set_line_zero(tol);
	m_surface.remove_knot();
	MGTolerance::set_line_zero(save);
}

//複数カーブの共通で削除できるノットを削除する。
//ただし、入力カーブは同じノットベクトルを持つものとする。
void remove_knot_curves(
	MGPvector<MGLBRep>& brepList,		//曲線列
	MGLBRep**		tp,	//接続面	input and output.
		//if tp[i] for crvl[i] was not null, converted new tp will be output.
	double tp_length
){
	MGLBRep& lb0=*(brepList[0]);
	int nCrv = brepList.size();	//曲線数
	int k = lb0.order(),n = lb0.bdim();//次元数&B表現次元
	int km1=k-1;
	int mid=(km1+n)/2;
	if(mid<=km1) mid=k;
	const double line0=MGTolerance::line_zero();
	const double rc0=tp_length*MGTolerance::angle_zero();
	size_t num_remained;
	int preDel, nDel, nDel2;

	std::vector<double> totalTol(nCrv,0.);
	std::vector<double> totalTolTP(nCrv,0.);
	std::vector<MGLBRep> lb(nCrv), tpSave(nCrv);
	std::vector<bool> tp_specified(nCrv);

	int i;
	for(i=n-1; i >= mid;){
		int j=0;
		for(; j<nCrv; j++){
			tp_specified[j]=false;if(tp){ if(tp[j]) tp_specified[j]=true;};
			MGLBRep& lbj=*(brepList[j]);lb[j]=lbj; //曲線は保存しておく
			nDel2=nDel = lbj.remove_knot_one(line0,i, totalTol[j],num_remained);
			if(tp_specified[j]){
				tpSave[j]=*tp[j];
				nDel2 = tp[j]->remove_knot_one(rc0,i, totalTolTP[j],num_remained);
			}
			if(j == 0) preDel = nDel;
			if(!nDel || nDel!=nDel2 || preDel!=nDel){//ノット削除に失敗したときの処理
				std::fill(totalTol.begin(), totalTol.end(), 0.0);	//誤差合計をクリアする
				std::fill(totalTolTP.begin(), totalTolTP.end(), 0.0);	//誤差合計をクリアする
				int j2=j;
				while(j2>=0){
					MGLBRep& lbj2=*(brepList[j2]); lbj2=lb[j2];	//曲線を元に戻す
					if(tp_specified[j2]) *(tp[j2])=tpSave[j2];
					j2--;
				}
				break;
			}
		}
		if(j>=nCrv)
			i-=(nDel+num_remained);
		else i--;
	}

	std::fill(totalTol.begin(), totalTol.end(), 0.0);	//誤差合計をクリアする
	std::fill(totalTolTP.begin(), totalTolTP.end(), 0.0);//誤差合計をクリアする
	i=k;
	size_t nt=lb0.bdim()-mid;
	while(lb0.bdim()-i>nt){
		int j=0;
		for(; j<nCrv; j++){
			tp_specified[j]=false;if(tp){ if(tp[j]) tp_specified[j]=true;};
			MGLBRep& lbj=*(brepList[j]);lb[j]=lbj;//曲線は保存しておく
			nDel2=nDel = lbj.remove_knot_one(line0,i, totalTol[j],num_remained);
			if(tp_specified[j]){
				tpSave[j]=*tp[j];
				nDel2 = tp[j]->remove_knot_one(rc0,i, totalTolTP[j],num_remained);
			}
			if(j == 0) preDel = nDel;
			if(!nDel || nDel!=nDel2 || preDel!=nDel){//ノット削除に失敗したときの処理
				std::fill(totalTol.begin(), totalTol.end(), 0.0);	//誤差合計をクリアする
				std::fill(totalTolTP.begin(), totalTolTP.end(), 0.0);	//誤差合計をクリアする
				int j2=j;
				while(j2>=0){
					MGLBRep& lbj2=*(brepList[j2]); lbj2=lb[j2];	//曲線を元に戻す
					if(tp_specified[j2]) *(tp[j2])=tpSave[j2];
					j2--;
				}
				break;
			}
		}
		if(j>=nCrv)
			i+=num_remained;
		else i++;
	}

}
