/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/Box.h"
#include "mg/Position.h"
#include "mg/Curve.h"
#include "mg/CompositeCurve.h"
#include "mg/CurveParameter.h"
#include "mg/CCisect_list.h"
#include "mg/CSisect_list.h"
#include "mg/FSurface.h"
#include "mg/Tolerance.h"

//#if defined(_DEBUG)
//#define new DEBUG_NEW
//#undef THIS_FILE
//static char THIS_FILE[] = __FILE__;
//#endif

#define SIN8DEG 0.022 //SIN8DEG=squaer of sin(8 degree).
//Private class for MGCurve::common(),
//Stores the information if the point of a curve's parameter value para is
//on the other curve or not. When the point is start or end point,
//the point is registered as off the curve if the two curves directions are
//opposite and one curve is departing from the other curve.
class MGCommonON{
public:
	MGCommonON(const MGCurve& curve1,const MGCurve& curve2,
		bool is_curve1Para,	//indicate if para is curve1's parameter or curve2.
		double para,
		int para_is_SE=0
			//=-1, or 1 if para is the start(-1) or end(1) parameter of curve1 or 2.
			//=0 means para is a mid point.
	);

	MGCommonON(int ondata=0, double para1=0., double para2=0.)
		:m_on(ondata), m_param1(para1), m_param2(para2){;};

	//Test if the param1 of this is in between the param1 of prev and next.
	bool in_between(const MGCommonON& prev, const MGCommonON& next)const{
		return prev.param1()<m_param1 && m_param1<next.param1();
	};

	//Get next data of this, given curve1's next parameter t1.
	MGCommonON get_next(double t1, const MGCurve& curve1, const MGCurve& curve2)const;
	void set_on(int ondata=1){m_on=ondata;};
	int on()const{return m_on;};
	bool on_and_in_between(const MGCommonON& prev, const MGCommonON& next)const{
		return m_on && in_between(prev,next);};
	double param1()const{return m_param1;};
	double param2()const{return m_param2;};

private:
	int m_on;//=0: two points are different points.
			//=1: two points are the same points.
	double m_param1, m_param2;//parameter values of two curves when m_on is true.
};

MGCommonON::MGCommonON(
	const MGCurve& curve1,const MGCurve& curve2,
	bool is_curve1Para,	//indicate if para is curve1's parameter or curve2.
	double para,
	int para_is_SE	//=-1, or 1 if para is the start(-1) or end(1) parameter of curve1 or 2.
					//=0 means para is a mid point.
):m_on(0){
	const MGCurve* para_curve; const MGCurve* other_curve;
	double* para_param; double* other_param;
	if(is_curve1Para){
		para_curve=&curve1; other_curve=&curve2;
		para_param=&m_param1; other_param=&m_param2;
	}else{
		para_curve=&curve2; other_curve=&curve1;
		para_param=&m_param2; other_param=&m_param1;
	}
	*para_param=para;
	double& t=*other_param;
	if(other_curve->on(para_curve->eval(para),t)){
		double perror=other_curve->param_error();
		int pra_is_SE_other=0;
		double ts=other_curve->param_s(), te=other_curve->param_e();
		if((t-ts)<=perror){
			pra_is_SE_other=-1;
			t=ts;
		}else if((te-t)<=perror){
			pra_is_SE_other=1;
			t=te;
		}

		MGVector v1=para_curve->eval(para,1), v2=other_curve->eval(t,1);
		MGVector v12=v1*v2;
		if((v12%v12) <= (v1%v1)*(v2%v2)*SIN8DEG){
			//v1 and v2 are parallel.
			if(para_is_SE){
				if(para_is_SE*pra_is_SE_other>0 && v1%v2<0.)
					return;
				if(para_is_SE*pra_is_SE_other<0 && v1%v2>0.)
					return;
			}
			m_on=1;
		}
	}
}

MGCommonON MGCommonON::get_next(
	double t1,//next parameter of curve1 to get next MGCommonON.
	const MGCurve& curve1, const MGCurve& curve2
)const{
	MGPosition pos=curve1.eval(t1);
	//のっているかどうか調べるのにOnではなくperp_guessを使用している
	double t2=param2();
	int found = curve2.perp_guess(1.,0.,pos,t2,t2);
	MGCommonON next(0,t1,t2);
	if(found){	
		if(pos == curve2.eval(t2))
			next.set_on();//t2 is not on this curve.
	}
	return next;
}

#define m_nMaxCommonEdge 4	//共通エッジは4個までしか求めない
//関数名： common
//目的：与えられた曲線と自身の共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			curve2,		(I/ )	与えられる曲線
//		std::vector<double>&	vecComSpan	( /O)	共通部分のパラメータ範囲
//		 4nの配列で、vecComSpan(4*i+0),vecComSpan(4*i+1)が自身のパラメータ範囲
//					(vecComSpan(4*i+0) < vecComSpan(4*i+1))、
//				 vecComSpan(4*i+2),vecComSpan(4*i+3)がcurve2のパラメータ範囲
//戻り値：
//		共通部分の数:	共通部分が求まった
//		0:				共通部分がなかった
//		-1:				共通エッジの収束計算エラー
//		-2:				共通エッジが４個以上求まった(のっていないと見なす)
//追記：
//	曲線が共通かどうかの誤差にはline_zero()を、パラメータ範囲の収束計算の誤差には、
//  パラメータ範囲*rc_zero()を使用した
int MGCurve::common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan
)const{
	vecComSpan.clear();

	const MGCompositeCurve* compo=dynamic_cast<const MGCompositeCurve*>(&curve2);
	size_t nv2=0;
	std::vector<double> vec2;
	if(compo){
		compo->common(*this,vec2);
		nv2=vec2.size();
		if(!nv2)
			return 0;
	}else{
		const MGTrimmedCurve* trim=dynamic_cast<const MGTrimmedCurve*>(&curve2);
		if(trim){
			trim->common(*this,vec2);
			nv2=vec2.size();
			if(!nv2)
				return 0;
		}
	}
	if(nv2){
		//When curve2 was a MGCompositeCurve or MGTrimmedCurve and common ranges are
		//obtained, replace the output order and return.
		size_t nv24=nv2/4;
		size_t idnv2=0;
		for(size_t j=0; j<nv24; j++){
			size_t j4=j*4;
			vecComSpan.push_back(vec2[j4+2]);
			vecComSpan.push_back(vec2[j4+3]);
			vecComSpan.push_back(vec2[j4]);
			vecComSpan.push_back(vec2[j4+1]);
		}
		return nv24;
	}

	//２曲線のボックス枠に重複があるかどうか調べる
	if(!has_common(curve2))
		return 0;

	//std::cout<<(*this)<<std::endl<<curve2<<std::endl;////////////
	//曲線同士が乗っているかを調べるので、誤差をラインゼロにする
	double dWcTol = MGTolerance::set_wc_zero(MGTolerance::line_zero());
	MGCommonON SEon[2]={
		MGCommonON(*this,curve2,false,curve2.param_s(),-1),
		MGCommonON(*this,curve2,false,curve2.param_e(),1)
	};
				//Data if curve2's start or end point is on this curve.

	int ndiv = intersect_dnum();
	double wholeSpan=param_span();
	double span=wholeSpan/ndiv;		//チェックポイントのパラメータスパン
	double error=wholeSpan*MGTolerance::rc_zero()*2.;

	int rc, num = 0;
	MGCommonON sparam(*this,curve2,true,param_s(),-1);
	//共通部分のパラメータを求める
try{
	do{
		double oneSpan[4];	//共通パラメータsparam1, eparam1, sparam2, eparam2の順で入っている
		rc = common_one_span(curve2,SEon,sparam,span,oneSpan);
		if(rc){
			double plen=oneSpan[1]-oneSpan[0];
			if(plen<=error)
				continue;//This is a intersection point, is discarded.
			for(size_t i=0; i<4; i++)
				vecComSpan.push_back(oneSpan[i]);
			num++;
			if(num>m_nMaxCommonEdge){
				return -2;	//最大繰り返し数を越えたらエラー終了
			}
		}
	}while(rc<0);	//curve1の途中であれば繰り返す
	MGTolerance::set_wc_zero(dWcTol);	//誤差を元に戻す
	return vecComSpan.size()/4;
}
catch(int err){
	//err	 0:ボックスが重なってなかった
	//		-1:繰り返し計算エラー
	//		-2:共通エッジ数が４を越えた
	MGTolerance::set_wc_zero(dWcTol);	//誤差を元に戻す
	return err;
}
}

//thisのスタートパラメータを与え共通パラメータ範囲を1つだけ取得する。
//sparamはスタートパラメータであり、次回検索に使用するパラメータを返す
//戻り値は	1:範囲が求まった(thisの最後まで)
//			0:範囲が求まらなかった(thisの最後まで)
//			-1:範囲が求まった(thisの途中まで、まだ共通範囲があるかもしれない)
int MGCurve::common_one_span(
	const MGCurve& curve2,
	MGCommonON SEon[2],//data if curve2's start or end point is on this curve.
	MGCommonON& sparam, //Input start parameter of this curve to search starting point of the next
					//common part.
					//On return from common_one_span, the next startign parameter will be set.
	double span,	//parameter span to advance the next on point check of this curve.
	double comSpan[4]//共通パラメータsparam1, eparam1, sparam2, eparam2を返却する
)const{
	const MGCurve& curve1=*this;
	double t1e=param_e();
	double t1em=t1e-span/2.;
	int endID=-1;//When endID=0 or 1, it indicates that SEon[endID] is the end point of the common part.
	MGCommonON ONprev(sparam), ONnext;
	double t1=sparam.param1();
	if(sparam.on()){
		int endID=-1;
			//When endID=0 or 1, it indicates that SEon[endID] is the end point of the common part.
		comSpan[0]=t1; comSpan[2]=sparam.param2();
		for(t1+=span; t1<t1em; t1+=span){
			ONnext=ONprev.get_next(t1,curve1,curve2);
			if(!ONnext.on()){
				ONnext=MGCommonON(*this,curve2,true,t1);
				if(!ONnext.on())
					break;
			}
			if(SEon[0].on_and_in_between(ONprev,ONnext)){
				endID=0;
				break;
			}else if(SEon[1].on_and_in_between(ONprev,ONnext)){
				endID=1;
				break;
			}
			ONprev=ONnext;
		}
		//Here, three cases:
		//case1=ONnext.on() is true and endID=0 or 1:common part passed through the start or end
		//                                           point of curve2.
		//case2=this curve was common(to t1>=t1em) to curve2.
		//case3=ONnext.on() is false: a boundary point is in between ONprev & ONnext.

		//case 1:parameter of curve1 pass through start or end point of curve2.
		if(endID>=0){
			if(endID==0){
				if(SEon[1].on_and_in_between(ONprev,ONnext)){
					double t2SE0=SEon[0].param2(), t2SE1=SEon[1].param2(), t2prev=ONprev.param2();
					if(fabs(t2SE1-t2prev)<fabs(t2SE0-t2prev))
						endID=1;
					sparam=SEon[(endID+1)%2];
				}else{
					sparam=ONnext;
				}
			}else{
				sparam=ONnext;
			}
			comSpan[1]=SEon[endID].param1(); comSpan[3]=SEon[endID].param2();
			SEon[endID].set_on(0);//This is to indicates SEon[endID] is already comsumed.
			return -1;
		}

		//case 2:parameter of curve1 reached to the end.
		if(t1>=t1em){
			ONnext=ONprev.get_next(t1e,curve1,curve2);
			if(ONnext.on()){
				comSpan[1]=t1e; comSpan[3]=ONnext.param2();
				return 1;			
			}
		}

		//case 3:curve1 took off from the common part. Get precise parameter.
		common_boundary(curve2,ONprev,ONnext,comSpan[1],comSpan[3]);
		sparam=ONnext;

		if(t1>=t1em)
			return 1;
		else
			return -1;
	}else{
		for(t1+=span; t1<t1em; t1+=span){
			ONnext=MGCommonON(curve1,curve2,true,t1);
			if(ONnext.on()){
				if(SEon[0].on_and_in_between(ONprev,ONnext)){
					endID=0;
				}else if(SEon[1].on_and_in_between(ONprev,ONnext)){
					endID=1;
				}
				break;
			}
			ONprev=ONnext;
		}
		//Here, ONprev.on() is false and three cases:
		//case1=ONnext.on() is true and endID=0 or 1:boundary point of the common is
		//      the start or end point of curve2.
		//case2=ONnext.on() is false: this curve was not common(to t1>=t1em) to curve2
		//case3=ONnext.on() is true and endID=-1, a boundary point of the common is in between ONprev and ONnext.
		
		//case 1:parameter of curve1 pass through start or end point of curve2.
		if(endID>=0){
			if(endID==0){
				if(SEon[1].on_and_in_between(ONprev,ONnext)){
					double t2SE0=SEon[0].param2(), t2SE1=SEon[1].param2(), t2prev=ONprev.param2();
					if(fabs(t2SE1-t2prev)<fabs(t2SE0-t2prev))
						endID=1;
				}
			}
			sparam=SEon[endID];
			SEon[endID].set_on(0);//This is to indicates SEon[endID] is already comsumed.
			return common_one_span(curve2,SEon,sparam,span,comSpan);
		}

		//case 2:parameter of curve1 reached to the end(here, ONprev==ONnext, and ONnext.on()==false).
		if(t1>=t1em){
			ONnext=MGCommonON(curve1,curve2,true,t1e,1);
			if(ONnext.on()){
				common_boundary(curve2,ONnext,ONprev,comSpan[0],comSpan[2]);
				comSpan[1]=t1e; comSpan[3]=ONnext.param2();
				return 1;			
			}else
				return 0;
		}

		//case 3:curve1 reached to a common part. Get the precise parameter and get the span.
		double t1on,t2on;
		common_boundary(curve2,ONnext,ONprev,t1on,t2on);
		sparam=MGCommonON(1,t1on,t2on);
		return common_one_span(curve2,SEon,sparam,span,comSpan);
	}
}

#define MAX_LOOP_COUNT 100
//this上のcurve2に乗っているパラメータonparamと乗っていないパラメータoffparam
//を与え、その境界点点を求める。
//チェックポイントの移動量 < (パラメータ範囲*rc_zero())になれば終了。
//戻り値は　収束して求まったパラメータ値である。
void MGCurve::common_boundary(
	const MGCurve& curve2,
	const MGCommonON& onparam,
	const MGCommonON& offparam,
	double& param1, double& param2//this curve's and curve2's parameter of the boundary.
)const{
	//移動量が誤差範囲になるまで繰り返し計算
	double tolerance=param_span()*MGTolerance::rc_zero();
	const MGCurve& curve1=*this;
	MGCommonON onp=onparam, offp=offparam;
	for(int ncount=0; ncount<MAX_LOOP_COUNT; ncount++){
		double t1=onp.param1();
		double t1next=(t1+offp.param1())*.5;
		if(fabs(t1next-t1)<=tolerance){
			param1=t1; param2=onp.param2();
			return;
		}

		MGCommonON next=onp.get_next(t1next,curve1,curve2);
		if(next.on())
			onp=next;
		else
			offp=next;
	}
	throw -1;	//収束しないときエラーで返す
}

//関数名：common
//目的：与えられた曲線と自身の交点もしくは共通部分があるかどうか調べる。
//引数：
//		const MGCurve&			curve2,		(I/ )	与えられる曲線
//		std::vector<double>&	vec_param	( /O)	共通部分のパラメータ範囲
//		MGCCisect_list&			isect		( /O)	交点
//				 4nの配列で、t(4*i+0),t(4*i+1)が自身のパラメータ範囲(t(4*i+0) < t(4*i+1))、
//							 t(4*i+2),t(4*i+3)が与曲線のパラメータ範囲(f(t(4*i+0))=f(t(4*i+2))
//戻り値：
//		3:交点も共通部分も求まった
//		2:交点のみが求まった
//		1:共通部分のみが求まった
//		0:交点も共通部分もなかった
//		-1:共通エッジの収束計算エラー
//		-2:共通エッジが４個以上求まった(のっていないと見なす)
//追記：
//	曲線が共通かどうかの誤差にはline_zero()、をパラメータ範囲の収束計算の
//	誤差には、パラメータ範囲*rc_zero()を使用した
int MGCurve::common(
	const MGCurve& curve2,
	std::vector<double>& vecComSpan,
	MGCCisect_list& isect
)const{
	isect.clear();

	//２曲線のボックス枠に重複があるかどうか調べる
	if(!has_common(curve2))
		return 0;

	//std::cout<<(*this)<<curve2<<std::endl;
	//共通部分のパラメータを求める
	int rc = common(curve2, vecComSpan);
	if(rc < 0)
		return rc;
	//std::cout<<"vecComSpan.size()="<<vecComSpan.size()<<std::endl;

	//交点を求める
	MGCCisect_list tmp_isect = this->isect(curve2);
	//std::cout<<"tmp_isect.size()="<<tmp_isect.size()<<std::endl;

	//共通部分に重なっている交点を削除する(this基準)
	MGCCisect_list::CCiterator i, isave;
	for(i=tmp_isect.begin();i != tmp_isect.end(); i++){
		double param1 = i->param1();
		std::vector<double>::iterator iter;
		for(iter = vecComSpan.begin(); iter != vecComSpan.end(); iter += 4){
			MGInterval itvl1(*iter, *(iter+1));
			if(itvl1 >> MGEReal(param1)){
				isave = i; isave++;
				tmp_isect.removeAt(i);
				i = --isave;
				break;
			}
		}
	}

	if(!tmp_isect.isEmpty()){
		isect = tmp_isect;
		if(rc > 0)
			return 3;	//交点と共通部分両方求まった
		else
			return	2;	//交点のみが求まった
	}else if(rc > 0)
		return 1;	//共通部分のみが求まった
	return 0;	//どっちも求まらなかった
}

//Compute the intersections of two objects.
MGisects MGCurve::intersection(const MGObject& obj2)const{
	MGisects isects=obj2.intersection(*this);
	isects.exchange12();
	return isects;
}

MGisects MGCurve::intersection(const MGCurve& obj2)const{
	MGCCisect_list isects2=isect(obj2);
	return MGisects(isects2);
}

MGisects MGCurve::intersection(const MGFSurface& obj2)const{
	MGCSisect_list isects2=obj2.isect(*this);
	return MGisects(isects2);
}
MGisects MGCurve::intersection(const MGSurface& obj2)const{
	MGCSisect_list isects2=isect(obj2);
	return MGisects(isects2);
}

//Get the paramer value of f(t)=0 that is defined by operator()(double)=0.
//Function's return code is:
//0: when the solution is successfully obtaine in t,
//1: There is no solution.
//-2:system error(usually this does not occur. If occured, some bugs are included.)
int MGCurveParameter::getCurveParameter(
	double& t	//input the guess parameter, the exact solution
		//will be returned when function's return value is 0.
){
	double f_new=(*this)(t);
	double f_new2=(*this)(t+m_delta);
	if(f_new*f_new2>0. && fabs(f_new)<fabs(f_new2))
		m_delta*=-1.;//since incremental was the opposite direction
					//(1st guess).
	int rc=getCurveParameter2(t);
	if(rc){
		m_delta*=-1.;//Try another direction(2nd guess).
		rc=getCurveParameter2(t);
	}
	return rc;
}

#define MAX_LOOP 20
//Get the paramer value of f(t)=0 that is defined by operator()(double)=0.
//Function's return code is:
//0: when the solution is successfully obtaine in t,
//1: solution not obtained for the guess parameter t(try with other guess parameter).
//-2:system error(usually this does not occur. If occured, some bugs are included.)
int MGCurveParameter::getCurveParameter2(
	double& t	//input the guess parameter, the exact solution
		//will be returned when function's return value is 0.
)const{
	double tnew=t, told;
	double f_new=(*this)(tnew), f_old;
	int loop_counter=0;
	do{
		told=tnew; f_old=f_new;
		tnew+=m_delta;
		tnew=m_prange.round_into_interval(tnew);
		f_new=(*this)(tnew);
		if(f_old*f_new<=0.){
			int ierr;
			t=mgNlbit(*this,told,tnew,m_error,20,ierr);
			if(ierr){
				return -2;
			}
			return 0;
		}
		if(++loop_counter>MAX_LOOP)
			break;
	}while(tnew>m_prange[0] && tnew<m_prange[1]);
	return 1;
}
