/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//! @file
//!	@brief  Implementaion for class MGIgesOfstream.
//!	@author DG Technologies(http://www.dgtech.co.jp/)

#include "MGCLStdAfx.h"
#include "mgIges/IgesOfstream.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif
using namespace std;

#define LINE_LENGTH 72
//Write out start and terminate section into output file stream.
void MGIgesOfstream::write_out_start_section(){
	int i=0, n=m_StartSection.length();
		//i is the position in m_StartSection to output next line out of Start section.

	int lineID=1;
	while(n-i>=LINE_LENGTH){
		m_ofstream<<m_StartSection.substr(i,LINE_LENGTH);
		m_ofstream<<'S'<<setw(7)<<setfill('0')<<lineID<<endl;
		i+=LINE_LENGTH;
		lineID++;
	}
	int redundant=n-i;
	if(redundant>0){
		m_ofstream<<left<<setw(LINE_LENGTH)
			<<setfill(' ')<<m_StartSection.substr(i,redundant);
		m_ofstream<<'S'<<setw(7)<<right<<setfill('0')<<lineID<<endl;
	}
}

//Write out all the Paramete Data Lines in m_plines to output file stream.
//m_plines[i] is one parameter data lines of IGES file
//for 0<=i<m_plines.size().
void MGIgesOfstream::write_out_PD_plines(){
	int plsize=m_plines.size();//パラメータデータラインの列の数を代入
	for(int i=0;i<plsize;i++){//その数だけループ
		MGIgesParamLine& igespline=*(m_plines[i]);//i番目のm_plines[i]の参照先の値を代入
		auto_ptr<string>& strng=igespline.paramLine();
			//paramLine()というアクセッサを使って今代入したものが入ってるm_paramLineを参照
		string line=*(strng);
		m_ofstream<<setw(65)<<setfill(' ')<<left<<line;
		int de=igespline.DEpointer();
		m_ofstream<<setw(7)<<setfill('0')<<right<<MGIges::DEpointer_to_lnumber(de); 
		m_ofstream<<'P';
		int ip1=i+1;
		m_ofstream<<setw(7)<<setfill('0')<<right<<ip1<<endl;
	}
}

//Write out DEpointer's Paramete Data Lines into m_plines.
//Paramete Data Lines of m_DrectoryEntries[DEpointer] will be
//output into m_plines.
//Function's return value is the number of lines output.
int MGIgesOfstream::write_out_PD_pline(int DEpointer){
	MGIgesDirectoryEntry* de=directoryEntry(DEpointer);
	const auto_ptr<MGIgesPD>& pd=de->paramData();

	MGPvector<string> plines;
	MGIges::put_integer(pd->type_number(),m_GSection,plines);
	pd->write_out_into_string(m_GSection,plines);
	MGIges::append_record_delimeter(m_GSection.recordDelimeter(),plines);

	int n=plines.size();
	for(int i=0; i<n; i++){
		auto_ptr<string> one_line(plines.release(i));
		m_plines.push_back(new MGIgesParamLine(one_line,DEpointer));
	}
	return n;
}

//Write out DE lines to output file stream, together with PD lines.
//write_out_DE_PD_lines outputs all the DEs and PDs using 
//write_out_PD_pline for each DE.
void MGIgesOfstream::write_out_DE_PD_lines(){
	//ＤＥの配列のサイズを求め、すべてのＤＥをplineに出力している
	int sum_DE = m_DirectoryEntries.size();
	for(int i=1;i<sum_DE;i++){//i starts from 1 since 1st DE is dummy.
		int PDpointer=get_next_param_line_count();
		int line_count=write_out_PD_pline(i);
		MGIgesDirectoryEntry& deNum = *(m_DirectoryEntries[i]);
		string DElines[2];
		deNum.put_to_string(PDpointer, line_count, i,DElines);
		m_ofstream<<DElines[0]<<endl;
		m_ofstream<<DElines[1]<<endl;
	}
	write_out_PD_plines();//出力用のファイルを生成し、ファイルに書き込む
}

void MGIgesOfstream::write_out_terminate_section(){
	int slineSize=MGIgesFstream::get_line_number_of_SS();
	m_ofstream<<'S';
	m_ofstream<<setw(7)<<setfill('0')<<slineSize;

	int nlinegsec=MGIgesFstream::get_line_number_of_GS();
	m_ofstream<<'G';
	m_ofstream<<setw(7)<<setfill('0')<<nlinegsec;

	int dlineSize=MGIgesFstream::get_line_number_of_DE();
	m_ofstream<<'D';
	m_ofstream<<setw(7)<<setfill('0')<<dlineSize;

	int mlineSize=MGIgesOfstream::get_line_number_of_PD();
	m_ofstream<<'P';
	m_ofstream<<setw(7)<<setfill('0')<<mlineSize;

	m_ofstream<<setw(41)<<setfill(' ')<<('T');
	m_ofstream<<setw(7)<<setfill('0')<<'1'<<endl;
}
