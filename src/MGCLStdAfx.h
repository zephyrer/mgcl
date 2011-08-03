/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// MGCLStdAfx.h : Include File to define standard indespencible include files
#if !defined(AFX_STDAFX_H__3869FDFD_2FCB_4F23_A280_3EB784B557F3__INCLUDED_)
#define AFX_STDAFX_H__3869FDFD_2FCB_4F23_A280_3EB784B557F3__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// C/C++ standard library
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <algorithm>
#include <bitset>
#include <deque>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

// Windows specific library
#if !defined(MGCL_NO_MFC)
#define VC_EXTRALEAN // Windows ヘッダーから殆ど使用されないスタッフを除外します。

#include <afxwin.h> // MFC のコアおよび標準コンポーネント
#include <afxext.h> // MFC の拡張部分
#include <afxdisp.h> // MFC のオートメーション クラス
#include <afxdtctl.h> // MFC の Internet Explorer 4 コモン コントロール サポート

#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h> // MFC の Windows コモン コントロール サポート
#endif // _AFX_NO_AFXCMN_SUPPORT

#if !defined(ULONG_PTR)
#define ULONG_PTR  ULONG   // typedef?
#endif  // ULONG_PTR

#include <gdiplus.h>
#include <Gdiplusimaging.h>

// OpenGL
#include <gl/gl.h>
#include <gl/glu.h>

#endif // !defined(MGCL_NO_MFC)

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ は前行の直前に追加の宣言を挿入します。

#endif // !defined(AFX_STDAFX_H__3869FDFD_2FCB_4F23_A280_3EB784B557F3__INCLUDED_)
