/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "MGCLStdAfx.h"
#include "mg/tolerance.h"
#include "mg/OfStream.h"
#include "mg/IfStream.h"
#include "MGStl.h"
#include "Tl/TLPoints.h"
#include "Tl/TLRect.h"

#include "Tl/TLDataVector.h"
#include "Tl/TLTriangle.h"

// 三角形の各頂点である3つのMGPositionを計算し、三角形の面法線ベクトルを取得する
// position1 - position0、position2 - position0
// という計算で2つのベクトルを算出し、次にこれらのベクトルの外積を求める
// そしてこれらのベクトルを正規化することにより三角形の面法線ベクトルを求め
// これを戻り値として返却する
MGUnit_vector CalcNormal(
	const MGPosition& position0, // 三角形の頂点の座標
	const MGPosition& position1,
	const MGPosition& position2
){
	MGVector vector1(position1 - position0);
	MGVector vector2(position2 - position0);
	return vector1 * vector2;
}

// ファイルの種類を調べる
// in[in]：ファイルストリーム
// 戻り値：バイナリファイル(true)、Asciiファイル(false)
// 事前条件：すでにオープンされているファイルストリームを渡す
bool IsBinaryFile(
	 std::ifstream& in // [in/out]:既にオープンしているファイルストリーム
){
	char buff[85];
	in.read(buff, 85);
	int readByte = in.gcount();
	for(int i = 80; i < readByte; i++){
		if(buff[i] == 0){
			// シーク位置をファイルの先頭に戻す
			in.seekg(0, std::ios::beg);
			return true;
		}
	}
	// シーク位置をファイルの先頭に戻す
	in.seekg(0, std::ios::beg);
	return false;
}

// コピーコンストラクタ
MGStl::MGStl(const MGStl& stl)
:MGObject(stl),m_box(stl.m_box), m_vecPos(stl.m_vecPos),
m_vecNormlTriang(stl.m_vecNormlTriang), m_indices(stl.m_indices){
}

//constructor from triangle data, index+vertices 
MGStl::MGStl(
	int nTriang, // 三角形の数
	const int* triang, // 頂点座標のインデックスの配列
	const double* verts // 頂点の座標
){
	// トレランスの値を求め、設定する(wc_zero = box対角線 * rc_zero)
	//double boxLen = m_box.length();
	// トレランスを設定する
	//MGTolerance::set_wc_zero(boxLen*MGTolerance::rc_zero());

	// Keep vertex position and index of vertex
	triangleMap VertexMap;
	// index of vertex
	int i0, i1, i2;
	for(int j = 0; j < nTriang; j++){
		int j3 = 3 * j;
		// Get index of cordinate of vertex
		i0 = triang[j3++] - 1;
		i1 = triang[j3++] - 1;
		i2 = triang[j3] - 1;
		// Get cordinate of vertex
		const MGPosition pos1(3,verts+3*i0);
		const MGPosition pos2(3,verts+3*i1);
		const MGPosition pos3(3,verts+3*i2);
		// まず、入力された3点から面の法線を求め、m_vecNormlTriangへpush_backする
		// 次に、指定した頂点の座標が既に保存されているか調べ
		// 保存されている場合は、該当する頂点のインデックスを受け取る
		// これをm_indicesへpush_backする
		// 保存されていない場合は、頂点をm_vecPosへpush_backし
		// 次にインデックスを+1し、それを受け取りm_indicesにpush_backする
		push_back_triangle(pos1, pos2, pos3, VertexMap);
	}
}

MGStl::MGStl(const mgTLData& tlData){// 変換コンストラクタ
	//IdentifyPosition使用用
	//三角形の頂点の座標、頂点の番号を格納するmap
	triangleMap VertexMap;

	AddTLData(tlData, VertexMap);
}

/////////KUROKI added
// mgTLDataVectorを引数にとるコンストラクタ
// MGShellのメッシュ化対応用
MGStl::MGStl(const mgTLDataVector& tlDataVector){
	
	//IdentifyPosition使用用
	//三角形の頂点の座標、頂点の番号を格納するmap
	triangleMap VertexMap;

	for(size_t i=0; i<tlDataVector.m_datas.size(); i++){
		AddTLData(tlDataVector.tldata(i), VertexMap);
	}
}

// デフォルトデストラクタ
MGStl::~MGStl(void){}

// 全ての頂点とボックスの座標に指定したMGVectorの値を加算
MGStl& MGStl::operator+=(const MGVector& v){
	int nPos = m_vecPos.size();
	for(int i = 0; i < nPos; i++){
		m_vecPos[i] += v;
	}
	m_box += v;
	return *this;
}

// 全ての頂点とボックスの座標から指定したMGVectorの値を引く
MGStl& MGStl::operator-=(const MGVector& v){
	int nPos = m_vecPos.size();
	for(int i = 0; i < nPos; i++){
		m_vecPos[i] -= v;
	}
	m_box -= v;
	return *this;
}

// 全ての頂点とボックスの座標に指定した値を掛ける
MGStl& MGStl::operator*=(double scale){
	int nPos = m_vecPos.size();
	for(int i = 0; i < nPos; i++){
		m_vecPos[i] *= scale; 
	}
	m_box *= scale;
	return *this;
}

// 与えられた変換を行う
MGStl& MGStl::operator*=(const MGMatrix& mat){
	int nVecPos = m_vecPos.size();
	for(int i = 0; i < nVecPos; i++){
		m_vecPos[i] *= mat;
	}
	int nVecTriang = m_vecNormlTriang.size();
	for(int j = 0; j < nVecTriang; j++){
		m_vecNormlTriang[j] *= mat;
	}
	m_box *= mat;
	return *this;
}

// 与えられた変換によるトランスフォームを行う
MGStl& MGStl::operator*=(const MGTransf& tr){
	int nPos = m_vecPos.size();
	for(int i = 0; i < nPos; i++)	{
		m_vecPos[i] *= tr;
	}
	int nVecTriang = m_vecNormlTriang.size();
	for(int j = 0; j < nVecTriang; j++)	{
		m_vecNormlTriang[j] *= tr.affine();
	}
	m_box *= tr;
	return *this;
}

// 二つのMGStlオブジェクトが等しいか比較する
bool MGStl::operator==(const MGStl& stl){
	// ボックス枠の比較
	if(this->m_box != stl.m_box)
		return false;

	// 座標の比較
	if(m_vecPos != stl.m_vecPos){
			return false;
	}

	// 面の法線ベクトルの比較
	if(m_vecNormlTriang != stl.m_vecNormlTriang){
		return false;
	}

	if(m_indices != stl.m_indices){
		return false;
	}

	// 一致している場合はtrueを返却する
	return true;
}

//////////KUROKI added
// 2つのMGStlオブジェクトを代入する
// ただし異なるクラスの場合，代入は行わない
MGStl& MGStl::operator=(const MGStl& stl){
	if(this==&stl)
        return *this;

	MGObject::operator=(stl);	//親クラスの代入文の実行

	m_box=stl.m_box;
	m_vecPos=stl.m_vecPos;
	m_vecNormlTriang=stl.m_vecNormlTriang;
	m_indices=stl.m_indices;
	return *this;
}

//Construct new object by copying to newed area.
//User must delete this copied object by "delete".
MGStl* MGStl::clone()const{
	return new MGStl(*this);
}

//Draw 3D point(vertex) in world coordinates.
//The object is converted to point(s) and is drawn.
//This is valid only for topology objects or MGPoint.
void MGStl::draw3DVertex()const{;}

// 指定した三角形の頂点座標が入れられている配列のインデックスを取得する
// 事前条件:0から三角形の枚数-1を超えた範囲の値を第１引数に指定しないこと
// 事後条件:pos[3]に指定した三角形の各頂点座標の配列の添え字が格納される
void MGStl::GetVertIndices(
	 int i,// [in]：三角形のインデックス(0 <= i < GetTriangleCount)
	 size_t pos[3]// [out]：指定した三角形の各頂点座標の配列の添え字[i, j, k]
)const{
	// 指定した三角形の頂点番号の配列のインデックスを求めておく
	size_t i3 = i*3;
	assert(i3 + 2 < m_indices.size());

	// 頂点座標のインデックスを代入
	pos[0] = m_indices[i3++];
	pos[1] = m_indices[i3++];
	pos[2] = m_indices[i3];
}

// 引数で指定したパスのSTLファイルを読み込みメンバに値を設定する
// またMGTolerance::wc_zeroに値を設定する
// 戻り値: =0ファイルの読み込が成功
//		  !=0 読み込まなかった。または失敗した(std::ifstreamのエラーコード）
// 事後条件:m_vecPos, m_vecNorml, m_indices, m_boxに図形の情報が格納される
//			MGTolerance::wc_zeroに値が設定される
int MGStl::LoadFile(
	const char* strFilePath // [in]:読み込むSTLファイルへのパス
){
	// ストリームをオープン
	std::ifstream in(strFilePath, std::ios_base::binary);
	if(!in){
		int state = in.rdstate();
		// オープンが失敗した場合、エラーコードを返却
		return state;
	}

	int error; // ファイル読み込みの結果を保持する
	std::vector<MGPosition> vecPos;
	MGBox box;
	if(IsBinaryFile(in)){
		// Load Binary File
		error = LoadBinary(in, vecPos);
	}else{
		// Load Ascii File
		error = LoadAscii(in, vecPos);
	}
	
	if(error){ // 0:正常に読み込みが完了 それ以外の値:読み込みが失敗
		// 読み込みが失敗した場合、エラーコードを返却
		return error;
	}

	// トレランスを設定
	double wc_zero_old = MGTolerance::set_wc_zero(m_box.length()*MGTolerance::rc_zero());
	// 読み込んだ座標値からm_indices, m_vecNormlTriangを設定
	set_mesh_data(vecPos);
	MGTolerance::set_wc_zero(wc_zero_old);

	// デバッグコード
	//std::cout << "stlファイル読み込み" << std::endl;
	//std::cout << "三角形数 " << m_vecNormlTriang.size() << std::endl;
	//std::cout << "頂点数 " << m_vecPos.size() << std::endl;
	//std::cout << "ボックス枠座標"<<  std::endl;
	//std::cout << "最小値 " << std::endl;
	//std::cout << "X " << m_box.low()(0) << std::endl;
	//std::cout << "Y " << m_box.low()(1) << std::endl;
	//std::cout << "Z " << m_box.low()(2) << std::endl;
	//std::cout << "最大値 " << std::endl;
	//std::cout << "X " << m_box.high()(0) << std::endl;
	//std::cout << "Y " << m_box.high()(1) << std::endl;
	//std::cout << "Z " << m_box.high()(2) << std::endl << std::endl;
	// エラーコードを返却(正常に読み込みが完了)
	return error;
}

// ファイルから読み込んだ全ての頂点の座標値が格納されている
// vecPosから各三角形の法線を計算し、m_vecNormlTriangに追加
// vecPosの各要素の座標値の重複をトレランスを元に判断し
// 重複を取り除き、m_vecPosに追加
// vecPosの各要素のインデックスをm_indicesに追加
// 事後条件:m_vecNormlTriang, m_indicesに値が設定される
void MGStl::set_mesh_data(
	const std::vector<MGPosition>& vecPos // [in]:ファイルから読み込んだ座標値の配列
){
	// ファイルから読んだ頂点数から三角形の数を取得
	int nTriang = vecPos.size()/3;
	// 座標値と頂点のインデックスを一時的に保持しておくマップ
	triangleMap VertexMap;
	for(int i = 0; i < nTriang; i++){
		// 三角形の頂点番号の配列のインデックスを求めておく
		int i3=i*3;
		const MGPosition& p0=vecPos[i3++];
		const MGPosition& p1=vecPos[i3++];
		const MGPosition& p2=vecPos[i3];
		// まず、入力された3点から面の法線を求め、m_vecNormlTriangへpush_backする
		// 次に、指定した頂点の座標が既に保存されているか調べ
		// 保存されている場合は、該当する頂点のインデックスを受け取る
		// これをm_indicesへpush_backする
		// 保存されていない場合は、頂点をm_vecPosへpush_backし
		// 次にインデックスを+1し、それを受け取りm_indicesにpush_backする
		push_back_triangle(p0, p1, p2, VertexMap);
	}
}

// Ascii形式のファイルを読み込み全ての座標値を取得する
// また図形のボックス枠を設定する
// 戻り値: =0ファイルの読み込が成功
//		  !=0 読み込みに失敗した(std::ifstreamのエラーコード)
// 事前条件:既にオープンされたファイルストリームを渡す
// 事後条件:vecPosに全ての座標値が格納され、 boxにボックス枠の座標値が設定される
// ファイルストリームが進む
int MGStl::LoadAscii(
	std::ifstream& in, // [in/out]:既にオープンされたファイルストリーム
	std::vector<MGPosition>& vecPos // [out]:ファイルから読み込んだ座標値の配列
){
	// メンバの値をリセット
	Initialize();
	// Keep one line data read from file
	std::string strLine;
	// Keep cordinate of vertex
	MGPosition position1(3), position2(3), position3(3);

	if(in.eof()){ // ファイルが空の場合
		return std::ios_base::eofbit;
	}

	// ファイル内にfacet normalがあるかを示す変数
	bool bStl = false;
	// エラーチェックの結果を格納する変数
	int readState;
	// Read file and set date to member
	while(getline(in, strLine)){
		// 読み込みエラーチェック
		readState=in.rdstate();
		if(readState == std::ios::failbit || readState == std::ios::badbit){
			return readState; // エラーコードを返却(読み込みが失敗)
		}

		if(strLine.find("facet normal") == -1){// 行内にfacet normalが存在するか
			continue;
		}
		// ファイル内にfacet normalを発見
		bStl = true;
		// Triangle data
		std::string dummy;
		// outer loop行を読み飛ばす
		in.ignore(10000, '\n');
		// Get cordinate of vertex
		in >> dummy >> position1(0) >> position1(1) >> position1(2);
		in >> dummy >> position2(0) >> position2(1) >> position2(2);
		in >> dummy >> position3(0) >> position3(1) >> position3(2);

		// 読み込みエラーチェック
		readState = in.rdstate();
		if(readState){
			return readState; // エラーコードを返却(読み込みが失敗)
		}

		// ファイルから読み込んだ座標値を引数のvectorにpush_backする
		vecPos.push_back(position1);
		vecPos.push_back(position2);
		vecPos.push_back(position3);
		// ボックス枠を3次元とし、expand
		m_box.expand(position1);
		m_box.expand(position2);
		m_box.expand(position3);
		// 改行する
		in.seekg(1, std::ios::cur);
		// endloop endfacet行を読み飛ばす
		in.ignore(10000, '\n');
		in.ignore(10000, '\n');
	}

	if(!bStl){
		// facet normalを一度も読まずに末行まで読み込んだ場合
		// つまり、別の種類のファイルを読んだ場合の処理
		return std::ios_base::eofbit;
	}

	// 正常終了を示す値を返却
	return 0;
}

// Binary形式のファイルを読み込み全ての座標値を取得する
// また図形のボックス枠を設定する
// 戻り値: =0ファイルの読み込が成功
//		  !=0 読み込みに失敗した(std::ifstreamのエラーコード)
// 事前条件:既にオープンされたファイルストリームを渡す
// 事後条件:vecPosに全ての座標値が格納され、 m_boxにボックス枠の座標値が設定される
//			ファイルストリームが進む
int MGStl::LoadBinary(
	std::ifstream& in, // [in/out]:既にオープンされたファイルストリーム
	std::vector<MGPosition>& vecPos //[out]:ファイルから読み込んだ座標値の配列
){
	// メンバの値をリセット
	Initialize();
	// Keep vertex and normal of triangle temporary 
	float fVertex1[3], fVertex2[3], fVertex3[3];
	// Keep the point of vertex of triangle
	MGPosition position1(3), position2(3), position3(3);	

	// Skip the top of 80 byte
	in.seekg(80, std::ios::beg);
	// Get count of triangle
	unsigned int nTriang;
	in.read((char*)&nTriang, 4);

	// 書き込みエラーチェック
	int readState = in.rdstate();
	if(readState){
		return readState; // エラーコードを返却(読み込みが失敗)
	}

	// Read and set cordinate of vertex and normal of triangle
	for(unsigned int i = 0; i < nTriang; i++){
		// normalが格納されている部分を読み飛ばす
		in.seekg(12, std::ios::cur);
		// Read cordinate of vertex from file
		for(int j = 0; j < 3; j++){
			in.read((char*)&fVertex1[j], 4);
			position1(j) = fVertex1[j];
		}
		for(int j = 0; j < 3; j++){
			in.read((char*)&fVertex2[j], 4);
			position2(j) = fVertex2[j];
		}
		for(int j = 0; j < 3; j++){
			in.read((char*)&fVertex3[j], 4);
			position3(j) = fVertex3[j];
		}

		// ファイルから読み込んだ座標値を引数のvectorにpush_backする
		vecPos.push_back(position1);
		vecPos.push_back(position2);
		vecPos.push_back(position3);
		// ボックス枠をexpand
		m_box.expand(position1);
		m_box.expand(position2);
		m_box.expand(position3);

		if(i == nTriang - 1){
			// 最後の三角形の場合は読み込み終了
			break;
		}
		// 次の三角形がある場合、2バイトスキップする
		in.seekg(2, SEEK_CUR);

		// 読み込みエラーチェック
		readState = in.rdstate();
		if(readState){
			return readState; // エラーコードを返却(読み込みが失敗)
		}
	}

	// 書き込みエラーチェック
	readState = in.rdstate();
	if(readState == std::ios::failbit || readState == std::ios::badbit){
		return readState; // エラーコードを返却(読み込みが失敗)
	}

	// 正常終了を示す値を返却
	return 0;
}

// 指定されたパスにAscii形式のSTLファイルを保存する
// 戻り値: =0ファイルの書き込みが成功
//		  !=0 書き込まなかった。または失敗した(std::ofstreamのエラーコード)
// 事前条件:m_vecPosとm_vecNormlとm_indicesに正しい情報が入れられている
// 事後条件:rSTLFilePathに指定したパスにAscii形式のSTLファイルが保存される
//仕様変更：引数にストリームを渡されるように変更
int MGStl::SaveAscii(
	std::ofstream& fout
	//const char* rSTLFilePath // [in]:ファイルの保存先のパス
)const{
	//std::ofstream fout(rSTLFilePath);

	// 書き込みエラーチェック
	int writeState = fout.rdstate();
	if(writeState){
		return writeState; // エラーコードを返却(読み込みが失敗)
	}

	fout.setf(std::ios::scientific);
	// solid
	fout << "solid ascii" << std::endl;
	// 三角形の数を取得(三角形の数 == 面法線ベクトルの個数)
	int nTriang = m_vecNormlTriang.size();
	// output normal of triangle and cordinate of vertex to file
	for(int i = 0; i < nTriang; i++){
		// Get triangle normal
		const MGUnit_vector& normal = m_vecNormlTriang[i];
		// facet normal
		fout << "    facet normal  " << normal(0) << " " << normal(1)<< " " << normal(2) << " " << std::endl;
		// outerloop
		fout << "        outer loop" << std::endl;
		
		// 頂点の番号の配列の添え字を計算しておく
		int index=i*3;
		// output cordinate of vertex to file
		const MGPosition& position1 = m_vecPos[m_indices[index++]];
		fout << "        vertex  " << position1(0) << " " << position1(1) << " " << position1(2) << " " << std::endl;
		const MGPosition& position2 = m_vecPos[m_indices[index++]];
		fout << "        vertex  " << position2(0) << " " << position2(1) << " " << position2(2) << " " << std::endl;
		const MGPosition& position3 = m_vecPos[m_indices[index]];
		fout << "        vertex  " << position3(0) << " " << position3(1) << " " << position3(2) << " " << std::endl;
		// outerloop
		fout << "        endloop" << std::endl;
		fout << "    endfacet" << std::endl;

		// 書き込みエラーチェック
		writeState = fout.rdstate();
		if(writeState){
			return writeState; // エラーコードを返却(読み込みが失敗)
		}
	}
	// endsolid
	fout << "endsolid" << std::endl;
	fout.unsetf(std::ios::scientific);

	writeState = fout.rdstate();
	if(writeState == std::ios::failbit || writeState == std::ios::badbit){
		return writeState; // エラーコードを返却(書き込みが失敗)
	}

	//fout.close();
	return 0;
}

// 指定されたパスにBinary形式のSTLファイルを保存する
// 戻り値: =0ファイルの書き込みが成功
//		  !=0 書き込まなかった。または失敗した(std::ofstreamのエラーコード)
// 事前条件:m_vecPosとm_vecNormlとm_indicesに正しい情報が入れられている
// 事後条件:rSTLFilePathに指定したパスにBinary形式のSTLファイルが保存される
int MGStl::SaveBinary(
	const char* rSTLFilePath  // [in]:ファイルの保存先のパス
)const{
	// File open
	std::ofstream fout(rSTLFilePath, std::ios::binary);
	// ストリームのエラーチェック
	int writeState = fout.rdstate();
	if(writeState)
		return writeState;

	// 三角形の数を取得(三角形の数 == 面法線ベクトルの個数)
	unsigned int nTriang = m_vecNormlTriang.size();
	// ファイル先頭80バイトをスキップ
	fout.seekp(80, std::ios::beg);
	// ファイルに三角形の枚数を書き込む
	fout.write((char*)&nTriang, 4);

	// 書き込みエラーチェック
	writeState = fout.rdstate();
	if(writeState){
		return writeState; // エラーコードを返却(書き込みに失敗)
	}

	// Keep vertex and normal of triangle temporary 
	float fNormal, fVertex;
	// 各三角形の情報を書き込む
	for(unsigned int i = 0; i < nTriang; i++){
		// 三角形の面法線ベクトルを書き込む
		const MGUnit_vector& normal = m_vecNormlTriang[i];
		for(int j = 0; j < 3; j++){
			fNormal = float(normal[j]);
			fout.write((char*)&fNormal, 4);
		}

		// 三角形の頂点番号の配列のインデックスを求めておく
		int i3 = i*3;

		// ベクターから頂点のインデックスが指している座標値を取り出し
		// 三角形の頂点の座標を書き込む
		const MGPosition& position1 = m_vecPos[m_indices[i3]];
		for(int j = 0; j < 3; j++){
			fVertex = float(position1[j]);
			fout.write((char*)&fVertex, 4);
		}
		
		const MGPosition& position2 = m_vecPos[m_indices[i3+1]];
		for(int j = 0; j < 3; j++)		{
			fVertex = float(position2[j]);
			fout.write((char*)&fVertex, 4);
		}
		
		const MGPosition& position3 = m_vecPos[m_indices[i3+2]];
		for(int j = 0; j < 3; j++){
			fVertex = float(position3[j]);
			fout.write((char*)&fVertex, 4);
		}

		// 最後の三角形の場合は書き込みを終了
		if(i == nTriang - 1){
			break;
		}

		// 三角形ごとに2バイト間隔を空ける
		fout.seekp(2, std::ios::cur);

		// 書き込みエラーチェック
		writeState = fout.rdstate();
		if(writeState){
			return writeState; // エラーコードを返却(書き込みに失敗)
		}
	}
	// 書き込みエラーチェック
	writeState = fout.rdstate();
	if(writeState == std::ios::failbit || writeState == std::ios::badbit){
		return writeState;// エラーコードを返却(書き込みが失敗)
	}

	fout.close();
	return 0;
}

// 入力であるpositionがVertexMapにすでに登録されていればそのm_vecPos添え字を返し
// 未登録の場合、新たにm_vecPosに格納し、positionとm_vecPos添え字のmapをVertexMapに
// 登録する
int MGStl::IdentifyPosition(
	const MGPosition& position, // 頂点の座標
	triangleMap& VertexMap // 重複のない頂点の座標、頂点のインデックスを保持する
){
	int nVertexIndex = m_vecPos.size();
	//iteratorをpositionの座標値以上となる要素がはじめて現れる場所を指すようにする
	triangleMap::iterator itAft = VertexMap.lower_bound(position), itPre;
	if(itAft != VertexMap.end()){	//iteratorが終わりを指していない場合
		if(itAft->first==position)
			return itAft->second;	//指している要素と一致した場合，そのindexを返す
	}
	if(itAft != VertexMap.begin()){	//iteratorが始めを指していない場合
		itPre=itAft;itPre--;	//一つ前の要素とpositionを比較する
		if(itPre->first==position){
			return itPre->second;	//一つ前の要素と一致した場合，そのindexを返す
		}
	}

	VertexMap.insert(itAft, std::make_pair(position, nVertexIndex));
	m_vecPos.push_back(position);
	return nVertexIndex;

	//original
	//triangleMap::iterator it = VertexMap.lower_bound(position);
	//if(it == VertexMap.end() || !(it->first==position)){//If not found.
	//	int nVertexIndex = m_vecPos.size();
	//	it = VertexMap.insert(it, std::make_pair(position, nVertexIndex));
	//	m_vecPos.push_back(it->first);
	//	
	//	std::cout<<"m_vecPos[i]"<<std::endl;
	//	for(size_t i=0; i<m_vecPos.size(); i++){
	//		std::cout<<m_vecPos[i]<<std::endl;
	//	}

	//	return nVertexIndex;
	//}else{//if found
	//	return it->second;
	//}
}

// メンバデータを書き込む関数
void MGStl::WriteMembers(MGOfstream& buf)const{
	MGObject::WriteMembers(buf);
	// 三角形の法線ベクトルを書き込む
	int nTriangle = m_vecNormlTriang.size();
	buf << nTriangle;
	for(int i = 0; i < nTriangle; i++){
		m_vecNormlTriang[i].dump(buf);
	}

	// 三角形の重複がない座標の頂点を書き込む
	int nPos = m_vecPos.size();
	buf << nPos;
	for(int k = 0; k < nPos; k++){
		m_vecPos[k].dump(buf);
	}

	// 三角形の各頂点のインデックスを書き込む
	int nIndex = m_indices.size();
	buf << nIndex;
	for(int l = 0; l < nIndex; l++){
		buf << m_indices[l];
	}
}

// メンバデータを読み込む関数
void MGStl::ReadMembers(MGIfstream& buf){
	MGObject::ReadMembers(buf);
	// 三角形の法線ベクトルを読み込む
	int nTriangle;
	buf >> nTriangle;
	m_vecNormlTriang.resize(nTriangle);
	for(int i = 0; i < nTriangle; i++){	
		m_vecNormlTriang[i].restore(buf);
	}

	// 三角形の各頂点の座標を重複せず読み込む
	int nPos;
	buf >> nPos;
	m_vecPos.resize(nPos);
	for(int k = 0; k < nPos; k++){
		m_vecPos[k].restore(buf);
	}

	// 三角形の各頂点のインデックスを読み込む
	int nIndex;
	buf >> nIndex;
	m_indices.resize(nIndex);
	for(int l = 0; l < nIndex; l++){
		buf >> m_indices[l];
	}
}

// 一度ファイルを読み、メンバに値が設定された状態で再びファイルを読むと
// メンバの値が上書きされる。これを避けるため、メンバの値を全てリセットする
void MGStl::Initialize(){
	m_box.set_null();
	m_vecPos.clear();
	m_vecNormlTriang.clear();
	m_indices.clear();
}

// 3点から作成した三角形の情報をメンバ変数に値を設定する
// 事前条件:入力するpos1, pos2, pos3は反時計回りになっていること
// 事後条件:m_vecPos, m_vecNormlTriang, m_indicesに三角形の情報が設定される
void MGStl::push_back_triangle(
	const MGPosition& pos1,
	const MGPosition& pos2,
	const MGPosition& pos3,
	triangleMap& VertexMap
){
	// 3点から面の法線を求め、メンバ変数にpush_backする
	m_vecNormlTriang.push_back(CalcNormal(pos1, pos2, pos3));
	// m_indicesに頂点のインデックスをpush_backし
	// m_vecPosに重複なくMGPositionをpush_backする
	m_indices.push_back(IdentifyPosition(pos1, VertexMap));
	m_indices.push_back(IdentifyPosition(pos2, VertexMap));
	m_indices.push_back(IdentifyPosition(pos3, VertexMap));
}

//TLDataをMGStlに追加する
void MGStl::AddTLData(const mgTLData& tlData, triangleMap& VertexMap){
	// サーフェイスを取得
	const MGSurface& surf = tlData.surface();
	const mgTLparameter& tlp = tlData.tlparam();
	double errorSave=MGTolerance::set_wc_zero(tlp.get_surface_error());

	// 全ての頂点のuv値を取得、mgTLDataの全頂点のuv値をMGPositonに変換し、m_vecPosに格納
	//対応表indexVecを作成する。頂点IDの対応を収めるvector:indexVec[旧ID]=新IDで対応をとる．
	const mgTLPoints& tlpoints = tlData.tlpoints();
	size_t npointOld=tlpoints.size();
	std::vector<int> indexVec(npointOld);
	m_box.set_null();
	mgTLPoints::const_iterator ip = tlpoints.begin();
	for(size_t oldID=0;oldID<npointOld;++oldID, ip++){
		MGPosition pos = surf.eval(*ip);
		indexVec[oldID]=IdentifyPosition(pos, VertexMap);
		m_box.expand(pos);
	}

	// 全てのポリゴンを取得(ファンもしくはストリップ)
	const mgTLTriangles& triangles = tlData.triangles();
	// ポリゴン配列のbeginとendを取得
	mgTLTriangles::const_triIterator j = triangles.begin(), je = triangles.end();
	for(;j!=je;++j){
	// ポリゴンごとにループ(indicesを取得する)
		const mgTLTriangle& tri=**j;// ポリゴンを取得
		size_t nVert = tri.size();// 頂点の個数を取得する
		if(nVert<3)
			continue;

		
		mgTESTRIANG geoType = tri.getGeometryType();//ポリゴンの形状(ファンもしくはストリップ)
		size_t nVm2 = nVert-2;
		size_t id0,id1,id2,id0uni,id1uni,id2uni;
		for(size_t i=0; i<nVm2; i++){
			if(geoType==mgTESTRIANG_FAN){
				id0=tri[0];
				id1=tri[i+1];
				id2=tri[i+2];
				id0uni=indexVec[id0];
				id1uni=indexVec[id1];
				id2uni=indexVec[id2];
			}else if(geoType==mgTESTRIANG_STRIP){
				id0=tri[i];
				id1=tri[i+1+(i%2)];
				id2=tri[i+2-(i%2)];
				id0uni=indexVec[id0];
				id1uni=indexVec[id1];
				id2uni=indexVec[id2];
			}
			if(m_vecPos[id0uni]==m_vecPos[id1uni])
				continue;
			if(m_vecPos[id0uni]==m_vecPos[id2uni])
				continue;
			if(m_vecPos[id1uni]==m_vecPos[id2uni])
				continue;
			m_indices.push_back(id0uni);
			m_indices.push_back(id1uni);
			m_indices.push_back(id2uni);

			// 三角形の中点uv値から法線ベクトルを求める
			MGPosition uv=(tlpoints[id0]+tlpoints[id1]+tlpoints[id2])/3.;
			m_vecNormlTriang.push_back(surf.unit_normal(uv));
		}
	}
	MGTolerance::set_wc_zero(errorSave);
}

// メンバーデータを直接セットする
void MGStl::set_all_data(
	const std::vector<MGPosition>& vertices,
	const std::vector<int>& indices)
{
	assert(indices.size() % 3 == 0);

	m_vecPos = vertices;
	m_indices = indices;

	// update normals and box
	update_bounds();
	update_normals();
}

// Update bounds.
void MGStl::update_bounds()
{
	MGBox work(3);
	const int nVert = m_vecPos.size();
	for(int i = 0; i < nVert; ++i){
		work.expand(m_vecPos[i]);
	}
	m_box = work;
}

// Update the all normals of the triangles.
void MGStl::update_normals()
{
	const int nTri = m_indices.size() / 3;
	std::vector<MGUnit_vector> work(nTri);

	for(int i = 0; i < nTri; ++i){
		size_t vid[3];
		GetVertIndices(i, vid);
		assert(vid[0] < m_vecPos.size());
		assert(vid[1] < m_vecPos.size());
		assert(vid[2] < m_vecPos.size());

		work[i] = CalcNormal(
			m_vecPos[vid[0]],
			m_vecPos[vid[1]],
			m_vecPos[vid[2]]);
	}

	work.swap(m_vecNormlTriang);
}
