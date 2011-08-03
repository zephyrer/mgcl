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

// �O�p�`�̊e���_�ł���3��MGPosition���v�Z���A�O�p�`�̖ʖ@���x�N�g�����擾����
// position1 - position0�Aposition2 - position0
// �Ƃ����v�Z��2�̃x�N�g�����Z�o���A���ɂ����̃x�N�g���̊O�ς����߂�
// �����Ă����̃x�N�g���𐳋K�����邱�Ƃɂ��O�p�`�̖ʖ@���x�N�g��������
// �����߂�l�Ƃ��ĕԋp����
MGUnit_vector CalcNormal(
	const MGPosition& position0, // �O�p�`�̒��_�̍��W
	const MGPosition& position1,
	const MGPosition& position2
){
	MGVector vector1(position1 - position0);
	MGVector vector2(position2 - position0);
	return vector1 * vector2;
}

// �t�@�C���̎�ނ𒲂ׂ�
// in[in]�F�t�@�C���X�g���[��
// �߂�l�F�o�C�i���t�@�C��(true)�AAscii�t�@�C��(false)
// ���O�����F���łɃI�[�v������Ă���t�@�C���X�g���[����n��
bool IsBinaryFile(
	 std::ifstream& in // [in/out]:���ɃI�[�v�����Ă���t�@�C���X�g���[��
){
	char buff[85];
	in.read(buff, 85);
	int readByte = in.gcount();
	for(int i = 80; i < readByte; i++){
		if(buff[i] == 0){
			// �V�[�N�ʒu���t�@�C���̐擪�ɖ߂�
			in.seekg(0, std::ios::beg);
			return true;
		}
	}
	// �V�[�N�ʒu���t�@�C���̐擪�ɖ߂�
	in.seekg(0, std::ios::beg);
	return false;
}

// �R�s�[�R���X�g���N�^
MGStl::MGStl(const MGStl& stl)
:MGObject(stl),m_box(stl.m_box), m_vecPos(stl.m_vecPos),
m_vecNormlTriang(stl.m_vecNormlTriang), m_indices(stl.m_indices){
}

//constructor from triangle data, index+vertices 
MGStl::MGStl(
	int nTriang, // �O�p�`�̐�
	const int* triang, // ���_���W�̃C���f�b�N�X�̔z��
	const double* verts // ���_�̍��W
){
	// �g�������X�̒l�����߁A�ݒ肷��(wc_zero = box�Ίp�� * rc_zero)
	//double boxLen = m_box.length();
	// �g�������X��ݒ肷��
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
		// �܂��A���͂��ꂽ3�_����ʂ̖@�������߁Am_vecNormlTriang��push_back����
		// ���ɁA�w�肵�����_�̍��W�����ɕۑ�����Ă��邩����
		// �ۑ�����Ă���ꍇ�́A�Y�����钸�_�̃C���f�b�N�X���󂯎��
		// �����m_indices��push_back����
		// �ۑ�����Ă��Ȃ��ꍇ�́A���_��m_vecPos��push_back��
		// ���ɃC���f�b�N�X��+1���A������󂯎��m_indices��push_back����
		push_back_triangle(pos1, pos2, pos3, VertexMap);
	}
}

MGStl::MGStl(const mgTLData& tlData){// �ϊ��R���X�g���N�^
	//IdentifyPosition�g�p�p
	//�O�p�`�̒��_�̍��W�A���_�̔ԍ����i�[����map
	triangleMap VertexMap;

	AddTLData(tlData, VertexMap);
}

/////////KUROKI added
// mgTLDataVector�������ɂƂ�R���X�g���N�^
// MGShell�̃��b�V�����Ή��p
MGStl::MGStl(const mgTLDataVector& tlDataVector){
	
	//IdentifyPosition�g�p�p
	//�O�p�`�̒��_�̍��W�A���_�̔ԍ����i�[����map
	triangleMap VertexMap;

	for(size_t i=0; i<tlDataVector.m_datas.size(); i++){
		AddTLData(tlDataVector.tldata(i), VertexMap);
	}
}

// �f�t�H���g�f�X�g���N�^
MGStl::~MGStl(void){}

// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵��MGVector�̒l�����Z
MGStl& MGStl::operator+=(const MGVector& v){
	int nPos = m_vecPos.size();
	for(int i = 0; i < nPos; i++){
		m_vecPos[i] += v;
	}
	m_box += v;
	return *this;
}

// �S�Ă̒��_�ƃ{�b�N�X�̍��W����w�肵��MGVector�̒l������
MGStl& MGStl::operator-=(const MGVector& v){
	int nPos = m_vecPos.size();
	for(int i = 0; i < nPos; i++){
		m_vecPos[i] -= v;
	}
	m_box -= v;
	return *this;
}

// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵���l���|����
MGStl& MGStl::operator*=(double scale){
	int nPos = m_vecPos.size();
	for(int i = 0; i < nPos; i++){
		m_vecPos[i] *= scale; 
	}
	m_box *= scale;
	return *this;
}

// �^����ꂽ�ϊ����s��
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

// �^����ꂽ�ϊ��ɂ��g�����X�t�H�[�����s��
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

// ���MGStl�I�u�W�F�N�g������������r����
bool MGStl::operator==(const MGStl& stl){
	// �{�b�N�X�g�̔�r
	if(this->m_box != stl.m_box)
		return false;

	// ���W�̔�r
	if(m_vecPos != stl.m_vecPos){
			return false;
	}

	// �ʂ̖@���x�N�g���̔�r
	if(m_vecNormlTriang != stl.m_vecNormlTriang){
		return false;
	}

	if(m_indices != stl.m_indices){
		return false;
	}

	// ��v���Ă���ꍇ��true��ԋp����
	return true;
}

//////////KUROKI added
// 2��MGStl�I�u�W�F�N�g��������
// �������قȂ�N���X�̏ꍇ�C����͍s��Ȃ�
MGStl& MGStl::operator=(const MGStl& stl){
	if(this==&stl)
        return *this;

	MGObject::operator=(stl);	//�e�N���X�̑�����̎��s

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

// �w�肵���O�p�`�̒��_���W��������Ă���z��̃C���f�b�N�X���擾����
// ���O����:0����O�p�`�̖���-1�𒴂����͈͂̒l���P�����Ɏw�肵�Ȃ�����
// �������:pos[3]�Ɏw�肵���O�p�`�̊e���_���W�̔z��̓Y�������i�[�����
void MGStl::GetVertIndices(
	 int i,// [in]�F�O�p�`�̃C���f�b�N�X(0 <= i < GetTriangleCount)
	 size_t pos[3]// [out]�F�w�肵���O�p�`�̊e���_���W�̔z��̓Y����[i, j, k]
)const{
	// �w�肵���O�p�`�̒��_�ԍ��̔z��̃C���f�b�N�X�����߂Ă���
	size_t i3 = i*3;
	assert(i3 + 2 < m_indices.size());

	// ���_���W�̃C���f�b�N�X����
	pos[0] = m_indices[i3++];
	pos[1] = m_indices[i3++];
	pos[2] = m_indices[i3];
}

// �����Ŏw�肵���p�X��STL�t�@�C����ǂݍ��݃����o�ɒl��ݒ肷��
// �܂�MGTolerance::wc_zero�ɒl��ݒ肷��
// �߂�l: =0�t�@�C���̓ǂݍ�������
//		  !=0 �ǂݍ��܂Ȃ������B�܂��͎��s����(std::ifstream�̃G���[�R�[�h�j
// �������:m_vecPos, m_vecNorml, m_indices, m_box�ɐ}�`�̏�񂪊i�[�����
//			MGTolerance::wc_zero�ɒl���ݒ肳���
int MGStl::LoadFile(
	const char* strFilePath // [in]:�ǂݍ���STL�t�@�C���ւ̃p�X
){
	// �X�g���[�����I�[�v��
	std::ifstream in(strFilePath, std::ios_base::binary);
	if(!in){
		int state = in.rdstate();
		// �I�[�v�������s�����ꍇ�A�G���[�R�[�h��ԋp
		return state;
	}

	int error; // �t�@�C���ǂݍ��݂̌��ʂ�ێ�����
	std::vector<MGPosition> vecPos;
	MGBox box;
	if(IsBinaryFile(in)){
		// Load Binary File
		error = LoadBinary(in, vecPos);
	}else{
		// Load Ascii File
		error = LoadAscii(in, vecPos);
	}
	
	if(error){ // 0:����ɓǂݍ��݂����� ����ȊO�̒l:�ǂݍ��݂����s
		// �ǂݍ��݂����s�����ꍇ�A�G���[�R�[�h��ԋp
		return error;
	}

	// �g�������X��ݒ�
	double wc_zero_old = MGTolerance::set_wc_zero(m_box.length()*MGTolerance::rc_zero());
	// �ǂݍ��񂾍��W�l����m_indices, m_vecNormlTriang��ݒ�
	set_mesh_data(vecPos);
	MGTolerance::set_wc_zero(wc_zero_old);

	// �f�o�b�O�R�[�h
	//std::cout << "stl�t�@�C���ǂݍ���" << std::endl;
	//std::cout << "�O�p�`�� " << m_vecNormlTriang.size() << std::endl;
	//std::cout << "���_�� " << m_vecPos.size() << std::endl;
	//std::cout << "�{�b�N�X�g���W"<<  std::endl;
	//std::cout << "�ŏ��l " << std::endl;
	//std::cout << "X " << m_box.low()(0) << std::endl;
	//std::cout << "Y " << m_box.low()(1) << std::endl;
	//std::cout << "Z " << m_box.low()(2) << std::endl;
	//std::cout << "�ő�l " << std::endl;
	//std::cout << "X " << m_box.high()(0) << std::endl;
	//std::cout << "Y " << m_box.high()(1) << std::endl;
	//std::cout << "Z " << m_box.high()(2) << std::endl << std::endl;
	// �G���[�R�[�h��ԋp(����ɓǂݍ��݂�����)
	return error;
}

// �t�@�C������ǂݍ��񂾑S�Ă̒��_�̍��W�l���i�[����Ă���
// vecPos����e�O�p�`�̖@�����v�Z���Am_vecNormlTriang�ɒǉ�
// vecPos�̊e�v�f�̍��W�l�̏d�����g�������X�����ɔ��f��
// �d������菜���Am_vecPos�ɒǉ�
// vecPos�̊e�v�f�̃C���f�b�N�X��m_indices�ɒǉ�
// �������:m_vecNormlTriang, m_indices�ɒl���ݒ肳���
void MGStl::set_mesh_data(
	const std::vector<MGPosition>& vecPos // [in]:�t�@�C������ǂݍ��񂾍��W�l�̔z��
){
	// �t�@�C������ǂ񂾒��_������O�p�`�̐����擾
	int nTriang = vecPos.size()/3;
	// ���W�l�ƒ��_�̃C���f�b�N�X���ꎞ�I�ɕێ����Ă����}�b�v
	triangleMap VertexMap;
	for(int i = 0; i < nTriang; i++){
		// �O�p�`�̒��_�ԍ��̔z��̃C���f�b�N�X�����߂Ă���
		int i3=i*3;
		const MGPosition& p0=vecPos[i3++];
		const MGPosition& p1=vecPos[i3++];
		const MGPosition& p2=vecPos[i3];
		// �܂��A���͂��ꂽ3�_����ʂ̖@�������߁Am_vecNormlTriang��push_back����
		// ���ɁA�w�肵�����_�̍��W�����ɕۑ�����Ă��邩����
		// �ۑ�����Ă���ꍇ�́A�Y�����钸�_�̃C���f�b�N�X���󂯎��
		// �����m_indices��push_back����
		// �ۑ�����Ă��Ȃ��ꍇ�́A���_��m_vecPos��push_back��
		// ���ɃC���f�b�N�X��+1���A������󂯎��m_indices��push_back����
		push_back_triangle(p0, p1, p2, VertexMap);
	}
}

// Ascii�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����
// �܂��}�`�̃{�b�N�X�g��ݒ肷��
// �߂�l: =0�t�@�C���̓ǂݍ�������
//		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��
// �������:vecPos�ɑS�Ă̍��W�l���i�[����A box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���
// �t�@�C���X�g���[�����i��
int MGStl::LoadAscii(
	std::ifstream& in, // [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[��
	std::vector<MGPosition>& vecPos // [out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��
){
	// �����o�̒l�����Z�b�g
	Initialize();
	// Keep one line data read from file
	std::string strLine;
	// Keep cordinate of vertex
	MGPosition position1(3), position2(3), position3(3);

	if(in.eof()){ // �t�@�C������̏ꍇ
		return std::ios_base::eofbit;
	}

	// �t�@�C������facet normal�����邩�������ϐ�
	bool bStl = false;
	// �G���[�`�F�b�N�̌��ʂ��i�[����ϐ�
	int readState;
	// Read file and set date to member
	while(getline(in, strLine)){
		// �ǂݍ��݃G���[�`�F�b�N
		readState=in.rdstate();
		if(readState == std::ios::failbit || readState == std::ios::badbit){
			return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}

		if(strLine.find("facet normal") == -1){// �s����facet normal�����݂��邩
			continue;
		}
		// �t�@�C������facet normal�𔭌�
		bStl = true;
		// Triangle data
		std::string dummy;
		// outer loop�s��ǂݔ�΂�
		in.ignore(10000, '\n');
		// Get cordinate of vertex
		in >> dummy >> position1(0) >> position1(1) >> position1(2);
		in >> dummy >> position2(0) >> position2(1) >> position2(2);
		in >> dummy >> position3(0) >> position3(1) >> position3(2);

		// �ǂݍ��݃G���[�`�F�b�N
		readState = in.rdstate();
		if(readState){
			return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}

		// �t�@�C������ǂݍ��񂾍��W�l��������vector��push_back����
		vecPos.push_back(position1);
		vecPos.push_back(position2);
		vecPos.push_back(position3);
		// �{�b�N�X�g��3�����Ƃ��Aexpand
		m_box.expand(position1);
		m_box.expand(position2);
		m_box.expand(position3);
		// ���s����
		in.seekg(1, std::ios::cur);
		// endloop endfacet�s��ǂݔ�΂�
		in.ignore(10000, '\n');
		in.ignore(10000, '\n');
	}

	if(!bStl){
		// facet normal����x���ǂ܂��ɖ��s�܂œǂݍ��񂾏ꍇ
		// �܂�A�ʂ̎�ނ̃t�@�C����ǂ񂾏ꍇ�̏���
		return std::ios_base::eofbit;
	}

	// ����I���������l��ԋp
	return 0;
}

// Binary�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����
// �܂��}�`�̃{�b�N�X�g��ݒ肷��
// �߂�l: =0�t�@�C���̓ǂݍ�������
//		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��
// �������:vecPos�ɑS�Ă̍��W�l���i�[����A m_box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���
//			�t�@�C���X�g���[�����i��
int MGStl::LoadBinary(
	std::ifstream& in, // [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[��
	std::vector<MGPosition>& vecPos //[out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��
){
	// �����o�̒l�����Z�b�g
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

	// �������݃G���[�`�F�b�N
	int readState = in.rdstate();
	if(readState){
		return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
	}

	// Read and set cordinate of vertex and normal of triangle
	for(unsigned int i = 0; i < nTriang; i++){
		// normal���i�[����Ă��镔����ǂݔ�΂�
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

		// �t�@�C������ǂݍ��񂾍��W�l��������vector��push_back����
		vecPos.push_back(position1);
		vecPos.push_back(position2);
		vecPos.push_back(position3);
		// �{�b�N�X�g��expand
		m_box.expand(position1);
		m_box.expand(position2);
		m_box.expand(position3);

		if(i == nTriang - 1){
			// �Ō�̎O�p�`�̏ꍇ�͓ǂݍ��ݏI��
			break;
		}
		// ���̎O�p�`������ꍇ�A2�o�C�g�X�L�b�v����
		in.seekg(2, SEEK_CUR);

		// �ǂݍ��݃G���[�`�F�b�N
		readState = in.rdstate();
		if(readState){
			return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}
	}

	// �������݃G���[�`�F�b�N
	readState = in.rdstate();
	if(readState == std::ios::failbit || readState == std::ios::badbit){
		return readState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
	}

	// ����I���������l��ԋp
	return 0;
}

// �w�肳�ꂽ�p�X��Ascii�`����STL�t�@�C����ۑ�����
// �߂�l: =0�t�@�C���̏������݂�����
//		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
// ���O����:m_vecPos��m_vecNorml��m_indices�ɐ�������񂪓�����Ă���
// �������:rSTLFilePath�Ɏw�肵���p�X��Ascii�`����STL�t�@�C�����ۑ������
//�d�l�ύX�F�����ɃX�g���[����n�����悤�ɕύX
int MGStl::SaveAscii(
	std::ofstream& fout
	//const char* rSTLFilePath // [in]:�t�@�C���̕ۑ���̃p�X
)const{
	//std::ofstream fout(rSTLFilePath);

	// �������݃G���[�`�F�b�N
	int writeState = fout.rdstate();
	if(writeState){
		return writeState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
	}

	fout.setf(std::ios::scientific);
	// solid
	fout << "solid ascii" << std::endl;
	// �O�p�`�̐����擾(�O�p�`�̐� == �ʖ@���x�N�g���̌�)
	int nTriang = m_vecNormlTriang.size();
	// output normal of triangle and cordinate of vertex to file
	for(int i = 0; i < nTriang; i++){
		// Get triangle normal
		const MGUnit_vector& normal = m_vecNormlTriang[i];
		// facet normal
		fout << "    facet normal  " << normal(0) << " " << normal(1)<< " " << normal(2) << " " << std::endl;
		// outerloop
		fout << "        outer loop" << std::endl;
		
		// ���_�̔ԍ��̔z��̓Y�������v�Z���Ă���
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

		// �������݃G���[�`�F�b�N
		writeState = fout.rdstate();
		if(writeState){
			return writeState; // �G���[�R�[�h��ԋp(�ǂݍ��݂����s)
		}
	}
	// endsolid
	fout << "endsolid" << std::endl;
	fout.unsetf(std::ios::scientific);

	writeState = fout.rdstate();
	if(writeState == std::ios::failbit || writeState == std::ios::badbit){
		return writeState; // �G���[�R�[�h��ԋp(�������݂����s)
	}

	//fout.close();
	return 0;
}

// �w�肳�ꂽ�p�X��Binary�`����STL�t�@�C����ۑ�����
// �߂�l: =0�t�@�C���̏������݂�����
//		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
// ���O����:m_vecPos��m_vecNorml��m_indices�ɐ�������񂪓�����Ă���
// �������:rSTLFilePath�Ɏw�肵���p�X��Binary�`����STL�t�@�C�����ۑ������
int MGStl::SaveBinary(
	const char* rSTLFilePath  // [in]:�t�@�C���̕ۑ���̃p�X
)const{
	// File open
	std::ofstream fout(rSTLFilePath, std::ios::binary);
	// �X�g���[���̃G���[�`�F�b�N
	int writeState = fout.rdstate();
	if(writeState)
		return writeState;

	// �O�p�`�̐����擾(�O�p�`�̐� == �ʖ@���x�N�g���̌�)
	unsigned int nTriang = m_vecNormlTriang.size();
	// �t�@�C���擪80�o�C�g���X�L�b�v
	fout.seekp(80, std::ios::beg);
	// �t�@�C���ɎO�p�`�̖�������������
	fout.write((char*)&nTriang, 4);

	// �������݃G���[�`�F�b�N
	writeState = fout.rdstate();
	if(writeState){
		return writeState; // �G���[�R�[�h��ԋp(�������݂Ɏ��s)
	}

	// Keep vertex and normal of triangle temporary 
	float fNormal, fVertex;
	// �e�O�p�`�̏�����������
	for(unsigned int i = 0; i < nTriang; i++){
		// �O�p�`�̖ʖ@���x�N�g������������
		const MGUnit_vector& normal = m_vecNormlTriang[i];
		for(int j = 0; j < 3; j++){
			fNormal = float(normal[j]);
			fout.write((char*)&fNormal, 4);
		}

		// �O�p�`�̒��_�ԍ��̔z��̃C���f�b�N�X�����߂Ă���
		int i3 = i*3;

		// �x�N�^�[���璸�_�̃C���f�b�N�X���w���Ă�����W�l�����o��
		// �O�p�`�̒��_�̍��W����������
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

		// �Ō�̎O�p�`�̏ꍇ�͏������݂��I��
		if(i == nTriang - 1){
			break;
		}

		// �O�p�`���Ƃ�2�o�C�g�Ԋu���󂯂�
		fout.seekp(2, std::ios::cur);

		// �������݃G���[�`�F�b�N
		writeState = fout.rdstate();
		if(writeState){
			return writeState; // �G���[�R�[�h��ԋp(�������݂Ɏ��s)
		}
	}
	// �������݃G���[�`�F�b�N
	writeState = fout.rdstate();
	if(writeState == std::ios::failbit || writeState == std::ios::badbit){
		return writeState;// �G���[�R�[�h��ԋp(�������݂����s)
	}

	fout.close();
	return 0;
}

// ���͂ł���position��VertexMap�ɂ��łɓo�^����Ă���΂���m_vecPos�Y������Ԃ�
// ���o�^�̏ꍇ�A�V����m_vecPos�Ɋi�[���Aposition��m_vecPos�Y������map��VertexMap��
// �o�^����
int MGStl::IdentifyPosition(
	const MGPosition& position, // ���_�̍��W
	triangleMap& VertexMap // �d���̂Ȃ����_�̍��W�A���_�̃C���f�b�N�X��ێ�����
){
	int nVertexIndex = m_vecPos.size();
	//iterator��position�̍��W�l�ȏ�ƂȂ�v�f���͂��߂Č����ꏊ���w���悤�ɂ���
	triangleMap::iterator itAft = VertexMap.lower_bound(position), itPre;
	if(itAft != VertexMap.end()){	//iterator���I�����w���Ă��Ȃ��ꍇ
		if(itAft->first==position)
			return itAft->second;	//�w���Ă���v�f�ƈ�v�����ꍇ�C����index��Ԃ�
	}
	if(itAft != VertexMap.begin()){	//iterator���n�߂��w���Ă��Ȃ��ꍇ
		itPre=itAft;itPre--;	//��O�̗v�f��position���r����
		if(itPre->first==position){
			return itPre->second;	//��O�̗v�f�ƈ�v�����ꍇ�C����index��Ԃ�
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

// �����o�f�[�^���������ފ֐�
void MGStl::WriteMembers(MGOfstream& buf)const{
	MGObject::WriteMembers(buf);
	// �O�p�`�̖@���x�N�g������������
	int nTriangle = m_vecNormlTriang.size();
	buf << nTriangle;
	for(int i = 0; i < nTriangle; i++){
		m_vecNormlTriang[i].dump(buf);
	}

	// �O�p�`�̏d�����Ȃ����W�̒��_����������
	int nPos = m_vecPos.size();
	buf << nPos;
	for(int k = 0; k < nPos; k++){
		m_vecPos[k].dump(buf);
	}

	// �O�p�`�̊e���_�̃C���f�b�N�X����������
	int nIndex = m_indices.size();
	buf << nIndex;
	for(int l = 0; l < nIndex; l++){
		buf << m_indices[l];
	}
}

// �����o�f�[�^��ǂݍ��ފ֐�
void MGStl::ReadMembers(MGIfstream& buf){
	MGObject::ReadMembers(buf);
	// �O�p�`�̖@���x�N�g����ǂݍ���
	int nTriangle;
	buf >> nTriangle;
	m_vecNormlTriang.resize(nTriangle);
	for(int i = 0; i < nTriangle; i++){	
		m_vecNormlTriang[i].restore(buf);
	}

	// �O�p�`�̊e���_�̍��W���d�������ǂݍ���
	int nPos;
	buf >> nPos;
	m_vecPos.resize(nPos);
	for(int k = 0; k < nPos; k++){
		m_vecPos[k].restore(buf);
	}

	// �O�p�`�̊e���_�̃C���f�b�N�X��ǂݍ���
	int nIndex;
	buf >> nIndex;
	m_indices.resize(nIndex);
	for(int l = 0; l < nIndex; l++){
		buf >> m_indices[l];
	}
}

// ��x�t�@�C����ǂ݁A�����o�ɒl���ݒ肳�ꂽ��ԂōĂуt�@�C����ǂނ�
// �����o�̒l���㏑�������B���������邽�߁A�����o�̒l��S�ă��Z�b�g����
void MGStl::Initialize(){
	m_box.set_null();
	m_vecPos.clear();
	m_vecNormlTriang.clear();
	m_indices.clear();
}

// 3�_����쐬�����O�p�`�̏��������o�ϐ��ɒl��ݒ肷��
// ���O����:���͂���pos1, pos2, pos3�͔����v���ɂȂ��Ă��邱��
// �������:m_vecPos, m_vecNormlTriang, m_indices�ɎO�p�`�̏�񂪐ݒ肳���
void MGStl::push_back_triangle(
	const MGPosition& pos1,
	const MGPosition& pos2,
	const MGPosition& pos3,
	triangleMap& VertexMap
){
	// 3�_����ʂ̖@�������߁A�����o�ϐ���push_back����
	m_vecNormlTriang.push_back(CalcNormal(pos1, pos2, pos3));
	// m_indices�ɒ��_�̃C���f�b�N�X��push_back��
	// m_vecPos�ɏd���Ȃ�MGPosition��push_back����
	m_indices.push_back(IdentifyPosition(pos1, VertexMap));
	m_indices.push_back(IdentifyPosition(pos2, VertexMap));
	m_indices.push_back(IdentifyPosition(pos3, VertexMap));
}

//TLData��MGStl�ɒǉ�����
void MGStl::AddTLData(const mgTLData& tlData, triangleMap& VertexMap){
	// �T�[�t�F�C�X���擾
	const MGSurface& surf = tlData.surface();
	const mgTLparameter& tlp = tlData.tlparam();
	double errorSave=MGTolerance::set_wc_zero(tlp.get_surface_error());

	// �S�Ă̒��_��uv�l���擾�AmgTLData�̑S���_��uv�l��MGPositon�ɕϊ����Am_vecPos�Ɋi�[
	//�Ή��\indexVec���쐬����B���_ID�̑Ή������߂�vector:indexVec[��ID]=�VID�őΉ����Ƃ�D
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

	// �S�Ẵ|���S�����擾(�t�@���������̓X�g���b�v)
	const mgTLTriangles& triangles = tlData.triangles();
	// �|���S���z���begin��end���擾
	mgTLTriangles::const_triIterator j = triangles.begin(), je = triangles.end();
	for(;j!=je;++j){
	// �|���S�����ƂɃ��[�v(indices���擾����)
		const mgTLTriangle& tri=**j;// �|���S�����擾
		size_t nVert = tri.size();// ���_�̌����擾����
		if(nVert<3)
			continue;

		
		mgTESTRIANG geoType = tri.getGeometryType();//�|���S���̌`��(�t�@���������̓X�g���b�v)
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

			// �O�p�`�̒��_uv�l����@���x�N�g�������߂�
			MGPosition uv=(tlpoints[id0]+tlpoints[id1]+tlpoints[id2])/3.;
			m_vecNormlTriang.push_back(surf.unit_normal(uv));
		}
	}
	MGTolerance::set_wc_zero(errorSave);
}

// �����o�[�f�[�^�𒼐ڃZ�b�g����
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
