/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#ifndef _MGStl_HH_
#define _MGStl_HH_

#include <map>
#include "mg/object.h"
#include "mg/Curve.h"
#include "mg/Surface.h"
#include "mg/FSurface.h"
#include "topo/Face.h"
#include "topo/Shell.h"
#include "Tl/TLData.h"

class mgTLData;

/** @addtogroup MGObjectRelated
 *  @{
 */

//Define MGStl Class.

///MGStl is a concrete class which represents an STL information.
///���̃N���X�͕����̎O�p�`��z��Ƃ��ĕێ�����.
///�S�Ă̎O�p�`�̒��_�̍��W�̓t�@�C������ǂݍ��񂾏��ŁA�d������菜���Ĕz��Ɋi�[�����.
///m_indices�ɂ̓t�@�C������ǂݍ��񂾏��ԂŊe�O�p�`�̒��_�̕��т��i�[����Ă���
///���̊e�v�f�ɂ͊Y�����钸�_���W�̔z��̓Y�������i�[�����.
///�Ⴆ�΁Ai�Ԗڂ̎O�p�`�̊e���_�̍��W�̓Y������indices[i*3]�Aindices[i*3+1]�Aindices[i*3+2]
///�Ƃ����ȏ�̏����Ŋi�[����Ă���.
class MGCLASS MGStl : public MGObject{

	/// map�ɂ�1�̃L�[�ɑ΂��ĕ����̒l�͎��߂��Ȃ�����
	/// �O�p�`���\�����钸�_��ID�����߂��\���̂𗘗p����
	class vertId{
	public:
		int id1;
		int id2;
		int id3;
	};

	/// �֐��I�u�W�F�N�g���i�[����\����
	class positionComp{
	public:
		bool operator()(const MGPosition& p1, const MGPosition& p2) const{
			/// lexicographical_compare style
			if(p1(0) < p2(0)){
				return true;
			}else if(p2(0) < p1(0)){
				return false;
			}
			if(p1(1) < p2(1)){
				return true;
			}else if(p2(1) < p1(1)){
				return false;
			}
			if(p1(2) < p2(2)){
				return true;
			}else if(p2(2) < p1(2)){
				return false;
			}
			/// p1 == p2
			return false;
		}		
	};

	/// �O�p�`�̒��_�̍��W�A���_�̔ԍ����i�[����map�̕ʖ����`.
	typedef std::map<MGPosition, int, positionComp> triangleMap;

	/// ���_�̃C���f�b�N�X�A���_�̍��W���i�[����map�̕ʖ����`.
	typedef std::map<int, MGPosition> IndexPosMap;

	/// �O�p�`�̃C���f�b�N�X�A3�̒��_�̃C���f�b�N�X���i�[����map�̕ʖ����`.
	typedef std::map<int, vertId> TriVertMap;

public:
	/// �f�t�H���g�R���X�g���N�^
	MGStl(void){;};

	/// �R�s�[�R���X�g���N�^
	MGStl(const MGStl& stl);

	///conversion constructor from tessellation data.
	MGStl(const mgTLDataVector& tlDataVector);
	MGStl(const mgTLData& tlData);

	///Constructor from triangle data, index+vertices.
	///This constructor uses the wc_zero to identify
	///different two positions as the same input position.
	MGStl(
		int nTriang, /// �O�p�`�̐�
		const int* triang, /// trang[i], [i+1], [i+2] make a triangle for i=0,...,nTriang-1
		const double* verts /// ���_�̍��W
	);

	/// ���z�f�X�g���N�^.
	~MGStl(void);

	/// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵��MGVector�̒l�����Z.
	MGStl& operator+=(const MGVector& v);

	/// �S�Ă̒��_�ƃ{�b�N�X�̍��W����w�肵��MGVector�̒l������.
	MGStl& operator-=(const MGVector& v);

	/// �S�Ă̒��_�ƃ{�b�N�X�̍��W�Ɏw�肵���l���|����.
	MGStl& operator*=(double scale);

	/// �^����ꂽ�ϊ����s��.
	MGStl& operator*=(const MGMatrix& mat);

	/// �^����ꂽ�ϊ��ɂ��g�����X�t�H�[�����s��.
	MGStl& operator*=(const MGTransf& tr);

	/// 2��MGStl�I�u�W�F�N�g�������������肷��.
	bool operator==(const MGStl& stl);

	/// 2��MGStl�I�u�W�F�N�g��������.
	MGStl& operator=(const MGStl& stl);

	std::ostream& out(std::ostream& ostrm)const;

	/// STL�t�@�C���̐}�`��bounding box���擾����.
	/// �߂�l�FSTL�t�@�C���̐}�`��bounding box
	const MGBox& box()const{return m_box;};

	///Construct new object by copying to newed area.
	///User must delete this copied object by "delete".
	MGStl* clone()const;

	/// Return This object's typeID.
	long identify_type()const;

	///Draw 3D curve in world coordinates.
	///The object is converted to curve(s) and is drawn.
	void drawWire(
		double span_length,	///<Line segment span length.
		int line_density=1	///<line density to draw a surface in wire mode.
	)const;

	///Draw 3D point(vertex) in world coordinates.
	///The object is converted to point(s) and is drawn.
	///This is valid only for topology objects or MGPoint.
	void draw3DVertex()const;

	///Shade the object in world coordinates.
	void shade(
		double span_length	///<Line segment span length.
	)const;

	///Compute the intersections of two objects.
	///***********MGStl does not support the intersection functions.******
	MGisects intersection(const MGObject& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGCurve& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGFSurface& obj2)const{return MGisects(this,obj2.object_pointer());};
	MGisects intersection(const MGSurface& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGFace& obj2)const{return MGisects(this,&obj2);};
	MGisects intersection(const MGShell& obj2)const{return MGisects(this,&obj2);};

	///Get manifold dimension.
	unsigned manifold_dimension()const{return 2;};

	///Write all member data.
	void WriteMembers(MGOfstream& buf)const;

	///Read all member data.
	void ReadMembers(MGIfstream& buf);

	///return the ref to m_vecPos.
	const std::vector<MGPosition>& positions()const{return m_vecPos;};
	std::vector<MGPosition>& positions(){return m_vecPos;};

	///return the ref to m_vecPos.
	const std::vector<MGUnit_vector>& normals()const{return m_vecNormlTriang;};
	std::vector<MGUnit_vector>& normals(){return m_vecNormlTriang;};

	/// �O�p�`�̌����擾����.
	/// �߂�l�F�O�p�`�̌�.
	int GetTriangleCount()const{return m_vecNormlTriang.size();};

	/// �w�肵���O�p�`�̒��_���W��������Ă���z��̃C���f�b�N�X���擾����.
	void GetVertIndices(
		int i, ///< [in]�F�O�p�`�̃C���f�b�N�X(0 <= i < GetTriangleCount)
		size_t pos[3] ///< [out]�F�w�肵���O�p�`�̊e���_���W�̔z��̓Y����[i, j, k]
	)const;

	/// �����Ŏw�肵���p�X��STL�t�@�C����ǂݍ��݃����o�ɒl��ݒ肷��.
	/// �܂�MGTolerance::wc_zero�ɒl��ݒ肷��.
	/// �߂�l: =0�t�@�C���̓ǂݍ�������
	///		  !=0 �ǂݍ��܂Ȃ������B�܂��͎��s����(std::ifstream�̃G���[�R�[�h�j
	/// �������:m_vecPos, m_vecNorml, m_indices, m_box�ɐ}�`�̏�񂪊i�[�����.
	///			MGTolerance::wc_zero�ɒl���ݒ肳���.
	int LoadFile(
		const char* strFilePath ///< [in]:�ǂݍ���STL�t�@�C���ւ̃p�X
	);

	/// �w�肳�ꂽ�p�X��Ascii�`����STL�t�@�C����ۑ�����.
	/// �߂�l: =0�t�@�C���̏������݂�����
	///		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
	/// �������:rSTLFilePath�Ɏw�肵���p�X��Ascii�`����STL�t�@�C�����ۑ������.
	int SaveAscii(
		std::ofstream& fout
	)const;

	/// �w�肳�ꂽ�p�X��Binary�`����STL�t�@�C����ۑ�����.
	/// �߂�l: =0�t�@�C���̏������݂�����
	///		  !=0 �������܂Ȃ������B�܂��͎��s����(std::ofstream�̃G���[�R�[�h)
	/// �������:rSTLFilePath�Ɏw�肵���p�X��Binary�`����STL�t�@�C�����ۑ������.
	///�d�l�ύX�F�����ɃX�g���[����n�����悤�ɕύX09/09/17.
	int SaveBinary(
		const char* rSTLFilePath ///< [in]:�t�@�C���̕ۑ���̃p�X
	)const;

	/// 3�_����쐬�����O�p�`�̏��������o�ϐ��ɒl��ݒ肷��.
	/// ���O����:���͂���pos1, pos2, pos3�͔����v���ɂȂ��Ă��邱��.
	/// �������:m_vecPos, m_vecNormlTriang, m_indices�ɎO�p�`�̏�񂪐ݒ肳���.
	void push_back_triangle(
		const MGPosition& pos1,
		const MGPosition& pos2,
		const MGPosition& pos3,
		triangleMap& VertexMap
	);

	/// �O�p�`���Ƃɖ@����\������.
	void display_arrows()const;

	/// �����o�[�f�[�^�𒼐ڃZ�b�g����.
	void set_all_data(
		const std::vector<MGPosition>& vertices,
		const std::vector<int>& indices);

	/// Update the bounds.
	void update_bounds();

	/// Update the all normals of the triangles.
	void update_normals();

	/// �I�u�W�F�N�g����ԋp����.
	std::string whoami()const{return "Stl";};

private:
	/// ��x�t�@�C����ǂ݁A�����o�ɒl���ݒ肳�ꂽ��ԂōĂуt�@�C����ǂނ�.
	/// �����o�̒l���㏑�������B���������邽�߁A�����o�̒l��S�ă��Z�b�g����.
	void Initialize();

	/// ���W�l��m_vecPos��push_back���Am_box�l��expand()����.
	/// �������:m_vecPos�ɗv�f��1�ǉ������.
	void push_back_Position(
		const MGPosition& pos ///< [in]:m_vecPos�ɒǉ����钸�_�̍��W.
		);

	/// �O�p�`���\�����钸�_ID��m_indices��push_back��
	/// �w�肵���O�p�`�ɂ��Ėʖ@���x�N�g�����v�Z��m_vecNormlTriang��push_back����.
	/// ���O����:m_vecPos�ɒl���i�[����Ă���.
	/// �������:�P�̎O�p�`�̏��m_vecNormlTriang�Am_indeces�Ɋi�[�����.
	void add_indices_and_calc_normal(
	 const int vertId[3], ///< [in]:�P�̎O�p�`�̒��_ID�̔z��.
	 int triIndex ///< [in]:�@�������߂�Ώۂ̎O�p�`�̃C���f�b�N�X
				  ///< (0 <= triIndex < GetTriangleCount).
	 );

	/// Ascii�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����.
	/// �܂��}�`�̃{�b�N�X�g��ݒ肷��.
	/// �߂�l: =0�t�@�C���̓ǂݍ�������
	///		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
	/// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��.
	/// �������:vecPos�ɑS�Ă̍��W�l���i�[����A box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���.
	/// �t�@�C���X�g���[�����i��.
	int LoadAscii(
		std::ifstream& in, ///< [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[�� .
		std::vector<MGPosition>& vecPos ///<[out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��.
		);

	/// Binary�`���̃t�@�C����ǂݍ��ݑS�Ă̍��W�l���擾����.
	/// �܂��}�`�̃{�b�N�X�g��ݒ肷��.
	/// �߂�l: =0�t�@�C���̓ǂݍ�������
	///		  !=0 �ǂݍ��݂Ɏ��s����(std::ifstream�̃G���[�R�[�h)
	/// ���O����:���ɃI�[�v�����ꂽ�t�@�C���X�g���[����n��.
	/// �������:vecPos�ɑS�Ă̍��W�l���i�[����A m_box�Ƀ{�b�N�X�g�̍��W�l���ݒ肳���.
	///			�t�@�C���X�g���[�����i��
	int LoadBinary(
		std::ifstream& in, ///< [in/out]:���ɃI�[�v�����ꂽ�t�@�C���X�g���[��.
		std::vector<MGPosition>& vecPos ///<[out]:�t�@�C������ǂݍ��񂾍��W�l�̔z��.
		);
	
	/// �t�@�C������ǂݍ��񂾑S�Ă̒��_�̍��W�l���i�[����Ă���.
	/// vecPos����e�O�p�`�̖@�����v�Z���Am_vecNormlTriang�ɒǉ�.
	/// vecPos�̊e�v�f�̍��W�l�̏d�����g�������X�����ɔ��f��
	/// �d������菜���Am_vecPos�ɒǉ�.
	/// vecPos�̊e�v�f�̃C���f�b�N�X��m_indices�ɒǉ�.
	/// �������:m_vecNormlTriang, m_indices�ɒl���ݒ肳���.
	void set_mesh_data(
		const std::vector<MGPosition>& vecPos ///< [in]:�t�@�C������ǂݍ��񂾍��W�l�̔z��.
	);

	/// ���͂ł���position��VertexMap�ɂ��łɓo�^����Ă��邩�`�F�b�N����.
	/// �o�^����Ă���΂���m_vecPos�Y������Ԃ�,
	/// ���o�^�̏ꍇ�A�V����m_vecPos�Ɋi�[���Aposition��m_vecPos�Y������map��VertexMap��
	/// �o�^����.
	int IdentifyPosition(
		const MGPosition& position, ///< [in]:���_�̍��W.
		triangleMap& VertexMap ///< [in/out]:���_�̍��W�A���_�̃C���f�b�N�X��ێ�����.
	);

	///TLData�̓ǂݍ��݂��s��.
	void AddTLData(const mgTLData& tlData, triangleMap& VertexMap);

	//////////

	/// STL�t�@�C���̐}�`��bounding box���i�[����.
	MGBox m_box;

	/// STL�t�@�C���̐}�`���\������e�O�p�`�́A���W�̏d�����Ȃ����_���W�̔z��.
	/// �t�@�C������ǂݍ��񂾏��ŁA���W�̏d������菜���A���W�l���i�[����Ă���.
	std::vector<MGPosition> m_vecPos;

	/// STL�t�@�C���̐}�`���\������e�O�p�`�̖@���x�N�g���̔z��.
	/// �t�@�C������ǂݍ��܂ꂽ���ŎO�p�`�̖@���x�N�g�����i�[����Ă���.
	///m_vecNormlTriang.size()*3=m_indices.size().
	///m_vecNormlTriang[i] is the normal of the triangle  m_indices[i], [i+1], [i+2] for
	///i=0, ..., m_indices.size()/3.
	std::vector<MGUnit_vector> m_vecNormlTriang;

	/// STL�t�@�C���̐}�`���\������e�O�p�`�̊e���_�ɑΉ�����
	/// ���W�̔z��̃C���f�b�N�X���i�[����z��.
	/// �t�@�C������ǂݍ��񂾏��ԂŊe�O�p�`�̒��_�̕��т��i�[����Ă���.
	/// �e�v�f�ɂ͊Y�����钸�_���W�̔z��̓Y�������i�[����Ă���.
	/// ��Fi�Ԗڂ̎O�p�`�̊e���_�̍��W�̔z��̓Y������
	/// m_indices[i*3]�A[i*3+1]�A[i*3+2]�Ƃ����ȏ�̏����Ŏ擾�ł���.
	std::vector<int> m_indices;
};

/** @} */ // end of MGObjectRelated group

#endif
