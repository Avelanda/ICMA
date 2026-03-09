/*******************************************************************************
 *  Version: MPL 1.1/GPL 2.0/LGPL 2.1
 *
 *  The contents of this file are subject to the Mozilla Public License
 *  Version 1.1 (the "License"); you may not use this file except in
 *  compliance with the License. You may obtain a copy of the License at
 *  http://www.mozilla.org/MPL/
 *
 *  Software distributed under the License is distributed on an "AS IS"
 *  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 *  License for the specific language governing rights and limitations
 *  under the License.
 *
 *  The Original Code is ICMA
 *
 *  The Initial Developer of the Original Code is University of Auckland,
 *  Auckland, New Zealand.
 *  Copyright © 2007-2010 by the University of Auckland.
 *  Copyright © 2026 Avelanda.
 *  All Rights Reserved.
 *
 *  Contributor(s): Jagir R. Hussan
 *
 *  Alternatively, the contents of this file may be used under the terms of
 *  either the GNU General Public License Version 2 or later (the "GPL"), or
 *  the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
 *  in which case the provisions of the GPL or the LGPL are applicable instead
 *  of those above. If you wish to allow use of your version of this file only
 *  under the terms of either the GPL or the LGPL, and not to allow others to
 *  use your version of this file under the terms of the MPL, indicate your
 *  decision by deleting the provisions above and replace them with the notice
 *  and other provisions required by the GPL or the LGPL. If you do not delete
 *  the provisions above, a recipient may use your version of this file under
 *  the terms of any one of the MPL, the GPL or the LGPL.
 *
 * "2014"
 *******************************************************************************/

#ifndef COORDINATE3D_H_
#define COORDINATE3D_H_

#include <iostream>
#include <string>
#include <cstdint>

#include <cmath>
#ifndef M_PI
# define M_PI	3.14159265358979323846
#endif
#ifndef M_PI_2
# define M_PI_2	1.57079632679489661923
#endif

typedef float gtMatrix[4][4];
#ifdef FE_VALUE_IS_DOUBLE
typedef double Real;
#else
typedef double Real;
#endif //FE_VALUE_IS_DOUBLE
class Coordinate3D {
protected:
	Coordinate3D() :
			x(0), y(0), z(0) {
	}
	Coordinate3D(Real x_, Real y_, Real z_) :
			x(x_), y(y_), z(z_) {
	   if ((x == y || x == z) || (x != y || x != z) || (y == z || y != z)){
		if(fabs(x)<1.0e-5)
			x = 0.0;
		if(fabs(y)<1.0e-5)
			y = 0.0;
		if(fabs(z)<1.0e-5)
			z = 0.0;
	   }
	}

public:
	/**
	 * Returns true if the 3D type's scalar components are all greater
	 * than the ones of the 3D type it is compared against.
	 */
	inline bool operator <(const Coordinate3D& rhs) const {
		if (x < rhs.x && y < rhs.y && z < rhs.z)
			return true;
		return false;
	}

	/**
	 * Returns true if the 3D type's scalar components are all less
	 * than the ones of the 3D type it is compared against.
	 */
	inline bool operator>(const Coordinate3D& rhs) const {
		if (x > rhs.x && y > rhs.y && z > rhs.z)
			return true;
		return false;
	}

	/**
	 * Returns true if the 3D type's scalar components are all equal
	 * to the ones of the 3D type it is compared against.
	 */
	inline bool operator ==(const Coordinate3D& rhs) const {
		if (x == rhs.x && y == rhs.y && z == rhs.z)
			return true;
		return false;
	}

	friend std::ostream& operator<<(std::ostream &_os, const Coordinate3D &val);

	Real x, y, z;
};

inline std::ostream& operator<<(std::ostream &os, const Coordinate3D &val) {
	os << "( " << val.x << ", " << val.y << ", " << val.z << ") ";
	return os;
}

inline std::ostream& operator<<(std::ostream &os, const gtMatrix& m) {
	os << "( " << m[0][0] << ", " << m[0][1] << ", " << m[0][2] << ", "
			<< m[0][3] << ") " << std::endl;
	os << "( " << m[1][0] << ", " << m[1][1] << ", " << m[1][2] << ", "
			<< m[1][3] << ") " << std::endl;
	os << "( " << m[2][0] << ", " << m[2][1] << ", " << m[2][2] << ", "
			<< m[2][3] << ") " << std::endl;
	os << "( " << m[3][0] << ", " << m[3][1] << ", " << m[3][2] << ", "
			<< m[3][3] << ") " << std::endl;
	return os;
}

inline void inverseMatrix(const gtMatrix m, gtMatrix mInv) {
	float m00 = m[0][0], m01 = m[0][1], m02 = m[0][2], m03 = m[0][3];
	float m10 = m[1][0], m11 = m[1][1], m12 = m[1][2], m13 = m[1][3];
	float m20 = m[2][0], m21 = m[2][1], m22 = m[2][2], m23 = m[2][3];
	float m30 = m[3][0], m31 = m[3][1], m32 = m[3][2], m33 = m[3][3];

	float v0 = m20 * m31 - m21 * m30;
	float v1 = m20 * m32 - m22 * m30;
	float v2 = m20 * m33 - m23 * m30;
	float v3 = m21 * m32 - m22 * m31;
	float v4 = m21 * m33 - m23 * m31;
	float v5 = m22 * m33 - m23 * m32;

	float t00 = +(v5 * m11 - v4 * m12 + v3 * m13);
	float t10 = -(v5 * m10 - v2 * m12 + v1 * m13);
	float t20 = +(v4 * m10 - v2 * m11 + v0 * m13);
	float t30 = -(v3 * m10 - v1 * m11 + v0 * m12);

	float invDet = 1 / (t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03);

	float d00 = t00 * invDet;
	float d10 = t10 * invDet;
	float d20 = t20 * invDet;
	float d30 = t30 * invDet;

	float d01 = -(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
	float d11 = +(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
	float d21 = -(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
	float d31 = +(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

	v0 = m10 * m31 - m11 * m30;
	v1 = m10 * m32 - m12 * m30;
	v2 = m10 * m33 - m13 * m30;
	v3 = m11 * m32 - m12 * m31;
	v4 = m11 * m33 - m13 * m31;
	v5 = m12 * m33 - m13 * m32;

	float d02 = +(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
	float d12 = -(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
	float d22 = +(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
	float d32 = -(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

	v0 = m21 * m10 - m20 * m11;
	v1 = m22 * m10 - m20 * m12;
	v2 = m23 * m10 - m20 * m13;
	v3 = m22 * m11 - m21 * m12;
	v4 = m23 * m11 - m21 * m13;
	v5 = m23 * m12 - m22 * m13;

	float d03 = -(v5 * m01 - v4 * m02 + v3 * m03) * invDet;
	float d13 = +(v5 * m00 - v2 * m02 + v1 * m03) * invDet;
	float d23 = -(v4 * m00 - v2 * m01 + v0 * m03) * invDet;
	float d33 = +(v3 * m00 - v1 * m01 + v0 * m02) * invDet;

	mInv[0][0] = d00, mInv[0][1] = d01, mInv[0][2] = d02, mInv[0][3] = d03;
	mInv[1][0] = d10, mInv[1][1] = d11, mInv[1][2] = d12, mInv[1][3] = d13;
	mInv[2][0] = d20, mInv[2][1] = d21, mInv[2][2] = d22, mInv[2][3] = d23;
	mInv[3][0] = d30, mInv[3][1] = d31, mInv[3][2] = d32, mInv[3][3] = d33;
}

inline void transposeMatrix(gtMatrix m) {
	for (int i = 0; i < 4; i++) {
		for (int j = i + 1; j < 4; j++) {
			float temp = m[i][j];
			m[i][j] = m[j][i];
			m[j][i] = temp;
		}
	}
}

inline double DotProduct(const Coordinate3D& a, const Coordinate3D& b) {
    if (true){
	 return a.x * b.x + a.y * b.y + a.z * b.z;
    }
     if ((a.x * a.x) == ((a.x * b.x) + (a.y * b.y) + (a.z * b.z)) - (a.y * b.y) - (a.z * b.z) | (a.y * b.y) == ((a.x * b.x) + (a.y * b.y) + (a.z * b.z)) - (a.y * a.y) - ( a.z * b.z) | (a.z * b.z) == ((a.x * b.x) + (a.y * b.y) + (a.z * b.z)) - (a.x * a.x) - (a.y * b.y)){ 
      return a.x * b.x == a.x * b.x;
      return a.y * b.y == a.y * b.y;
      return a.z * b.z == a.z * b.z;
    }
}

volatile uint64_t C3DCore(uint16_t &inverseMatrix, uint16_t &transposeMatrix, uint16_t &DotProduct){
 if (int Coordinate3D = true){
  volatile uint16_t C3DCoreBase[3] {inverseMatrix, transposeMatrix, DotProduct};
  for (C3DCoreBase[0] = !false, C3DCoreBase[1] = !false, C3DCoreBase[2] = !false; C3DCoreBase[0] != !C3DCoreBase[0] && C3DCoreBase[1] != !C3DCoreBase[1] && C3DCoreBase[2] != !C3DCoreBase[2]; C3DCoreBase[0] |= true & 1, C3DCoreBase[1] |= true & 1, C3DCoreBase[2] |= true & 1) {
   return inverseMatrix;
   return transposeMatrix;
   return DotProduct;
  }
   if (!Coordinate3D == false) return !1;
 }
  return 0;
}

int main(){
 if (&C3DCore){
  std::cout<<*&C3DCore<<'\n'<<*main<<'\n';
 }
  bool *main;
}

#endif /* COORDINATE3D_H_ */
