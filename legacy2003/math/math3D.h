// ============================================================================
// 
// Definition mathematischer Funktionen
//
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 9 $
// 
// ============================================================================
#ifndef MATH3D_H
#define MATH3D_H

// C++ Includes

#include <stdio.h>
#include <math.h>
#include <set>
#include <vector>
#include <map>

// Local Includes
#include "vec.h"
#include "../def.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define EPSILON  1e-8
#define EPSILONF 1e-8f

bool ismatherr();
void setmatherrflag(bool flag);



/*****************************************************************************
 * Hilfsfunktionen
 *****************************************************************************/
#define MAX(x,y) ((x)>(y) ? (x):(y))
#define MIN(x,y) ((x)<(y) ? (x):(y))
#define SGN(x) ((x) > 0 ? 1 : ((x) < 0 ? -1 : 0))
#define SQR(x) ((x)*(x))

/*****************************************************************************
 swap two floating point values
 *****************************************************************************/
inline void swap(float **a, float **b)
{
	register float *tmp;
	tmp = *a;
	*a = *b;
	*b = tmp;
}

template <class T> inline void swap(T &a, T&b)
{
   T t;
   t = a;
   a = b;
   b = t;
}
/*
template <class T> inline void swap(T *a, T *b)
{
   T t;
   t = *a;
   *a = *b;
   *b = t;
}
*/
// define a complex number
class Complex
{
public:
   float r;
   float i;

   /***************************************************************************
    * construct a default Complex number
    **************************************************************************/
   Complex() : r(0.0f), i(0.0f)
   {
   }

   /***************************************************************************
    * construct a Complex number
    **************************************************************************/
   Complex(float r, float i)
   {
      this->r = r;
      this->i = i;
   }

   /***************************************************************************
    * define division by a real number
    **************************************************************************/
	inline Complex operator /(const float s)
	{
		return Complex((r / s), (i / s));
	}

   /***************************************************************************
    * define division by a real number
    **************************************************************************/
	inline void operator /=(const float s)
	{
		r = r / s;
      i = i / s;
	}

   /***************************************************************************
    * return the absolut value of the complex number
    **************************************************************************/
   float absolut()
   {
      return sqrtf(r*r + i*i);
   }
};

// define a texture coordinate
typedef struct
{
	float u;
	float v;
}
Tex2;


// define a edge
typedef struct Edge
{
	int a;
	int b;
} 
Edge;



// define a color value
typedef struct rgb
{
	unsigned char r, g, b;

   /***************************************************************************
    * initialization
    **************************************************************************/
   inline void init(unsigned char r, unsigned char g, unsigned char b)
   {
      this->r = r;
      this->g = g;
      this->b = b;
   }

   /***************************************************************************
    * basic math operations for rgb values
    **************************************************************************/
   inline void operator +=(rgb &c)
   {
      r += c.r;
      g += c.g;
      b += c.b;
   }

   /***************************************************************************
    * equal
    **************************************************************************/
   inline bool operator ==(rgb &c)
   {
      return r == c.r && g == c.g && b == c.b;
   }

   /***************************************************************************
    * equal
    **************************************************************************/
   inline bool operator !=(rgb &c)
   {
      return r != c.r || g != c.g || b != c.b;
   }
} rgb;

typedef struct rgbf
{
	float r, g, b;
} rgbf;

typedef struct rgba
{
	unsigned char r, g, b, a;
} rgba;

typedef struct rgbaf
{
	float r, g, b, a;
} rgbaf;

// define a face
typedef struct Face
{
	int a;
	int b;
	int c;
} 
Face;

/******************************************************************************
 *
 * Typedefs für Container und Iteratoren
 *
 *****************************************************************************/

// integers
typedef std::set<int> int_set;
typedef std::vector<int> int_array;
typedef std::set<int>::iterator int_iterator;
typedef std::vector<int> int_vector;


/*****************************************************************************
 3x3 Matrix definitions

 access elements with mat.m[row][column]
 *****************************************************************************/
struct mat3
{
	float m[3][3];

	/***************************************************************************
	 * matrix multiplication
	 **************************************************************************/
	inline mat3 operator *(const mat3 m2) const
	{
		mat3  result;
		float accu;

		for (int i=0; i<3; i++) {
			for (int j=0; j<3; j++) {
				accu = m[i][0] * m2.m[0][j];
				accu+= m[i][1] * m2.m[1][j];
				accu+= m[i][2] * m2.m[2][j];
			
            result.m[i][j] = accu;
         }
		}
		return result;
	}

	/***************************************************************************
	 * multiplication of matrix and vector
	 **************************************************************************/
	inline Vec3 operator *(const Vec3 v) const
	{
		Vec3 result;
		//
		result.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
		result.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
		result.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z;
		//
		return result;
	}

	/***************************************************************************
	 * multiplication of transposed matrix and vector
	 **************************************************************************/
	inline Vec3 trmul(const Vec3 v) const
	{
		Vec3 result;
		//
		result.x = m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z;
		result.y = m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z;
		result.z = m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z;
		//
		return result;
	}

	/**
	 * define matrix aggregation
	 */
	inline mat3 operator +(const mat3 m2)
	{
		mat3 result;
		int i;
		//
		for (i=0; i<3; i++) {
			result.m[i][0] = m[i][0] + m2.m[i][0];
			result.m[i][1] = m[i][1] + m2.m[i][1];
			result.m[i][2] = m[i][2] + m2.m[i][2];
		}
		return result;
	}

  	inline mat3 operator *(const float s)
	{
		mat3 result;
		//
		for (int i=0; i<3; i++) {
			result.m[i][0] = s * m[i][0];
			result.m[i][1] = s * m[i][1];
			result.m[i][2] = s * m[i][2];
		}
		return result;
	}

	/**
	 * define matrix aggregation
	 */
	inline void operator +=(const mat3 m2)
	{
		for (int i=0; i<3; i++) {
			m[i][0] += m2.m[i][0];
			m[i][1] += m2.m[i][1];
			m[i][2] += m2.m[i][2];
		}
	}

   inline void normv()
   {
      for (int i=0; i<3; i++) {
         float len = m[0][i]*m[0][i] + m[1][i]*m[1][i] + m[2][i]*m[2][i];
         len = sqrtf(len);

         // normalize
         if (!IsFZero(len)) {
            m[0][i] /= len;
            m[1][i] /= len;
            m[2][i] /= len;
         }
      }
   }

   /**
	 * inverse of this matrix
	 *
	 * Gauss - Jordan Algorithm
	 */
	bool inverse(mat3 &result)
	{
		// define a new matrix  (this | identity)
		float gj[3][6];
		float s;
		// 
		// define row - pointers
		float *row0 = &gj[0][0];
		float *row1 = &gj[1][0];
		float *row2 = &gj[2][0];

		// initalize gauss jordan matrix
		row0[0] = m[0][0]; row0[1] = m[0][1]; row0[2] = m[0][2];
		row1[0] = m[1][0]; row1[1] = m[1][1]; row1[2] = m[1][2];
		row2[0] = m[2][0]; row2[1] = m[2][1]; row2[2] = m[2][2];
		//
		row0[3] = 1.0f; row0[4] = 0.0f; row0[5] = 0.0f;
		row1[3] = 0.0f; row1[4] = 1.0f; row1[5] = 0.0f;
		row2[3] = 0.0f; row2[4] = 0.0f; row2[5] = 1.0f;


		// search pivot, if there is no element m[0][x] != 0, matrix is not regular
		if (fabs(row0[0]) < fabs(row1[0])) swap(&row0, &row1);
		if (fabs(row0[0]) < fabs(row2[0])) swap(&row0, &row2);
		if (row0[0] == 0.0f) return false;
		//
		// eliminate first element
		s = row1[0] / row0[0];
		row1[0]-= s*row0[0]; row1[1]-= s*row0[1]; row1[2]-= s*row0[2];
		row1[3]-= s*row0[3]; row1[4]-= s*row0[4]; row1[5]-= s*row0[5];
		s = row2[0] / row0[0];
		row2[0]-= s*row0[0]; row2[1]-= s*row0[1]; row2[2]-= s*row0[2];
		row2[3]-= s*row0[3]; row2[4]-= s*row0[4]; row2[5]-= s*row0[5];
		s = 1.0f / row0[0];
		row0[0]*= s; row0[1] *= s; row0[2] *= s;
		row0[3]*= s; row0[4] *= s; row0[5] *= s;


		// search pivot, if there is no element m[1][x] != 0, matrix is not regular
		if (fabs(row1[1]) < fabs(row2[1])) swap(&row1, &row2);
		if (row1[1] == 0.0f) return false;
		//
		// eliminate second element
		s = 1.0f / row1[1];
		row1[0]*= s; row1[1] *= s; row1[2] *= s;
		row1[3]*= s; row1[4] *= s; row1[5] *= s;
		s = row2[1] / row1[1];
		row2[0]-= s*row1[0]; row2[1]-= s*row1[1]; row2[2]-= s*row1[2];
		row2[3]-= s*row1[3]; row2[4]-= s*row1[4]; row2[5]-= s*row1[5];


		// check if matrix is regular
		if (row2[2] == 0.0f) return false;
		//
		// eliminate last element
		s = 1.0f / row2[2];
		row2[3] *= s; row2[4] *= s; row2[5] *= s;

		// substitute back
		s = row1[2];
		row1[3]-= s*row2[3]; row1[4]-= s*row2[4]; row1[5]-= s*row2[5];
		s = row0[2];
		row0[3]-= s*row2[3]; row0[4]-= s*row2[4]; row0[5]-= s*row2[5];
		//
		s = row0[1];
		row0[3]-= s*row1[3]; row0[4]-= s*row1[4]; row0[5]-= s*row1[5];

		result.m[0][0]=row0[3]; result.m[0][1]=row0[4]; result.m[0][2]=row0[5];
		result.m[1][0]=row1[3]; result.m[1][1]=row1[4]; result.m[1][2]=row1[5];
		result.m[2][0]=row2[3]; result.m[2][1]=row2[4]; result.m[2][2]=row2[5];

		return true;
	}


	void set_rotate_x(float degrees)
	{
		float rad = degrees * 3.14159f / 180.0f;
		float cr = (float)cos(rad);
		float sr = (float)sin(rad);

		m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f;
		m[1][0] = 0.0f; m[1][1] = cr;   m[1][2] = -sr;
		m[2][0] = 0.0f; m[2][1] = sr; 	m[2][2] = cr;
	}

	void set_rotate_y(float degrees)
	{
		float rad = degrees * 3.14159f / 180.0f;
		float cr = (float)cos(rad);
		float sr = -(float)sin(rad);

		m[0][0] = cr;   m[0][1] = 0.0f; m[0][2] = sr;
		m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f;
		m[2][0] = -sr;  m[2][1] = 0.0f; m[2][2] = cr;
	}

	void set_rotate_z(float degrees)
	{
		float rad = degrees * 3.14159f / 180.0f;
		float cr = (float)cos(rad);
		float sr = (float)sin(rad);

		m[0][0] = cr;   m[0][1] = -sr;  m[0][2] = 0.0f;
		m[1][0] = sr;   m[1][1] = cr;   m[1][2] = 0.0f;
		m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f;
	}
	
   /*
   * A function for creating a rotation matrix that rotates a vector called
   * "from" into another vector called "to".
   * Input : from[3], to[3] which both must be *normalized* non-zero vectors
   * Output: mtx[3][3] -- a 3x3 matrix in colum-major form
   *
   * Authors: Tomas Möller, John Hughes
   *          "Efficiently Building a Matrix to Rotate One Vector to Another"
   *          Journal of Graphics Tools, 4(4):1-4, 1999
   *
   * Original sourcecode can be found at
   * http://www.acm.org/jgt/papers/MollerHughes99/code.html
   */
   void fromToRotation(Vec3 &from, Vec3 &to) 
   {
      Vec3 v  = from ^ to;
      float e = from * to;
      float f = (e < 0) ? -e : e;

      /* "from" and "to"-vector almost parallel */
      if (f > 1.0f - EPSILONF) {
         // vector most nearly orthogonal to "from"
         Vec3 x;
         x.x = (from.x > 0.0) ? from.x : -from.x;
         x.y = (from.y > 0.0) ? from.y : -from.y;
         x.z = (from.z > 0.0) ? from.z : -from.z;

         if (x.x < x.y) {
            if (x.x < x.z) x = vec3_e_x;
            else           x = vec3_e_z;
         }
         else {
            if (x.y < x.z) x = vec3_e_y;
            else           x = vec3_e_z;
         }
         Vec3 u = x - from;
         Vec3 v = x - to;

         float c1 = 2.0f / (u * u);
         float c2 = 2.0f / (v * v);
         float c3 = c1 * c2  * (u * v);

         m[0][0] = - c1 * u.x * u.x - c2 * v.x * v.x + c3 * v.x * u.x + 1.0f;
         m[0][1] = - c1 * u.x * u.y - c2 * v.x * v.y + c3 * v.x * u.y;
         m[0][2] = - c1 * u.x * u.z - c2 * v.x * v.z + c3 * v.x * u.z;

         m[1][0] = - c1 * u.y * u.x - c2 * v.y * v.x + c3 * v.y * u.x;
         m[1][1] = - c1 * u.y * u.y - c2 * v.y * v.y + c3 * v.y * u.y + 1.0f;
         m[1][2] = - c1 * u.y * u.z - c2 * v.y * v.z + c3 * v.y * u.z;

         m[2][0] = - c1 * u.z * u.x - c2 * v.z * v.x + c3 * v.z * u.x;
         m[2][1] = - c1 * u.z * u.y - c2 * v.z * v.y + c3 * v.z * u.y;
         m[2][2] = - c1 * u.z * u.z - c2 * v.z * v.z + c3 * v.z * u.z + 1.0f;
      }
      // the most common case, unless "from"="to", or "from"=-"to"
      else {
         // hand optimized version by Gottfried Chen (9 mults less)
         float h = 1.0f / (1.0f + e);
         float hvx = h * v.x;
         float hvz = h * v.z;
         float hvxy = hvx * v.y;
         float hvxz = hvx * v.z;
         float hvyz = hvz * v.y;

         m[0][0] = e + hvx * v.x;
         m[0][1] = hvxy - v.z;
         m[0][2] = hvxz + v.y;

         m[1][0] = hvxy + v.z;
         m[1][1] = e + h * v.y * v.y;
         m[1][2] = hvyz - v.x;

         m[2][0] = hvxz - v.y;
         m[2][1] = hvyz + v.x;
         m[2][2] = e + hvz * v.z;
      }
   }




	/**
	 * clears all elements of the matrix
	 *
	 * because float(0) = %00000000 00000000 00000000 00000000
	 */
	inline void clear()
	{
		memset(m, 0, sizeof(m));
	}

	/**
	 * sets the identity
	 */
	inline void identity()
	{
		memset(m, 0, sizeof(m));
		//
		m[0][0] = m[1][1] = m[2][2] = 1.0f;
	}

	/**
	 * prints the matrix contents
	 */
	void print()
	{
		for (int i=0; i<3; i++)
			printf("[%3.3f, %3.3f, %3.3f]\n", m[i][0], m[i][1], m[i][2]);
	}
};



/*****************************************************************************
 4x4 Matrix definitions

 access elements with mat.m[row][column]
 *****************************************************************************/
class mat4
{
public:

	float m[4][4];

	/**
	 * define matrix multiplication
	 */
	inline mat4 operator *(const mat4 m2)
	{
		mat4 result;
		int i, j;
		float accu;
		//
		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				accu = m[i][0] * m2.m[0][j];
				accu+= m[i][1] * m2.m[1][j];
				accu+= m[i][2] * m2.m[2][j];
				accu+= m[i][3] * m2.m[3][j];
			
            result.m[i][j] = accu;
         }
		}
		return result;
	}


	/**
	 * define matrix aggregation
	 */
	inline mat4 operator +(const mat4 m2)
	{
		mat4 result;
		int i;
		//
		for (i=0; i<3; i++) {
			result.m[i][0] = m[i][0] + m2.m[i][0];
			result.m[i][1] = m[i][1] + m2.m[i][1];
			result.m[i][2] = m[i][2] + m2.m[i][2];
			result.m[i][3] = m[i][3] + m2.m[i][3];
		}
		return result;
	}





	/**
	 * Inverser einer Matrix berechnen
	 *
	 * Gauss - Jordan Algorithm
    *
    * die 4*4 matrix wird zu einer 4*8 matrix 
	 */
	bool inverse(mat4 &result)
	{
		// define a new matrix  (this | identity)
		float gj[4][8];
		float s;
		// 
		// define row - pointers
		float *row0 = &gj[0][0];
		float *row1 = &gj[1][0];
		float *row2 = &gj[2][0];
		float *row3 = &gj[3][0];

		// initalize gauss jordan matrix
		row0[0] = m[0][0]; row0[1] = m[0][1]; row0[2] = m[0][2]; row0[3] = m[0][3];
		row1[0] = m[1][0]; row1[1] = m[1][1]; row1[2] = m[1][2]; row1[3] = m[1][3];
		row2[0] = m[2][0]; row2[1] = m[2][1]; row2[2] = m[2][2]; row2[3] = m[2][3];
		row3[0] = m[3][0]; row3[1] = m[3][1]; row3[2] = m[3][2]; row3[3] = m[3][3];
		//
		row0[4] = 1.0f; row0[5] = 0.0f; row0[6] = 0.0f; row0[7] = 0.0f;
		row1[4] = 0.0f; row1[5] = 1.0f; row1[6] = 0.0f; row1[7] = 0.0f;
		row2[4] = 0.0f; row2[5] = 0.0f; row2[6] = 1.0f; row2[7] = 0.0f;
		row3[4] = 0.0f; row3[5] = 0.0f; row3[6] = 0.0f; row3[7] = 1.0f;


		// search pivot, if there is no element m[0][x] != 0, matrix is not regular
		if (fabs(row0[0]) < fabs(row1[0])) swap(&row0, &row1);
		if (fabs(row0[0]) < fabs(row2[0])) swap(&row0, &row2);
		if (fabs(row0[0]) < fabs(row3[0])) swap(&row0, &row3);
		if (row0[0] == 0.0f) return false;
		//
		// eliminate first element
		s = row1[0] / row0[0];
		row1[0]-= s*row0[0]; row1[1]-= s*row0[1]; row1[2]-= s*row0[2]; row1[3]-= s*row0[3]; 
		row1[4]-= s*row0[4]; row1[5]-= s*row0[5]; row1[6]-= s*row0[6]; row1[7]-= s*row0[7];
		s = row2[0] / row0[0];
		row2[0]-= s*row0[0]; row2[1]-= s*row0[1]; row2[2]-= s*row0[2]; row2[3]-= s*row0[3]; 
		row2[4]-= s*row0[4]; row2[5]-= s*row0[5]; row2[6]-= s*row0[6]; row2[7]-= s*row0[7]; 
		s = row3[0] / row0[0];
		row3[0]-= s*row0[0]; row3[1]-= s*row0[1]; row3[2]-= s*row0[2]; row3[3]-= s*row0[3]; 
		row3[4]-= s*row0[4]; row3[5]-= s*row0[5]; row3[6]-= s*row0[6]; row3[7]-= s*row0[7]; 
		s = 1.0f / row0[0];
		row0[0] *= s; row0[1] *= s; row0[2] *= s; row0[3] *= s;
		row0[4] *= s; row0[5] *= s; row0[6] *= s; row0[7] *= s;


		// search pivot, if there is no element m[1][x] != 0, matrix is not regular
		if (fabs(row1[1]) < fabs(row2[1])) swap(&row1, &row2);
		if (fabs(row1[1]) < fabs(row3[1])) swap(&row1, &row3);
		if (row1[1] == 0.0f) return false;
		//
		// eliminate second element
		s = row2[1] / row1[1];
		row2[0]-= s*row1[0]; row2[1]-= s*row1[1]; row2[2]-= s*row1[2]; row2[3]-= s*row1[3]; 
		row2[4]-= s*row1[4]; row2[5]-= s*row1[5]; row2[6]-= s*row1[6]; row2[7]-= s*row1[7]; 
		s = row3[1] / row1[1];
		row3[0]-= s*row1[0]; row3[1]-= s*row1[1]; row3[2]-= s*row1[2]; row3[3]-= s*row1[3]; 
		row3[4]-= s*row1[4]; row3[5]-= s*row1[5]; row3[6]-= s*row1[6]; row3[7]-= s*row1[7]; 
		s = 1.0f / row1[1];
		row1[0] *= s; row1[1] *= s; row1[2] *= s; row1[3] *= s; 
		row1[4] *= s; row1[5] *= s; row1[6] *= s; row1[7] *= s;

		// search pivot, if there is no element m[1][x] != 0, matrix is not regular
		if (fabs(row2[2]) < fabs(row3[2])) swap(&row2, &row3);
		if (row2[2] == 0.0f) return false;
		//
		// eliminate third element
		s = row3[2] / row2[2];
		row3[0]-= s*row0[0]; row3[1]-= s*row0[1]; row3[2]-= s*row0[2]; row3[3]-= s*row0[3]; 
		row3[4]-= s*row0[4]; row3[5]-= s*row0[5]; row3[6]-= s*row0[6]; row3[7]-= s*row0[7]; 
		s = 1.0f / row2[2];
		row2[0] *= s; row2[1] *= s; row2[2] *= s; row2[3] *= s;
		row2[4] *= s; row2[5] *= s; row2[6] *= s; row2[7] *= s;


		// check if matrix is regular
		if (row3[3] == 0.0f) return false;
		//
		// eliminate last element
		s = 1.0f / row3[3];
		row3[3]=1.0f; row3[4]*=s; row3[5]*=s; row3[6]*=s; row3[7]*=s;


		// substitute back
		s = row2[3];
		row2[4]-=s*row3[4]; row2[5]-=s*row3[5]; row2[6]-=s*row3[6]; row2[7]-=s*row3[7];
		s = row1[3];
		row1[4]-=s*row3[4]; row1[5]-=s*row3[5]; row1[6]-=s*row1[6]; row1[7]-=s*row3[7];
		s = row0[3];
		row0[4]-=s*row3[4]; row0[5]-=s*row3[5]; row0[6]-=s*row1[6]; row0[7]-=s*row3[7];
		//
		s = row1[2];
		row1[4]-=s*row2[4]; row1[5]-=s*row2[5]; row1[6]-=s*row2[6]; row1[7]-=s*row2[7];
		s = row0[2];
		row0[4]-=s*row2[4]; row0[5]-=s*row2[5]; row0[6]-=s*row2[6]; row0[7]-=s*row2[7];
		//
		s = row0[1];
		row0[4]-=s*row1[4]; row0[5]-=s*row1[5]; row0[6]-=s*row1[6]; row0[7]-=s*row1[7];


		result.m[0][0]=row0[4]; result.m[0][1]=row0[5]; result.m[0][2]=row0[6]; result.m[0][3]=row0[7];
		result.m[1][0]=row1[4]; result.m[1][1]=row1[5]; result.m[1][2]=row1[6]; result.m[1][3]=row1[7];
		result.m[2][0]=row2[4]; result.m[2][1]=row2[5]; result.m[2][2]=row2[6]; result.m[2][3]=row2[7];
		result.m[3][0]=row3[4]; result.m[3][1]=row3[5]; result.m[3][2]=row3[6]; result.m[3][3]=row3[7];

		return true;
	}


	/**
	 * evaluate the plane equation 
	 *
	 * plane is defined by the three given vectors
	 * plane equation is stored within a matrix
	 * in the followinf form:
	 *
	 * plane equation: ax + by + cz + d
	 * where [a, b, c] is the unit-length normal to the plane
	 *
	 * | a2 ab ac ad |
	 * | ab b2 bc bd |
	 * | ac bc c2 cd |
	 * | ad bd cd d2 |
	 *
	 * a2 meaning a^2
	 */
	mat4 getPlaneEquation(Vec3 &v1, Vec3 &v2, Vec3 &v3)
	{
		mat4 m;
		float d;

		// evaluate plane normal
		Vec3 av = v2 - v1;
		Vec3 bv = v3 - v1;
		Vec3 n = av ^ bv;
		// 
		n.norm();
		
		// evaluate distance factor
		d = -(n.x * v1.x) - (n.y * v1.y) - (n.z * v1.z);

		// set up the quadric
		m.m[0][0] = n.x * n.x;
		m.m[1][1] = n.y * n.y;
		m.m[2][2] = n.z * n.z;
		m.m[3][3] = d * d;
		//
		m.m[1][0] = m.m[0][1]=n.x * n.y;
		m.m[2][0] = m.m[0][2]=n.x * n.z;
		m.m[3][0] = m.m[0][3]=n.x * d;
		//
		m.m[2][1] = m.m[1][2]=n.y * n.z;
		m.m[3][1] = m.m[1][3]=n.y * d;
		//
		m.m[3][2] = m.m[2][3]=n.z * d;

		return m;
	}


	/**
	 * clears all elements of the matrix
	 *
	 * because float(0) = %00000000 00000000 00000000 00000000
	 */
	inline void clear()
	{
		memset(m, 0, sizeof(m));
	}


	/**
	 * sets the identity
	 */
	inline void identity()
	{
		memset(m, 0, sizeof(m));
		//
		m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
	}


	void set_translate_x(float x)
	{
		m[0][3] = x;
	}

	void set_translate_y(float y)
	{
		m[1][3] = y;
	}

	void set_translate_z(float z)
	{
		m[2][3] = z;
	}


	mat4 &set_rotate_x(float degrees)
	{
		float rad = degrees * 3.14159f / 180.0f;
		float cr = (float)cos(rad);
		float sr = (float)sin(rad);

		m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = 0.0f;
		m[1][0] = 0.0f; m[1][1] = cr;   m[1][2] = -sr;  m[1][3] = 0.0f;
		m[2][0] = 0.0f; m[2][1] = sr; 	m[2][2] = cr;   m[2][3] = 0.0f;
		m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;

		return *this;
	}

	mat4 &set_rotate_y(float degrees)
	{
		float rad = degrees * 3.14159f / 180.0f;
		float cr = (float)cos(rad);
		float sr = (float)sin(rad);

		m[0][0] = cr;   m[0][1] = 0.0f; m[0][2] = sr;   m[0][3] = 0.0f;
		m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = 0.0f;
		m[2][0] = -sr;  m[2][1] = 0.0f; m[2][2] = cr;   m[2][3] = 0.0f;
		m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;

		return *this;
	}

	mat4 &rotate_z(float degrees)
	{
		float rad = degrees * 3.14159f / 180.0f;
		float cr = (float)cos(rad);
		float sr = (float)sin(rad);

		m[0][0] = cr;   m[0][1] = -sr;  m[0][2] = 0.0f; m[0][3] = 0.0f;
		m[1][0] = sr;   m[1][1] = cr;   m[1][2] = 0.0f; m[1][3] = 0.0f;
		m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = 0.0f;
		m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;

		return *this;
	}

	/**
	 * prints the matrix contents
	 */
	void print()
	{
		for (int i=0; i<4; i++)
			printf("[%3.3f, %3.3f, %3.3f, %3.3f]\n", m[i][0], m[i][1], m[i][2], m[i][3]);
	}

   	/**
	 * prints the matrix contents
	 */
	void print(float *m)
	{
		printf("[%3.3f, %3.3f, %3.3f, %3.3f]\n", m[0], m[1], m[2], m[3]);
	}
};






/*
 * Compute inverse of 4x4 transformation SINGLE-PRECISION matrix.
 * Code contributed by Jacques Leroy <jle@star.be>
 * Code lifted from Brian Paul's Mesa freeware OpenGL implementation.
 * Return GL_TRUE for success, GL_FALSE for failure (singular matrix)
 */
inline bool invertMatrixf(float *out, const float *m)
{
/* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { float *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(m,r,c) (m)[(c)*4+(r)]

  float wtmp[4][8];
  float m0, m1, m2, m3, s;
  float *r0, *r1, *r2, *r3;

  r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

  r0[0] = MAT(m,0,0), r0[1] = MAT(m,0,1),
  r0[2] = MAT(m,0,2), r0[3] = MAT(m,0,3),
  r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0,

  r1[0] = MAT(m,1,0), r1[1] = MAT(m,1,1),
  r1[2] = MAT(m,1,2), r1[3] = MAT(m,1,3),
  r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0,

  r2[0] = MAT(m,2,0), r2[1] = MAT(m,2,1),
  r2[2] = MAT(m,2,2), r2[3] = MAT(m,2,3),
  r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0,

  r3[0] = MAT(m,3,0), r3[1] = MAT(m,3,1),
  r3[2] = MAT(m,3,2), r3[3] = MAT(m,3,3),
  r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;

  /* choose pivot - or die */
  if (fabs(r3[0])>fabs(r2[0])) SWAP_ROWS(r3, r2);
  if (fabs(r2[0])>fabs(r1[0])) SWAP_ROWS(r2, r1);
  if (fabs(r1[0])>fabs(r0[0])) SWAP_ROWS(r1, r0);
  if (0.0 == r0[0]) {
    return false;
  }

  /* eliminate first variable     */
  m1 = r1[0]/r0[0]; m2 = r2[0]/r0[0]; m3 = r3[0]/r0[0];
  s = r0[1]; r1[1] -= m1 * s; r2[1] -= m2 * s; r3[1] -= m3 * s;
  s = r0[2]; r1[2] -= m1 * s; r2[2] -= m2 * s; r3[2] -= m3 * s;
  s = r0[3]; r1[3] -= m1 * s; r2[3] -= m2 * s; r3[3] -= m3 * s;
  s = r0[4];
  if (s != 0.0) { r1[4] -= m1 * s; r2[4] -= m2 * s; r3[4] -= m3 * s; }
  s = r0[5];
  if (s != 0.0) { r1[5] -= m1 * s; r2[5] -= m2 * s; r3[5] -= m3 * s; }
  s = r0[6];
  if (s != 0.0) { r1[6] -= m1 * s; r2[6] -= m2 * s; r3[6] -= m3 * s; }
  s = r0[7];
  if (s != 0.0) { r1[7] -= m1 * s; r2[7] -= m2 * s; r3[7] -= m3 * s; }

  /* choose pivot - or die */
  if (fabs(r3[1])>fabs(r2[1])) SWAP_ROWS(r3, r2);
  if (fabs(r2[1])>fabs(r1[1])) SWAP_ROWS(r2, r1);
  if (0.0 == r1[1]) {
    return false;
  }

  /* eliminate second variable */
  m2 = r2[1]/r1[1]; m3 = r3[1]/r1[1];
  r2[2] -= m2 * r1[2]; r3[2] -= m3 * r1[2];
  r2[3] -= m2 * r1[3]; r3[3] -= m3 * r1[3];
  s = r1[4]; if (0.0 != s) { r2[4] -= m2 * s; r3[4] -= m3 * s; }
  s = r1[5]; if (0.0 != s) { r2[5] -= m2 * s; r3[5] -= m3 * s; }
  s = r1[6]; if (0.0 != s) { r2[6] -= m2 * s; r3[6] -= m3 * s; }
  s = r1[7]; if (0.0 != s) { r2[7] -= m2 * s; r3[7] -= m3 * s; }

  /* choose pivot - or die */
  if (fabs(r3[2])>fabs(r2[2])) SWAP_ROWS(r3, r2);
  if (0.0 == r2[2]) {
    return false;
  }

  /* eliminate third variable */
  m3 = r3[2]/r2[2];
  r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
  r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6],
  r3[7] -= m3 * r2[7];

  /* last check */
  if (0.0 == r3[3]) {
    return false;
  }

  s = 1.0f/r3[3];              /* now back substitute row 3 */
  r3[4] *= s; r3[5] *= s; r3[6] *= s; r3[7] *= s;

  m2 = r2[3];                 /* now back substitute row 2 */
  s  = 1.0f/r2[2];
  r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
  r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
  m1 = r1[3];
  r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
  r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
  m0 = r0[3];
  r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
  r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

  m1 = r1[2];                 /* now back substitute row 1 */
  s  = 1.0f/r1[1];
  r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
  r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
  m0 = r0[2];
  r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
  r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

  m0 = r0[1];                 /* now back substitute row 0 */
  s  = 1.0f/r0[0];
  r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
  r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

  MAT(out,0,0) = r0[4]; MAT(out,0,1) = r0[5],
  MAT(out,0,2) = r0[6]; MAT(out,0,3) = r0[7],
  MAT(out,1,0) = r1[4]; MAT(out,1,1) = r1[5],
  MAT(out,1,2) = r1[6]; MAT(out,1,3) = r1[7],
  MAT(out,2,0) = r2[4]; MAT(out,2,1) = r2[5],
  MAT(out,2,2) = r2[6]; MAT(out,2,3) = r2[7],
  MAT(out,3,0) = r3[4]; MAT(out,3,1) = r3[5],
  MAT(out,3,2) = r3[6]; MAT(out,3,3) = r3[7]; 

  return true;

#undef MAT
#undef SWAP_ROWS
}


/*
 * Compute inverse of 4x4 transformation SINGLE-PRECISION matrix.
 * Code contributed by Jacques Leroy <jle@star.be>
 * Code lifted from Brian Paul's Mesa freeware OpenGL implementation.
 * Return GL_TRUE for success, GL_FALSE for failure (singular matrix)
 */
inline bool invertMatrixd(double *out, const double *m)
{
/* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { double *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(m,r,c) (m)[(c)*4+(r)]

  double wtmp[4][8];
  double m0, m1, m2, m3, s;
  double *r0, *r1, *r2, *r3;

  r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

  r0[0] = MAT(m,0,0); r0[1] = MAT(m,0,1);
  r0[2] = MAT(m,0,2); r0[3] = MAT(m,0,3);
  r0[4] = 1.0; r0[5] = r0[6] = r0[7] = 0.0,

  r1[0] = MAT(m,1,0); r1[1] = MAT(m,1,1);
  r1[2] = MAT(m,1,2); r1[3] = MAT(m,1,3);
  r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0;

  r2[0] = MAT(m,2,0); r2[1] = MAT(m,2,1);
  r2[2] = MAT(m,2,2); r2[3] = MAT(m,2,3);
  r2[6] = 1.0; r2[4] = r2[5] = r2[7] = 0.0;

  r3[0] = MAT(m,3,0); r3[1] = MAT(m,3,1);
  r3[2] = MAT(m,3,2); r3[3] = MAT(m,3,3);
  r3[7] = 1.0; r3[4] = r3[5] = r3[6] = 0.0;

  /* choose pivot - or die */
  if (fabs(r3[0])>fabs(r2[0])) SWAP_ROWS(r3, r2);
  if (fabs(r2[0])>fabs(r1[0])) SWAP_ROWS(r2, r1);
  if (fabs(r1[0])>fabs(r0[0])) SWAP_ROWS(r1, r0);
  if (0.0 == r0[0]) {
    return false;
  }

  /* eliminate first variable     */
  m1 = r1[0]/r0[0]; m2 = r2[0]/r0[0]; m3 = r3[0]/r0[0];
  s = r0[1]; r1[1] -= m1 * s; r2[1] -= m2 * s; r3[1] -= m3 * s;
  s = r0[2]; r1[2] -= m1 * s; r2[2] -= m2 * s; r3[2] -= m3 * s;
  s = r0[3]; r1[3] -= m1 * s; r2[3] -= m2 * s; r3[3] -= m3 * s;
  s = r0[4];
  if (s != 0.0) { r1[4] -= m1 * s; r2[4] -= m2 * s; r3[4] -= m3 * s; }
  s = r0[5];
  if (s != 0.0) { r1[5] -= m1 * s; r2[5] -= m2 * s; r3[5] -= m3 * s; }
  s = r0[6];
  if (s != 0.0) { r1[6] -= m1 * s; r2[6] -= m2 * s; r3[6] -= m3 * s; }
  s = r0[7];
  if (s != 0.0) { r1[7] -= m1 * s; r2[7] -= m2 * s; r3[7] -= m3 * s; }

  /* choose pivot - or die */
  if (fabs(r3[1])>fabs(r2[1])) SWAP_ROWS(r3, r2);
  if (fabs(r2[1])>fabs(r1[1])) SWAP_ROWS(r2, r1);
  if (0.0 == r1[1]) {
    return false;
  }

  /* eliminate second variable */
  m2 = r2[1]/r1[1]; m3 = r3[1]/r1[1];
  r2[2] -= m2 * r1[2]; r3[2] -= m3 * r1[2];
  r2[3] -= m2 * r1[3]; r3[3] -= m3 * r1[3];
  s = r1[4]; if (0.0 != s) { r2[4] -= m2 * s; r3[4] -= m3 * s; }
  s = r1[5]; if (0.0 != s) { r2[5] -= m2 * s; r3[5] -= m3 * s; }
  s = r1[6]; if (0.0 != s) { r2[6] -= m2 * s; r3[6] -= m3 * s; }
  s = r1[7]; if (0.0 != s) { r2[7] -= m2 * s; r3[7] -= m3 * s; }

  /* choose pivot - or die */
  if (fabs(r3[2])>fabs(r2[2])) SWAP_ROWS(r3, r2);
  if (0.0 == r2[2]) {
    return false;
  }

  /* eliminate third variable */
  m3 = r3[2]/r2[2];
  r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
  r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6],
  r3[7] -= m3 * r2[7];

  /* last check */
  if (0.0 == r3[3]) {
    return false;
  }

  s = 1.0/r3[3];              /* now back substitute row 3 */
  r3[4] *= s; r3[5] *= s; r3[6] *= s; r3[7] *= s;

  m2 = r2[3];                 /* now back substitute row 2 */
  s  = 1.0/r2[2];
  r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
  r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
  m1 = r1[3];
  r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
  r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
  m0 = r0[3];
  r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
  r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

  m1 = r1[2];                 /* now back substitute row 1 */
  s  = 1.0/r1[1];
  r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
  r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
  m0 = r0[2];
  r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
  r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

  m0 = r0[1];                 /* now back substitute row 0 */
  s  = 1.0/r0[0];
  r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
  r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

  MAT(out,0,0) = r0[4]; MAT(out,0,1) = r0[5],
  MAT(out,0,2) = r0[6]; MAT(out,0,3) = r0[7],
  MAT(out,1,0) = r1[4]; MAT(out,1,1) = r1[5],
  MAT(out,1,2) = r1[6]; MAT(out,1,3) = r1[7],
  MAT(out,2,0) = r2[4]; MAT(out,2,1) = r2[5],
  MAT(out,2,2) = r2[6]; MAT(out,2,3) = r2[7],
  MAT(out,3,0) = r3[4]; MAT(out,3,1) = r3[5],
  MAT(out,3,2) = r3[6]; MAT(out,3,3) = r3[7]; 

  return true;

#undef MAT
#undef SWAP_ROWS
}

#endif