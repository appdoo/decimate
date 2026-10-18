//=============================================================================
//
// Vec.h 
//
// Definition von Vektoren
//
// $Author: Andor $
// $Date: 17.12.03 0:30 $
// $Revision: 7 $
// 
//=============================================================================
#ifndef VEC_H
#define VEC_H

// local includes
#include "../def.h"

// C++ Includes
#include <math.h>

// ----------------------------------------------------------------------------
// type forwarding
// ----------------------------------------------------------------------------
typedef struct Vec2 Vec2;
typedef struct Vec3 Vec3;
typedef struct Vertex vertex;

// ----------------------------------------------------------------------------
// extern vector constants
// ----------------------------------------------------------------------------
extern Vec3 vec3_e_x;
extern Vec3 vec3_e_y;
extern Vec3 vec3_e_z;


// ----------------------------------------------------------------------------
// vec2
//
// twodimensional vector type
//
// ----------------------------------------------------------------------------
typedef struct Vec2
{
   float x;
   float y;
} 
Vec2;



// ----------------------------------------------------------------------------
// Vec3
//
// threedimensional vector type
//
// ----------------------------------------------------------------------------
typedef struct Vec3
{
	float x, y, z;

   // -------------------------------------------------------------------------
   // 3D vector addition
   // -------------------------------------------------------------------------
	inline Vec3 operator +(const Vec3 &v) const
   {
      Vec3 result;

	   result.x = x + v.x;
	   result.y = y + v.y;
	   result.z = z + v.z;

		return result;
	}

   // -------------------------------------------------------------------------
   // 3D vector addition
   // -------------------------------------------------------------------------
	inline void operator +=(const Vec3 &v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
	}

   // -------------------------------------------------------------------------
   // 3D vector subtraktion
   // -------------------------------------------------------------------------
	inline Vec3 operator -(const Vec3 &v) const
	{
		Vec3 result;

      result.x = x - v.x;
      result.y = y - v.y;
		result.z = z - v.z;

		return result;
	}

   // -------------------------------------------------------------------------
   // 3D vector subtraktion
   // -------------------------------------------------------------------------
	inline void operator -=(const Vec3 &v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

   // -------------------------------------------------------------------------
   // vector cross product as ^ operator
   // -------------------------------------------------------------------------
	inline Vec3 operator ^(const Vec3 &v) const
	{
		Vec3 result;

      result.x = y * v.z - z * v.y;
      result.y = z * v.x - x * v.z;
      result.z = x * v.y - y * v.x;

		return result;
	}

   // -------------------------------------------------------------------------
   // scalar multiply
   // -------------------------------------------------------------------------
	inline Vec3 operator *(const float s) const
	{
		Vec3 result;

		result.x = x * s;
		result.y = y * s;
		result.z = z * s;

		return result;
	}

   // -------------------------------------------------------------------------
   // scalar divide
   // -------------------------------------------------------------------------
	inline Vec3 operator /(const float s) const
	{
		Vec3 result;

		result.x = x / s;
		result.y = y / s;
		result.z = z / s;

		return result;
	}

   // -------------------------------------------------------------------------
   // vector dot product
   // -------------------------------------------------------------------------
	inline float dotproduct(const Vec3 &v)
	{
		return x*v.x + y*v.y + z*v.z;
	}

   // -------------------------------------------------------------------------
   // vector dot product
   // -------------------------------------------------------------------------
	inline float operator *(const Vec3 &v) const
	{
		return x*v.x + y*v.y + z*v.z;
	}
   
   // -------------------------------------------------------------------------
   // vector normalization
   // -------------------------------------------------------------------------
	inline Vec3 &norm()
	{
      register float norm;
      
      norm  = sqrtf(x*x + y*y + z*z);

      if (!IsFZero(norm)) {
		   x /= norm;
		   y /= norm;
		   z /= norm;
      }
      return *this;
	}

   // -------------------------------------------------------------------------
   // vector length
   // -------------------------------------------------------------------------
	inline float length()
   {
      float len = sqrtf(x*x + y*y + z*z);
      return IsFZero(len) ? 0.0f : len;
	}

   // -------------------------------------------------------------------------
   // squared vector length
   // -------------------------------------------------------------------------
	inline float length2()
	{
		return x*x + y*y + z*z;
	}
} Vec3;


// ----------------------------------------------------------------------------
//
// vertex
//
// a weighted vertex definition
// (TODO: add weight count)
// ----------------------------------------------------------------------------

typedef struct Vertex
{
   // texture coordinates
   float u, v;
   // point index
   int   index;

} Vertex;



#endif