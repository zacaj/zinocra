#pragma once

#include "main.h"
#include "vec3.h"
//#include "mat4f.h"
class mat4f;
#include <assimp.h>        // Plain-C interface
#include "vec4.h"
#include "mat4f.h"
#include <stdio.h>
#ifdef USE_BULLET
class btQuaternion;
#endif
/// Quaternion.  Look it up.

class Quaternion
{
public:
	float x,y,z,w;
	Quaternion(void);
	/**
		Makes a quaternion from a YXZ Euler rotation (degrees)

		*/
	Quaternion(vec3f rot);
	/**
		\param a Angle in degrees

		*/
	Quaternion(vec3f axis,float angle);
	Quaternion(float X,float Y,float Z,float W);
	Quaternion(vec4f v)
	{
		x=v.x;
		y=v.y;
		z=v.z;
		w=v.w;
	}
#ifdef USE_BULLET
	Quaternion(const btQuaternion &q);
#endif
	StaticSerializer(Quaternion)
	{
		SF(x)
		SF(y)
		SF(z)
		SF(w)
	}
	Quaternion normalized();
	void normalize();
	Quaternion conjugate();
	/**
		\param a Angle in degrees

		*/
	void axisAngle(vec3f v,float a);
	void YXZ(vec3f a);
	Quaternion& operator = (const Quaternion &q);
	Quaternion& operator = (const aiQuaternion &q);
	operator vec4f();
	operator String();
	operator mat4f();
#ifdef USE_BULLET
	operator btQuaternion();
#endif
	Quaternion operator * (const Quaternion &q)
	{
		Quaternion C;

		C.x = w*q.x + x*q.w + y*q.z - z*q.y;
		C.y = w*q.y - x*q.z + y*q.w + z*q.x;
		C.z = w*q.z + x*q.y - y*q.x + z*q.w;
		C.w = w*q.w - x*q.x - y*q.y - z*q.z;
		return C;
	}
	Quaternion operator *= (const Quaternion &q)
	{
		(*this)=(*this) * q;
		return *this;
	}
};
