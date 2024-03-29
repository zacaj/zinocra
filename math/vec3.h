#pragma once

#include "main.h"
#include <assimp.h>        // Plain-C interface
template <class T> class vec4;//#include "vec4.h"

#ifdef USE_BULLET
#include <btBulletDynamicsCommon.h>
#endif
class mat4f;
//#include "Serialization.h"
class SerializationHandler;
class Buffer;
enum SerializationMode;

template <class T>
class vec3
{
public:
	T x,y,z;
	vec3()
	{
		x=y=z=0;
	}
	vec3(T t)
	{
		x=y=z=t;
	}
	vec3(T _x,T _y)
	{
		x=_x;
		y=_y;
		z=0;
	}
	vec3(T _x,T _y,T _z)
	{
		x=_x;
		y=_y;
		z=_z;
	}
	void operator() (T _x,T _y,T _z)
	{
		x=_x;
		y=_y;
		z=_z;
	}
	vec3(vec4<T> v)
	{
		x=v.x;
		y=v.y;
		z=v.z;
	}
	StaticSerializer(vec3)
	{
		SF(x)
		SF(y)
		SF(z)
	}
	float dot(vec3<T> v)
	{
		return x*v.x+y*v.y+z*v.z;
	}
	inline float magnitude()
	{
		return sqrt(x*x+y*y+z*z);
	}
	inline T magnitudesq()
	{
		return x*x+y*y+z*z;
	}
	vec3<T> normalized()
	{
		float mag=magnitude();
		return vec3(x/mag,y/mag,z/mag);
	}
	void normalize()
	{
		float mag=magnitude();
		x/=mag;
		y/=mag;
		z/=mag;
	}
	float dist(vec3<T> b)
	{
		vec3<T> v=*this-b;
		return v.magnitude();
	}
	vec3<T> operator * (const float &f) const
	{
		return vec3<T>(x*f,y*f,z*f);
	}
	vec3<T> operator / (const float &f) const
	{
		return vec3<T>(x/f,y/f,z/f);
	}
	vec3<T> operator + (const vec3<T> &v) const
	{
		return vec3<T>(x+v.x,y+v.y,z+v.z);
	}
	vec3<T> operator - (const vec3<T> &v) const
	{
		return vec3<T>(x-v.x,y-v.y,z-v.z);
	}
	vec3<T> operator - (const vec4<T> &v) const
	{
		return vec3<T>(x-v.x,y-v.y,z-v.z);
	}
	vec3<T> operator * (const vec3<T> &v) const
	{
		return vec3<T>(x*v.x,y*v.y,z*v.z);
	}
	vec3<T> operator += (const T &f)
	{
		x+=f;
		y+=f;
		z+=f;
		return *this;
	}
	vec3<T> operator -= (const T &f)
	{
		x-=f;
		y-=f;
		z-=f;
		return *this;
	}
	vec3<T> operator *= (const T &f)
	{
		x*=f;
		y*=f;
		z*=f;
		return *this;
	}
	vec3<T> operator /= (const T &f)
	{
		x/=f;
		y/=f;
		z/=f;
		return *this;
	}
	vec3<T> operator - () const
	{
		return vec3<T>(-x,-y,-z);
	}
	vec3<T> operator += (const vec3<T> &v)
	{
		x+=v.x;
		y+=v.y;
		z+=v.z;
		return *this;
	}
	vec3<T> operator -= (const vec3<T> &v)
	{
		x-=v.x;
		y-=v.y;
		z-=v.z;
		return *this;
	}
	vec3<T> operator /= (const vec3<T> &v)
	{
		x/=v.x;
		y/=v.y;
		z/=v.z;
		return *this;
	}
	vec3<T> operator * (const mat4f &m) const//note probably not correct for axis oriantation used
	{
		/*vec3<T> r;
		r.x=x * m.m[0] +
			y * m.m[1] +
			z * m.m[2] +
			m.m[3];

		r.y=x * m.m[4] +
			y * m.m[5] +
		//	z * m.m[6] +//assimp uses left handed y up
			//m.m[7];
			z*m.m[6] - 
			m.m[11];

		r.z=x * m.m[8] +
			y * m.m[9] +
		//	z * m.m[10] +
			//m.m[11];
			z * m.m[10] +
			m.m[7];*/
		BREAK;//IF THIS IS CALLED MAKE THE CALLER NOT CALL IT SWITCH THE ORDER BUGS BE HERE
		return vec3<T>(0);
	}
	vec3<T> operator * (const aiMatrix4x4 &m) const
	{
		vec3<T> r;
		r.x=x * m.a1 +
			y * m.a2 +
			z * m.a3 +
			m.a4;

		r.y=x * m.b1 +
			y * m.b2 +
			z * m.b3 +//assimp uses left handed y up
			m.b4;

		r.z=x * m.c1 +
			y * m.c2 +
			z * m.c3 +
			m.c4;
		return r;
	}
	vec3<T> operator *= (const aiMatrix4x4 &m)
	{
		(*this)=(*this) * m;
		return *this;
	}
	vec3<T> operator *= (const mat4f &m)
	{
		(*this)=m*(*this);
		return *this;
	}
	operator vec4<T>()
	{
		return vec4<T>(x,y,z,1);
	}
	bool operator >= (T n)
	{
		return x>=n && y>=n && z>=n;
	}
	bool operator < (T n)
	{
		return x<n && y<n && z<n;
	}
#ifdef USE_BULLET
	operator btVector3()
	{
		return btVector3(x,y,z);
	}
	vec3<T>(btVector3 &v)
	{
		x=v.getX();
		y=v.getY();
		z=v.getZ();
	}
#endif
};

typedef vec3<float> vec3f;
typedef vec3<double> vec3d;
typedef vec3<int> vec3i;
typedef vec3<unsigned int> vec3u;
typedef vec3<unsigned char> vec3b;
typedef vec3<char> vec3c;

