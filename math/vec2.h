#pragma once
#ifdef NLIB
#include "main.h"
#include "vec3.h"
#endif
template <class T>
class vec2
{
public:
	T x,y;
	vec2()
	{
		x=y=0;
	}
	vec2<T> operator () (T _n)
	{
		x=_n;
		y=_n;
		return *this;
	}
	vec2(T t)
	{
		x=y=t;
	}
	vec2(T _x,T _y)
	{
		x=_x;
		y=_y;
	}
#ifdef NLIB
	vec2(vec3<T> v)
	{
		x=v.x;
		y=v.y;
	}
}
#endif
#ifdef SERIALIZE
	StaticSerializer(vec2)
	{
		SF(x)
		SF(y)
	}
#endif
	vec2<T> operator + (const vec2<T> &v) const
	{
		return vec2<T>(x+v.x,y+v.y);
	}
	vec2<T> operator - (const vec2<T> &v) const
	{
		return vec2<T>(x-v.x,y-v.y);
	}
	vec2<T> operator * (const float &f) const
	{
		return vec2<T>(x*f,y*f);
	}
	vec2<T> operator * (const vec2<T> &f) const
	{
		return vec2<T>(x*f.x,y*f.y);
	}
	vec2<T> operator / (const float &f) const
	{
		return vec2<T>(x/f,y/f);
	}
	vec2<T>& operator += (const vec2<T> &s)
	{
		x+=s.x;
		y+=s.y;
		return *this;
	}
	vec2<T>& operator /= (const float &s)
	{
		x/=s;
		y/=s;
		return *this;
	}
	vec2<T>& operator *= (const float &s)
	{
		x*=s;
		y*=s;
		return *this;
	}
	vec2<T>& operator -= (const float &s)
	{
		x-=s;
		y-=s;
		return *this;
	}
	bool operator == (const vec2<T> &v)
	{
		return v.x==x && v.y==y;
	}
	bool operator != (const vec2<T> &v)
	{
		return v.x!=x || v.y!=y;
	}
	vec2<T> operator() (T _x,T _y)
	{
		x=_x;
		y=_y;
		return *this;
	}
	float dot(vec2<T> v)
	{
		return x*v.x+y*v.y;
	}
	inline float magnitude()
	{
		return sqrt(x*x+y*y);
	}
	inline T magnitudesq()
	{
		return x*x+y*y;
	}
	vec2<T> normalized()
	{
		float mag=magnitude();
		return vec2(x/mag,y/mag);
	}
	void normalize()
	{
		float mag=magnitude();
		x/=mag;
		y/=mag;
	}
};

typedef vec2<float> vec2f;
typedef vec2<double> vec2d;
typedef vec2<int> vec2i;
#ifdef NLIB
typedef vec2<uint> vec2u;
typedef vec2<uchar> vec2b;
#endif
typedef vec2<char> vec2c;