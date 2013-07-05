#pragma once

#include "main.h"
#include "vec3.h"
#include "vec4.h"
#ifdef Model
struct BoneTransform;
#endif
class Quaternion;
#ifdef USE_ASSIMP
#include <assimp.h>        // Plain-C interface
#endif


/**
	4x4 float matrix

	0  1  2  3
	4  5  6  7
	8  9  10 11
	12 13 14 15

|   0,0	0,1	0,2	0,3	|
|	1,0	1,1	1,2	1,3	|
|	2,0	2,1	2,2	2,3	|
|	3,0	3,1	3,2	3,3	|
	*/
SerializableClass(mat4f)
{
public:
	float m[16];
	mat4f(void);///< Initializes the mat4f to an identity matrix
	~mat4f(void);
	mat4f(bool t);
	mat4f(float s);///< Initializes the mat4f to all \a s
	mat4f(const float *data);///< Initializes the mat4f with 16 floats from \a data
#ifdef USE_ASSIMP
	mat4f(aiMatrix4x4 &mat);
#endif
	mat4f(float x,float y,float z=0);
	//void serialize(SerializationHandler *s,Buffer *b,bool d,SerializationMode mode=COMPLETE)
	Serializer(mat4f)
	{
		for(int i=0;i<16;i++)
		{
			float t=m[i];
			SF(t)
			m[i]=t;
		}
	}
	/**
		Creates an axis angle orientation matrix

		\param axis The axis to rotate around
		\param angle The angle to rotate around \a axis
		*/
	void axisAngle(vec3f axis,float angle);
	/**
		Creates a perspective projection matrix

		\param fov Field of View of the projection (45-90 recommended)
		\param zNear The distance to the near plane of the viewing frustum
		\param zFar The distance to the far plane of the viewing frustum
		\param aspect The aspect ratio for the projection (w/h)
	*/
	void perspective(float fov,float znear,float zfar,float aspect);
	void setAxes(vec3f right,vec3f up,vec3f backward);
	/**
		Creates an orthograpic projection matrix

		If you really cant figure this out, google glOrtho
		*/
	void orthographic( float left,float right,float bottom,float top,float zNear,float zFar );
	void zero(); ///< Zeros the matrix;
	/**
		Creates a matrix that would scale a 2x2x2 cube centered on origin to the extents given

		*/
	void scaleUnitToExtents(vec3f min,vec3f max);
	void identity();///< Creates an identity matrix;
	/**
		\returns A transposed version of the matrix

		*/
	mat4f transposed() const;
	void transpose();///< Transposes the matrix;
	mat4f inversed() const;
	void invert();///< Transposes the matrix;
	void lookAt(float centerx, float centery, float centerz,float upx, float upy, float upz);

	void upload(); ///< Uploads the matrix to the current Shader's transform matrix
	void multAndUpload(); ///< Multiplies the current Shader's transorfm matrix
	void uploadProjection();
	/**
		Sets the first 3 diagonals to \a s
		*/
	void scaleTo(vec3f s);
	/**
		Makes a new mat4f that only applies the orientation from this mat4f

		\returns That matrix
		*/
	mat4f mat3f();
	/// Extracts position from matrix
	vec3f getPosition();
	void setPosition(vec3f p);
	mat4f& operator = (const mat4f &mat);
	operator float*()
	{
		return m;
	}
	char* c_str();
	operator String();
	void printf();
	/*float& operator[] (unsigned i)
	{
		return m[i];
	}*///?fix
#ifdef Model
	mat4f operator = (const BoneTransform &p);
#endif
#ifdef USE_ASSIMP
	mat4f operator = (const aiQuaternion &q);
	mat4f & operator = (const Quaternion &q);
#endif

	mat4f operator + (const mat4f &mat) const;
	mat4f& operator += (const mat4f &mat);
	mat4f operator - (const mat4f &mat) const;
	mat4f& operator -= (const mat4f &mat);
	mat4f& operator - ();
	vec3f operator * (const vec3f &v) const;
	vec4f operator * (const vec4f &v) const;
	mat4f operator * (const mat4f &mat) const;
	mat4f operator *= (const aiMatrix4x4 &mat);
	mat4f operator * (const aiMatrix4x4 &mat);
	mat4f operator *= (const mat4f &mat);
	mat4f operator * (const float &s) const;
	mat4f& operator *= (const float &s);
	mat4f operator / (const float &s) const;
	mat4f& operator /= (const float &s);
};

extern mat4f &currentTransform;