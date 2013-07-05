#include "StdAfx.h"
#include "mat4f.h"
#include "Model.h"
#include <math.h>
#include "gl.h"
#include "utility.h"
mat4f tempm4f;
mat4f &currentTransform=tempm4f;

mat4f::mat4f(void)
{
	identity();
}

mat4f::mat4f(float s)
{
	for(int i=0;i<16;i++)
		m[i]=s;
}

mat4f::mat4f( const float *data )
{
	//for(int i=0;i<16;i++)
	//	m[i]=data[i];
	memcpy(m,data,sizeof(float)*16);
}

mat4f::mat4f( bool t )
{

}

mat4f::mat4f( aiMatrix4x4 &mat )
{
	m[0]=mat.a1;
	m[1]=mat.a2;
	m[2]=mat.a3;
	m[3]=mat.a4;

	m[4]=mat.b1;
	m[5]=mat.b2;
	m[6]=mat.b3;
	m[7]=mat.b4;

	m[8]=mat.c1;
	m[9]=mat.c2;
	m[10]=mat.c3;
	m[11]=mat.c4;

	m[12]=mat.d1;
	m[13]=mat.d2;
	m[14]=mat.d3;
	m[15]=mat.d4;
}

mat4f::mat4f( float x,float y,float z/*=0*/ )
{
	identity();
	m[3]=x;
	m[7]=y;
	m[11]=z;
}

mat4f::~mat4f(void)
{
}

mat4f mat4f::transposed() const
{
	float r[16];
	r[0]=m[0];
	r[1]=m[4];
	r[2]=m[8];
	r[3]=m[12];
	r[4]=m[1];
	r[5]=m[5];
	r[6]=m[9];
	r[7]=m[13];
	r[8]=m[2];
	r[9]=m[6];
	r[10]=m[10];
	r[11]=m[14];
	r[12]=m[3];
	r[13]=m[7];
	r[14]=m[11];
	r[15]=m[15];
	return mat4f(r);
}

void mat4f::perspective(float fov,float znear,float zfar,float aspect)
{
    float xymax = znear * (float)tan(fov * (M_PI/360.f));
    float ymin = -xymax;
    float xmin = -xymax;

    float width = xymax - xmin;
    float height = xymax - ymin;

    float depth = zfar - znear;
    float q = -(zfar + znear) / depth;
    float qn = -2 * (zfar * znear) / depth;

    float w = 2 * znear / width;
    w = w / aspect;
    float h = 2 * znear / height;
	m[0]  = w;m[1]  = 0;m[2]  = 0;m[3] = 0;

	m[4]  = 0;m[5]  = h;m[6]  = 0;m[7] = 0;

	m[8]  = 0;m[9]  = 0;m[10] = q;m[11] =qn;

	m[12]  = 0;m[13]  = 0;m[14] = -1;m[15] = 0;/**/

}

void mat4f::axisAngle(vec3f axis,float angle)
{
	float c = (float)cos(angle*DEG2RAD);
	float s = (float)sin(angle*DEG2RAD);
	float t = 1.f - c;
	//  if axis is not already normalized then uncomment this
	// double magnitude = Math.sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
	// if (magnitude==0) throw error;
	// a1.x /= magnitude;
	// a1.y /= magnitude;
	// a1.z /= magnitude;
	m[0] = c + axis.x*axis.x*t;
	m[5] = c + axis.y*axis.y*t;
	m[10] = c + axis.z*axis.z*t;

	float tmp1 = axis.x*axis.y*t;
	float tmp2 = axis.z*s;
	m[4] = tmp1 + tmp2;
	m[1] = tmp1 - tmp2;
	tmp1 = axis.x*axis.z*t;
	tmp2 = axis.y*s;
	m[8] = tmp1 - tmp2;
	m[2] = tmp1 + tmp2;    tmp1 = axis.y*axis.z*t;
	tmp2 = axis.x*s;
	m[9] = tmp1 + tmp2;
	m[6] = tmp1 - tmp2;
}
void mat4f::printf()
{
	print("\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n"
		,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12],m[13],m[14],m[15]);
}
void mat4f::lookAt(float centerx, float centery, float centerz,
	float upx, float upy, float upz)
{
    float x[3], y[3], z[3];
    float mag;  
    /* Make rotation matrix */
    
    /* Z vector */
    z[0] = - centerx;
    z[1] = - centery;
    z[2] = - centerz;
    mag = sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
    if (mag) {          /* mpichler, 19950515 */
        z[0] /= mag;
        z[1] /= mag;
        z[2] /= mag;
    }
    
    /* Y vector */
    y[0] = upx;
    y[1] = upy;
    y[2] = upz;
    
    /* X vector = Y cross Z */
    x[0] = y[1] * z[2] - y[2] * z[1];
    x[1] = -y[0] * z[2] + y[2] * z[0];
    x[2] = y[0] * z[1] - y[1] * z[0];
    
    /* Recompute Y = Z cross X */
    y[0] = z[1] * x[2] - z[2] * x[1];
    y[1] = -z[0] * x[2] + z[2] * x[0];
    y[2] = z[0] * x[1] - z[1] * x[0];
    
    /* mpichler, 19950515 */
    /* cross product gives area of parallelogram, which is < 1.0 for
     * non-perpendicular unit-length vectors; so normalize x, y here
     */
    
    mag = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    if (mag) {
        x[0] /= mag;
        x[1] /= mag;
        x[2] /= mag;
    }
    
    mag = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    if (mag) {
        y[0] /= mag;
        y[1] /= mag;
        y[2] /= mag;
    }
    
#define M(row,col)  m[col+row*4]
    M(0, 0) = x[0];
    M(0, 1) = x[1];
    M(0, 2) = x[2];
    M(0, 3) = 0.0;
    M(1, 0) = y[0];
    M(1, 1) = y[1];
    M(1, 2) = y[2];
    M(1, 3) = 0.0;
    M(2, 0) = z[0];
    M(2, 1) = z[1];
    M(2, 2) = z[2];
    M(2, 3) = 0.0;
    M(3, 0) = 0.0;
    M(3, 1) = 0.0;
    M(3, 2) = 0.0;
    M(3, 3) = 1.0;
#undef M    
}

void mat4f::upload()
{
	//currentTransform=*this;
#ifdef GL1
	transpose();
	glLoadMatrixf(m);
	transpose();
#elif defined GL234
	glUniformMatrix4fv(currentShader->transformID,1,GL_TRUE,m);
#elif defined SW
	currentTransformMatrix=*this;
#endif
}

void mat4f::uploadProjection()
{
#ifdef GL1
	glMatrixMode(GL_PROJECTION);
	upload();
	glMatrixMode(GL_MODELVIEW);
#elif defined GL234
	glUniformMatrix4fv(currentShader->projectID,1,GL_TRUE,m);
#elif defined SW
	currentProjectionMatrix=*this;
#endif
}


void mat4f::multAndUpload()
{
	error("Not implemented\n",0);

	//*this=currentTransform* *this;
	//glUniformMatrix4fv(currentShader->transformID,1,GL_TRUE,m);
}

void mat4f::scaleTo( vec3f s )
{
	m[0]=s.x;
	m[5]=s.y;
	m[10]=s.z;
}

mat4f mat4f::mat3f()
{
	float r[16];
	r[0]=m[0];
	r[1]=m[1];
	r[2]=m[2];
	r[3]=0;
	r[4]=m[4];
	r[5]=m[5];
	r[6]=m[6];
	r[7]=0;
	r[8]=m[8];
	r[9]=m[9];
	r[10]=m[10];
	r[11]=0;
	r[12]=0;
	r[13]=0;
	r[14]=0;
	r[15]=1;
	return mat4f(r);
}

void mat4f::orthographic( float left,float right,float bottom,float top,float zNear,float zFar )
{
	m[0]=2.f/(right-left);
	m[4]=m[8]=0;
	m[3]=-(right+left)/(right-left);
	m[1]=0;
	m[5]=2.f/(top-bottom);
	m[9]=0;
	m[7]=-(top+bottom)/(top-bottom);
	m[2]=m[6]=0;
	m[10]=-2.f/(zFar-zNear);
	m[11]=-(zFar+zNear)/(zFar-zNear);
	m[12]=m[13]=m[14]=0;
	m[15]=1;
}
mat4f mat4f::operator = (const BoneTransform &p)
{
	*this=p.q;
	m[3]=p.p.x;
	m[7]=p.p.y;
	m[11]=p.p.z;
	return *this;
}
mat4f & mat4f::operator = (const Quaternion &q)
{
	float sqw = q.w*q.w;
	float sqx = q.x*q.x;
	float sqy = q.y*q.y;
	float sqz = q.z*q.z;

	// invs (inverse square length) is only required if quaternion is not already normalised
	float invs = 1 ;/// (sqx + sqy + sqz + sqw);
	m[0] = ( sqx - sqy - sqz + sqw)*invs ; // since sqw + sqx + sqy + sqz =1/invs*invs
	m[5] = (-sqx + sqy - sqz + sqw)*invs ;
	m[10] = (-sqx - sqy + sqz + sqw)*invs ;

	float tmp1 = q.x*q.y;
	float tmp2 = q.z*q.w;
	m[4] = (float)2.0 * (tmp1 + tmp2)*invs ;
	m[1] = (float)2.0 * (tmp1 - tmp2)*invs ;

	tmp1 = q.x*q.z;
	tmp2 = q.y*q.w;
	m[8] = (float)2.0 * (tmp1 - tmp2)*invs ;
	m[2] = (float)2.0 * (tmp1 + tmp2)*invs ;
	tmp1 = q.y*q.z;
	tmp2 = q.x*q.w;
	m[9] = (float)2.0 * (tmp1 + tmp2)*invs ;
	m[6] = (float)2.0 * (tmp1 - tmp2)*invs ;
	m[3]=m[7]=m[11]=0.0;
	m[15]=1;
	m[3]=m[7]=m[11]=m[12]=m[13]=m[14]=0;
	return *this;
}

mat4f& mat4f::operator=( const mat4f &mat )
{
	for(int i=0;i<16;i++)
		m[i]=mat.m[i];
	return *this;
}

mat4f mat4f::operator=( const aiQuaternion &q )
{
	float sqw = q.w*q.w;
	float sqx = q.x*q.x;
	float sqy = q.y*q.y;
	float sqz = q.z*q.z;

	// invs (inverse square length) is only required if quaternion is not already normalised
	float invs = 1 ;/// (sqx + sqy + sqz + sqw);
	m[0] = ( sqx - sqy - sqz + sqw)*invs ; // since sqw + sqx + sqy + sqz =1/invs*invs
	m[5] = (-sqx + sqy - sqz + sqw)*invs ;
	m[10] = (-sqx - sqy + sqz + sqw)*invs ;

	float tmp1 = q.x*q.y;
	float tmp2 = q.z*q.w;
	m[4] = (float)2.0 * (tmp1 + tmp2)*invs ;
	m[1] = (float)2.0 * (tmp1 - tmp2)*invs ;

	tmp1 = q.x*q.z;
	tmp2 = q.y*q.w;
	m[8] = (float)2.0 * (tmp1 - tmp2)*invs ;
	m[2] = (float)2.0 * (tmp1 + tmp2)*invs ;
	tmp1 = q.y*q.z;
	tmp2 = q.x*q.w;
	m[9] = (float)2.0 * (tmp1 + tmp2)*invs ;
	m[6] = (float)2.0 * (tmp1 - tmp2)*invs ;
	m[3]=m[7]=m[11]=0.0;
	m[15]=1;
	m[3]=m[7]=m[11]=m[12]=m[13]=m[14]=0;
	return *this;
}

mat4f::operator String()
{
	char temp[256];
	sprintf(temp,"\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n"
		,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12],m[13],m[14],m[15]);
	String ret=temp;
	return ret;
}

mat4f mat4f::inversed() const
/*{
	mat4f inv(0.f);
	mat4f invOut(0.f);
	inv[0] =   m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15]
	+ m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
	inv[4] =  -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15]
	- m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
	inv[8] =   m[4]*m[9]*m[15] - m[4]*m[11]*m[13] - m[8]*m[5]*m[15]
	+ m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
	inv[12] = -m[4]*m[9]*m[14] + m[4]*m[10]*m[13] + m[8]*m[5]*m[14]
	- m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
	inv[1] =  -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15]
	- m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
	inv[5] =   m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15]
	+ m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
	inv[9] =  -m[0]*m[9]*m[15] + m[0]*m[11]*m[13] + m[8]*m[1]*m[15]
	- m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
	inv[13] =  m[0]*m[9]*m[14] - m[0]*m[10]*m[13] - m[8]*m[1]*m[14]
	+ m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
	inv[2] =   m[1]*m[6]*m[15] - m[1]*m[7]*m[14] - m[5]*m[2]*m[15]
	+ m[5]*m[3]*m[14] + m[13]*m[2]*m[7] - m[13]*m[3]*m[6];
	inv[6] =  -m[0]*m[6]*m[15] + m[0]*m[7]*m[14] + m[4]*m[2]*m[15]
	- m[4]*m[3]*m[14] - m[12]*m[2]*m[7] + m[12]*m[3]*m[6];
	inv[10] =  m[0]*m[5]*m[15] - m[0]*m[7]*m[13] - m[4]*m[1]*m[15]
	+ m[4]*m[3]*m[13] + m[12]*m[1]*m[7] - m[12]*m[3]*m[5];
	inv[14] = -m[0]*m[5]*m[14] + m[0]*m[6]*m[13] + m[4]*m[1]*m[14]
	- m[4]*m[2]*m[13] - m[12]*m[1]*m[6] + m[12]*m[2]*m[5];
	inv[3] =  -m[1]*m[6]*m[11] + m[1]*m[7]*m[10] + m[5]*m[2]*m[11]
	- m[5]*m[3]*m[10] - m[9]*m[2]*m[7] + m[9]*m[3]*m[6];
	inv[7] =   m[0]*m[6]*m[11] - m[0]*m[7]*m[10] - m[4]*m[2]*m[11]
	+ m[4]*m[3]*m[10] + m[8]*m[2]*m[7] - m[8]*m[3]*m[6];
	inv[11] = -m[0]*m[5]*m[11] + m[0]*m[7]*m[9] + m[4]*m[1]*m[11]
	- m[4]*m[3]*m[9] - m[8]*m[1]*m[7] + m[8]*m[3]*m[5];
	inv[15] =  m[0]*m[5]*m[10] - m[0]*m[6]*m[9] - m[4]*m[1]*m[10]
	+ m[4]*m[2]*m[9] + m[8]*m[1]*m[6] - m[8]*m[2]*m[5];

	float det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
	if (det == 0)
		return false;

	det = 1.0 / det;

	for (int i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return invOut;
}*/
{
	float a0 = m[ 0]*m[ 5] - m[ 1]*m[ 4];
	float a1 = m[ 0]*m[ 6] - m[ 2]*m[ 4];
	float a2 = m[ 0]*m[ 7] - m[ 3]*m[ 4];
	float a3 = m[ 1]*m[ 6] - m[ 2]*m[ 5];
	float a4 = m[ 1]*m[ 7] - m[ 3]*m[ 5];
	float a5 = m[ 2]*m[ 7] - m[ 3]*m[ 6];
	float b0 = m[ 8]*m[13] - m[ 9]*m[12];
	float b1 = m[ 8]*m[14] - m[10]*m[12];
	float b2 = m[ 8]*m[15] - m[11]*m[12];
	float b3 = m[ 9]*m[14] - m[10]*m[13];
	float b4 = m[ 9]*m[15] - m[11]*m[13];
	float b5 = m[10]*m[15] - m[11]*m[14];

	mat4f inverse;

	float det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;
	//if (Math<float>::FAbs(det) > epsilon)
	{
		//Matrix4 inverse;
		inverse.m[ 0] = + m[ 5]*b5 - m[ 6]*b4 + m[ 7]*b3;
		inverse.m[ 4] = - m[ 4]*b5 + m[ 6]*b2 - m[ 7]*b1;
		inverse.m[ 8] = + m[ 4]*b4 - m[ 5]*b2 + m[ 7]*b0;
		inverse.m[12] = - m[ 4]*b3 + m[ 5]*b1 - m[ 6]*b0;
		inverse.m[ 1] = - m[ 1]*b5 + m[ 2]*b4 - m[ 3]*b3;
		inverse.m[ 5] = + m[ 0]*b5 - m[ 2]*b2 + m[ 3]*b1;
		inverse.m[ 9] = - m[ 0]*b4 + m[ 1]*b2 - m[ 3]*b0;
		inverse.m[13] = + m[ 0]*b3 - m[ 1]*b1 + m[ 2]*b0;
		inverse.m[ 2] = + m[13]*a5 - m[14]*a4 + m[15]*a3;
		inverse.m[ 6] = - m[12]*a5 + m[14]*a2 - m[15]*a1;
		inverse.m[10] = + m[12]*a4 - m[13]*a2 + m[15]*a0;
		inverse.m[14] = - m[12]*a3 + m[13]*a1 - m[14]*a0;
		inverse.m[ 3] = - m[ 9]*a5 + m[10]*a4 - m[11]*a3;
		inverse.m[ 7] = + m[ 8]*a5 - m[10]*a2 + m[11]*a1;
		inverse.m[11] = - m[ 8]*a4 + m[ 9]*a2 - m[11]*a0;
		inverse.m[15] = + m[ 8]*a3 - m[ 9]*a1 + m[10]*a0;

		float invDet = ((float)1)/det;
		inverse.m[ 0] *= invDet;
		inverse.m[ 1] *= invDet;
		inverse.m[ 2] *= invDet;
		inverse.m[ 3] *= invDet;
		inverse.m[ 4] *= invDet;
		inverse.m[ 5] *= invDet;
		inverse.m[ 6] *= invDet;
		inverse.m[ 7] *= invDet;
		inverse.m[ 8] *= invDet;
		inverse.m[ 9] *= invDet;
		inverse.m[10] *= invDet;
		inverse.m[11] *= invDet;
		inverse.m[12] *= invDet;
		inverse.m[13] *= invDet;
		inverse.m[14] *= invDet;
		inverse.m[15] *= invDet;

		return inverse;
	}

//	return ZERO;
}

void mat4f::setAxes( vec3f right,vec3f up,vec3f backward )
{
	m[0]=right.x;
	m[4]=right.y;
	m[8]=right.z;
	m[1]=up.x;
	m[5]=up.y;
	m[9]=up.z;
	m[2]=backward.x;
	m[6]=backward.y;
	m[10]=backward.z;
}

void mat4f::scaleUnitToExtents( vec3f min,vec3f max )
{
	m[0]=(fabs(max.x)+fabs(min.x))/2;
	m[5]=(fabs(max.y)+fabs(min.y))/2;
	m[10]=(fabs(max.z)+fabs(min.z))/2;
	m[3]=-(fabs(min.x)+fabs(max.x))/2+max.x;
	m[7]=-(fabs(min.y)+fabs(max.y))/2+max.y;
	m[11]=-(fabs(min.z)+fabs(max.z))/2+max.z;
}

char* mat4f::c_str()
{
	char *temp=new char[256];
	sprintf(temp,"\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n%2.3f %2.3f %2.3f %2.3f\n"
		,m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12],m[13],m[14],m[15]);
	return temp;
}

void mat4f::zero() /*/< Zeros the matrix */
{
	//for(int i=0;i<16;i++)
	//	m[i]=0;
	memset(m,0,sizeof(float)*16);
}

void mat4f::identity() /*/< Creates an identity matrix */
{
	zero();
	m[0]=m[5]=m[10]=m[15]=1;
}

void mat4f::transpose() /*/< Transposes the matrix */
{
	*this=transposed();
}

void mat4f::invert() /*/< Transposes the matrix */
{
	*this=inversed();
}

vec3f mat4f::getPosition()
{
	return vec3f(m[3],m[7],m[11]);
}

void mat4f::setPosition( vec3f p )
{
	m[3]=p.x;
	m[7]=p.y;
	m[11]=p.z;
}

mat4f& mat4f::operator-()
{
	for(int i=0;i<16;i++)
		m[i]=-m[i];
	return *this;
}

mat4f mat4f::operator + (const mat4f &mat) const
{
	mat4f r;
	for(int i=0;i<16;i++)
		r.m[i]=m[i]+mat.m[i];
	return r;
}

mat4f& mat4f::operator += (const mat4f &mat)
{
	for(int i=0;i<16;i++)
		m[i]+=mat.m[i];
	return *this;
}

mat4f mat4f::operator - (const mat4f &mat) const
{
	mat4f r;
	for(int i=0;i<16;i++)
		r.m[i]=m[i]-mat.m[i];
	return r;
}
mat4f& mat4f::operator -= (const mat4f &mat)
{
	for(int i=0;i<16;i++)
		m[i]-=mat.m[i];
	return *this;
}


	vec3f mat4f::operator * (const vec3f &v) const
	{
		return vec3f(v.x * m[0] +
			v.y * m[1] +
			v.z * m[2] +
			m[3],

		v.x * m[4] +
			v.y * m[5] +
			v.z * m[6] +
			m[7],

		v.x * m[8] +
			v.y * m[9] +
			v.z * m[10] +
			m[11]);
		/*r.w=v.x*m[12]+
			v.y*m[13]+
			v.z*m[14]+
			m[15]*/
	}
	/*vec4f operator * (const vec3f &v) const
	{
		vec4f r;
		r.x=v.x * m[0] +
			v.y * m[1] +
			v.z * m[2] +
			m[3];

		r.y=v.x * m[4] +
			v.y * m[5] +
			v.z * m[6] +
			m[7];

		r.z=v.x * m[8] +
			v.y * m[9] +
			v.z * m[10] +
			m[11];
		r.w=m[12]+
			m[13]+
			m[14]+
			m[15];
		return r;
	}*/
	vec4f mat4f::operator * (const vec4f &v) const
	{
		return vec4f(v.x * m[0] +
			v.y * m[1] +
			v.z * m[2] +
			v.w*m[3],

		v.x * m[4] +
			v.y * m[5] +
			v.z * m[6] +
			v.w*m[7],

		v.x * m[8] +
			v.y * m[9] +
			v.z * m[10] +
			v.w*m[11],
		v.x*m[12]+
			v.y*m[13]+
			v.z*m[14]+
			v.w*m[15]);
	}
	mat4f mat4f::operator * (const mat4f &mat) const
	{
		//mat4f r;
		float r[16];
		r[0]=m[0]*mat.m[0]+m[1]*mat.m[4]+m[2]*mat.m[8]+m[3]*mat.m[12];
		r[1]=m[0]*mat.m[1]+m[1]*mat.m[5]+m[2]*mat.m[9]+m[3]*mat.m[13];
		r[2]=m[0]*mat.m[2]+m[1]*mat.m[6]+m[2]*mat.m[10]+m[3]*mat.m[14];
		r[3]=m[0]*mat.m[3]+m[1]*mat.m[7]+m[2]*mat.m[11]+m[3]*mat.m[15];

		r[4]=m[4]*mat.m[0]+m[5]*mat.m[4]+m[6]*mat.m[8]+m[7]*mat.m[12];
		r[5]=m[4]*mat.m[1]+m[5]*mat.m[5]+m[6]*mat.m[9]+m[7]*mat.m[13];
		r[6]=m[4]*mat.m[2]+m[5]*mat.m[6]+m[6]*mat.m[10]+m[7]*mat.m[14];
		r[7]=m[4]*mat.m[3]+m[5]*mat.m[7]+m[6]*mat.m[11]+m[7]*mat.m[15];

		r[8]=m[8]*mat.m[0]+m[9]*mat.m[4]+m[10]*mat.m[8]+m[11]*mat.m[12];
		r[9]=m[8]*mat.m[1]+m[9]*mat.m[5]+m[10]*mat.m[9]+m[11]*mat.m[13];
		r[10]=m[8]*mat.m[2]+m[9]*mat.m[6]+m[10]*mat.m[10]+m[11]*mat.m[14];
		r[11]=m[8]*mat.m[3]+m[9]*mat.m[7]+m[10]*mat.m[11]+m[11]*mat.m[15];

		r[12]=m[12]*mat.m[0]+m[13]*mat.m[4]+m[14]*mat.m[8]+m[15]*mat.m[12];
		r[13]=m[12]*mat.m[1]+m[13]*mat.m[5]+m[14]*mat.m[9]+m[15]*mat.m[13];
		r[14]=m[12]*mat.m[2]+m[13]*mat.m[6]+m[14]*mat.m[10]+m[15]*mat.m[14];
		r[15]=m[12]*mat.m[3]+m[13]*mat.m[7]+m[14]*mat.m[11]+m[15]*mat.m[15];
		return mat4f(r);
	}
	mat4f mat4f::operator *= (const aiMatrix4x4 &mat)
	{
		(*this)=(*this) * mat;
		return *this;
	}
	mat4f mat4f::operator * (const aiMatrix4x4 &mat)
	{
		//mat4f r;
		float r[16];
		r[0]=m[0]*mat.a1+m[1]*mat.b1+m[2]*mat.c1+m[3]*mat.d1;
		r[1]=m[0]*mat.a2+m[1]*mat.b2+m[2]*mat.c2+m[3]*mat.d2;
		r[2]=m[0]*mat.a3+m[1]*mat.b3+m[2]*mat.c3+m[3]*mat.d3;
		r[3]=m[0]*mat.a4+m[1]*mat.b4+m[2]*mat.c4+m[3]*mat.d4;

		r[4]=m[4]*mat.a1+m[5]*mat.b1+m[6]*mat.c1+m[7]*mat.d1;
		r[5]=m[4]*mat.a2+m[5]*mat.b2+m[6]*mat.c2+m[7]*mat.d2;
		r[6]=m[4]*mat.a3+m[5]*mat.b3+m[6]*mat.c3+m[7]*mat.d3;
		r[7]=m[4]*mat.a4+m[5]*mat.b4+m[6]*mat.c4+m[7]*mat.d4;

		r[8]=m[8]*mat.a1+m[9]*mat.b1+m[10]*mat.c1+m[11]*mat.d1;
		r[9]=m[8]*mat.a2+m[9]*mat.b2+m[10]*mat.c2+m[11]*mat.d2;
		r[10]=m[8]*mat.a3+m[9]*mat.b3+m[10]*mat.c3+m[11]*mat.d3;
		r[11]=m[8]*mat.a4+m[9]*mat.b4+m[10]*mat.c4+m[11]*mat.d4;

		r[12]=m[12]*mat.a1+m[13]*mat.b1+m[14]*mat.c1+m[15]*mat.d1;
		r[13]=m[12]*mat.a2+m[13]*mat.b2+m[14]*mat.c2+m[15]*mat.d2;
		r[14]=m[12]*mat.a3+m[13]*mat.b3+m[14]*mat.c3+m[15]*mat.d3;
		r[15]=m[12]*mat.a4+m[13]*mat.b4+m[14]*mat.c4+m[15]*mat.d4;
		return mat4f(r);
	}
	mat4f mat4f::operator *= (const mat4f &mat)
	{
		(*this)=(*this) * mat;
		return *this;
	}
	mat4f mat4f::operator * (const float &s) const
	{
		mat4f r;
		for(int i=0;i<16;i++)
			r.m[i]=m[i]*s;
		return r;
	}
	mat4f& mat4f::operator *= (const float &s)
	{
		for(int i=0;i<16;i++)
			m[i]*=s;
		return *this;
	}
	mat4f mat4f::operator / (const float &s) const
	{
		mat4f r;
		float inverse=1.0f/s;
		for(int i=0;i<16;i++)
			r.m[i]=m[i]*inverse;
		return r;
	}
	mat4f& mat4f::operator /= (const float &s)
	{
		float inverse=1.0f/s;
		for(int i=0;i<16;i++)
			m[i]*=inverse;
		return *this;
	}