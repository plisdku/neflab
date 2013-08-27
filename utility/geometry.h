/*
 *  geometry.h
 *  vmlibtest
 *
 *  Created by Paul Hansen on 2/8/08.
 *  Copyright 2008 Stanford University. All rights reserved.
 *
 *  $Rev:: 184                           $:  Revision of last commit
 *  $Author:: pch                        $:  Author of last commit
 *
 *  $Date: 2009-10-16 16:37:28 -0700 (Fri, 16 Oct 2009) $:
 *  $Id: geometry.h 184 2009-10-16 23:37:28Z pch $:
 *
 */

#ifndef _GEOMETRY_
#define _GEOMETRY_

#include "VectorMatrix.h"
#include <vector>

typedef Matrix3i Mat3i;
typedef Matrix3f Mat3f;
typedef Matrix3d Mat3d;
typedef Matrix3b Mat3b;

template<typename T>
class Triangle3;


template<typename T>
class ClosedInterval
{
public:
    ClosedInterval()
    {
        mEndpoints[0] = T();
        mEndpoints[1] = T();
    }
    
    ClosedInterval(const T & t1, const T & t2)
    {
        mEndpoints[0] = t1;
        mEndpoints[1] = t2;
    }
    
    const T & operator[](unsigned int nn) const { return mEndpoints[nn]; }
    T & operator[](unsigned int nn) { return mEndpoints[nn]; }
    
    T length() const { return mEndpoints[1] - mEndpoints[0]; }
    
    bool intersects(const ClosedInterval<T> & i2) const
    {
        return !(mEndpoints[1] < i2[0] || mEndpoints[0] > i2[1]);
    }
    
private:
    T mEndpoints[2];
};

template<typename T>
class Rect
{
public:
	Rect();
	Rect(T i0, T j0, T k0, T i1, T j1, T k1);
	Rect(const Vector3<T> & p0, const Vector3<T> & p1);
	//Rect(const Rect<T> & copyMe);
    
    template<class S>
    Rect(const Rect<S> & copyMe);
	
    // size is calculated as p2-p1
	T size(unsigned int dim) const;
	Vector3<T> size() const;
    
    /**
     *  Return the signed volume of this rectangle.  The volume is the product
     *  of size() in each dimension.  If p1 and p2 represent the corners of
     *  the rectangle then this is the volume (for instance if this rectangle
     *  is floating-point-valued).  If p1 and p2 are inclusive indices of cells
     *  at corners, consider using count() instead.
     */
    T volume() const;
    
    // num is calculated as p2-p1+1 (so, "inclusive" size suitable for ints)
    T num(unsigned int dim) const;
    Vector3<T> num() const;
    
    /**
     *  Return the signed volume of this rectangle.  The volume is the product
     *  of num() in each dimension.  If p1 and p2 are inclusive indices of cells
     *  at corners, this will return the total number of cells in the rectangle.
     *  If p1 and p2 represent corners, consider using volume() instead.
     */
    T count() const;
	
	bool
	encloses(const Rect<T> & inRect) const;
	
	bool
	encloses(const Vector3<T> & inPt) const;
	
	bool
	encloses(T x, T y, T z) const;
    
    bool
    encloses(const Triangle3<T> & tri) const;
	
	bool
	intersects(const Rect<T> & inRect) const;
    
    bool
    intersects(const Triangle3<T> & tri) const;
	
	int
	numNonSingularDims() const;
    
    bool operator==(const Rect<T> & rhs) const
    {
        return (p1 == rhs.p1 && p2 == rhs.p2);
    }
    bool operator!=(const Rect<T> & rhs) const
    {
        return (p1 != rhs.p1 || p2 != rhs.p2);
    }
	
    ClosedInterval<T> projectionAlong(const Vector3<T> & v) const;
    
	Vector3<T> p1;
	Vector3<T> p2; // must be indexwise < p1.
};

typedef Rect<long> Rect3i;
typedef Rect<double> Rect3d;

template<typename T>
Rect<T>
operator * (const Rect<T> & lhs, const T rhs);

template<typename T>
Rect<T>
operator / (const Rect<T> & lhs, const T rhs);

template<typename T>
Rect<T>
operator * (const T scalar, const Rect<T> & rhs);

template<typename T>
Rect<T> &
operator *= (Rect<T> & lhs, const T scalar);

template<typename T>
Rect<T> &
operator /= (Rect<T> & lhs, const T scalar);

template<typename T>
Rect<T>
operator * (const Matrix3<T> & lhs, const Rect<T> & rhs);

template<typename T>
Rect<T>
operator + (const Rect<T> & lhs, const Vector3<T> & rhs);

template<typename T>
Rect<T>
operator - (const Rect<T> & lhs, const Vector3<T> & rhs);

template<typename T>
Rect<T>
operator + (const Vector3<T> & lhs, const Rect<T> & rhs);

template<typename T>
Rect<T>
operator - (const Vector3<T> & lhs, const Rect<T> & rhs);

template<typename T, typename S>
Rect<T>
inset( const Rect<T> & inRect, S dx0, S dy0, S dz0, S dx1, S dy1, S dz1);

template<typename T, typename S>
Rect<T>
inset( const Rect<T> & inRect, S distance);

template<typename T>
Vector3<T>
clip( const Rect<T> & clipRect, const Vector3<T> & v );

template<typename T>
Rect<T>
clip( const Rect<T> & clipRect, const Rect<T> & rectToClip );

template<typename T>
Rect<T>
intersection( const Rect<T> & r1, const Rect<T> & r2 );

template<typename T>
Rect<T>
cyclicPermute(const Rect<T> & r, unsigned int nn);

template<typename T>
bool operator<(const Rect<T> & lhs, const Rect<T> & rhs);

template<typename T>
bool operator>(const Rect<T> & lhs, const Rect<T> & rhs);

#pragma mark *** OrientedRect ***

template<typename T>
struct OrientedRect
{
	OrientedRect();
	OrientedRect(const Rect<T> & r, const Vector3<T> & n) :
		rect(r),
		normal(n)
	{}
	
	Rect<T> rect;
	Vector3<T> normal;
};
typedef OrientedRect<long> OrientedRect3i;
typedef OrientedRect<double> OrientedRect3d;

template<typename T>
OrientedRect<T>
cyclicPermute(const OrientedRect<T> & r, unsigned int nn);

#pragma mark *** Plane ***

// dot(plane.normal(), v) + plane.constant() = 0 for v on the plane.
template<typename T>
class Plane3
{
public:
    Plane3();
    Plane3(const Vector3<T> & normal, const T & constant);
    
    const Vector3<T> & normal() const { return mNormal; }
    void normal(const Vector3<T> & normal) { mNormal = normal; }
    const T & constant() const { return mConstant; }
    void constant(const T & constant) { mConstant = constant; }
    
protected:
    Vector3<T> mNormal;
    T mConstant;
};
typedef Plane3<long> Plane3i;
typedef Plane3<double> Plane3d;

template<class T>
Vector3<T>
intersection(const Plane3<T> & plane, const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection);

// return time t such that rayOrigin + t*rayDirection lies on the plane.
template<class T>
T intersectionTime(const Plane3<T> & plane, const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection);

// return the derivative of the intersection time t (see above) with respect
// to the plane normal and plane constant.  The components of the derivative
// are returned in a Plane3 structure.
template<class T>
Plane3<T>
intersectionTimeJacobian(const Plane3<T> & plane, const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection);

template<typename T>
std::ostream & operator<<(std::ostream & str, const Plane3<T> & plane);

template<typename T>
std::istream & operator>>(std::istream & str, Plane3<T> & plane);



#pragma mark *** Polygon ***

template<typename T>
class Polygon3
{
public:
    Polygon3();
    Polygon3(const std::vector<Vector3<T>*> & outerBoundary,
        const std::vector<std::vector<Vector3<T>*> > & holes);
    Polygon3(const std::vector<std::vector<Vector3<T>*> > & contours);
    
    const std::vector<Vector3<T>*> & outerBoundary() const
    {
        return mOuterBoundary;
    }
    
    const std::vector<std::vector<Vector3<T>*> > & holes() const
    {
        return mHoles;
    }
    
private:
    std::vector<Vector3<T>*> mOuterBoundary;
    std::vector<std::vector<Vector3<T>*> > mHoles;
};
typedef Polygon3<double> Polygon3d;

#pragma mark *** Triangle ***

template<typename T>
class Triangle3
{
public:
    Triangle3()
    {
    }
    
    Triangle3(const Vector3<T> & p0, const Vector3<T> & p1,
        const Vector3<T> & p2)
    {
        mVertices[0] = p0;
        mVertices[1] = p1;
        mVertices[2] = p2;
    }
    
    Rect<T> bounds() const
    {
        return Rect<T>(
            vec_min(vec_min(mVertices[0], mVertices[1]), mVertices[2]),
            vec_max(vec_max(mVertices[0], mVertices[1]), mVertices[2]));
    }
    
    Vector3<T> & operator[](int vertex) { return mVertices[vertex]; }
    const Vector3<T> & operator[](int vertex) const
        { return mVertices[vertex]; }
    
    Vector3<T> normal() const;
    Plane3<T> plane() const;
    T area() const { return norm(cross(mVertices[1]-mVertices[0],
        mVertices[2]-mVertices[0]))/2; }
    
    Vector3<T> edgeDirection(int nn) const;
    
    bool enclosesProjection(const Vector3<T> & point) const;
    
    ClosedInterval<T> projectionAlong(const Vector3<T> & v) const;
    
    // Return the derivative of the normal vector with respect to perturbation
    // of a given vertex.
    Matrix3<T> dndv(unsigned int vertex) const;
    
    // Return the derivative of the plane constant with respect to perturbation
    // of a given vertex.
    Vector3<T> dcdv(unsigned int vertex) const;
    
    // Return the derivative of the UNIT normal vector with respect to
    // perturbation of a given vertex.  Convenience function.
    Matrix3<T> dudv(unsigned int vertex) const;
    
    /**
     * Express a ray-triangle intersection in the basis
     *   v1-v0
     *   v2-v0
     *   rayDirection
     * such that the ray-triangle intersection can be written
     *   v0 + a*(v1-v0) + b*(v2-v0)
     * and the distance t from the ray origin is defined so that
     *   intersection point = rayOrigin + t*rayDirection.
     * The return value is the coefficient vector (a, b, t).
    **/
    Vector3<T> intersectionCoefficients(const Vector3<T> & rayOrigin,
        const Vector3<T> & rayDirection) const;
    
    Matrix3<T> intersectionCoefficientsJacobian(const Vector3<T> & rayOrigin,
        const Vector3<T> & rayDirection, int vertex) const;
    
    Matrix3<T> intersectionJacobian(const Vector3<T> & rayOrigin,
        const Vector3<T> & rayDirection, int vertex) const;
    
    bool operator==(const Triangle3<T> & rhs) const;
private:
    Vector3<T> mVertices[3];
};

template<typename T>
bool sameVertices(const Triangle3<T> & t1, const Triangle3<T> & t2)
{
    return ( (t1[0] == t2[0] || t1[0] == t2[1] || t1[0] == t2[2]) &&
        (t1[1] == t2[0] || t1[1] == t2[1] || t1[1] == t2[2]) &&
        (t1[2] == t2[0] || t1[2] == t2[1] || t1[2] == t2[2]) );
}

template<typename T>
bool parallelUnitNormals(const Triangle3<T> & t1, const Triangle3<T> & t2,
    const T & tolerance)
{
    return dot(unit(t1.normal()), unit(t2.normal())) > (T(1) - tolerance);
}

template<typename T>
bool oppositeUnitNormals(const Triangle3<T> & t1, const Triangle3<T> & t2,
    const T & tolerance)
{
    return dot(unit(t1.normal()), -unit(t2.normal())) > (T(1) - tolerance);
}

template<typename T>
Triangle3<T>
operator * (const Matrix3<T> & lhs, const Triangle3<T> & rhs);

template<typename T>
Triangle3<T>
operator + (const Triangle3<T> & tri, const Vector3<T> & v);

template<typename T>
Triangle3<T>
operator + (const Vector3<T> & v, const Triangle3<T> & tri);

typedef Triangle3<double> Triangle3d;


#pragma mark *** Point and line stuff ***

//template<typename T>
//T
//pointLineDistance(const Vector3<T> & point,
//    const Vector3<T> & pointOnLine0, const Vector3<T> & pointOnLine1)
//{
//    return sqrt(pointLineDistanceSquared(point, pointOnLine0, pointOnLine1));
//}

template<typename T>
T
pointLineDistanceSquared(const Vector3<T> & point,
    const Vector3<T> & pointOnLine0, const Vector3<T> & pointOnLine1)
{
    // Area = 0.5*base*height
    // Area = 0.5*cross product
    // The point-line distance formula follows, without any square roots.
    
    return sumSquares(cross(point - pointOnLine0, point - pointOnLine1)) /
        sumSquares(pointOnLine1 - pointOnLine0);
}


template<typename T>
T
pointSegmentDistanceSquared(const Vector3<T> & point,
    const Vector3<T> & segmentBegin, const Vector3<T> & segmentEnd)
{
    Vector3<T> direction = segmentEnd - segmentBegin;
    
    double distAlongLine = dot(point, direction);
    
    if (distAlongLine > dot(segmentEnd, direction))
        return sumSquares(point - segmentEnd);
    else if (distAlongLine < dot(segmentBegin, direction))
        return sumSquares(point - segmentBegin);
    else
        return pointLineDistanceSquared(point, segmentBegin, segmentEnd);
}


#pragma mark *** MeshFace ***

// TODO: Get this out of here!!!
struct MeshFace
{
    MeshFace()
    {
        vertices[0] = vertices[1] = vertices[2] = vertices[3] = -1;
    }
    
    MeshFace(int v1, int v2, int v3)
    {
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;
        vertices[3] = -1;
    }
    
    MeshFace(int v1, int v2, int v3, int v4)
    {
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;
        vertices[3] = v4;
    }
    
    int numVertices() const
    {
        for (int nn = 0; nn < 4; nn++)
        if (vertices[nn] == -1)
            return nn;
        return 4;
    }
    
    int vertices[4];
};

inline std::istream & operator>>(std::istream & str, MeshFace & face)
{
    int v[4];
    str >> v[0] >> v[1] >> v[2];
//    face = MeshFace(v[0], v[1], v[2]);
    if (str >> v[3])
        face = MeshFace(v[0], v[1], v[2], v[3]);
    else
    {
        str.clear(); // clear the error flags from failing str >> v[3]
        assert(str);
        face = MeshFace(v[0], v[1], v[2]);
    }
    return str;
}

#pragma mark *** Input/Output ***

template<typename T>
std::ostream & operator<<(std::ostream & str, const Rect<T> & rect);

template<typename T>
std::istream & operator>>(std::istream & str, Rect<T> & rect);

template<typename T>
std::ostream & operator<<(std::ostream & str, const Triangle3<T> & triangle);

template<typename T>
std::istream & operator>>(std::istream & str, Triangle3<T> & triangle);

template<typename T>
std::ostream & operator<<(std::ostream & str, const OrientedRect<T> & orect);

#include "geometry-inl.h"

#endif
