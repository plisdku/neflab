/*
 *  geometry.cpp
 *  vmlibtest
 *
 *  Created by Paul Hansen on 2/8/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef _GEOMETRY_

#pragma mark *** Rect ***


template<typename T>
Rect<T>::Rect() :
    p1(0,0,0),
    p2(0,0,0)
{
}

template<typename T>
Rect<T>::Rect(T i0, T j0, T k0, T i1, T j1, T k1) :
    p1(i0, j0, k0),
    p2(i1, j1, k1)
{
}

template<typename T>
Rect<T>::Rect( const Vector3<T> & inP1,
    const Vector3<T> & inP2 ) :
	p1(inP1), p2(inP2)
{
}

//template<typename T>
//Rect<T>::Rect(const Rect<T> & copyMe) :
//    p1(copyMe.p1),
//    p2(copyMe.p2)
//{
//}

template<typename T>
template<typename S>
Rect<T>::Rect(const Rect<S> & copyMe) :
    p1(copyMe.p1),
    p2(copyMe.p2)
{
}

template<typename T>
T Rect<T>::size(unsigned int dim) const
{
	return p2[dim] - p1[dim];
}

template<typename T>
Vector3<T> Rect<T>::size() const
{
	return p2 - p1 ;
}

template<typename T>
T Rect<T>::volume() const
{
    return size(0)*size(1)*size(2);
}

template<typename T>
T Rect<T>::num(unsigned int dim) const
{
	return p2[dim] - p1[dim] + 1;
}

template<typename T>
Vector3<T> Rect<T>::num() const
{
	return p2 - p1 + 1;
}

template<typename T>
T Rect<T>::count() const
{
    return num(0)*num(1)*num(2);
}


template<typename T>
bool Rect<T>::encloses(const Rect<T> & inRect) const
{
    return (encloses(inRect.p1) && encloses(inRect.p2));
}

template<typename T>
bool Rect<T>::encloses(const Vector3<T> & inPt) const
{
    if ( inPt[0] >= p1[0] && inPt[0] <= p2[0] &&
         inPt[1] >= p1[1] && inPt[1] <= p2[1] &&
         inPt[2] >= p1[2] && inPt[2] <= p2[2] )
         return 1;
    return 0;
}

template<typename T>
bool Rect<T>::encloses(T x, T y, T z) const
{
    if ( x >= p1[0] && x <= p2[0] &&
         y >= p1[1] && y <= p2[1] &&
         z >= p1[2] && z <= p2[2] )
         return 1;
    return 0;
}

template<typename T>
bool Rect<T>::encloses(const Triangle3<T> & triangle) const
{
    return encloses(triangle[0]) &&
        encloses(triangle[1]) &&
        encloses(triangle[2]);
}

template<typename T>
bool Rect<T>::intersects(const Rect<T> & inRect) const
{
	return (vec_ge(p2, inRect.p1) && vec_le(p1, inRect.p2));
}

template<typename T>
bool Rect<T>::intersects(const Triangle3<T> & tri) const
{
    // Use the Separating Axis Theorem (SAT) to test for intersection.
    // Test thirteen candidate separation axes:
    //  the x, y and z directions (surface normals of the rectangle)
    //  the triangle surface normal
    //  all the normal vectors of supporting planes of a rectangle edge and a
    //   triangle edge; three rectangle edge directions and three triangle edges
    //   multiply to a total of nine more candidate axes.
    
    // First three axes: the rectangle's face normals, a.k.a. the XYZ axes.
    for (int xyz = 0; xyz < 3; xyz++)
    {
        if (!projectionAlong(Vector3d::unit(0)).intersects(
            tri.projectionAlong(Vector3d::unit(0))))
        {
            return false;
        }
    }
    
    // Fourth axis: the triangle normal.
    Vector3d triNormal = tri.normal();
    if (!projectionAlong(triNormal).intersects(tri.projectionAlong(triNormal)))
        return false;
    
    // Nine more axes!
    for (int tt = 0; tt < 3; tt++) // edges of the triangle
    for (int rr = 0; rr < 3; rr++) // directions of rectangle edges
    {
        Vector3d axis(cross(tri.edgeDirection(tt), Vector3d::unit(rr)));
        
        if (!projectionAlong(axis).intersects(tri.projectionAlong(axis)))
            return false;
    }
    
    return true;
}

template<typename T>
int Rect<T>::numNonSingularDims() const
{
    int ndims = 0;
    if (p1[0] != p2[0])
        ndims++;
    if (p1[1] != p2[1])
        ndims++;
    if (p1[2] != p2[2])
        ndims++;
    
    return ndims;
}

template<typename T>
ClosedInterval<T> Rect<T>::
projectionAlong(const Vector3<T> & v) const
{
    Vector3<T> corner1, corner2;
    
    for (int nn = 0; nn < 3; nn++)
    {
        if (v[nn] > 0)
        {
            corner1[nn] = std::min(p1[nn], p2[nn]);
            corner2[nn] = std::max(p1[nn], p2[nn]);
        }
        else
        {
            corner1[nn] = std::max(p2[nn], p1[nn]);
            corner2[nn] = std::min(p2[nn], p1[nn]);
        }   
    }
    
    return ClosedInterval<T>(dot(corner1, v), dot(corner2, v));
}


#pragma mark *** Other functions ***


template<typename T>
Rect<T>
operator * (const Rect<T> & lhs, const T rhs)
{
    return Rect<T>(lhs.p1*rhs, lhs.p2*rhs);
}

template<typename T>
Rect<T>
operator / (const Rect<T> & lhs, const T rhs)
{
    return Rect<T>(lhs.p1/rhs, lhs.p2/rhs);
}

template<typename T>
Rect<T>
operator * (const T scalar, const Rect<T> & rhs)
{
    return Rect<T>(rhs.p1*scalar, rhs.p2*scalar);
}

template<typename T>
Rect<T> &
operator *= (Rect<T> & lhs, const T scalar)
{
    lhs.p1 *= scalar;
    lhs.p2 *= scalar;
    return lhs;
}

template<typename T>
Rect<T> &
operator /= (Rect<T> & lhs, const T scalar)
{
    lhs.p1 /= scalar;
    lhs.p2 /= scalar;
    return lhs;
}


template<typename T>
Rect<T>
operator * (const Matrix3<T> & lhs, const Rect<T> & rhs)
{
    return Rect<T>(lhs*rhs.p1, lhs*rhs.p2);
}


template<typename T>
Rect<T>
operator + (const Rect<T> & lhs, const Vector3<T> & rhs)
{
	return Rect<T>(lhs.p1 + rhs, lhs.p2 + rhs);
}

template<typename T>
Rect<T>
operator - (const Rect<T> & lhs, const Vector3<T> & rhs)
{
	return Rect<T>(lhs.p1 - rhs, lhs.p2 - rhs);
}

template<typename T>
Rect<T>
operator + (const Vector3<T> & lhs, const Rect<T> & rhs)
{
	return Rect<T>(rhs.p1 + lhs, rhs.p2 + lhs);
}

template<typename T>
Rect<T>
operator - (const Vector3<T> & lhs, const Rect<T> & rhs)
{
	return Rect<T>(lhs - rhs.p1, lhs - rhs.p2);
}



template<typename T, typename S>
Rect<T>
inset( const Rect<T> & inRect, S dx0, S dy0, S dz0, S dx1, S dy1, S dz1 )
{
    Rect<T> out(inRect);
    out.p1[0] += dx0;
    out.p1[1] += dy0;
    out.p1[2] += dz0;
    out.p2[0] -= dx1;
    out.p2[1] -= dy1;
    out.p2[2] -= dz1;
    return out;
}

template<typename T, typename S>
Rect<T>
inset( const Rect<T> & inRect, S distance )
{
    return inset(inRect, distance, distance, distance, distance, distance,
        distance);
}
                   
template<typename T>
Vector3<T>  clip(const Rect<T> & clipRect, const Vector3<T> & v)
{
    Vector3<T> out(v);
    if (out[0] < clipRect.p1[0])
        out[0] = clipRect.p1[0];
    else if (out[0] > clipRect.p2[0])
        out[0] = clipRect.p2[0];
    
    if (out[1] < clipRect.p1[1])
        out[1] = clipRect.p1[1];
    else if (out[1] > clipRect.p2[1])
        out[1] = clipRect.p2[1];
    
    if (out[2] < clipRect.p1[2])
        out[2] = clipRect.p1[2];
    else if (out[2] > clipRect.p2[2])
        out[2] = clipRect.p2[2];
    
    return out;
}

template<typename T>
Rect<T> clip( const Rect<T> & clipRect, const Rect<T> & rectToClip)
{
	Rect<T> out( clip(clipRect, rectToClip.p1), clip(clipRect, rectToClip.p2) );
	return out;
}

template<typename T>
Rect<T> intersection( const Rect<T> & r1, const Rect<T> & r2)
{
    return Rect<T>(vec_max(r1.p1, r2.p1), vec_min(r1.p2, r2.p2));
}

template<typename T>
Rect<T> cyclicPermute(const Rect<T> & r, unsigned int nn)
{
    return Rect<T>(cyclicPermute(r.p1,nn), cyclicPermute(r.p2,nn));
}

template<typename T>
bool operator<(const Rect<T> & lhs, const Rect<T> & rhs)
{
    if (lhs.p1 < rhs.p1)
        return 1;
    else if (lhs.p1 > rhs.p1)
        return 0;
    else if (lhs.p2 < rhs.p2)
        return 1;
    return 0;
}


template<typename T>
bool operator>(const Rect<T> & lhs, const Rect<T> & rhs)
{
    if (lhs.p1 > rhs.p1)
        return 1;
    else if (lhs.p1 < rhs.p1)
        return 0;
    else if (lhs.p2 > rhs.p2)
        return 1;
    return 0;
}

#pragma mark *** Plane ***

template<typename T>
Plane3<T>::
Plane3()
{
}

template<typename T>
Plane3<T>::
Plane3(const Vector3<T> & normal, const T & constant) :
    mNormal(normal),
    mConstant(constant)
{
}

template<typename T>
Vector3<T>
intersection(const Plane3<T> & plane, const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection)
{
    double t = (-plane.constant() - dot(rayOrigin, plane.normal())) / 
        dot(rayDirection, plane.normal());
    return t*rayDirection + rayOrigin;
}

template<typename T>
T
intersectionTime(const Plane3<T> & plane, const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection)
{
    double t = (-plane.constant() - dot(rayOrigin, plane.normal())) / 
        dot(rayDirection, plane.normal());
    return t;
}

template<typename T>
Plane3<T>
intersectionTimeJacobian(const Plane3<T> & plane, const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection)
{
    double nDotDir = dot(plane.normal(), rayDirection);
    double nDotOrigin = dot(plane.normal(), rayOrigin);
    return Plane3<T>(
        ( rayDirection*(plane.constant() + nDotOrigin) - rayOrigin*nDotDir )
            / (nDotDir*nDotDir),
        -1/nDotDir);
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const Plane3<T> & plane)
{
    str << "(" << plane.normal() << ", " << plane.constant() << ")";
    return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Plane3<T> & plane)
{
    Vector3<T> v;
    T c;
    str >> v >> c;
    plane = Plane3<T>(v,c);
    return str;
}

#pragma mark *** Polygon ***

template<typename T>
Polygon3<T>::
Polygon3()
{
}

template<typename T>
Polygon3<T>::
Polygon3(const std::vector<Vector3<T>*> & outerBoundary,
    const std::vector<std::vector<Vector3<T>*> > & holes) :
    mOuterBoundary(outerBoundary),
    mHoles(holes)
{
}

template<typename T>
Polygon3<T>::
Polygon3(const std::vector<std::vector<Vector3<T>*> > & contours)
{
    if (contours.size() > 0)
        mOuterBoundary = contours[0];
    if (contours.size() > 1)
    {
        mHoles.resize(contours.size()-1);
        std::copy(contours.begin() + 1, contours.end(), mHoles.begin());
    }
}


#pragma mark *** Triangle ***

template<typename T>
Vector3<T> Triangle3<T>::
normal() const
{
    return cross(mVertices[1] - mVertices[0], mVertices[2] - mVertices[0]);
}

template<typename T>
Plane3<T> Triangle3<T>::
plane() const
{
    Vector3<T> n = normal();
    return Plane3<T>(n, -dot(n, mVertices[0]));
}

template<typename T>
Vector3<T> Triangle3<T>::
edgeDirection(int nn) const
{
    return Vector3<T>(mVertices[(nn+1)%3] - mVertices[nn]);
}

template<typename T>
bool Triangle3<T>::
enclosesProjection(const Vector3<T> & point) const
{
    Vector3<T> ab = cross(mVertices[1]-mVertices[0], point-mVertices[0]);
    Vector3<T> bc = cross(mVertices[2]-mVertices[1], point-mVertices[1]);
    Vector3<T> ca = cross(mVertices[0]-mVertices[2], point-mVertices[2]);
    
    return (dot(ab,bc) >= 0 && dot(bc, ca) >= 0);
}

template<typename T>
ClosedInterval<T> Triangle3<T>::
projectionAlong(const Vector3<T> & v) const
{
    T dist0(dot(v, mVertices[0]));
    T dist1(dot(v, mVertices[1]));
    T dist2(dot(v, mVertices[2]));
    
    return ClosedInterval<T>(std::min(std::min(dist0, dist1), dist2),
        std::max(std::max(dist0, dist1), dist2));
}

template<typename T>
Matrix3<T> Triangle3<T>::
dndv(unsigned int vertex) const
{
    const Vector3<T> u1 = mVertices[1] - mVertices[0];
    const Vector3<T> u2 = mVertices[2] - mVertices[0];
    
    const Vector3<T> x(1,0,0), y(0,1,0), z(0,0,1);
    
    if (vertex == 0)
    {
        return Matrix3<T>::withColumns(
            -cross(x, u2) - cross(u1, x),
            -cross(y, u2) - cross(u1, y),
            -cross(z, u2) - cross(u1, z) );
    }
    else if (vertex == 1)
    {
        return Matrix3<T>::withColumns(
            cross(x, u2), cross(y, u2), cross(z, u2));
    }
    else
    {
        return Matrix3<T>::withColumns(
            cross(u1, x), cross(u1, y), cross(u1, z));
    }
}

template<typename T>
Vector3<T> Triangle3<T>::
dcdv(unsigned int vertex) const
{
    if (vertex == 0)
    {
        return -mVertices[0]*dndv(vertex) - normal();
    }
    else if (vertex == 1)
    {
        return -mVertices[0]*dndv(vertex);
    }
    else
    {
        return -mVertices[0]*dndv(vertex);
    }
}

template<typename T>
Matrix3<T> Triangle3<T>::
dudv(unsigned int vertex) const
{
    Vector3d normalVector(normal());
    double norm_n = norm(normalVector);
    Matrix3d frontMatrix = Matrix3d::diagonal(1.0/norm_n) -
        (1.0/norm_n/norm_n/norm_n)*
        outerProduct(normalVector, normalVector);
    
    return frontMatrix*dndv(vertex);
//    jacobian[0] = frontMatrix*dndv(0);
//    jacobian[1] = frontMatrix*dndv(1);
//    jacobian[2] = frontMatrix*dndv(2);
}

template<typename T>
Vector3<T> Triangle3<T>::
intersectionCoefficients(const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection) const
{
    Matrix3<T> inverseBasis(inverse(Matrix3<T>::withColumns(
        mVertices[1] - mVertices[0],
        mVertices[2] - mVertices[0],
        -rayDirection)));
    
    return inverseBasis*(rayOrigin - mVertices[0]);
}

template<typename T>
Matrix3<T> Triangle3<T>::
intersectionCoefficientsJacobian(const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection, int vertex) const
{
    Matrix3<T> inverseBasis(inverse(Matrix3<T>::withColumns(
        mVertices[1] - mVertices[0],
        mVertices[2] - mVertices[0],
        -rayDirection)));
    
    Vector3<T> coefficients = inverseBasis*(rayOrigin - mVertices[0]);
    
    if (vertex == 0)
    {
        return (coefficients[0]+coefficients[1]-1)*inverseBasis;
    }
    else if (vertex == 1)
    {
        return -coefficients[0]*inverseBasis;
    }
    else
    {
        return -coefficients[1]*inverseBasis;
    }
}

template<typename T>
Matrix3<T> Triangle3<T>::
intersectionJacobian(const Vector3<T> & rayOrigin,
    const Vector3<T> & rayDirection, int vertex) const
{
    return outerProduct(rayDirection,
        intersectionCoefficientsJacobian(rayOrigin, rayDirection, vertex)
        .row(2));
}

template<typename T>
bool Triangle3<T>::
operator==(const Triangle3<T> & rhs) const
{
    return mVertices[0] == rhs.mVertices[0] &&
        mVertices[1] == rhs.mVertices[1] &&
        mVertices[2] == rhs.mVertices[2];
}


template<typename T>
Triangle3<T>
operator * (const Matrix3<T> & lhs, const Triangle3<T> & rhs)
{
    return Triangle3<T>(lhs*rhs[0], lhs*rhs[1], lhs*rhs[2]);
}

template<typename T>
Triangle3<T>
operator + (const Triangle3<T> & lhs, const Vector3<T> & rhs)
{
    return rhs + lhs;
}

template<typename T>
Triangle3<T>
operator + (const Vector3<T> & lhs, const Triangle3<T> & rhs)
{
    return Triangle3<T>(lhs + rhs[0], lhs + rhs[1], lhs + rhs[2]);
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const Triangle3<T> & triangle)
{
	str << "[" << triangle[0] << "; " << triangle[1] << "; "
        << triangle[2] << "]";
	return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Triangle3<T> & triangle)
{
	str >> triangle[0] >> triangle[1] >> triangle[2];
    
	return str;
}

#pragma mark *** OrientedRect ***

template<typename T>
OrientedRect<T> cyclicPermute(const OrientedRect<T> & r, unsigned int nn)
{
    return OrientedRect<T>(cyclicPermute(r.rect,nn),
        cyclicPermute(r.normal,nn));
}

#pragma mark *** Input/Output ***

template<typename T>
std::ostream & operator<<(std::ostream & str, const Rect<T> & rect)
{
	str << "[" << rect.p1 << ", " << rect.p2 << "]";
	return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Rect<T> & rect)
{
	str >> rect.p1[0] >> rect.p1[1] >> rect.p1[2]
		>> rect.p2[0] >> rect.p2[1] >> rect.p2[2];
	return str;
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const OrientedRect<T> & orect)
{
	str << orect.rect << ", " << orect.normal;
	return str;
}

#endif
