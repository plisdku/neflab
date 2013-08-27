/*
 *  VectorMatrix.cpp
 *  MyVectorMatrix
 *
 *  Created by Paul Hansen on 5/18/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 */

#ifdef _VECTORMATRIX_

#include <cassert>

#ifndef _SMAX_ETCETERA_
#define _SMAX_ETCETERA_
template<typename T>
static T smax(T lhs, T rhs)
{
    return rhs > lhs ? rhs : lhs;
}

template<typename T>
static T smin(T lhs, T rhs)
{
    return lhs < rhs ? lhs : rhs;
}

template<typename T>
static T sabs(T val)
{
    return val < 0 ? -val : val;
}
#endif

#pragma mark *** VECTOR ***

template<typename T>
Vector3<T>::
Vector3()
{
    v[0] = v[1] = v[2] = T();
}

template<typename T>
Vector3<T>::
Vector3(const Vector3<T> & copyMe)
{
    v[0] = copyMe[0];
    v[1] = copyMe[1];
    v[2] = copyMe[2];
}

template<typename T>
template<typename T2>
Vector3<T>::
Vector3(const Vector3<T2> & copyMe )
{
    v[0] = copyMe[0];
    v[1] = copyMe[1];
    v[2] = copyMe[2];
}

template<typename T>
Vector3<T>::
Vector3(T v0, T v1, T v2)
{
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
}

template<typename T>
template<typename T2>
Vector3<T>::
Vector3( T2 v0, T2 v1, T2 v2 )
{
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
}

template<typename T>
template<typename T2>
Vector3<T>::
Vector3( const T2 * vv )
{
    v[0] = vv[0];
    v[1] = vv[1];
    v[2] = vv[2];
}

template<typename T>
Vector3<T> Vector3<T>::
unit(unsigned int xyz)
{
    return Vector3<T>(xyz == 0, xyz == 1, xyz == 2);
}

template<typename T, typename S>
Vector3<T> operator+(const Vector3<T> & lhs, const Vector3<S> & rhs)
{
    return Vector3<T>(lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator-(const Vector3<T> & lhs, const Vector3<S> & rhs)
{
    return Vector3<T>(lhs[0]-rhs[0], lhs[1]-rhs[1], lhs[2]-rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator*(const Vector3<T> & lhs, const Vector3<S> & rhs)
{
    return Vector3<T>(lhs[0]*rhs[0], lhs[1]*rhs[1], lhs[2]*rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator/(const Vector3<T> & lhs, const Vector3<S> & rhs)
{
    return Vector3<T>(lhs[0]/rhs[0], lhs[1]/rhs[1], lhs[2]/rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator%(const Vector3<T> & lhs, const Vector3<S> & rhs)
{
    return Vector3<T>(lhs[0]%rhs[0], lhs[1]%rhs[1], lhs[2]%rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator+(const Vector3<T> & lhs, S rhs)
{
    return Vector3<T>(lhs[0]+rhs, lhs[1]+rhs, lhs[2]+rhs);
}

template<typename T, typename S>
Vector3<T> operator-(const Vector3<T> & lhs, S rhs)
{
    return Vector3<T>(lhs[0]-rhs, lhs[1]-rhs, lhs[2]-rhs);
}

template<typename T, typename S>
Vector3<T> operator*(const Vector3<T> & lhs, S rhs)
{
    return Vector3<T>(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs);
}

template<typename T, typename S>
Vector3<T> operator/(const Vector3<T> & lhs, S rhs)
{
    return Vector3<T>(lhs[0]/rhs, lhs[1]/rhs, lhs[2]/rhs);
}

template<typename T, typename S>
Vector3<T> operator%(const Vector3<T> & lhs, S rhs)
{
    return Vector3<T>(lhs[0]%rhs, lhs[1]%rhs, lhs[2]%rhs);
}


template<typename T, typename S>
Vector3<T> operator+(S lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(lhs+rhs[0], lhs+rhs[1], lhs+rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator-(S lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(lhs-rhs[0], lhs-rhs[1], lhs-rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator*(S lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(lhs*rhs[0], lhs*rhs[1], lhs*rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator/(S lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(lhs/rhs[0], lhs/rhs[1], lhs/rhs[2]);
}

template<typename T, typename S>
Vector3<T> operator%(S lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(lhs%rhs[0], lhs%rhs[1], lhs%rhs[2]);
}


template<typename T>
Vector3<T> & Vector3<T>::
operator+=(const Vector3<T> & rhs)
{
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    v[2] += rhs.v[2];
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator-=(const Vector3<T> & rhs)
{
    v[0] -= rhs.v[0];
    v[1] -= rhs.v[1];
    v[2] -= rhs.v[2];
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator*=(const Vector3<T> & rhs)
{
    v[0] *= rhs.v[0];
    v[1] *= rhs.v[1];
    v[2] *= rhs.v[2];
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator/=(const Vector3<T> & rhs)
{
    v[0] /= rhs.v[0];
    v[1] /= rhs.v[1];
    v[2] /= rhs.v[2];
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator%=(const Vector3<T> & rhs)
{
    v[0] %= rhs.v[0];
    v[1] %= rhs.v[1];
    v[2] %= rhs.v[2];
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator+=(const T & rhs)
{
    v[0] += rhs;
    v[1] += rhs;
    v[2] += rhs;
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator-=(const T & rhs)
{
    v[0] -= rhs;
    v[1] -= rhs;
    v[2] -= rhs;
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator*=(const T & rhs)
{
    v[0] *= rhs;
    v[1] *= rhs;
    v[2] *= rhs;
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator/=(const T & rhs)
{
    v[0] /= rhs;
    v[1] /= rhs;
    v[2] /= rhs;
    return *this;
}

template<typename T>
Vector3<T> & Vector3<T>::
operator%=(const T & rhs)
{
    v[0] %= rhs;
    v[1] %= rhs;
    v[2] %= rhs;
    return *this;
}

template<typename T>
Vector3<T> operator+(const Vector3<T> & rhs)
{
    return rhs;
}

template<typename T>
Vector3<T> operator-(const Vector3<T> & rhs)
{
    return Vector3<T>(-rhs[0], -rhs[1], -rhs[2]);
}

template<typename T>
bool operator==(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2]);
}

template<typename T>
bool operator!=(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return !(lhs == rhs);
}

template<typename T>
bool operator<(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    if (lhs[0] < rhs[0]) return 1;
    else if (lhs[0] > rhs[0]) return 0;
    if (lhs[1] < rhs[1]) return 1;
    else if (lhs[1] > rhs[1]) return 0;
    if (lhs[2] < rhs[2]) return 1;
    
    return 0;
    //else if (lhs[2] > rhs[2]) return 0;
}

template<typename T>
bool operator>(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    if (lhs[0] > rhs[0]) return 1;
    else if (lhs[0] < rhs[0]) return 0;
    if (lhs[1] > rhs[1]) return 1;
    else if (lhs[1] < rhs[1]) return 0;
    if (lhs[2] > rhs[2]) return 1;
    
    return 0;
    //else if (lhs[2] < rhs[2]) return 0;
}

template<typename T>
bool operator<=(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return !(lhs > rhs);
}

template<typename T>
bool operator>=(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return !(lhs < rhs);
}

template<typename T>
Vector3<T> operator!(const Vector3<T> & v)
{
    return Vector3<T>(!v[0], !v[1], !v[2]);
}


template<typename T>
T dot(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
}

template<typename T>
Vector3<T> cross(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(lhs[1]*rhs[2] - lhs[2]*rhs[1],
        lhs[2]*rhs[0] - lhs[0]*rhs[2],
        lhs[0]*rhs[1] - lhs[1]*rhs[0]);
}

template<typename T>
T norm(const Vector3<T> & v)
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

template<typename T>
T norm1(const Vector3<T> & v)
{
    return sabs(v[0]) + sabs(v[1]) + sabs(v[2]);
}

template<typename T>
T normInf(const Vector3<T> & v)
{
    return smax(smax(sabs(v[0]), sabs(v[1])), sabs(v[2]));
}

template<typename T>
T sumSquares(const Vector3<T> & v)
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


template<typename T>
Vector3<T> dominantComponent(const Vector3<T> & rhs)
{
    T a2, b2, c2;
    a2 = rhs[0]*rhs[0];
    b2 = rhs[1]*rhs[1];
    c2 = rhs[2]*rhs[2];
    
    if (a2 > b2 && a2 > c2)
        return Vector3<T>(rhs[0], 0, 0);
    else if (b2 > c2)
        return Vector3<T>(0, rhs[1], 0);
    return Vector3<T>(0, 0, rhs[2]);
}

template<typename T>
int dominantDirection(const Vector3<T> & rhs)
{
    T a2, b2, c2;
    a2 = rhs[0]*rhs[0];
    b2 = rhs[1]*rhs[1];
    c2 = rhs[2]*rhs[2];
    
    if (a2 > b2 && a2 > c2)
        return 0;
    else if (b2 > c2)
        return 1;
    else
        return 2;
}

template<typename T>
Vector3<T> cyclicPermute(const Vector3<T> & rhs, unsigned int nn)
{
    return Vector3<T>(rhs[(3-nn)%3], rhs[(4-nn)%3], rhs[(5-nn)%3]);
}

template<typename T>
Vector3<T> unit(const Vector3<T> & v)
{
    if (v[0] || v[1] || v[2])
        return Vector3<T>(v/norm(v));
    else
        return v;
}

template<typename T>
Vector3<T> floor(const Vector3<T> & v)
{
    return Vector3<T>(floor(v[0]), floor(v[1]), floor(v[2]));
}

template<typename T>
Vector3<T> ceil(const Vector3<T> & v)
{
    return Vector3<T>(floor(v[0]), floor(v[1]), floor(v[2]));
}

//template<typename T>
//Vector3<T> projection(const Vector3<T> & v, const Vector3<T> & direction)
//{
//    return unit(v)*dot(v, unit(direction));
//}


template <typename T, typename S>
bool vec_eq(const Vector3<T> & lhs, const S & rhs)
{
    return (lhs[0] == rhs && lhs[1] == rhs && lhs[2] == rhs);
}

template <typename T, typename S>
bool vec_lt(const Vector3<T>& lhs, const S & rhs)
{
    return (lhs[0] < rhs && lhs[1] < rhs && lhs[2] < rhs);
}

template <typename T, typename S>
bool vec_gt(const Vector3<T>& lhs, const S & rhs)
{
    return (lhs[0] > rhs && lhs[1] > rhs && lhs[2] > rhs);
}

template <typename T, typename S>
bool vec_le(const Vector3<T>& lhs, const S & rhs)
{
    return (lhs[0] <= rhs && lhs[1] <= rhs && lhs[2] <= rhs);
}

template <typename T, typename S>
bool vec_ge(const Vector3<T>& lhs, const S & rhs)
{
    return (lhs[0] >= rhs && lhs[1] >= rhs && lhs[2] >= rhs);
}

template <typename T, typename S>
bool vec_lt(const Vector3<T>& lhs, const Vector3<S>& rhs)
{
    return (lhs[0] < rhs[0] && lhs[1] < rhs[1] && lhs[2] < rhs[2]);
}

template <typename T, typename S>
bool vec_gt(const Vector3<T>& lhs, const Vector3<S>& rhs)
{
    return (lhs[0] > rhs[0] && lhs[1] > rhs[1] && lhs[2] > rhs[2]);
}

template <typename T, typename S>
bool vec_le(const Vector3<T>& lhs, const Vector3<S>& rhs)
{
    return (lhs[0] <= rhs[0] && lhs[1] <= rhs[1] && lhs[2] <= rhs[2]);
}

template <typename T, typename S>
bool vec_ge(const Vector3<T>& lhs, const Vector3<S>& rhs)
{
    return (lhs[0] >= rhs[0] && lhs[1] >= rhs[1] && lhs[2] >= rhs[2]);
}


template<typename T>
Vector3<T> vec_max(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(smax(lhs[0], rhs[0]), smax(lhs[1], rhs[1]),
        smax(lhs[2], rhs[2]));
}

template<typename T>
Vector3<T> vec_min(const Vector3<T> & lhs, const Vector3<T> & rhs)
{
    return Vector3<T>(smin(lhs[0], rhs[0]), smin(lhs[1], rhs[1]),
        smin(lhs[2], rhs[2]));
}

template<typename T>
Vector3<T> vec_max(const Vector3<T> & lhs, T rhs)
{
    return Vector3<T>(smax(lhs[0], rhs), smax(lhs[1], rhs), smax(lhs[2], rhs));
}

template<typename T>
Vector3<T> vec_min(const Vector3<T> & lhs, T rhs)
{
    return Vector3<T>(smin(lhs[0], rhs), smin(lhs[1], rhs), smin(lhs[2], rhs));
}

template<typename T>
Vector3<T> vec_floor(const Vector3<T> & lhs)
{
    return Vector3<T>(floor(lhs[0]), floor(lhs[1]), floor(lhs[2]));
}

template<typename T>
Vector3<T> vec_abs(const Vector3<T> & lhs)
{
    return Vector3<T>(sabs(lhs[0]), sabs(lhs[1]), sabs(lhs[2]));
}

template<class T>
T pointLineDistance(const Vector3<T> & p0, const Vector3<T> & p1,
    const Vector3<T> & test)
{
    Vector3<T> d0 = test - p0;
    Vector3<T> d1 = test - p1;
    
    // We need to be careful here.  Geometry dictates that our calculation will
    // give a distance >= 0, but floating point reality may interfere.
    
    T distAlongLine = dot(test-p0, unit(p1-p0));
    T difference = sumSquares(test-p0) - distAlongLine*distAlongLine;
    
    if (difference < 0)
        return 0;
    else
        return sqrt(sumSquares(test-p0) - distAlongLine*distAlongLine);
    
    //return distPerpToLine;
}



template<typename T>
std::ostream & operator<<(std::ostream & str, const Vector3<T> & rhs)
{
    str << "[" << rhs[0] << ", " << rhs[1] << ", " << rhs[2] << "]";
    return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Vector3<T> & rhs)
{
    str >> rhs[0] >> rhs[1] >> rhs[2];
    return str;
}




#pragma mark *** MATRIX ***



template<typename T>
Matrix3<T>::
Matrix3()
{
    m[0] = m[1] = m[2] = m[3] = m[4] = m[5] = m[6] = m[7] = m[8] = T();
}

template<typename T>
Matrix3<T>::
Matrix3(const Matrix3<T> & copyMe)
{
    m[0] = copyMe.m[0];
    m[1] = copyMe.m[1];
    m[2] = copyMe.m[2];
    m[3] = copyMe.m[3];
    m[4] = copyMe.m[4];
    m[5] = copyMe.m[5];
    m[6] = copyMe.m[6];
    m[7] = copyMe.m[7];
    m[8] = copyMe.m[8];
}

template<typename T>
template<typename T2>
Matrix3<T>::
Matrix3(const Matrix3<T2> & copyMe )
{
    m[0] = copyMe[0];
    m[1] = copyMe[1];
    m[2] = copyMe[2];
    m[3] = copyMe[3];
    m[4] = copyMe[4];
    m[5] = copyMe[5];
    m[6] = copyMe[6];
    m[7] = copyMe[7];
    m[8] = copyMe[8];
}

template<typename T>
template<typename S>
Matrix3<T>::
Matrix3(S m00, S m01, S m02, S m10, S m11, S m12, S m20, S m21, S m22)
{
    m[0] = m00;
    m[1] = m01;
    m[2] = m02;
    m[3] = m10;
    m[4] = m11;
    m[5] = m12;
    m[6] = m20;
    m[7] = m21;
    m[8] = m22;
}

template<typename T>
Matrix3<T> Matrix3<T>::
zero()
{
    return Matrix3<T>(0,0,0,0,0,0,0,0,0);
}

template<typename T>
Matrix3<T> Matrix3<T>::
eye()
{
    return Matrix3<T>(1,0,0,0,1,0,0,0,1);
}

template<typename T>
Matrix3<T> Matrix3<T>::
all(const T & val)
{
    return Matrix3<T>(val, val, val, val, val, val, val, val, val);
}

template<typename T>
template<typename T2>
Matrix3<T> Matrix3<T>::
withColumns(const Vector3<T2> & c1, const Vector3<T2> & c2,
    const Vector3<T> & c3)
{
    return Matrix3<T>(
        c1[0],c2[0],c3[0],
        c1[1],c2[1],c3[1],
        c1[2],c2[2],c3[2]);
}

template<typename T>
template<typename T2>
Matrix3<T> Matrix3<T>::
withRows(const Vector3<T2> & c1, const Vector3<T2> & c2,
    const Vector3<T> & c3)
{
    return Matrix3<T>(
        c1[0],c1[1],c1[2],
        c2[0],c2[1],c2[2],
        c3[0],c3[1],c3[2]);
}

template<typename T>
template<typename T2>
Matrix3<T> Matrix3<T>::
diagonal(const Vector3<T2> & d)
{
    return Matrix3<T>(
        d[0], T(0),    T(0),
        T(0),    d[1], T(0),
        T(0),    T(0),    d[2]);
}

template<typename T>
template<typename T2>
Matrix3<T> Matrix3<T>::
diagonal(T2 d)
{
    return Matrix3<T>(
        d, T(0), T(0),
        T(0), d, T(0),
        T(0), T(0), d);
}

template<typename T>
Matrix3<T> Matrix3<T>::
cyclicPermutation()
{
    return Matrix3<T>(
        T(0), T(0), T(1),
        T(1), T(0), T(0),
        T(0), T(1), T(0));
}

template<typename T>
Matrix3<T> Matrix3<T>::
cyclicPermutation(unsigned int n)
{
    return Matrix3<T>(
        T(n%3 == 0), T(n%3 == 2), T(n%3 == 1),
        T(n%3 == 1), T(n%3 == 0), T(n%3 == 2),
        T(n%3 == 2), T(n%3 == 1), T(n%3 == 0));
}

template<typename T>
Matrix3<T> Matrix3<T>::
rotation(T radians, const Vector3<T> & axis)
{
    Vector3<T> u0 = unit(axis);
    Vector3<T> u1 = cyclicPermute(axis, 1);
    Vector3<T> u2 = cyclicPermute(axis, 2);
    
    Matrix3<T> basis = Matrix3<T>::withColumns(u0, u1, u2);
    Matrix3<T> invBasis = transpose(basis);
    
    Matrix3<T> rotation2d(T(1.0), T(0.0), T(0.0),
        T(0.0), T(cos(radians)), T(-sin(radians)),
        T(0.0), T(sin(radians)), T(cos(radians)));
    
    return basis*rotation2d*invBasis;
}

template<typename T>
Vector3<T> Matrix3<T>::
row(unsigned int rr) const
{
    return Vector3<T>(m[3*rr], m[3*rr+1], m[3*rr+2]);
}

template<typename T>
Vector3<T> Matrix3<T>::
column(unsigned int cc) const
{
    return Vector3<T>(m[cc], m[3+cc], m[6+cc]);
}

template<typename T>
Vector3<T> Matrix3<T>::
diagonal() const
{
    return Vector3<T>(m[0], m[4], m[8]);
}

template<typename T>
void Matrix3<T>::
row(unsigned int rr, const Vector3<T> & newRow)
{
    (*this)(rr,0) = newRow[0];
    (*this)(rr,1) = newRow[1];
    (*this)(rr,2) = newRow[2];
}

template<typename T>
void Matrix3<T>::
column(unsigned int cc, const Vector3<T> & newCol)
{
    (*this)(0, cc) = newCol[0];
    (*this)(1, cc) = newCol[1];
    (*this)(2, cc) = newCol[2];
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator+=(const Matrix3<T> & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] += rhs.m[nn];
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator-=(const Matrix3<T> & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] -= rhs.m[nn];
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator*=(const Matrix3<T> & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] *= rhs.m[nn];
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator/=(const Matrix3<T> & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] /= rhs.m[nn];
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator%=(const Matrix3<T> & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] %= rhs.m[nn];
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator+=(const T & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] += rhs;
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator-=(const T & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] -= rhs;
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator*=(const T & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] *= rhs;
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator/=(const T & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] /= rhs;
    return *this;
}

template<typename T>
Matrix3<T> & Matrix3<T>::
operator%=(const T & rhs)
{
    for (int nn = 0; nn < 9; nn++)
        m[nn] %= rhs;
    return *this;
}

template<typename T>
Matrix3<T> transpose(const Matrix3<T> & rhs)
{
    return Matrix3<T>(
        rhs[0], rhs[3], rhs[6],
        rhs[1], rhs[4], rhs[7],
        rhs[2], rhs[5], rhs[8]);
}

template<typename T>
T determinant(const Matrix3<T> & rhs)
{
    return rhs[0]*(rhs[4]*rhs[8] - rhs[7]*rhs[5])
        - rhs[1]*(rhs[3]*rhs[8] - rhs[6]*rhs[5])
        + rhs[2]*(rhs[3]*rhs[7] - rhs[6]*rhs[4]);
}

template<typename T>
T trace(const Matrix3<T> & rhs)
{
    return rhs[0] + rhs[4] + rhs[8];
}

template<typename T>
Matrix3<T> inverse(const Matrix3<T> & rhs)
{
    // http://www.dr-lex.be/random/matrix_inv.html;
    return Matrix3<T>(
        rhs(2,2)*rhs(1,1)-rhs(2,1)*rhs(1,2),
        -(rhs(2,2)*rhs(0,1)-rhs(2,1)*rhs(0,2)),
        rhs(1,2)*rhs(0,1)-rhs(1,1)*rhs(0,2),
        -(rhs(2,2)*rhs(1,0)-rhs(2,0)*rhs(1,2)),
        rhs(2,2)*rhs(0,0)-rhs(2,0)*rhs(0,2),
        -(rhs(1,2)*rhs(0,0)-rhs(1,0)*rhs(0,2)),
        rhs(2,1)*rhs(1,0)-rhs(2,0)*rhs(1,1),
        -(rhs(2,1)*rhs(0,0)-rhs(2,0)*rhs(0,1)),
        rhs(1,1)*rhs(0,0)-rhs(1,0)*rhs(0,1)
        ) / determinant(rhs);
}

template<typename T>
Matrix3<T> outerProduct(const Vector3<T> & u, const Vector3<T> & v)
{
    return Matrix3<T>(
        u[0]*v[0], u[0]*v[1], u[0]*v[2],
        u[1]*v[0], u[1]*v[1], u[1]*v[2],
        u[2]*v[0], u[2]*v[1], u[2]*v[2]);
}

template<typename T>
T norm1(const Matrix3<T> & matrix)
{
    // maximum absolute column sum of the matrix
    return smax(smax(norm1(matrix.column(0)), norm1(matrix.column(1))),
        norm1(matrix.column(2)));
}

template<typename T>
T normInf(const Matrix3<T> & matrix)
{
    // maximum absolute row sum of the matrix
    return smax(smax(norm1(matrix.row(0)), norm1(matrix.row(1))),
        norm1(matrix.row(2)));
}

template<typename T>
T vectorNorm(const Matrix3<T> & matrix)
{
    T sum = 0.0;
    
    for (int nn = 0; nn < 9; nn++)
    {
        sum += matrix[nn]*matrix[nn];
    }
    
    return sqrt(sum);
}

template<typename T>
Matrix3<T> unitTrace(const Matrix3<T> & matrix)
{
    return matrix/trace(matrix);
}

template<typename T>
Matrix3<T> diagonalPartOf(const Matrix3<T> & matrix)
{
    return Matrix3<T>(
        matrix[0], T(0), T(0),
        T(0), matrix[4], T(0),
        T(0), T(0), matrix[8]);
}

template<typename T>
Matrix3<T> offDiagonalPartOf(const Matrix3<T> & matrix)
{
    return Matrix3<T>(
        T(0), matrix[1], matrix[2],
        matrix[3], T(0), matrix[5],
        matrix[6], matrix[7], T(0));
}

template<typename T>
Matrix3<T> operator+(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    return Matrix3<T>(lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2],
        lhs[3]+rhs[3], lhs[4]+rhs[4], lhs[5]+rhs[5], lhs[6]+rhs[6],
        lhs[7]+rhs[7], lhs[8]+rhs[8]);
}

template<typename T>
Matrix3<T> operator-(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    return Matrix3<T>(lhs[0]-rhs[0], lhs[1]-rhs[1], lhs[2]-rhs[2],
        lhs[3]-rhs[3], lhs[4]-rhs[4], lhs[5]-rhs[5], lhs[6]-rhs[6],
        lhs[7]-rhs[7], lhs[8]-rhs[8]);
}

template<typename T>
Matrix3<T> operator*(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    return Matrix3<T>(
        lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[6],
        lhs[0]*rhs[1]+lhs[1]*rhs[4]+lhs[2]*rhs[7],
        lhs[0]*rhs[2]+lhs[1]*rhs[5]+lhs[2]*rhs[8],
        lhs[3]*rhs[0]+lhs[4]*rhs[3]+lhs[5]*rhs[6],
        lhs[3]*rhs[1]+lhs[4]*rhs[4]+lhs[5]*rhs[7],
        lhs[3]*rhs[2]+lhs[4]*rhs[5]+lhs[5]*rhs[8],
        lhs[6]*rhs[0]+lhs[7]*rhs[3]+lhs[8]*rhs[6],
        lhs[6]*rhs[1]+lhs[7]*rhs[4]+lhs[8]*rhs[7],
        lhs[6]*rhs[2]+lhs[7]*rhs[5]+lhs[8]*rhs[8]);
        
}

template<typename T, typename S>
Matrix3<T> operator%(const Matrix3<T> & lhs, const Matrix3<S> & rhs)
{
    return Matrix3<T>(lhs[0]%rhs[0], lhs[1]%rhs[1], lhs[2]%rhs[2],
        lhs[3]%rhs[3], lhs[4]%rhs[4], lhs[5]%rhs[5], lhs[6]%rhs[6],
        lhs[7]%rhs[7], lhs[8]%rhs[8]);
}

template<typename T, typename S>
Vector3<S> operator*(const Matrix3<T> & lhs, const Vector3<S> & rhs)
{
    return Vector3<S>(
        lhs[0]*rhs[0]+lhs[1]*rhs[1]+lhs[2]*rhs[2],
        lhs[3]*rhs[0]+lhs[4]*rhs[1]+lhs[5]*rhs[2],
        lhs[6]*rhs[0]+lhs[7]*rhs[1]+lhs[8]*rhs[2]);
}

template<typename T, typename S>
Vector3<S> operator*(const Vector3<S> & lhs, const Matrix3<T> & rhs)
{
    return Vector3<S>(
        lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[6],
        lhs[0]*rhs[1]+lhs[1]*rhs[4]+lhs[2]*rhs[7],
        lhs[0]*rhs[2]+lhs[1]*rhs[5]+lhs[2]*rhs[8]);
}


template<typename T>
Matrix3<T> operator+(const Matrix3<T> & lhs, const T & rhs)
{
    return Matrix3<T>(lhs[0]+rhs, lhs[1]+rhs, lhs[2]+rhs, lhs[3]+rhs,
        lhs[4]+rhs, lhs[5]+rhs, lhs[6]+rhs, lhs[7]+rhs, lhs[8]+rhs);
}

template<typename T>
Matrix3<T> operator-(const Matrix3<T> & lhs, const T & rhs)
{
    return Matrix3<T>(lhs[0]-rhs, lhs[1]-rhs, lhs[2]-rhs, lhs[3]-rhs,
        lhs[4]-rhs, lhs[5]-rhs, lhs[6]-rhs, lhs[7]-rhs, lhs[8]-rhs);
}

template<typename T>
Matrix3<T> operator*(const Matrix3<T> & lhs, const T & rhs)
{
    return Matrix3<T>(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs, lhs[3]*rhs,
        lhs[4]*rhs, lhs[5]*rhs, lhs[6]*rhs, lhs[7]*rhs, lhs[8]*rhs);
}

template<typename T>
Matrix3<T> operator/(const Matrix3<T> & lhs, const T & rhs)
{
    return Matrix3<T>(lhs[0]/rhs, lhs[1]/rhs, lhs[2]/rhs, lhs[3]/rhs,
        lhs[4]/rhs, lhs[5]/rhs, lhs[6]/rhs, lhs[7]/rhs, lhs[8]/rhs);
}
template<typename T, typename S>
Matrix3<T> operator/(const Matrix3<T> & lhs, const S & rhs)
{
    return Matrix3<T>(lhs[0]/rhs, lhs[1]/rhs, lhs[2]/rhs, lhs[3]/rhs,
        lhs[4]/rhs, lhs[5]/rhs, lhs[6]/rhs, lhs[7]/rhs, lhs[8]/rhs);
}



template<typename T>
Matrix3<T> operator+(const T & lhs, const Matrix3<T> rhs)
{
    return Matrix3<T>(rhs[0]+lhs, rhs[1]+lhs, rhs[2]+lhs, rhs[3]+lhs,
        rhs[4]+lhs, rhs[5]+lhs, rhs[6]+lhs, rhs[7]+lhs, rhs[8]+lhs);
}

template<typename T>
Matrix3<T> operator-(const T & lhs, const Matrix3<T> rhs)
{
    return Matrix3<T>(rhs[0]-lhs, rhs[1]-lhs, rhs[2]-lhs, rhs[3]-lhs,
        rhs[4]-lhs, rhs[5]-lhs, rhs[6]-lhs, rhs[7]-lhs, rhs[8]-lhs);
}

template<typename T>
Matrix3<T> operator*(const T & lhs, const Matrix3<T> rhs)
{
    return Matrix3<T>(rhs[0]*lhs, rhs[1]*lhs, rhs[2]*lhs, rhs[3]*lhs,
        rhs[4]*lhs, rhs[5]*lhs, rhs[6]*lhs, rhs[7]*lhs, rhs[8]*lhs);
}

template<typename T>
Matrix3<T> operator/(const T & lhs, const Matrix3<T> rhs)
{
    return Matrix3<T>(rhs[0]/lhs, rhs[1]/lhs, rhs[2]/lhs, rhs[3]/lhs,
        rhs[4]/lhs, rhs[5]/lhs, rhs[6]/lhs, rhs[7]/lhs, rhs[8]/lhs);
}




template<typename T>
Matrix3<T> operator+(const Matrix3<T> & rhs)
{
    return rhs;
}

template<typename T>
Matrix3<T> operator-(const Matrix3<T> & rhs)
{
    return Matrix3<T>(-rhs[0], -rhs[1], -rhs[2], -rhs[3], -rhs[4], -rhs[5],
        -rhs[6], -rhs[7], -rhs[8]);
}

template<typename T>
bool operator==(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2] &&
        lhs[3] == rhs[3] && lhs[4] == rhs[4] && lhs[5] == rhs[5] &&
        lhs[6] == rhs[6] && lhs[7] == rhs[7] && lhs[8] == rhs[8]);
}

template<typename T>
bool operator!=(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    return !(lhs == rhs);
}

template<typename T>
bool operator<(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    if (lhs[0] < rhs[0]) return 1;
    else if (lhs[0] > rhs[0]) return 0;
    if (lhs[1] < rhs[1]) return 1;
    else if (lhs[1] > rhs[1]) return 0;
    if (lhs[2] < rhs[2]) return 1;
    else if (lhs[2] > rhs[2]) return 0;
    if (lhs[3] < rhs[3]) return 1;
    else if (lhs[3] > rhs[3]) return 0;
    if (lhs[4] < rhs[4]) return 1;
    else if (lhs[4] > rhs[4]) return 0;
    if (lhs[5] < rhs[5]) return 1;
    else if (lhs[5] > rhs[5]) return 0;
    if (lhs[6] < rhs[6]) return 1;
    else if (lhs[6] > rhs[6]) return 0;
    if (lhs[7] < rhs[7]) return 1;
    else if (lhs[7] > rhs[7]) return 0;
    if (lhs[8] < rhs[8]) return 1;
    return 0;
    //else if (lhs[8] > rhs[8]) return 0;
}

template<typename T>
bool operator>(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    if (lhs[0] > rhs[0]) return 1;
    else if (lhs[0] < rhs[0]) return 0;
    if (lhs[1] > rhs[1]) return 1;
    else if (lhs[1] < rhs[1]) return 0;
    if (lhs[2] > rhs[2]) return 1;
    else if (lhs[2] < rhs[2]) return 0;
    if (lhs[3] > rhs[3]) return 1;
    else if (lhs[3] < rhs[3]) return 0;
    if (lhs[4] > rhs[4]) return 1;
    else if (lhs[4] < rhs[4]) return 0;
    if (lhs[5] > rhs[5]) return 1;
    else if (lhs[5] < rhs[5]) return 0;
    if (lhs[6] > rhs[6]) return 1;
    else if (lhs[6] < rhs[6]) return 0;
    if (lhs[7] > rhs[7]) return 1;
    else if (lhs[7] < rhs[7]) return 0;
    if (lhs[8] > rhs[8]) return 1;
    return 0;
}

template<typename T>
bool operator<=(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    return !(lhs > rhs);
}
template<typename T>
bool operator>=(const Matrix3<T> & lhs, const Matrix3<T> & rhs)
{
    return !(lhs < rhs);
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const Matrix3<T> & rhs)
{
    str << "[ " << rhs[0] << ", " << rhs[1] << ", " << rhs[2] << ",\n";
    str << "  " << rhs[3] << ", " << rhs[4] << ", " << rhs[5] << ",\n";
    str << "  " << rhs[6] << ", " << rhs[7] << ", " << rhs[8] << "]";
    return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Matrix3<T> & rhs)
{
    str >> rhs[0] >> rhs[1] >> rhs[2] >> rhs[3] >> rhs[4] >> rhs[5]
        >> rhs[6] >> rhs[7] >> rhs[8];
    return str;
}


#endif