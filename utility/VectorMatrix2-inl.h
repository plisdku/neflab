/*
 *  VectorMatrix.cpp
 *  MyVectorMatrix
 *
 *  Created by Paul Hansen on 5/18/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 */

#ifdef _VECTORMATRIX2_

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
Vector2<T>::
Vector2()
{
    v[0] = v[1] = 0;
}

template<typename T>
Vector2<T>::
Vector2(const Vector2<T> & copyMe)
{
    v[0] = copyMe[0];
    v[1] = copyMe[1];
}

template<typename T>
template<typename T2>
Vector2<T>::
Vector2(const Vector2<T2> & copyMe )
{
    v[0] = copyMe[0];
    v[1] = copyMe[1];
}

template<typename T>
Vector2<T>::
Vector2(T v0, T v1)
{
    v[0] = v0;
    v[1] = v1;
}

template<typename T>
template<typename T2>
Vector2<T>::
Vector2( T2 v0, T2 v1)
{
    v[0] = v0;
    v[1] = v1;
}

template<typename T>
Vector2<T> Vector2<T>::
unit(unsigned int xyz)
{
    return Vector2<T>(xyz == 0, xyz == 1);
}

template<typename T, typename S>
Vector2<T> operator+(const Vector2<T> & lhs, const Vector2<S> & rhs)
{
    return Vector2<T>(lhs[0]+rhs[0], lhs[1]+rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator-(const Vector2<T> & lhs, const Vector2<S> & rhs)
{
    return Vector2<T>(lhs[0]-rhs[0], lhs[1]-rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator*(const Vector2<T> & lhs, const Vector2<S> & rhs)
{
    return Vector2<T>(lhs[0]*rhs[0], lhs[1]*rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator/(const Vector2<T> & lhs, const Vector2<S> & rhs)
{
    return Vector2<T>(lhs[0]/rhs[0], lhs[1]/rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator%(const Vector2<T> & lhs, const Vector2<S> & rhs)
{
    return Vector2<T>(lhs[0]%rhs[0], lhs[1]%rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator+(const Vector2<T> & lhs, S rhs)
{
    return Vector2<T>(lhs[0]+rhs, lhs[1]+rhs);
}

template<typename T, typename S>
Vector2<T> operator-(const Vector2<T> & lhs, S rhs)
{
    return Vector2<T>(lhs[0]-rhs, lhs[1]-rhs);
}

template<typename T, typename S>
Vector2<T> operator*(const Vector2<T> & lhs, S rhs)
{
    return Vector2<T>(lhs[0]*rhs, lhs[1]*rhs);
}

template<typename T, typename S>
Vector2<T> operator/(const Vector2<T> & lhs, S rhs)
{
    return Vector2<T>(lhs[0]/rhs, lhs[1]/rhs);
}

template<typename T, typename S>
Vector2<T> operator%(const Vector2<T> & lhs, S rhs)
{
    return Vector2<T>(lhs[0]%rhs, lhs[1]%rhs);
}


template<typename T, typename S>
Vector2<T> operator+(S lhs, const Vector2<T> & rhs)
{
    return Vector2<T>(lhs+rhs[0], lhs+rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator-(S lhs, const Vector2<T> & rhs)
{
    return Vector2<T>(lhs-rhs[0], lhs-rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator*(S lhs, const Vector2<T> & rhs)
{
    return Vector2<T>(lhs*rhs[0], lhs*rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator/(S lhs, const Vector2<T> & rhs)
{
    return Vector2<T>(lhs/rhs[0], lhs/rhs[1]);
}

template<typename T, typename S>
Vector2<T> operator%(S lhs, const Vector2<T> & rhs)
{
    return Vector2<T>(lhs%rhs[0], lhs%rhs[1]);
}


template<typename T>
Vector2<T> & Vector2<T>::
operator+=(const Vector2<T> & rhs)
{
    v[0] += rhs.v[0];
    v[1] += rhs.v[1];
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator-=(const Vector2<T> & rhs)
{
    v[0] -= rhs.v[0];
    v[1] -= rhs.v[1];
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator*=(const Vector2<T> & rhs)
{
    v[0] *= rhs.v[0];
    v[1] *= rhs.v[1];
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator/=(const Vector2<T> & rhs)
{
    v[0] /= rhs.v[0];
    v[1] /= rhs.v[1];
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator%=(const Vector2<T> & rhs)
{
    v[0] %= rhs.v[0];
    v[1] %= rhs.v[1];
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator+=(const T & rhs)
{
    v[0] += rhs;
    v[1] += rhs;
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator-=(const T & rhs)
{
    v[0] -= rhs;
    v[1] -= rhs;
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator*=(const T & rhs)
{
    v[0] *= rhs;
    v[1] *= rhs;
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator/=(const T & rhs)
{
    v[0] /= rhs;
    v[1] /= rhs;
    return *this;
}

template<typename T>
Vector2<T> & Vector2<T>::
operator%=(const T & rhs)
{
    v[0] %= rhs;
    v[1] %= rhs;
    return *this;
}

template<typename T>
Vector2<T> operator+(const Vector2<T> & rhs)
{
    return rhs;
}

template<typename T>
Vector2<T> operator-(const Vector2<T> & rhs)
{
    return Vector2<T>(-rhs[0], -rhs[1]);
}

template<typename T>
bool operator==(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1]);
}

template<typename T>
bool operator!=(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return !(lhs == rhs);
}

template<typename T>
bool operator<(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    if (lhs[0] < rhs[0]) return 1;
    else if (lhs[0] > rhs[0]) return 0;
    if (lhs[1] < rhs[1]) return 1;
    return 0;
}

template<typename T>
bool operator>(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    if (lhs[0] > rhs[0]) return 1;
    else if (lhs[0] < rhs[0]) return 0;
    if (lhs[1] > rhs[1]) return 1;
    return 0;
}

template<typename T>
bool operator<=(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return !(lhs > rhs);
}

template<typename T>
bool operator>=(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return !(lhs < rhs);
}

template<typename T>
Vector2<T> operator!(const Vector2<T> & v)
{
    return Vector2<T>(!v[0], !v[1]);
}


template<typename T>
T dot(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return lhs[0]*rhs[0] + lhs[1]*rhs[1];
}

template<typename T>
T cross(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return lhs[0]*rhs[1] - lhs[1]*rhs[0];
}

template<typename T>
T norm(const Vector2<T> & v)
{
    return sqrt(v[0]*v[0]+v[1]*v[1]);
}

template<typename T>
T norm1(const Vector2<T> & v)
{
    return smax(smax(sabs(v[0]), sabs(v[1])));
}

template<typename T>
T sumSquares(const Vector2<T> & v)
{
    return v[0]*v[0] + v[1]*v[1];
}

template<typename T>
Vector2<T> dominantComponent(const Vector2<T> & rhs)
{
    T a2, b2;
    a2 = rhs[0]*rhs[0];
    b2 = rhs[1]*rhs[1];
    
    if (a2 > b2)
        return Vector2<T>(rhs[0], 0);
    return Vector2<T>(0, rhs[1]);
}

template<typename T>
Vector2<T> cyclicPermute(const Vector2<T> & rhs, unsigned int nn)
{
    return Vector2<T>(rhs[(2-nn)%2], rhs[(3-nn)%2]);
}


template <typename T, typename S>
bool vec_eq(const Vector2<T> & lhs, const S & rhs)
{
    return (lhs[0] == rhs && lhs[1] == rhs);
}

template <typename T, typename S>
bool vec_lt(const Vector2<T>& lhs, const S & rhs)
{
    return (lhs[0] < rhs && lhs[1] < rhs);
}

template <typename T, typename S>
bool vec_gt(const Vector2<T>& lhs, const S & rhs)
{
    return (lhs[0] > rhs && lhs[1] > rhs);
}

template <typename T, typename S>
bool vec_le(const Vector2<T>& lhs, const S & rhs)
{
    return (lhs[0] <= rhs && lhs[1] <= rhs);
}

template <typename T, typename S>
bool vec_ge(const Vector2<T>& lhs, const S & rhs)
{
    return (lhs[0] >= rhs && lhs[1] >= rhs);
}

template <typename T, typename S>
bool vec_lt(const Vector2<T>& lhs, const Vector2<S>& rhs)
{
    return (lhs[0] < rhs[0] && lhs[1] < rhs[1]);
}

template <typename T, typename S>
bool vec_gt(const Vector2<T>& lhs, const Vector2<S>& rhs)
{
    return (lhs[0] > rhs[0] && lhs[1] > rhs[1]);
}

template <typename T, typename S>
bool vec_le(const Vector2<T>& lhs, const Vector2<S>& rhs)
{
    return (lhs[0] <= rhs[0] && lhs[1] <= rhs[1]);
}

template <typename T, typename S>
bool vec_ge(const Vector2<T>& lhs, const Vector2<S>& rhs)
{
    return (lhs[0] >= rhs[0] && lhs[1] >= rhs[1]);
}


template<typename T>
Vector2<T> vec_max(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return Vector2<T>(smax(lhs[0], rhs[0]), smax(lhs[1], rhs[1]));
}

template<typename T>
Vector2<T> vec_min(const Vector2<T> & lhs, const Vector2<T> & rhs)
{
    return Vector2<T>(smin(lhs[0], rhs[0]), smin(lhs[1], rhs[1]));
}

template<typename T>
Vector2<T> vec_max(const Vector2<T> & lhs, T rhs)
{
    return Vector2<T>(smax(lhs[0], rhs), smax(lhs[1], rhs));
}

template<typename T>
Vector2<T> vec_min(const Vector2<T> & lhs, T rhs)
{
    return Vector2<T>(smin(lhs[0], rhs), smin(lhs[1], rhs));
}

template<typename T>
Vector2<T> vec_floor(const Vector2<T> & lhs)
{
    return Vector2<T>(floor(lhs[0]), floor(lhs[1]));
}

template<typename T>
Vector2<T> vec_abs(const Vector2<T> & lhs)
{
    return Vector2<T>(sabs(lhs[0]), sabs(lhs[1]));
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const Vector2<T> & rhs)
{
    str << "(" << rhs[0] << ", " << rhs[1] << ")";
    return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Vector2<T> & rhs)
{
    str >> rhs[0] >> rhs[1];
    return str;
}




#pragma mark *** MATRIX ***



template<typename T>
Matrix2<T>::
Matrix2()
{
    m[0] = m[1] = m[2] = m[3] = 0;
}

template<typename T>
Matrix2<T>::
Matrix2(const Matrix2<T> & copyMe)
{
    m[0] = copyMe.m[0];
    m[1] = copyMe.m[1];
    m[2] = copyMe.m[2];
    m[3] = copyMe.m[3];
}

template<typename T>
template<typename T2>
Matrix2<T>::
Matrix2(const Matrix2<T2> & copyMe )
{
    m[0] = copyMe[0];
    m[1] = copyMe[1];
    m[2] = copyMe[2];
    m[3] = copyMe[3];
}

template<typename T>
template<typename S>
Matrix2<T>::
Matrix2(S m00, S m01, S m10, S m11)
{
    m[0] = m00;
    m[1] = m01;
    m[2] = m10;
    m[3] = m11;
}

template<typename T>
Matrix2<T> Matrix2<T>::
eye()
{
    return Matrix2<T>(1,0,0,1);
}

template<typename T>
template<typename T2>
Matrix2<T> Matrix2<T>::
withColumns(const Vector2<T2> & c1, const Vector2<T2> & c2)
{
    return Matrix2<T>(
        c1[0],c2[0],
        c1[1],c2[1]);
}

template<typename T>
template<typename T2>
Matrix2<T> Matrix2<T>::
withRows(const Vector2<T2> & c1, const Vector2<T2> & c2)
{
    return Matrix2<T>(
        c1[0],c1[1],
        c2[0],c2[1]);
}

template<typename T>
template<typename T2>
Matrix2<T> Matrix2<T>::
diagonal(const Vector2<T2> & d)
{
    return Matrix2<T>(
        d[0], 0L,
        0L, d[1]);
}

template<typename T>
template<typename T2>
Matrix2<T> Matrix2<T>::
diagonal(T2 d)
{
    return Matrix2<T>(
        d, 0,
        0, d);
}

template<typename T>
Matrix2<T> Matrix2<T>::
cyclicPermutation()
{
    return Matrix2<T>(
        0, 1,
        1, 0);
}

template<typename T>
Matrix2<T> transpose(const Matrix2<T> & rhs)
{
    return Matrix2<T>(
        rhs[0], rhs[2],
        rhs[1], rhs[3]);
}

template<typename T>
T determinant(const Matrix2<T> & rhs)
{
    return rhs[0]*rhs[3] - rhs[1]*rhs[2];
}

template<typename T>
Matrix2<T> inverse(const Matrix2<T> & rhs)
{
    return Matrix2<T>(
        rhs[3], -rhs[1],
        -rhs[2], rhs[0]) / determinant(rhs);
}

template<typename T, typename S>
Matrix2<T> operator+(const Matrix2<T> & lhs, const Matrix2<S> & rhs)
{
    return Matrix2<T>(lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2],
        lhs[3]+rhs[3]);
}

template<typename T, typename S>
Matrix2<T> operator-(const Matrix2<T> & lhs, const Matrix2<S> & rhs)
{
    return Matrix2<T>(lhs[0]-rhs[0], lhs[1]-rhs[1], lhs[2]-rhs[2],
        lhs[3]-rhs[3]);
}

template<typename T, typename S>
Matrix2<T> operator*(const Matrix2<T> & lhs, const Matrix2<S> & rhs)
{
    return Matrix2<T>(
        lhs[0]*rhs[0] + lhs[1]*rhs[2], lhs[0]*rhs[1] + lhs[1]*rhs[3],
        lhs[2]*rhs[0] + lhs[3]*rhs[2], lhs[2]*rhs[1] + lhs[3]*rhs[3]);
}

template<typename T, typename S>
Matrix2<T> operator%(const Matrix2<T> & lhs, const Matrix2<S> & rhs)
{
    return Matrix2<T>(lhs[0]%rhs[0], lhs[1]%rhs[1], lhs[2]%rhs[2],
        lhs[3]%rhs[3]);
}

template<typename T, typename S>
Vector2<S> operator*(const Matrix2<T> & lhs, const Vector2<S> & rhs)
{
    return Vector2<S>(
        lhs[0]*rhs[0] + lhs[1]*rhs[1],
        lhs[2]*rhs[0] + lhs[3]*rhs[1]);
}

template<typename T, typename S>
Vector2<S> operator*(const Vector2<S> & lhs, const Matrix2<T> & rhs)
{
    return Vector2<S>(
        lhs[0]*rhs[0] + lhs[1]*rhs[2],
        lhs[0]*rhs[1] + lhs[1]*rhs[3]);
}


template<typename T, typename S>
Matrix2<T> operator+(const Matrix2<T> & lhs, S & rhs)
{
    return Matrix2<T>(lhs[0]+rhs, lhs[1]+rhs, lhs[2]+rhs, lhs[3]+rhs);
}

template<typename T, typename S>
Matrix2<T> operator-(const Matrix2<T> & lhs, S & rhs)
{
    return Matrix2<T>(lhs[0]-rhs, lhs[1]-rhs, lhs[2]-rhs, lhs[3]-rhs);
}

template<typename T, typename S>
Matrix2<T> operator*(const Matrix2<T> & lhs, S rhs)
{
    return Matrix2<T>(lhs[0]*rhs, lhs[1]*rhs, lhs[2]*rhs, lhs[3]*rhs);
}

template<typename T, typename S>
Matrix2<T> operator/(const Matrix2<T> & lhs, S rhs)
{
    return Matrix2<T>(lhs[0]/rhs, lhs[1]/rhs, lhs[2]/rhs, lhs[3]/rhs);
}



template<typename T, typename S>
Matrix2<T> operator+(S lhs, const Matrix2<T> rhs)
{
    return Matrix2<T>(rhs[0]+lhs, rhs[1]+lhs, rhs[2]+lhs, rhs[3]+lhs);
}

template<typename T, typename S>
Matrix2<T> operator-(S lhs, const Matrix2<T> rhs)
{
    return Matrix2<T>(rhs[0]-lhs, rhs[1]-lhs, rhs[2]-lhs, rhs[3]-lhs);
}

template<typename T, typename S>
Matrix2<T> operator*(S lhs, const Matrix2<T> rhs)
{
    return Matrix2<T>(rhs[0]*lhs, rhs[1]*lhs, rhs[2]*lhs, rhs[3]*lhs);
}

template<typename T, typename S>
Matrix2<T> operator/(S lhs, const Matrix2<T> rhs)
{
    return Matrix2<T>(rhs[0]/lhs, rhs[1]/lhs, rhs[2]/lhs, rhs[3]/lhs);
}




template<typename T>
Matrix2<T> operator+(const Matrix2<T> & rhs)
{
    return rhs;
}

template<typename T>
Matrix2<T> operator-(const Matrix2<T> & rhs)
{
    return Matrix2<T>(-rhs[0], -rhs[1], -rhs[2], -rhs[3]);
}

template<typename T>
bool operator==(const Matrix2<T> & lhs, const Matrix2<T> & rhs)
{
    return (lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2] &&
        lhs[3] == rhs[3]);
}

template<typename T>
bool operator!=(const Matrix2<T> & lhs, const Matrix2<T> & rhs)
{
    return !(lhs == rhs);
}

template<typename T>
bool operator<(const Matrix2<T> & lhs, const Matrix2<T> & rhs)
{
    if (lhs[0] < rhs[0]) return 1;
    else if (lhs[0] > rhs[0]) return 0;
    if (lhs[1] < rhs[1]) return 1;
    else if (lhs[1] > rhs[1]) return 0;
    if (lhs[2] < rhs[2]) return 1;
    else if (lhs[2] > rhs[2]) return 0;
    if (lhs[3] < rhs[3]) return 1;
    return 0;
}

template<typename T>
bool operator>(const Matrix2<T> & lhs, const Matrix2<T> & rhs)
{
    if (lhs[0] > rhs[0]) return 1;
    else if (lhs[0] < rhs[0]) return 0;
    if (lhs[1] > rhs[1]) return 1;
    else if (lhs[1] < rhs[1]) return 0;
    if (lhs[2] > rhs[2]) return 1;
    else if (lhs[2] < rhs[2]) return 0;
    if (lhs[3] > rhs[3]) return 1;
    return 0;
}

template<typename T>
bool operator<=(const Matrix2<T> & lhs, const Matrix2<T> & rhs)
{
    return !(lhs > rhs);
}
template<typename T>
bool operator>=(const Matrix2<T> & lhs, const Matrix2<T> & rhs)
{
    return !(lhs < rhs);
}

template<typename T>
std::ostream & operator<<(std::ostream & str, const Matrix2<T> & rhs)
{
    str << "[ " << rhs[0] << ", " << rhs[1] << ",\n";
    str << "  " << rhs[2] << ", " << rhs[3] << "]";
    return str;
}

template<typename T>
std::istream & operator>>(std::istream & str, Matrix2<T> & rhs)
{
    str >> rhs[0] >> rhs[1] >> rhs[2] >> rhs[3];
    return str;
}


#endif