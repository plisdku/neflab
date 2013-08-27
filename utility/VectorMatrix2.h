/*
 *  VectorMatrix.h
 *  MyVectorMatrix
 *
 *  Created by Paul Hansen on 5/18/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  THESE VECTORS ETC. DO NOT BEHAVE CORRECTLY FOR COMPLEX NUMBERS.
 *  There are no complex conjugations anywhere, e.g. in dot products and in
 *  absolute values/norms.
 */

#ifndef _VECTORMATRIX2_
#define _VECTORMATRIX2_

#include <iostream>
#include <cmath>

template<typename T>
class Vector2;
typedef Vector2<long> Vector2i;
typedef Vector2<double> Vector2d;
typedef Vector2<float> Vector2f;
typedef Vector2<bool> Vector2b;

template<typename T>
class Matrix2;
typedef Matrix2<long> Matrix2i;
typedef Matrix2<double> Matrix2d;
typedef Matrix2<float> Matrix2f;
typedef Matrix2<bool> Matrix2b;


template<typename T>
class Vector2
{
public:
    Vector2();
    Vector2(const Vector2<T> & copyMe);
    
    template<typename T2>
    explicit Vector2( const Vector2<T2> & copyMe );
    
    Vector2(T v0, T v1);
    
    static Vector2<T> unit(unsigned int xyz);
    
    template<typename T2>
    explicit Vector2( T2 v0, T2 v1 );
    
    T & operator[](unsigned int n)
        { return v[n]; }
    const T & operator[](unsigned int n) const
        { return v[n]; }
        
    const T* asArray() const { return v; }
    T* asArray() { return v; }
    
    Vector2<T> & operator+=(const Vector2<T> & rhs);
    Vector2<T> & operator-=(const Vector2<T> & rhs);
    Vector2<T> & operator*=(const Vector2<T> & rhs);
    Vector2<T> & operator/=(const Vector2<T> & rhs);
    Vector2<T> & operator%=(const Vector2<T> & rhs);
    
    Vector2<T> & operator+=(const T & rhs);
    Vector2<T> & operator-=(const T & rhs);
    Vector2<T> & operator*=(const T & rhs);
    Vector2<T> & operator/=(const T & rhs);
    Vector2<T> & operator%=(const T & rhs);
private:
    T v[2];
};


template<typename T, typename S>
Vector2<T> operator+(const Vector2<T> & lhs, const Vector2<S> & rhs);
template<typename T, typename S>
Vector2<T> operator-(const Vector2<T> & lhs, const Vector2<S> & rhs);
template<typename T, typename S>
Vector2<T> operator*(const Vector2<T> & lhs, const Vector2<S> & rhs);
template<typename T, typename S>
Vector2<T> operator/(const Vector2<T> & lhs, const Vector2<S> & rhs);
template<typename T, typename S>
Vector2<T> operator%(const Vector2<T> & lhs, const Vector2<S> & rhs);

template<typename T, typename S>
Vector2<T> operator+(const Vector2<T> & lhs, S rhs);
template<typename T, typename S>
Vector2<T> operator-(const Vector2<T> & lhs, S rhs);
template<typename T, typename S>
Vector2<T> operator*(const Vector2<T> & lhs, S rhs);
template<typename T, typename S>
Vector2<T> operator/(const Vector2<T> & lhs, S rhs);
template<typename T, typename S>
Vector2<T> operator%(const Vector2<T> & lhs, S rhs);

template<typename T, typename S>
Vector2<T> operator+(S lhs, const Vector2<T> & rhs);
template<typename T, typename S>
Vector2<T> operator-(S lhs, const Vector2<T> & rhs);
template<typename T, typename S>
Vector2<T> operator*(S lhs, const Vector2<T> & rhs);
template<typename T, typename S>
Vector2<T> operator/(S lhs, const Vector2<T> & rhs);
template<typename T, typename S>
Vector2<T> operator%(S lhs, const Vector2<T> & rhs);

template<typename T>
Vector2<T> operator+(const Vector2<T> & rhs);
template<typename T>
Vector2<T> operator-(const Vector2<T> & rhs);
template<typename T>
bool operator==(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
bool operator!=(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
bool operator<(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
bool operator>(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
bool operator<=(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
bool operator>=(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
Vector2<T> operator!(const Vector2<T> & v);

template<typename T>
T dot(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
T cross(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
T norm(const Vector2<T> & v);
template<typename T>
T norm1(const Vector2<T> & v);
template<typename T>
T sumSquares(const Vector2<T> & v);
template<typename T>
Vector2<T> dominantComponent(const Vector2<T> & rhs);
template<typename T>
Vector2<T> cyclicPermute(const Vector2<T> & rhs, unsigned int nn);

template <typename T, typename S>
bool vec_eq(const Vector2<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_lt(const Vector2<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_gt(const Vector2<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_le(const Vector2<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_ge(const Vector2<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_lt(const Vector2<T>& lhs, const Vector2<S>& rhs);
template <typename T, typename S>
bool vec_gt(const Vector2<T>& lhs, const Vector2<S>& rhs);
template <typename T, typename S>
bool vec_le(const Vector2<T>& lhs, const Vector2<S>& rhs);
template <typename T, typename S>
bool vec_ge(const Vector2<T>& lhs, const Vector2<S>& rhs);


template<typename T>
Vector2<T> vec_max(const Vector2<T> & lhs, const Vector2<T> & rhs);
template<typename T>
Vector2<T> vec_min(const Vector2<T> & lhs, const Vector2<T> & rhs); 
template<typename T>
Vector2<T> vec_max(const Vector2<T> & lhs, T rhs);
template<typename T>
Vector2<T> vec_min(const Vector2<T> & lhs, T rhs); 

template<typename T>
Vector2<T> vec_floor(const Vector2<T> & lhs);
template<typename T>
Vector2<T> vec_abs(const Vector2<T> & lhs);

template<typename T>
std::ostream & operator<<(std::ostream & str, const Vector2<T> & rhs);

template<typename T>
std::istream & operator>>(std::istream & str, Vector2<T> & rhs);


template<typename T>
class Matrix2
{
public:
    Matrix2();
    Matrix2(const Matrix2<T> & copyMe);
    
    template<typename T2>
    explicit Matrix2(const Matrix2<T2> & copyMe );
    
    template<typename S>
    Matrix2(S m00, S m01, S m10, S m11);
    
    static Matrix2<T> eye();
    
    template<typename T2>
    static Matrix2<T> withColumns(const Vector2<T2> & c1,
        const Vector2<T2> & c2);
    
    template<typename T2>
    static Matrix2<T> withRows(const Vector2<T2> & c1, const Vector2<T2> & c2);
    
    template<typename T2>
    static Matrix2<T> diagonal(const Vector2<T2> & d);
    
    template<typename T2>
    static Matrix2<T> diagonal(T2 d);
    
    static Matrix2<T> cyclicPermutation();
    
    
    T & operator[](unsigned int n)
        { return m[n]; }
    const T & operator[](unsigned int n) const
        { return m[n]; }
    T & operator()(unsigned int mm, unsigned int nn)
        { return m[2*mm+nn]; }
    const T & operator()(unsigned int mm, unsigned int nn) const
        { return m[2*mm+nn]; }
    
private:
    T m[4];
};

template<typename T>
Matrix2<T> transpose(const Matrix2<T> & lhs);
template<typename T>
T determinant(const Matrix2<T> & lhs);
template<typename T>
Matrix2<T> inverse(const Matrix2<T> & lhs);

template<typename T, typename S>
Matrix2<T> operator+(const Matrix2<T> & lhs, const Matrix2<S> & rhs);
template<typename T, typename S>
Matrix2<T> operator-(const Matrix2<T> & lhs, const Matrix2<S> & rhs);
template<typename T, typename S>
Matrix2<T> operator*(const Matrix2<T> & lhs, const Matrix2<S> & rhs);
template<typename T, typename S>
Matrix2<T> operator%(const Matrix2<T> & lhs, const Matrix2<S> & rhs);

template<typename T, typename S>
Vector2<S> operator*(const Matrix2<T> & lhs, const Vector2<S> & rhs);

// this is like trans(v)*M, but no conjugation or nuthin' happens here...
template<typename T, typename S>
Vector2<S> operator*(const Vector2<S> & lhs, const Matrix2<T> & rhs);

template<typename T, typename S>
Matrix2<T> operator+(const Matrix2<T> & lhs, S & rhs);
template<typename T, typename S>
Matrix2<T> operator-(const Matrix2<T> & lhs, S & rhs);
template<typename T, typename S>
Matrix2<T> operator*(const Matrix2<T> & lhs, S rhs);
template<typename T, typename S>
Matrix2<T> operator/(const Matrix2<T> & lhs, S rhs);

template<typename T, typename S>
Matrix2<S> operator+(S lhs, const Matrix2<T> rhs);
template<typename T, typename S>
Matrix2<S> operator-(S lhs, const Matrix2<T> rhs);
template<typename T, typename S>
Matrix2<S> operator*(S lhs, const Matrix2<T> rhs);
template<typename T, typename S>
Matrix2<S> operator/(S lhs, const Matrix2<T> rhs);

template<typename T>
Matrix2<T> operator+(const Matrix2<T> & rhs);
template<typename T>
Matrix2<T> operator-(const Matrix2<T> & rhs);
template<typename T>
bool operator==(const Matrix2<T> & lhs, const Matrix2<T> & rhs);
template<typename T>
bool operator!=(const Matrix2<T> & lhs, const Matrix2<T> & rhs);
template<typename T>
bool operator<(const Matrix2<T> & lhs, const Matrix2<T> & rhs);
template<typename T>
bool operator>(const Matrix2<T> & lhs, const Matrix2<T> & rhs);
template<typename T>
bool operator<=(const Matrix2<T> & lhs, const Matrix2<T> & rhs);
template<typename T>
bool operator>=(const Matrix2<T> & lhs, const Matrix2<T> & rhs);

template<typename T>
std::ostream & operator<<(std::ostream & str, const Matrix2<T> & rhs);

template<typename T>
std::istream & operator>>(std::istream & str, Matrix2<T> & rhs);



#include "VectorMatrix2-inl.h"
#endif
