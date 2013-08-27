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

#ifndef _VECTORMATRIX_
#define _VECTORMATRIX_


#include <iostream>
#include <cmath>

template<typename T>
class Vector3;
typedef Vector3<long> Vector3i;
typedef Vector3<double> Vector3d;
typedef Vector3<float> Vector3f;
typedef Vector3<bool> Vector3b;

template<typename T>
class Matrix3;
typedef Matrix3<long> Matrix3i;
typedef Matrix3<double> Matrix3d;
typedef Matrix3<float> Matrix3f;
typedef Matrix3<bool> Matrix3b;


template<typename T>
class Vector3
{
public:
    Vector3();
    Vector3(const Vector3<T> & copyMe);
    
    template<typename T2>
    explicit Vector3( const Vector3<T2> & copyMe );
    
    Vector3(T v0, T v1, T v2);
    
    static Vector3<T> unit(unsigned int xyz);
    
    template<typename T2>
    explicit Vector3( T2 v0, T2 v1, T2 v2 );
    
    template<typename T2>
    explicit Vector3(const T2 * vv);
    
    T & operator[](unsigned int n)
        { return v[n]; }
    const T & operator[](unsigned int n) const
        { return v[n]; }
    
    const T* asArray() const { return v; }
    T* asArray() { return v; }
    
    Vector3<T> & operator+=(const Vector3<T> & rhs);
    Vector3<T> & operator-=(const Vector3<T> & rhs);
    Vector3<T> & operator*=(const Vector3<T> & rhs);
    Vector3<T> & operator/=(const Vector3<T> & rhs);
    Vector3<T> & operator%=(const Vector3<T> & rhs);
    
    Vector3<T> & operator+=(const T & rhs);
    Vector3<T> & operator-=(const T & rhs);
    Vector3<T> & operator*=(const T & rhs);
    Vector3<T> & operator/=(const T & rhs);
    Vector3<T> & operator%=(const T & rhs);
    
    struct SortAlong
    {
    public:
        SortAlong() : mAxis(0) {}
        SortAlong(unsigned int axis) : mAxis(axis) {}
        
        bool operator()(const Vector3<T> & lhs, const Vector3<T> & rhs) const
        {
            return lhs[mAxis] < rhs[mAxis];
        }
        
        unsigned int mAxis;
    };
    
    struct IsNonzero
    {
    public:
        bool operator()(const Vector3<T> & v) const
        {
            return v[0] != 0 || v[1] != 0 || v[2] != 0;
        }
    };
    
    struct GetElement
    {
        GetElement(int i): m_i(i) {}
        
        const T & operator()(const Vector3<T> & vec) const
        {
            return vec[m_i];
        }
        
        int m_i;
    };
    
private:
    T v[3];
};


template<typename T, typename S>
Vector3<T> operator+(const Vector3<T> & lhs, const Vector3<S> & rhs);
template<typename T, typename S>
Vector3<T> operator-(const Vector3<T> & lhs, const Vector3<S> & rhs);
template<typename T, typename S>
Vector3<T> operator*(const Vector3<T> & lhs, const Vector3<S> & rhs);
template<typename T, typename S>
Vector3<T> operator/(const Vector3<T> & lhs, const Vector3<S> & rhs);
template<typename T, typename S>
Vector3<T> operator%(const Vector3<T> & lhs, const Vector3<S> & rhs);

template<typename T, typename S>
Vector3<T> operator+(const Vector3<T> & lhs, S rhs);
template<typename T, typename S>
Vector3<T> operator-(const Vector3<T> & lhs, S rhs);
template<typename T, typename S>
Vector3<T> operator*(const Vector3<T> & lhs, S rhs);
template<typename T, typename S>
Vector3<T> operator/(const Vector3<T> & lhs, S rhs);
template<typename T, typename S>
Vector3<T> operator%(const Vector3<T> & lhs, S rhs);

template<typename T, typename S>
Vector3<T> operator+(S lhs, const Vector3<T> & rhs);
template<typename T, typename S>
Vector3<T> operator-(S lhs, const Vector3<T> & rhs);
template<typename T, typename S>
Vector3<T> operator*(S lhs, const Vector3<T> & rhs);
template<typename T, typename S>
Vector3<T> operator/(S lhs, const Vector3<T> & rhs);
template<typename T, typename S>
Vector3<T> operator%(S lhs, const Vector3<T> & rhs);

template<typename T>
Vector3<T> operator+(const Vector3<T> & rhs);
template<typename T>
Vector3<T> operator-(const Vector3<T> & rhs);
template<typename T>
bool operator==(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
bool operator!=(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
bool operator<(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
bool operator>(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
bool operator<=(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
bool operator>=(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
Vector3<T> operator!(const Vector3<T> & v);

template<typename T>
T dot(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
Vector3<T> cross(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
T norm(const Vector3<T> & v); // this is the 2-norm
template<typename T>
T norm1(const Vector3<T> & v); // this is the 1-norm, the sum of absolute vals
template<typename T>
T normInf(const Vector3<T> & v); // this is the inf-norm, the max absolute val
template<typename T>
T sumSquares(const Vector3<T> & v);
template<typename T>
Vector3<T> dominantComponent(const Vector3<T> & rhs);
template<typename T>
int dominantDirection(const Vector3<T> & rhs);
template<typename T>
Vector3<T> cyclicPermute(const Vector3<T> & rhs, unsigned int nn);
template<typename T>
Vector3<T> unit(const Vector3<T> & v);
//template<typename T>
//Vector3<T> projection(const Vector3<T> & v, const Vector3<T> direction);
template<typename T>
Vector3<T> floor(const Vector3<T> & v);
template<typename T>
Vector3<T> ceil(const Vector3<T> & v);

template <typename T, typename S>
bool vec_eq(const Vector3<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_lt(const Vector3<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_gt(const Vector3<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_le(const Vector3<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_ge(const Vector3<T>& lhs, const S & rhs);
template <typename T, typename S>
bool vec_lt(const Vector3<T>& lhs, const Vector3<S>& rhs);
template <typename T, typename S>
bool vec_gt(const Vector3<T>& lhs, const Vector3<S>& rhs);
template <typename T, typename S>
bool vec_le(const Vector3<T>& lhs, const Vector3<S>& rhs);
template <typename T, typename S>
bool vec_ge(const Vector3<T>& lhs, const Vector3<S>& rhs);


template<typename T>
Vector3<T> vec_max(const Vector3<T> & lhs, const Vector3<T> & rhs);
template<typename T>
Vector3<T> vec_min(const Vector3<T> & lhs, const Vector3<T> & rhs); 
template<typename T>
Vector3<T> vec_max(const Vector3<T> & lhs, T rhs);
template<typename T>
Vector3<T> vec_min(const Vector3<T> & lhs, T rhs); 

template<typename T>
Vector3<T> vec_floor(const Vector3<T> & lhs);
template<typename T>
Vector3<T> vec_abs(const Vector3<T> & lhs);

template<class T>
T pointLineDistance(const Vector3<T> & p0, const Vector3<T> & p1,
    const Vector3<T> & test);

template<typename T>
std::ostream & operator<<(std::ostream & str, const Vector3<T> & rhs);

template<typename T>
std::istream & operator>>(std::istream & str, Vector3<T> & rhs);


template<typename T>
class Matrix3
{
public:
    Matrix3();
    
    Matrix3(const Matrix3<T> & copyMe);
    
    template<typename T2>
    explicit Matrix3(const Matrix3<T2> & copyMe );
    
    template<typename S>
    Matrix3(S m00, S m01, S m02, S m10, S m11, S m12, S m20, S m21, S m22);
    
    static Matrix3<T> zero();
    static Matrix3<T> eye();
    static Matrix3<T> all(const T & val);
    
    template<typename T2>
    static Matrix3<T> withColumns(const Vector3<T2> & c1,
        const Vector3<T2> & c2, const Vector3<T> & c3);
    
    template<typename T2>
    static Matrix3<T> withRows(const Vector3<T2> & c1, const Vector3<T2> & c2,
        const Vector3<T> & c3);
    
    template<typename T2>
    static Matrix3<T> diagonal(const Vector3<T2> & d);
    
    template<typename T2>
    static Matrix3<T> diagonal(T2 d);
    
    static Matrix3<T> cyclicPermutation();
    static Matrix3<T> cyclicPermutation(unsigned int n);
    static Matrix3<T> rotation(T radians, const Vector3<T> & axis);
    
    T & operator[](unsigned int n)
        { return m[n]; }
    const T & operator[](unsigned int n) const
        { return m[n]; }
    T & operator()(unsigned int mm, unsigned int nn)
        { return m[3*mm+nn]; }
    const T & operator()(unsigned int mm, unsigned int nn) const
        { return m[3*mm+nn]; }
    
    struct GetElement
    {
        GetElement(int i, int j): m_i(i), m_j(j) {}
        
        const T & operator()(const Matrix3<T> & matrix) const
        {
            return matrix(m_i, m_j);
        }
        
        int m_i;
        int m_j;
    };
    
    Vector3<T> row(unsigned int rr) const;
    Vector3<T> column(unsigned int cc) const;
    Vector3<T> diagonal() const;
    
    void row(unsigned int rr, const Vector3<T> & newRow);
    void column(unsigned int rr, const Vector3<T> & newCol);
    
    Matrix3<T> & operator+=(const Matrix3<T> & rhs);
    Matrix3<T> & operator-=(const Matrix3<T> & rhs);
    Matrix3<T> & operator*=(const Matrix3<T> & rhs);
    Matrix3<T> & operator/=(const Matrix3<T> & rhs);
    Matrix3<T> & operator%=(const Matrix3<T> & rhs);
    
    Matrix3<T> & operator+=(const T & rhs);
    Matrix3<T> & operator-=(const T & rhs);
    Matrix3<T> & operator*=(const T & rhs);
    Matrix3<T> & operator/=(const T & rhs);
    Matrix3<T> & operator%=(const T & rhs);
private:
    T m[9];
};

template<typename T>
Matrix3<T> transpose(const Matrix3<T> & lhs);
template<typename T>
T determinant(const Matrix3<T> & lhs);
template<typename T>
T trace(const Matrix3<T> & lhs);
template<typename T>
Matrix3<T> inverse(const Matrix3<T> & lhs);
template<typename T>
Matrix3<T> outerProduct(const Vector3<T> & u, const Vector3<T> & v);

template<typename T>
T norm1(const Matrix3<T> & matrix);
template<typename T>
T normInf(const Matrix3<T> & matrix);
template<typename T>
T vectorNorm(const Matrix3<T> & matrix);

template<typename T>
Matrix3<T> unitTrace(const Matrix3<T> & matrix);

template<typename T>
Matrix3<T> diagonalPartOf(const Matrix3<T> & matrix);
template<typename T>
Matrix3<T> offDiagonalPartOf(const Matrix3<T> & matrix);

template<typename T>
Matrix3<T> operator+(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
Matrix3<T> operator-(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
Matrix3<T> operator*(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
Matrix3<T> operator%(const Matrix3<T> & lhs, const Matrix3<T> & rhs);

template<typename T, typename S>
Vector3<S> operator*(const Matrix3<T> & lhs, const Vector3<S> & rhs);

// this is like trans(v)*M, but no conjugation or nuthin' happens here...
template<typename T, typename S>
Vector3<S> operator*(const Vector3<S> & lhs, const Matrix3<T> & rhs);

template<typename T>
Matrix3<T> operator+(const Matrix3<T> & lhs, const T & rhs);
template<typename T>
Matrix3<T> operator-(const Matrix3<T> & lhs, const T & rhs);
template<typename T>
Matrix3<T> operator*(const Matrix3<T> & lhs, const T & rhs);
template<typename T>
Matrix3<T> operator/(const Matrix3<T> & lhs, const T & rhs);
template<typename T, typename S>
Matrix3<T> operator/(const Matrix3<T> & lhs, const S & rhs);

template<typename T>
Matrix3<T> operator+(const T & lhs, const Matrix3<T> rhs);
template<typename T>
Matrix3<T> operator-(const T & lhs, const Matrix3<T> rhs);
template<typename T>
Matrix3<T> operator*(const T & lhs, const Matrix3<T> rhs);
template<typename T>
Matrix3<T> operator/(const T & lhs, const Matrix3<T> rhs);

template<typename T>
Matrix3<T> operator+(const Matrix3<T> & rhs);
template<typename T>
Matrix3<T> operator-(const Matrix3<T> & rhs);
template<typename T>
bool operator==(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
bool operator!=(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
bool operator<(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
bool operator>(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
bool operator<=(const Matrix3<T> & lhs, const Matrix3<T> & rhs);
template<typename T>
bool operator>=(const Matrix3<T> & lhs, const Matrix3<T> & rhs);

template<typename T>
std::ostream & operator<<(std::ostream & str, const Matrix3<T> & rhs);

template<typename T>
std::istream & operator>>(std::istream & str, Matrix3<T> & rhs);



#include "VectorMatrix-inl.h"
#endif
