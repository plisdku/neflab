/*
 *  Polynomial.h
 *  Snapdragon
 *
 *  Created by Paul Hansen on 4/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifndef _POLYNOMIAL_
#define _POLYNOMIAL_


template<class T = float>
class Polynomial
{
public:
    typedef T Real;
    Polynomial(); // default polynomial = 0
    Polynomial(const Polynomial<T> & copyMe);
    Polynomial(Real constantValue);
    Polynomial(Real c0, Real c1);
    Polynomial(Real c0, Real c1, Real c2);
    Polynomial(Real c0, Real c1, Real c2, Real c3);
    
    Real operator[](unsigned int order) const { return coefficient(order); }
    Real operator()(Real x) const { return value(x); }
    
    Real coefficient(unsigned int whichOrder) const;
    void coefficient(unsigned int whichOrder, Real newValue);
    Real value(Real x) const;
    unsigned int order() const { return mOrder; }
    //void order(unsigned int newOrder); // does not change the value
    
    Polynomial<T> & operator=(const Polynomial<T> & rhs);
    Polynomial<T> & operator+=(const Polynomial<T> & rhs);
    Polynomial<T> & operator-=(const Polynomial<T> & rhs);
    Polynomial<T> & operator*=(const Polynomial<T> & rhs);
    
    static unsigned int maximumOrder() { return kMaximumPolynomialOrder; }
protected:
    enum { kMaximumPolynomialOrder = 20 };
    Real mCoefficients[kMaximumPolynomialOrder+1];
    unsigned int mOrder;
};

// test equality of order and all coefficients up to order()
template<class T>
bool operator==(const Polynomial<T> & p1, const Polynomial<T> & p2);

template<class T>
bool operator!=(const Polynomial<T> & p1, const Polynomial<T> & p2);


template<class T>
bool operator==(const Polynomial<T> & lhs, const T & rhs);
template<class T>
bool operator==(const Polynomial<T> & lhs, const T & rhs);

template<class T>
bool operator==(const T & lhs, const Polynomial<T> & rhs);
template<class T>
bool operator==(const T & lhs, const Polynomial<T> & rhs);

// lexicographic ordering
//bool operator<(const Polynomial & p1, const Polynomial & p2);
//bool operator<=(const Polynomial & p1, const Polynomial & p2);
//bool operator>(const Polynomial & p1, const Polynomial & p2);
//bool operator>=(const Polynomial & p1, const Polynomial & p2);

template<class T>
Polynomial<T> operator+(const Polynomial<T> & lhs, const Polynomial<T> & rhs);
template<class T>
Polynomial<T> operator-(const Polynomial<T> & lhs, const Polynomial<T> & rhs);
template<class T>
Polynomial<T> operator*(const Polynomial<T> & lhs, const Polynomial<T> & rhs);

template<class T>
Polynomial<T> operator+(const Polynomial<T> & lhs, const T & rhs);
template<class T>
Polynomial<T> operator+(const T & lhs, const Polynomial<T> & rhs);
template<class T>
Polynomial<T> operator-(const Polynomial<T> & lhs, const T & rhs);
template<class T>
Polynomial<T> operator-(const T & lhs, const Polynomial<T> & rhs);
template<class T>
Polynomial<T> operator*(const Polynomial<T> & lhs, const T & rhs);
template<class T>
Polynomial<T> operator*(const T & lhs, const Polynomial<T> & rhs);
template<class T>
Polynomial<T> operator/(const Polynomial<T> & lhs, const T & rhs);

template<class T>
Polynomial<T> operator-(const Polynomial<T> & p);

template<class T, class STREAM>
STREAM & operator<<(STREAM & str, const Polynomial<T> & rhs);

#include "Polynomial-inl.h"

#endif
