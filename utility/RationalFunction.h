/*
 *  RationalFunction.h
 *  Snapdragon
 *
 *  Created by Paul Hansen on 4/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifndef _RATIONALFUNCTION_
#define _RATIONALFUNCTION_

#include "Polynomial.h"

template<class T>
class RationalFunction
{
public:
    typedef T Real;
    typedef Polynomial<T> Poly;
    RationalFunction();
    RationalFunction(Real constantValue);
    RationalFunction(const Polynomial<T> & poly);
    RationalFunction(const Polynomial<T> & numer, const Polynomial<T> & denom);
    RationalFunction(const RationalFunction<T> & copyMe);
    
    static RationalFunction<T> irregular(const Polynomial<T> & numerator,
        const Polynomial<T> & denominator);
    
    Polynomial<T> & numerator() { return mNumerator; }
    Polynomial<T> & denominator() { return mDenominator; }
    const Polynomial<T> & numerator() const { return mNumerator; }
    const Polynomial<T> & denominator() const { return mDenominator; }
    
    Real value(Real x) const;
    Real operator()(Real x) const { return value(x); }
    unsigned int order() const;
    
    RationalFunction<T> & operator=(const RationalFunction<T> & rhs);
    RationalFunction<T> & operator+=(const RationalFunction<T> & rhs);
    RationalFunction<T> & operator-=(const RationalFunction<T> & rhs);
    RationalFunction<T> & operator*=(const RationalFunction<T> & rhs);
    RationalFunction<T> & operator/=(const RationalFunction<T> & rhs);
    
    static unsigned int maximumOrder() { return Polynomial<T>::maximumOrder(); }
private:
    Polynomial<T> mNumerator;
    Polynomial<T> mDenominator;
};

// OPERATIONS: RationalFunction x RationalFunction
template<class T>
bool operator==(const RationalFunction<T> & lhs, const RationalFunction<T> & rhs);
template<class T>
bool operator!=(const RationalFunction<T> & lhs, const RationalFunction<T> & rhs);

template<class T>
RationalFunction<T> operator+(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs);
template<class T>
RationalFunction<T> operator-(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs);
template<class T>
RationalFunction<T> operator*(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs);
template<class T>
RationalFunction<T> operator/(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs);


// RationalFunction x Scalar
template<class T>
bool operator==(const RationalFunction<T> & lhs, const T & rhs);
template<class T>
bool operator!=(const RationalFunction<T> & lhs, const T & rhs);

template<class T>
RationalFunction<T> operator+(const RationalFunction<T> & lhs, const T & rhs);
template<class T>
RationalFunction<T> operator-(const RationalFunction<T> & lhs, const T & rhs);
template<class T>
RationalFunction<T> operator*(const RationalFunction<T> & lhs, const T & rhs);
template<class T>
RationalFunction<T> operator/(const RationalFunction<T> & lhs, const T & rhs);

// Scalar x RationalFunction
template<class T>
bool operator==(const T & lhs, const RationalFunction<T> & rhs);
template<class T>
bool operator!=(const T & lhs, const RationalFunction<T> & rhs);

template<class T>
RationalFunction<T> operator+(const T & lhs, const RationalFunction<T> & rhs);
template<class T>
RationalFunction<T> operator-(const T & lhs, const RationalFunction<T> & rhs);
template<class T>
RationalFunction<T> operator*(const T & lhs, const RationalFunction<T> & rhs);
template<class T>
RationalFunction<T> operator/(const T & lhs, const RationalFunction<T> & rhs);


// Unary
template<class T>
RationalFunction<T> operator-(const RationalFunction<T> & r);

template<class T>
RationalFunction<T> reciprocal(const RationalFunction<T> & r);

template<class T>
const Polynomial<T> & numerator(const RationalFunction<T> & r)
{
    return r.numerator();
}

template<class T>
const Polynomial<T> & denominator(const RationalFunction<T> & r)
{
    return r.denominator();
}

struct NumeratorTerm
{
public:
    NumeratorTerm(int order) : mOrder(order) {}
    
    template<class T>
    T operator()(const RationalFunction<T> & r)
    {
        return r.numerator()[mOrder];
    }
    int mOrder;
};

struct DenominatorTerm
{
public:
    DenominatorTerm(int order) : mOrder(order) {}
    
    template<class T>
    T operator()(const RationalFunction<T> & r)
    {
        return r.denominator()[mOrder];
    }
    int mOrder;
};

template<class T>
struct EvaluateAt
{
public:
    EvaluateAt(T arg) : mArgument(arg) {}
    
    T operator()(const RationalFunction<T> & r)
    {
        return r(mArgument);
    }
    
    T mArgument;
};

struct IsOrder
{
public:
    IsOrder(int numerOrder, int denomOrder) :
        mNumer(numerOrder), mDenom(denomOrder)
    {
    }
    
    template<class T>
    bool operator()(const RationalFunction<T> & r) const
    {
        if ( (mNumer == -1 || mNumer == r.numerator().order()) &&
             (mDenom == -1 || mDenom == r.denominator().order()) )
        {
            return true;
        }
        return false;
    }
    
    int mNumer, mDenom;
};

struct OrderGreaterThan
{
public:
    OrderGreaterThan(int minNumerOrder, int minDenomOrder) :
        mNumer(minNumerOrder), mDenom(minDenomOrder)
    {
    }
    
    template<class T>
    bool operator()(const RationalFunction<T> & r) const
    {
        if (mNumer == -1)
            return (r.denominator().order() >= mDenom);
        else if (mDenom == -1)
            return (r.numerator().order() >= mNumer);
        else if (r.numerator().order() >= mNumer &&
            r.denominator().order() >= mDenom)
        {
            return true;
        }
        return false;
    }
    
    int mNumer, mDenom;
};

struct SelectOrder
{
public:
    SelectOrder(int numeratorOrder, int denominatorOrder) :
        mNumer(numeratorOrder), mDenom(denominatorOrder)
    {
    }
    
    template<class T>
    RationalFunction<T> operator()(const RationalFunction<T> & r) const
    {
        if (mNumer == -1 && r.denominator().order() == mDenom)
            return r;
        else if (mDenom == -1 && r.numerator().order() == mNumer)
            return r;
        else if (r.numerator().order() == mNumer &&
            r.denominator().order() == mDenom)
        {
            return r;
        }
        else
            return RationalFunction<T>(0.0);
    }

    int mNumer, mDenom;
};


struct IsNonzero
{
public:
    template<class T>
    bool operator()(const RationalFunction<T> & r) const
    {
        if (r.numerator() == T(0.0))
            return false;
        return true;
    }
};

template<class T, class STREAM>
STREAM & operator<<(STREAM & str, const RationalFunction<T> & rhs);

#include "RationalFunction-inl.h"
#endif
