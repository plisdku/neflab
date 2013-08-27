/*
 *  Polynomial-inl.h
 *  Snapdragon
 *
 *  Created by Paul Hansen on 4/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifdef _POLYNOMIAL_
#include <cassert>
#include <stdexcept>

template<class T>
Polynomial<T>::
Polynomial() :
    mOrder(0)
{
    mCoefficients[0] = T(0.0);
}

template<class T>
Polynomial<T>::
Polynomial(const Polynomial<T> & copyMe) :
    mOrder(copyMe.mOrder)
{
    assert(mOrder <= kMaximumPolynomialOrder);
    for (int mm = 0; mm <= mOrder; mm++)
        mCoefficients[mm] = copyMe.mCoefficients[mm];
}

template<class T>
Polynomial<T>::
Polynomial(Real constantValue) :
    mOrder(0)
{
    mCoefficients[0] = constantValue;
}

template<class T>
Polynomial<T>::
Polynomial(Real c0, Real c1) :
    mOrder(1)
{
    mCoefficients[0] = c0;
    mCoefficients[1] = c1;
}

template<class T>
Polynomial<T>::
Polynomial(Real c0, Real c1, Real c2) :
    mOrder(2)
{
    mCoefficients[0] = c0;
    mCoefficients[1] = c1;
    mCoefficients[2] = c2;
}

template<class T>
Polynomial<T>::
Polynomial(Real c0, Real c1, Real c2, Real c3) :
    mOrder(3)
{
    mCoefficients[0] = c0;
    mCoefficients[1] = c1;
    mCoefficients[2] = c2;
    mCoefficients[3] = c3;
}

template<class T>
typename Polynomial<T>::Real Polynomial<T>::
coefficient(unsigned int whichOrder) const
{
    if (whichOrder > kMaximumPolynomialOrder)
        throw(std::logic_error("Polynomial order out of bounds"));
    else if (whichOrder > order())
        return T(0.0);
    else
        return mCoefficients[whichOrder];
}

template<class T>
void Polynomial<T>::
coefficient(unsigned int whichOrder, Real newValue)
{
    if (whichOrder > kMaximumPolynomialOrder)
        throw(std::logic_error("Polynomial order out of bounds"));
    
    if (whichOrder > order())
    {
        for (int nn = order()+1; nn <= whichOrder; nn++)
            mCoefficients[nn] = T(0.0);
        mOrder = whichOrder;
    }
        
    mCoefficients[whichOrder] = newValue;
    
    if (whichOrder == order())
    {
        int ord = whichOrder;
        while (ord > 0 && mCoefficients[ord] == T(0.0))
            ord--;
        mOrder = ord;
    }
    
//    if (whichOrder > order())
//        throw(std::logic_error("Polynomial order out of bounds"));
//    mCoefficients[whichOrder] = newValue;
}

template<class T>
typename Polynomial<T>::Real Polynomial<T>::
value(Real x) const
{
    Real xToTheN = T(1.0); // begin with x^0
    Real sum = T(0.0);
    for (int mm = 0; mm <= mOrder; mm++)
    {
        sum += xToTheN*mCoefficients[mm];
        xToTheN *= x;
    }
    return sum;
}

template<class T>
Polynomial<T> & Polynomial<T>::
operator=(const Polynomial<T> & rhs)
{
    if (this == &rhs)
        return *this;
    
    mOrder = rhs.mOrder;
    assert(mOrder <= kMaximumPolynomialOrder);
    for (int mm = 0; mm <= mOrder; mm++)
        mCoefficients[mm] = rhs.mCoefficients[mm];
    
    return *this;
}

template<class T>
Polynomial<T> & Polynomial<T>::
operator+=(const Polynomial<T> & rhs)
{
//    if (order() < rhs.order())
//        order(rhs.order());
    
    for (int mm = 0; mm <= rhs.order(); mm++)
        coefficient(mm, coefficient(mm) + rhs.coefficient(mm));
    
    return *this;
}

template<class T>
Polynomial<T> & Polynomial<T>::
operator-=(const Polynomial<T> & rhs)
{
//    if (order() < rhs.order())
//        order(rhs.order());
    
    for (int mm = 0; mm <= rhs.order(); mm++)
        coefficient(mm, coefficient(mm) - rhs.coefficient(mm));
    
    return *this;
}

template<class T>
Polynomial<T> & Polynomial<T>::
operator*=(const Polynomial<T> & rhs)
{
    *this = (*this * rhs);
    return *this;
}

template<class T>
bool operator==(const Polynomial<T> & lhs, const Polynomial<T> & rhs)
{
    if (lhs.order() != rhs.order())
        return false;
    for (int mm = 0; mm <= lhs.order(); mm++)
    if (lhs[mm] != rhs[mm])
        return false;
    return true;
}

template<class T>
bool operator!=(const Polynomial<T> & lhs, const Polynomial<T> & rhs)
{
    return !(lhs == rhs);
}


template<class T>
bool operator==(const Polynomial<T> & lhs, const T & rhs)
{
    return (lhs.order() == 0 && lhs[0] == rhs);
}

template<class T>
bool operator!=(const Polynomial<T> & lhs, const T & rhs)
{
    return !(lhs == rhs);
}

template<class T>
bool operator==(const T & lhs, const Polynomial<T> & rhs)
{
    return (rhs.order() == 0 && rhs[0] == lhs);
}

template<class T>
bool operator!=(const T & lhs, const Polynomial<T> & rhs)
{
    return !(lhs == rhs);
}



template<class T>
Polynomial<T> operator+(const Polynomial<T> & lhs, const Polynomial<T> & rhs)
{
    if (lhs.order() >= rhs.order())
    {
        Polynomial<T> result;
//        result.order(lhs.order());
        for (int mm = 0; mm <= rhs.order(); mm++)
            result.coefficient(mm, lhs.coefficient(mm) + rhs.coefficient(mm));
        for (int mm = rhs.order()+1; mm <= lhs.order(); mm++)
            result.coefficient(mm, lhs.coefficient(mm));
        return result;
    }
    else
    {
        return (rhs + lhs);
    }
}

template<class T>
Polynomial<T> operator-(const Polynomial<T> & lhs, const Polynomial<T> & rhs)
{
    Polynomial<T> result;
    if (lhs.order() >= rhs.order())
    {
//        result.order(lhs.order());
        for (int mm = 0; mm <= rhs.order(); mm++)
            result.coefficient(mm, lhs.coefficient(mm) - rhs.coefficient(mm));
        for (int mm = rhs.order()+1; mm <= lhs.order(); mm++)
            result.coefficient(mm, lhs.coefficient(mm));
    }
    else
    {
//        result.order(rhs.order());
        for (int mm = 0; mm <= lhs.order(); mm++)
            result.coefficient(mm, lhs.coefficient(mm) - rhs.coefficient(mm));
        for (int mm = lhs.order()+1; mm <= rhs.order(); mm++)
            result.coefficient(mm, -1.0 * rhs.coefficient(mm));
    }
    return result;
}

template<class T>
Polynomial<T> operator*(const Polynomial<T> & lhs, const Polynomial<T> & rhs)
{
    Polynomial<T> result;
//    result.order(lhs.order() + rhs.order());
    
    for (int ll = 0; ll <= lhs.order(); ll++)
    for (int rr = 0; rr <= rhs.order(); rr++)
    {
        result.coefficient(ll+rr, result.coefficient(ll+rr) + lhs[ll]*rhs[rr]);
    }
    
    return result;
}

template<class T>
Polynomial<T> operator+(const Polynomial<T> & lhs, const T & rhs)
{
    Polynomial<T> result;
    
    for (int ll = 0; ll <= lhs.order(); ll++)
        result.coefficient(ll, lhs[ll]+rhs);
    
    return result;
}
template<class T>
Polynomial<T> operator-(const Polynomial<T> & lhs, const T & rhs)
{
    Polynomial<T> result;
    
    for (int ll = 0; ll <= lhs.order(); ll++)
        result.coefficient(ll, lhs[ll]-rhs);
    
    return result;
}
template<class T>
Polynomial<T> operator*(const Polynomial<T> & lhs, const T & rhs)
{
    Polynomial<T> result;
    
    for (int ll = 0; ll <= lhs.order(); ll++)
        result.coefficient(ll, lhs[ll]*rhs);
    
    return result;
}
template<class T>
Polynomial<T> operator/(const Polynomial<T> & lhs, const T & rhs)
{
    Polynomial<T> result;
    
    for (int ll = 0; ll <= lhs.order(); ll++)
        result.coefficient(ll, lhs[ll]/rhs);
    
    return result;
}

template<class T>
Polynomial<T> operator+(const T & lhs, const Polynomial<T> & rhs)
{
    Polynomial<T> result;
    
    for (int rr = 0; rr <= rhs.order(); rr++)
        result.coefficient(rr, lhs+rhs[rr]);
    
    return result;
}
template<class T>
Polynomial<T> operator-(const T & lhs, const Polynomial<T> & rhs)
{
    Polynomial<T> result;
    
    for (int rr = 0; rr <= rhs.order(); rr++)
        result.coefficient(rr, lhs-rhs[rr]);
    
    return result;
}
template<class T>
Polynomial<T> operator*(const T & lhs, const Polynomial<T> & rhs)
{
    Polynomial<T> result;
    
    for (int rr = 0; rr <= rhs.order(); rr++)
        result.coefficient(rr, lhs*rhs[rr]);
    
    return result;
}

template<class T>
Polynomial<T> operator-(const Polynomial<T> & p)
{
    return T(-1.0) * p;
}


template<class T, class STREAM>
STREAM & operator<<(STREAM & str, const Polynomial<T> & rhs)
{
    str << "(";
    for (int nn = 0; nn < rhs.order(); nn++)
    {
        str << rhs[nn] << "x^" << nn << " + ";
    }
    str << rhs[rhs.order()] << "x^" << rhs.order() << ")";
    return str;
}




#endif