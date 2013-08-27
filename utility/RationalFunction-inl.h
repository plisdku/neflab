/*
 *  RationalFunction-inl.h
 *  Snapdragon
 *
 *  Created by Paul Hansen on 4/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifdef _RATIONALFUNCTION_

template<class T>
RationalFunction<T>::
RationalFunction() :
    mNumerator(T(0.0)),
    mDenominator(T(1.0))
{
}

template<class T>
RationalFunction<T>::
RationalFunction(T constantValue) :
    mNumerator(constantValue),
    mDenominator(T(1.0))
{
}

template<class T>
RationalFunction<T>::
RationalFunction(const Polynomial<T> & poly) :
    mNumerator(poly),
    mDenominator(T(1.0))
{
}

template<class T>
RationalFunction<T>::
RationalFunction(const Polynomial<T> & numer, const Polynomial<T> & denom) :
    mNumerator(numer),
    mDenominator(denom)
{
    if (mNumerator == T(0.0))
        mDenominator = T(1.0);
}

template<class T>
RationalFunction<T>::
RationalFunction(const RationalFunction<T> & copyMe) :
    mNumerator(copyMe.numerator()),
    mDenominator(copyMe.denominator())
{
}

template<class T>
RationalFunction<T> RationalFunction<T>::
irregular(const Polynomial<T> & numer, const Polynomial<T> & denom)
{
    RationalFunction<T> r;
    r.numerator() = numer;
    r.denominator() = denom;
    return r;
}

template<class T>
T RationalFunction<T>::
value(T x) const
{
    return mNumerator(x) / mDenominator(x);
}

template<class T>
unsigned int RationalFunction<T>::
order() const
{
    return mNumerator.order() > mDenominator.order() ?
        mNumerator.order() : mDenominator.order();
}

template<class T>
RationalFunction<T> & RationalFunction<T>::
operator=(const RationalFunction<T> & rhs)
{
    if (this == &rhs)
        return *this;
    
    mNumerator = rhs.numerator();
    mDenominator = rhs.denominator();
    
    return *this;
}

template<class T>
RationalFunction<T> & RationalFunction<T>::
operator+=(const RationalFunction<T> & rhs)
{
    if (denominator() == rhs.denominator())
        numerator() += rhs.numerator();
    else
    {
        mNumerator = numerator()*rhs.denominator() +
            denominator()*rhs.numerator();
    }
    if (mNumerator == T(0.0))
        mDenominator = T(1.0);
    else
        mDenominator = denominator() * rhs.denominator();
    
    return *this;
}

template<class T>
RationalFunction<T> & RationalFunction<T>::
operator-=(const RationalFunction<T> & rhs)
{
    if (denominator() == rhs.denominator())
        numerator() -= rhs.numerator();
    else
    {
        mNumerator = numerator()*rhs.denominator() -
            denominator()*rhs.numerator();
    }
    if (mNumerator == T(0.0))
        mDenominator = T(1.0);
    else
        mDenominator = denominator() * rhs.denominator();
    
    return *this;
}

template<class T>
RationalFunction<T> & RationalFunction<T>::
operator*=(const RationalFunction<T> & rhs)
{
    mNumerator *= rhs.numerator();
    if (mNumerator == T(0.0))
        mDenominator = T(1.0);
    else
        mDenominator *= rhs.denominator();
    return *this;
}

template<class T>
RationalFunction<T> & RationalFunction<T>::
operator/=(const RationalFunction<T> & rhs)
{
    mNumerator *= rhs.denominator();
    mDenominator *= rhs.numerator();
    return *this;
}

template<class T>
bool operator==(const RationalFunction<T> & lhs, const RationalFunction<T> & rhs)
{
    return lhs.numerator() == rhs.numerator() &&
        lhs.denominator() == rhs.denominator();
}

template<class T>
bool operator!=(const RationalFunction<T> & lhs, const RationalFunction<T> & rhs)
{
    return !(lhs == rhs);
}

template<class T>
RationalFunction<T> operator+(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs)
{
    if (lhs.denominator() == rhs.denominator())
    {
        return RationalFunction<T>(lhs.numerator() + rhs.numerator(),
            lhs.denominator());
    }
    else
    {
        return RationalFunction<T>(lhs.numerator()*rhs.denominator() +
            lhs.denominator()*rhs.numerator(),
            lhs.denominator() * rhs.denominator());
    }
}

template<class T>
RationalFunction<T> operator-(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs)
{
    if (lhs.denominator() == rhs.denominator())
    {
        return RationalFunction<T>(lhs.numerator() - rhs.numerator(),
            lhs.denominator());
    }
    else
    {
        return RationalFunction<T>(lhs.numerator()*rhs.denominator() -
            lhs.denominator()*rhs.numerator(),
            lhs.denominator() * rhs.denominator());
    }
}

template<class T>
RationalFunction<T> operator*(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs)
{
    return RationalFunction<T>(lhs.numerator()*rhs.numerator(),
        lhs.denominator()*rhs.denominator());
}

template<class T>
RationalFunction<T> operator/(const RationalFunction<T> & lhs,
    const RationalFunction<T> & rhs)
{
    return RationalFunction<T>(lhs.numerator()*rhs.denominator(),
        lhs.denominator()*rhs.numerator());
}






template<class T>
bool operator==(const RationalFunction<T> & lhs, const T & rhs)
{
    return (lhs.numerator().order() == 0 && lhs.denominator().order() == 0 &&
        lhs.numerator()[0]/lhs.denominator()[0] == rhs);
}

template<class T>
bool operator!=(const RationalFunction<T> & lhs, const T & rhs)
{
    return !(lhs == rhs);
}

template<class T>
RationalFunction<T> operator+(const RationalFunction<T> & lhs,
    const T & rhs)
{
    return RationalFunction<T>(lhs.numerator() + rhs*lhs.denominator(),
        lhs.denominator());
}

template<class T>
RationalFunction<T> operator-(const RationalFunction<T> & lhs,
    const T & rhs)
{
    return RationalFunction<T>(lhs.numerator() - rhs*lhs.denominator(),
        lhs.denominator());
}

template<class T>
RationalFunction<T> operator*(const RationalFunction<T> & lhs,
    const T & rhs)
{
    return RationalFunction<T>(lhs.numerator()*rhs, lhs.denominator());
}

template<class T>
RationalFunction<T> operator/(const RationalFunction<T> & lhs,
    const T & rhs)
{
    return RationalFunction<T>(lhs.numerator() / rhs, lhs.denominator());
}




template<class T>
bool operator==(const T & lhs, const RationalFunction<T> & rhs)
{
    return (rhs == lhs);
}

template<class T>
bool operator!=(const T & lhs, const RationalFunction<T> & rhs)
{
    return !(lhs == rhs);
}

template<class T>
RationalFunction<T> operator+(const T & lhs, const RationalFunction<T> & rhs)
{
    return RationalFunction<T>(rhs.numerator() + lhs*rhs.denominator(),
        rhs.denominator());
}

template<class T>
RationalFunction<T> operator-(const T & lhs, const RationalFunction<T> & rhs)
{
    return RationalFunction<T>(lhs*rhs.denominator() - rhs.numerator(),
        rhs.denominator());
}

template<class T>
RationalFunction<T> operator*(const T & lhs, const RationalFunction<T> & rhs)
{
    return RationalFunction<T>(rhs.numerator()*lhs, rhs.denominator());
}

template<class T>
RationalFunction<T> operator/(const T & lhs, const RationalFunction<T> & rhs)
{
    return RationalFunction<T>(lhs*rhs.denominator(), rhs.numerator());
}

template<class T>
RationalFunction<T> operator-(const RationalFunction<T> & r)
{
    return RationalFunction<T>(-r.numerator(), r.denominator());
}







template<class T>
RationalFunction<T> reciprocal(const RationalFunction<T> & r)
{
    if (r.denominator() == T(0.0))
        throw(std::logic_error("Reciprocal of zero function is not defined"));
    return RationalFunction<T>(r.denominator(), r.numerator());
}


template<class T, class STREAM>
STREAM & operator<<(STREAM & str, const RationalFunction<T> & rhs)
{
    str << "[" << rhs.numerator() << " / " << rhs.denominator() << "]";
    return str;
}


#endif
