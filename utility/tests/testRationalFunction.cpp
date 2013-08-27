/*
 *  testRatFuncF.cpp
 *  Snapdragon
 *
 *  Created by Paul Hansen on 4/7/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test3D

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include <iostream>
#include <cstdlib>
#include "RationalFunction.h"

typedef RationalFunction<float> RatFuncF;

using namespace std;

BOOST_AUTO_TEST_CASE( TestConstructors )
{
    RatFuncF zed;
    BOOST_CHECK_EQUAL(zed.value(2342.13), 0.0);
    BOOST_CHECK_EQUAL(zed.order(), 0);
    BOOST_CHECK_EQUAL(zed.numerator()[0], 0.0);
    BOOST_CHECK_EQUAL(zed.denominator()[0], 1.0);
    BOOST_CHECK(zed == 0.0f);
    
    RatFuncF uno(1.0);
    BOOST_CHECK_EQUAL(uno.value(5.2), 1.0);
    BOOST_CHECK_EQUAL(uno.order(), 0);
    BOOST_CHECK_EQUAL(uno.numerator()[0], 1.0);
    BOOST_CHECK_EQUAL(uno.denominator()[0], 1.0);
    BOOST_CHECK(uno == 1.0f);
    
    RatFuncF secondOne(uno);
    RatFuncF secondZero(zed);
    BOOST_CHECK(secondOne == uno);
    BOOST_CHECK(secondZero == zed);
    BOOST_CHECK_EQUAL(secondOne(3.1), 1.0);
    BOOST_CHECK_EQUAL(secondZero(-12), 0.0);
    BOOST_CHECK_EQUAL(secondOne.order(), uno.order());
    BOOST_CHECK_EQUAL(secondZero.order(), zed.order());
    
    RatFuncF thirdOne, thirdZero;
    thirdOne = uno;
    thirdZero = zed;
    BOOST_CHECK(thirdOne == uno);
    BOOST_CHECK(thirdZero == zed);
    BOOST_CHECK(thirdOne != zed);
    BOOST_CHECK(thirdZero != uno);
    
    RatFuncF::Poly numer, denom;
    numer.coefficient(0, 1.0);
    numer.coefficient(1, 2.0);
    numer.coefficient(2, 3.0);
    denom.coefficient(0, 1.0);
    denom.coefficient(1, 2.0);
    
    RatFuncF r(numer, denom);
    BOOST_CHECK_EQUAL(r.order(), r.numerator().order());
    RatFuncF rInverse(denom, numer);
    BOOST_CHECK_EQUAL(rInverse.order(), rInverse.denominator().order());
    BOOST_CHECK(r.numerator() == numer);
    BOOST_CHECK(r.denominator() == denom);
    BOOST_CHECK(r.numerator() == rInverse.denominator());
    BOOST_CHECK(r.denominator() == rInverse.numerator());
}

BOOST_AUTO_TEST_CASE( Value )
{
    RatFuncF three(3.0, 1.0);
    RatFuncF three2(6.0, 2.0);
    BOOST_CHECK_EQUAL(three(0), 3.0);
    BOOST_CHECK_EQUAL(three2(0), 3.0);
    
    RatFuncF::Poly numer, denom;
    const int ORDER = 10;
    for (int nn = 0; nn <= ORDER; nn++)
    {
        numer.coefficient(nn, drand48());
        denom.coefficient(nn, drand48());
    }
    RatFuncF r(numer, denom);
    RatFuncF s(denom, numer);
    
    BOOST_CHECK_EQUAL(r(5.0), numer(5.0)/denom(5.0));
    BOOST_CHECK_CLOSE(s(5.0), 1.0/r(5.0), 1e-5);
}

BOOST_AUTO_TEST_CASE( TestOrder )
{
    RatFuncF::Poly o0(0), o1(0,1), o2(0,1,2), o3(0,1,2,3);
    
    RatFuncF z(o0, o1); // this is 0/x, which will be reset to 0/1.
    BOOST_CHECK_EQUAL(z.order(), 0);
    RatFuncF r(o1, o1);
    BOOST_CHECK_EQUAL(r.order(), 1);
    RatFuncF s(o3, o1);
    BOOST_CHECK_EQUAL(s.order(), 3);
    RatFuncF t(o2, o1);
    BOOST_CHECK_EQUAL(t.order(), 2);
}

BOOST_AUTO_TEST_CASE( AccessCoefficients )
{
    RatFuncF::Poly f;
    for (int nn = 0; nn <= 5; nn++)
        f.coefficient(nn, -1.0 * nn);
    
    BOOST_CHECK_EQUAL(f[0], 0.0);
    BOOST_CHECK_EQUAL(f[1], -1.0);
    BOOST_CHECK_EQUAL(f[2], -2.0);
    BOOST_CHECK_EQUAL(f[3], f.coefficient(3));
    BOOST_CHECK_EQUAL(f[4], f.coefficient(4));
    BOOST_CHECK_EQUAL(f[5], f.coefficient(5));
}

static RatFuncF::Poly polyOfOrderUpTo(unsigned int n)
{
    RatFuncF::Poly p;
    int order = lrand48()%n;
    for (int mm = 0; mm <= order; mm++)
        p.coefficient(mm, drand48());
    
    return p;
}

static RatFuncF ratFuncOfOrderUpTo(unsigned int n)
{
    return RatFuncF(polyOfOrderUpTo(n), polyOfOrderUpTo(n));
}

BOOST_AUTO_TEST_CASE( Arithmetic )
{
    int numTests = 50;
    for (int mm = 0; mm < numTests; mm++)
    {
        RatFuncF r(ratFuncOfOrderUpTo(10));
        RatFuncF s(ratFuncOfOrderUpTo(10));
        
        double x = drand48();
        
        BOOST_CHECK_CLOSE( (r+s)(x), r(x)+s(x), 1.0 ); // test to 1% tolerance
        BOOST_CHECK_CLOSE( (r-s)(x), r(x)-s(x), 1.0 );
        BOOST_CHECK_CLOSE( (r*s)(x), r(x)*s(x), 1.0 );
        BOOST_CHECK_CLOSE( (r/s)(x), r(x)/s(x), 1.0 );
        
        BOOST_CHECK_CLOSE( (r*s)(x), (s*r)(x), 1.0 );
        BOOST_CHECK_CLOSE( (r+s)(x), (s+r)(x), 1.0 );
        BOOST_CHECK_CLOSE( (r/s)(x), 1.0/(s/r)(x), 1.0 );
        BOOST_CHECK_CLOSE( (r-s)(x), -(s-r)(x), 1.0 );
    }
}

BOOST_AUTO_TEST_CASE( InPlaceAdd )
{
    RatFuncF::Poly numer;
    numer.coefficient(0, 1.0);
    numer.coefficient(1, 2.0);
    numer.coefficient(2, 3.0);
    
    RatFuncF::Poly denom;
    denom.coefficient(0, 2.0);
    denom.coefficient(1, -1.0);
    
    RatFuncF r(numer, denom);
    RatFuncF s(denom, numer);
    
    RatFuncF r2(r), r3(r), r4(r), r5(r);
    BOOST_CHECK(r2 == r);
    r2 += s;
    r3 -= s;
    r4 *= s;
    r5 /= s;
    
    BOOST_CHECK(r2 == r + s);
    BOOST_CHECK(r3 == r - s);
    BOOST_CHECK(r4 == r * s);
    BOOST_CHECK(r5 == r / s);
}


BOOST_AUTO_TEST_CASE( InPlaceArithmetic )
{
    RatFuncF r(ratFuncOfOrderUpTo(10)),
        s(ratFuncOfOrderUpTo(10));
    
    RatFuncF r2(r), r3(r), r4(r), r5(r);
    
    r2 += s;
    r3 -= s;
    r4 *= s;
    r5 /= s;
    
    int numTests = 50;
    for (int mm = 0; mm < numTests; mm++)
    {
        double x = 50 + 100*drand48();
        
        BOOST_CHECK_CLOSE(r2(x), (r+s)(x), 1.0);
        BOOST_CHECK_CLOSE(r3(x), (r-s)(x), 1.0);
        BOOST_CHECK_CLOSE(r4(x), (r*s)(x), 1.0);
        BOOST_CHECK_CLOSE(r5(x), (r/s)(x), 1.0);
    }
}

BOOST_AUTO_TEST_CASE( Print )
{
    BOOST_CHECK(true);
    
    RatFuncF defaultRF;
    BOOST_TEST_MESSAGE("Default RatFuncF " << defaultRF);
}
