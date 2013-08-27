/*
 *  testPolyF.cpp
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

#include <iostream>
#include <cstdlib>
#include "Polynomial.h"

using namespace std;

typedef Polynomial<float> PolyF;

BOOST_AUTO_TEST_CASE( TestConstructors )
{
    PolyF zed;
    BOOST_CHECK_EQUAL(zed.value(2342.13), 0.0);
    BOOST_CHECK_EQUAL(zed.order(), 0);
    BOOST_CHECK_EQUAL(zed[0], 0.0);
    BOOST_CHECK_EQUAL(zed[1], 0.0);
    BOOST_CHECK_EQUAL(zed[2], 0.0);
    
    PolyF uno(1.0);
    BOOST_CHECK_EQUAL(uno.value(5.2), 1.0);
    BOOST_CHECK_EQUAL(uno.order(), 0);
    BOOST_CHECK_EQUAL(uno[0], 1.0);
    BOOST_CHECK_EQUAL(uno[1], 0.0);
    BOOST_CHECK_EQUAL(uno[2], 0.0);
    
    PolyF secondOne(uno);
    PolyF secondZero(zed);
    BOOST_CHECK(secondOne == uno);
    BOOST_CHECK(secondZero == zed);
    BOOST_CHECK_EQUAL(secondOne(3.1), 1.0);
    BOOST_CHECK_EQUAL(secondZero(-12), 0.0);
    BOOST_CHECK_EQUAL(secondOne.order(), uno.order());
    BOOST_CHECK_EQUAL(secondZero.order(), zed.order());
    
    PolyF linear(1.0, 2.0);
    BOOST_CHECK_EQUAL(linear.order(), 1);
    BOOST_CHECK_EQUAL(linear[0], 1.0);
    BOOST_CHECK_EQUAL(linear[1], 2.0);
    BOOST_CHECK_EQUAL(linear[2], 0.0);
    BOOST_CHECK_EQUAL(linear[3], 0.0);
    
    PolyF quadratic(2.0, 3.0, -4.0);
    BOOST_CHECK_EQUAL(quadratic.order(), 2);
    BOOST_CHECK_EQUAL(quadratic[0], 2.0);
    BOOST_CHECK_EQUAL(quadratic[1], 3.0);
    BOOST_CHECK_EQUAL(quadratic[2], -4.0);
    BOOST_CHECK_EQUAL(quadratic[3], 0.0);
    BOOST_CHECK_EQUAL(quadratic[4], 0.0);
    
    PolyF cubic(-0.5, 0.5, 1.5, 2.5);
    BOOST_CHECK_EQUAL(cubic.order(), 3);
    BOOST_CHECK_EQUAL(cubic[0], -0.5);
    BOOST_CHECK_EQUAL(cubic[1], 0.5);
    BOOST_CHECK_EQUAL(cubic[2], 1.5);
    BOOST_CHECK_EQUAL(cubic[3], 2.5);
    BOOST_CHECK_EQUAL(cubic[4], 0.0);
    BOOST_CHECK_EQUAL(cubic[cubic.maximumOrder()], 0.0);
}

BOOST_AUTO_TEST_CASE( Value )
{
    PolyF p;
    p.coefficient(0, 1.0);
    p.coefficient(1, 2.0);
    p.coefficient(2, -2.0);
    p.coefficient(3, -1.0);
    p.coefficient(4, 1.0);
    
    BOOST_CHECK_EQUAL(p(0.0), p[0]);
    
    for (double xx = -4.0; xx <= 4.0; xx += 1.0)
    {
        double val = 1.0 + 2.0*xx - 2.0*xx*xx - 1.0*xx*xx*xx
            + 1.0*xx*xx*xx*xx;
        BOOST_CHECK_EQUAL(p(xx), val);
    }
}


BOOST_AUTO_TEST_CASE( ChangeOrder )
{
    PolyF p(1, 2, 3, 4);
    BOOST_CHECK_EQUAL(p.order(), 3);
    BOOST_CHECK_EQUAL(p[4], 0.0);
    p.coefficient(1,0);
    BOOST_CHECK_EQUAL(p.order(), 3);
    p.coefficient(2,0);
    BOOST_CHECK_EQUAL(p.order(), 3);
    p.coefficient(3,0);
    BOOST_CHECK_EQUAL(p.order(), 0);
}

BOOST_AUTO_TEST_CASE( AccessCoefficients )
{
    PolyF f;
    for (int nn = 0; nn <= 5; nn++)
        f.coefficient(nn, -1.0 * nn);
    
    BOOST_CHECK_EQUAL(f[0], 0.0);
    BOOST_CHECK_EQUAL(f[1], -1.0);
    BOOST_CHECK_EQUAL(f[2], -2.0);
    BOOST_CHECK_EQUAL(f[3], f.coefficient(3));
    BOOST_CHECK_EQUAL(f[4], f.coefficient(4));
    BOOST_CHECK_EQUAL(f[5], f.coefficient(5));
}

BOOST_AUTO_TEST_CASE( SumDifference )
{
    PolyF p, q, r;
    
    p.coefficient(0, 1.0);
    p.coefficient(1, -1.0);
    p.coefficient(2, 2.0);
    BOOST_CHECK_EQUAL(p.order(), 2);
    
    q.coefficient(5.0);
    BOOST_CHECK_EQUAL(q.order(), 0);
    
    r.coefficient(0, -1.0);
    r.coefficient(1, 2.0);
    BOOST_CHECK_EQUAL(r.order(), 1);
    
    PolyF sumPQ = p + q;
    BOOST_CHECK_EQUAL(sumPQ.order(), 2);
    BOOST_CHECK_EQUAL(sumPQ[0], p[0] + q[0]);
    BOOST_CHECK_EQUAL(sumPQ[1], p[1]);
    BOOST_CHECK_EQUAL(sumPQ[2], p[2]);
    BOOST_CHECK_EQUAL(sumPQ(3.0), p(3.0) + q(3.0));
    PolyF sumQP = q + p;
    BOOST_CHECK(sumQP == sumPQ);
    
    PolyF diffQR = q - r;
    BOOST_CHECK_EQUAL(diffQR.order(), 1);
    BOOST_CHECK_EQUAL(diffQR[0], q[0] - r[0]);
    BOOST_CHECK_EQUAL(diffQR[1], -r[1]);
    BOOST_CHECK_EQUAL(diffQR(-4.0), q(-4.0) - r(-4.0));
    PolyF diffRQ = r - q;
    BOOST_CHECK(diffRQ == -1.0f*diffQR);
    BOOST_CHECK(diffRQ == diffQR*(-1.0f));
    
    PolyF sumPQR = p + q + r;
    BOOST_CHECK_EQUAL(sumPQR(-10), p(-10) + q(-10) + r(-10));
    
    PolyF diffPQR = p - q - r;
    BOOST_CHECK_EQUAL(diffPQR(-5), p(-5) - q(-5) - r(-5));
    
    PolyF p2(p);
    p2 += q;
    BOOST_CHECK(p2 == sumPQ);
    
    PolyF q2(q);
    q2 -= r;
    BOOST_CHECK(q2 == diffQR);
}

BOOST_AUTO_TEST_CASE( Product )
{
    PolyF p, q, r;
    
    p.coefficient(0, 1.0);
    p.coefficient(1, -1.0);
    p.coefficient(2, 2.0);
    
    q.coefficient(0, 5.0);
    
    r.coefficient(0, -1.0);
    r.coefficient(1, 2.0);
    
    PolyF prodPQ = p*q;
    BOOST_CHECK_EQUAL(prodPQ(0.0), prodPQ[0]);
    BOOST_CHECK_EQUAL(prodPQ(0.0), p(0.0)*q(0.0));
    BOOST_CHECK_EQUAL(prodPQ.order(), 2);
    BOOST_CHECK_EQUAL(prodPQ[0], p[0]*q[0]);
    BOOST_CHECK_EQUAL(prodPQ[1], p[1]*q[0]);
    BOOST_CHECK_EQUAL(prodPQ[2], p[2]*q[0]);
    PolyF prodQP = q*p;
    BOOST_CHECK(prodPQ == prodQP);
    
    PolyF prodPR = p*r;
    BOOST_CHECK_EQUAL(prodPR.order(), p.order() + r.order());
    BOOST_CHECK_EQUAL(prodPR[0], p[0]*r[0]);
    BOOST_CHECK_EQUAL(prodPR[1], p[1]*r[0] + p[0]*r[1]);
    BOOST_CHECK_EQUAL(prodPR[2], p[2]*r[0] + p[1]*r[1]);
    
    BOOST_CHECK_EQUAL(prodPR(0.0), p(0.0)*r(0.0));
    BOOST_CHECK_EQUAL(prodPR(1.0), p(1.0)*r(1.0));
    BOOST_CHECK_EQUAL(prodPR(-2.0), p(-2.0)*r(-2.0));
    
    PolyF prodRP = r*p;
    BOOST_CHECK(prodPR == prodRP);
    
    PolyF prodPQR = p*q*r;
    BOOST_CHECK_EQUAL(prodPQR(22), p(22)*q(22)*r(22));
    
    // test in-place multiplication
    PolyF p2(p);
    p2 *= q;
    BOOST_CHECK(p2 == prodPQ);
    PolyF p3(p);
    p3 *= r;
    BOOST_CHECK(p3 == prodPR);
}

BOOST_AUTO_TEST_CASE( TestPrint )
{
    PolyF nil;
    PolyF zed(0.0);
    PolyF uno(1.0);
    
    BOOST_CHECK(true);
    BOOST_TEST_MESSAGE("Default zero polynomial: " << nil);
    BOOST_TEST_MESSAGE("Zero polynomial: " << zed);
    BOOST_TEST_MESSAGE("Unit polynomial: " << uno);
}

