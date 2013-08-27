/*
 *  testGeometry.cpp
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
#include "geometry.h"

using namespace std;

//static Vector3d randomVector()
//{
//    return Vector3d(drand48(), drand48(), drand48());
//}
//
//static Triangle3d randomTriangle()
//{
//    return Triangle3d(randomVector(), randomVector(), randomVector());
//}

BOOST_AUTO_TEST_CASE( TestPlane )
{
    Plane3d x_1(Vector3d(1,0,0), -1.0);
    Plane3d y_3(Vector3d(0,1,0), -3.0);
    Plane3d z_2(Vector3d(0,0,-1), 2.0);
    
    Vector3d rayX(1,0,0);
    Vector3d rayY(0,-13,0);
    Vector3d rayZ(0,0,200);
    
    BOOST_CHECK(intersection(x_1, Vector3d(0,0,0), rayX) == Vector3d(1,0,0));
    BOOST_CHECK(intersection(x_1, Vector3d(0,1,5), rayX) == Vector3d(1,1,5));
    
    BOOST_CHECK(intersection(y_3, Vector3d(0,0,0), rayY) == Vector3d(0,3,0));
    BOOST_CHECK(intersection(y_3, Vector3d(1,1,1), rayY) == Vector3d(1,3,1));
    
    BOOST_CHECK(intersection(z_2, Vector3d(0,0,0), rayZ) == Vector3d(0,0,2));
    BOOST_CHECK(intersection(z_2, Vector3d(-1,1,1), rayZ) == Vector3d(-1,1,2));
    
    // plane parallel to x=y.
//    Plane3d xy45(Vector3d(-1,1,0), -1.0); // should run through (-.707, .707)
//    Vector3d rayXY(unitVector(Vector3d(-1,1,0)));
}

BOOST_AUTO_TEST_CASE( TestTriPlane )
{
    Triangle3d xy(Vector3d(1,0,0), Vector3d(0,1,0), Vector3d(0,0,0));
    
    Vector3d normalVector(xy.normal());
    Plane3d plane(xy.plane());
    
    BOOST_CHECK(normalVector == plane.normal());
    BOOST_CHECK_EQUAL(normalVector[0], 0.0);
    BOOST_CHECK_EQUAL(normalVector[1], 0.0);
    BOOST_CHECK_EQUAL(normalVector[2], 1.0); // only if two sides are length 1
}

BOOST_AUTO_TEST_CASE( TestTriPlanePerturbation_Z )
{
    Triangle3d tri(Vector3d(0,0,1),
        Vector3d(1,0,1),
        Vector3d(0,1,1));
    
//    BOOST_TEST_MESSAGE("Plane = " << tri.plane());
    double epsilon = 1e-1;
    for (int nn = 0; nn < 3; nn++)
    {
//        BOOST_TEST_MESSAGE("dcdv" << nn << " = " << tri.dcdv(nn) );
//        BOOST_TEST_MESSAGE("dndv" << nn << " =\n" << tri.dndv(nn) );
//        BOOST_TEST_MESSAGE("dndv" << nn << " * v0 = " << tri.dndv(nn)*tri[0]);
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Vector3d perturbation = epsilon*Vector3d::unit(xyz);
            Triangle3d tri2(tri);
            tri2[nn] += perturbation;
            
            double diffc = tri2.plane().constant() - tri.plane().constant();
            double dc = tri.dcdv(nn)[xyz];
            
//            BOOST_TEST_MESSAGE("dcdv" << nn << "" << char('x'+xyz)
//                << " = " << diffc/epsilon
//                << " expected " << dc);
        }
    }
}

BOOST_AUTO_TEST_CASE( TestTriPlanePerturbation_Angled )
{
    Triangle3d tri(Vector3d(0,0,0), Vector3d(1,1,1), Vector3d(0,1,2));
    double epsilon = 1e-5;
    double tolerance = 1e-3;
    
    for (int vertex = 0; vertex < 3; vertex++)
    for (int direction = 0; direction < 3; direction++)
    {
        Vector3d dv[] = { Vector3d(0,0,0), Vector3d(0,0,0), Vector3d(0,0,0) };
        dv[vertex][direction] = epsilon;
        
        Triangle3d tri2(tri[0] + dv[0], tri[1]+dv[1], tri[2]+dv[2]);
        
        Vector3d dNormal = tri.dndv(0)*dv[0] +
            tri.dndv(1)*dv[1] +
            tri.dndv(2)*dv[2];
        Vector3d dUnitNormal = tri.dudv(0)*dv[0] +
            tri.dudv(1)*dv[1] +
            tri.dudv(2)*dv[2];
        
        // A plane is defined by dot(a,x) + b = 0.  b is the constant.  Here is
        // its sensitivity, dConstant.
        double dConstant = dot(tri.dcdv(0), dv[0]) +
            dot(tri.dcdv(1), dv[1]) +
            dot(tri.dcdv(2), dv[2]);
        
        // Relative error < tolerance:
        // |delta normal - dNormal| < tolerance * |dNormal|
        BOOST_CHECK_SMALL(
            norm(tri2.plane().normal() - tri.plane().normal() - dNormal),
            tolerance*norm(dNormal));
        
        // Error < tolerance*epsilon.
        // Loosely speaking, the constant is expected to be as big as the
        // perturbation for a triangle with O(1) dimensions.  If I change
        // the triangle dimensions I'll need to rethink the test tolerance.
        BOOST_CHECK_SMALL(
            tri2.plane().constant() - tri.plane().constant() - dConstant,
            tolerance*epsilon);
        
        // The (separate) calculation of dUnitNormal needs to be tested too.
        // The finite difference is O(epsilon), so try to keep the error down
        // to O(tolerance*epsilon).  In this case, using a larger triangle
        // we'd be able to push the tolerance lower.
        BOOST_CHECK_SMALL(
            norm(tri2.normal()/norm(tri2.normal()) -
                tri.normal()/norm(tri.normal()) -
                dUnitNormal),
            tolerance*epsilon);
    }
}


BOOST_AUTO_TEST_CASE( TestTriJacobianSimple )
{
    Triangle3d t(Vector3d(0,0,0), Vector3d::unit(0), Vector3d::unit(1));
    
    for (int vnum = 0; vnum < 3; vnum++)
    {
//        BOOST_TEST_MESSAGE("Sensitivity to " << vnum << ":\n"
//            << t.dndv(vnum));
        
        for (int xyz = 0; xyz < 3; xyz++)
        {
            Vector3d perturbation(Vector3d::unit(xyz));
            double epsilon = 1e-5;
            Triangle3d t2(t);
            t2[vnum] += epsilon*perturbation;
            
            Vector3d predictedNormal = t.normal() +
                t.dndv(vnum)*epsilon*perturbation;
            
            BOOST_CHECK_SMALL(norm(predictedNormal - t2.normal()), 1e-8);
//            BOOST_TEST_MESSAGE("Perturb v[" << vnum << "] along "
//                << perturbation << ": expected " << predictedNormal
//                << " got " << t2.normal() );
        }
    }
}

BOOST_AUTO_TEST_CASE( TestTriJacobian )
{
    Triangle3d tri(Vector3d(0,0,0), Vector3d(1,1,1), Vector3d(0,1,2));
    double epsilon = 1e-5;
    double tolerance = 1e-3;
        
    for (int vertex = 0; vertex < 3; vertex++)
    for (int direction = 0; direction < 3; direction++)
    {
        Vector3d perturbation(epsilon*Vector3d::unit(direction));
        Triangle3d tri2(tri);
        tri2[vertex] += perturbation;
        
        Vector3d predictedNormal = tri.normal() +
            tri.dndv(vertex)*perturbation;
        
        Vector3d diffN = tri2.normal() - tri.normal();
        Vector3d dn = tri.dndv(vertex)*perturbation;
//            BOOST_TEST_MESSAGE("Perturb v[" << vnum << "] along "
//                << perturbation << ": dn " << dn
//                << " got " << diffN );
        
        // For a triangle of roughly unit dimensions, the unit normal changes
        // by roughly epsilon.  Multiply by tolerance.
        BOOST_CHECK_SMALL(norm(predictedNormal - tri2.normal()),
            epsilon*tolerance*norm(tri.normal()));
//            BOOST_TEST_MESSAGE("Perturb v[" << vnum << "] along "
//                << perturbation << ": expected " << predictedNormal
//                << " got " << t2.normal() );
    }
}

BOOST_AUTO_TEST_CASE( TestRayTriangleDecomposition )
{
    Triangle3d t(Vector3d(0,0,0), Vector3d(1,1,1), Vector3d(0,1,2));
    double tolerance = 1e-5;
    
    Vector3d rayOrigin(-1,-2,-3);
    Vector3d rayDir(3, 2, -1);
    
    // Coefficient vector (a, b, t) for ray-triangle intersection
    Vector3d abt = t.intersectionCoefficients(rayOrigin, rayDir);
    Vector3d intersectPt = (t[1]-t[0])*abt[0] + (t[2]-t[0])*abt[1] + t[0];
    Vector3d intersectionCheck = intersection(t.plane(), rayOrigin, rayDir);
    
    BOOST_CHECK_SMALL(norm(intersectPt - intersectionCheck), tolerance);
}

// This test was written while debugging the test below.  It is not terribly
// enlightening.  (-:
BOOST_AUTO_TEST_CASE ( RayPlaneIntersectionPerturbationZ )
{
    Plane3d p(Vector3d(0,0,1), -1.0); // z = 1
    
    Vector3d rayOrigin(0,0,0);
    Vector3d rayDir(0,0,1);
    
    double t = intersectionTime(p, rayOrigin, rayDir);
    BOOST_CHECK_EQUAL(t, 1.0);
    
    double epsilon = 1e-1;
    Plane3d pPerturbedNormal(Vector3d(0,0,1+epsilon), -1.0);
    Plane3d pPerturbedConstant(Vector3d(0,0,1), -1.0+epsilon);
    
    Plane3d dtdp = intersectionTimeJacobian(p, rayOrigin, rayDir);
//    BOOST_TEST_MESSAGE(dtdp);
    
    double tN = intersectionTime(pPerturbedNormal, rayOrigin, rayDir);
    double tC = intersectionTime(pPerturbedConstant, rayOrigin, rayDir);
    
    double dtN = dtdp.normal()[2]*epsilon;
    double dtC = dtdp.constant()*epsilon;
    
//    BOOST_TEST_MESSAGE("t = " << t);
//    BOOST_TEST_MESSAGE("tN = " << tN << ", dtN = " << dtN);
//    BOOST_TEST_MESSAGE("tC = " << tC << ", dtC = " << dtC);
    
    double dtN2 = dot(Vector3d(0,0,epsilon), dtdp.normal());
    double dtC2 = epsilon*dtdp.constant();
    
    BOOST_CHECK_EQUAL(dtN, dtN2);
    BOOST_CHECK_EQUAL(dtC, dtC2);
}

BOOST_AUTO_TEST_CASE( RayPlaneIntersectionPerturbation )
{
    Plane3d p(Vector3d(1,2,3), 4);
    Vector3d rayOrigin(-1,-2,-3);
    Vector3d rayDir(3, 2, -1);
    
    double epsilon = 1e-5;
    double tolerance = 1e-5;

    double t = intersectionTime(p, rayOrigin, rayDir);
    
    // These should be identical values.
    BOOST_CHECK_SMALL(norm(intersection(p, rayOrigin, rayDir) -
        rayOrigin - t*rayDir), 1e-10);
    
    Plane3d dtdp = intersectionTimeJacobian(p, rayOrigin, rayDir);
    
    for (int xyz = 0; xyz < 4; xyz++) // change normal AND constant
    {
        Vector3d perturbNormal(0,0,0);
        double perturbConstant(0);
        
        if (xyz < 3)
            perturbNormal[xyz] = epsilon;
        else
            perturbConstant = epsilon;
        
        Plane3d perturbation(perturbNormal, perturbConstant);
        Plane3d p2(p.normal() + perturbation.normal(),
            p.constant() + perturbation.constant());
        double t2 = intersectionTime(p2, rayOrigin, rayDir);
        
        double dt = dot(perturbation.normal(), dtdp.normal())
            + perturbation.constant()*dtdp.constant();
        double tExpected = t + dt;
        
    //        BOOST_TEST_MESSAGE("delta t " << t2-t << " dt " << dt );
        // Calculate the relative error in dt.
        BOOST_CHECK_SMALL(t2 - tExpected, tolerance*fabs(dt));
    }
}

BOOST_AUTO_TEST_CASE( TriangleIntersectionPerturbation )
{
    vector<Triangle3d> testTriangles;
    vector<Vector3d> testRayOrigins;
    vector<Vector3d> testRayDirections;
    
    double epsilon = 1e-5;
    double tolerance = 1e-5;
    
    testTriangles.push_back(Triangle3d(
        Vector3d(0,0,0),
        Vector3d(1,0,0),
        Vector3d(0,1,0)));
    testTriangles.push_back(Triangle3d(
        Vector3d(0,0,0),
        Vector3d(0,1,0),
        Vector3d(0,0,1)));
    testTriangles.push_back(Triangle3d(
        Vector3d(10,0,0),
        Vector3d(20,0,40),
        Vector3d(0,0,40)));
    testRayOrigins.push_back(Vector3d(0.25,0.25,10));
    testRayOrigins.push_back(Vector3d(-20,0.3, 0.3));
    testRayOrigins.push_back(Vector3d(25,100,120));
    testRayDirections.push_back(Vector3d(0,0,-1));
    testRayDirections.push_back(Vector3d(100,0,0));
    testRayDirections.push_back(Vector3d(-1,-50,-10));
    
    for (unsigned int tt = 0; tt < testTriangles.size(); tt++)
    {
        Triangle3d & t = testTriangles.at(tt);
        Vector3d & rayOrigin = testRayOrigins.at(tt);
        Vector3d & rayDir = testRayDirections.at(tt);
        
//        Vector3d abt = t.intersectionCoefficients(rayOrigin, rayDir);
//        Vector3d intersectPt = (t[1]-t[0])*abt[0] + (t[2]-t[0])*abt[1] + t[0];
        Vector3d intersectPt = intersection(t.plane(), rayOrigin, rayDir);
        for (int vertex = 0; vertex < 3; vertex++)
        for (int direction = 0; direction < 3; direction++)
        {
            Vector3d perturbation(0,0,0);
            perturbation[direction] = epsilon;
            Triangle3d perturbed(t);
            perturbed[vertex] += perturbation;
            
            Vector3d perturbedIntersect = intersection(perturbed.plane(),
                rayOrigin, rayDir);
            Vector3d perturbedExpect = intersectPt +
                t.intersectionJacobian(rayOrigin, rayDir, vertex)*perturbation;
            
            // The intersection point should move by O(epsilon).
            // Check that it's accurate to tolerance. 
            BOOST_CHECK_SMALL(norm(perturbedIntersect - perturbedExpect),
                tolerance*epsilon*norm(intersectPt));
        }
    }
}

BOOST_AUTO_TEST_CASE( TestClosedInterval )
{
    ClosedInterval<double> i1(3, 4), iLeft, iLeftOverlap, iRight, iRightOverlap,
        iBigger;
    
    iLeft[0] = 1;
    iLeft[1] = 2;
    BOOST_CHECK(!i1.intersects(iLeft));
    BOOST_CHECK(!iLeft.intersects(i1));
    
    iRight[0] = 5;
    iRight[1] = 6;
    BOOST_CHECK(!i1.intersects(iRight));
    BOOST_CHECK(!iRight.intersects(i1));
    
    iLeftOverlap[0] = -1;
    iLeftOverlap[1] = 3.5;
    BOOST_CHECK(i1.intersects(iLeftOverlap));
    BOOST_CHECK(iLeftOverlap.intersects(i1));
    
    iRightOverlap[0] = 3.5;
    iRightOverlap[1] = 7;
    BOOST_CHECK(i1.intersects(iRightOverlap));
    BOOST_CHECK(iRightOverlap.intersects(i1));
    
    iBigger[0] = 0;
    iBigger[1] = 10;
    BOOST_CHECK(i1.intersects(iBigger));
    BOOST_CHECK(iBigger.intersects(i1));
    
    ClosedInterval<double> backwardsInterval(1, -1);
    ClosedInterval<double> onePoint(0,0);
    
    BOOST_CHECK(!backwardsInterval.intersects(onePoint));
    BOOST_CHECK(!onePoint.intersects(backwardsInterval));
    
    iBigger[0] = -2;
    iBigger[1] = 2;
    BOOST_CHECK(iBigger.intersects(backwardsInterval));
    BOOST_CHECK(backwardsInterval.intersects(iBigger));
}

BOOST_AUTO_TEST_CASE( TestRectProjection )
{
    Rect3d r(0,0,0,1,1,1);
    
    ClosedInterval<double> interval;
    
    interval = r.projectionAlong(Vector3d::unit(0));
    BOOST_CHECK_CLOSE(interval[0], 0, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 1, 1e-6);
    
    interval = r.projectionAlong(-Vector3d::unit(0));
    BOOST_CHECK_CLOSE(interval[0], -1, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 0, 1e-6);
    
    interval = r.projectionAlong(Vector3d(1,1,1));
    BOOST_CHECK_CLOSE(interval[0], 0, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 3.0, 1e-6);
    
    interval = r.projectionAlong(Vector3d(1,-1,0));
    BOOST_CHECK_CLOSE(interval[0], -1.0, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 1.0, 1e-6);
}

BOOST_AUTO_TEST_CASE( TestTriangleProjection )
{
    Triangle3d tri(Vector3d(0,0,0), Vector3d(1,0,0), Vector3d(0,1,1));
    
    ClosedInterval<double> interval;
    
    interval = tri.projectionAlong(Vector3d::unit(1));
    BOOST_CHECK_CLOSE(interval[0], 0, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 1, 1e-6);
    
    interval = tri.projectionAlong(-Vector3d::unit(1));
    BOOST_CHECK_CLOSE(interval[0], -1, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 0, 1e-6);
    
    interval = tri.projectionAlong(Vector3d(1,1,1));
    BOOST_CHECK_CLOSE(interval[0], 0, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 2, 1e-6);
    
    interval = tri.projectionAlong(Vector3d(-1,1,0));
    BOOST_CHECK_CLOSE(interval[0], -1, 1e-6);
    BOOST_CHECK_CLOSE(interval[1], 1, 1e-6);
}

BOOST_AUTO_TEST_CASE( TestRectTriangleIntersection )
{
    Triangle3d tri(Vector3d(0,0,0), Vector3d(1,0,0), Vector3d(0,1,1));
    
    Rect3d r(0,0,0,1,1,1);
    BOOST_CHECK(r.intersects(tri));
    
    r = Rect3d(0,0,1.001,1,1,2); // too high
    BOOST_CHECK(!r.intersects(tri));
    
    r = Rect3d(1.001, 0, 0, 2, 1, 1); // too far right
    BOOST_CHECK(!r.intersects(tri));
    
    r = Rect3d(0.8, 0.8, -1, 1.5, 1.5, 2);
    BOOST_CHECK(!r.intersects(tri));
    
    r = Rect3d(-10, -10, -10, 10, 10, 10);
    BOOST_CHECK(r.intersects(tri));
    
    r = Rect3d(-0.5, -0.5, 0.5, 0.5, 0.5, 1);
    BOOST_CHECK(r.intersects(tri));
    
    // this is the only case, I think, which will depend on the last nine axis
    // checks.
    tri = Triangle3d(Vector3d(0,0,0), Vector3d(1,1,0.5), Vector3d(0,1,0.5));
    r = Rect3d(0.5, -0.5, -1, 1, 0.49, 1);
    BOOST_CHECK(!r.intersects(tri));
    
    r = Rect3d(0.5, -0.5, -1, 1, 0.51, 1);
    BOOST_CHECK(r.intersects(tri));
}



