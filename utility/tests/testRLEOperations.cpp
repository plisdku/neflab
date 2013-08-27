/*
 *  testRLEOperations.cpp
 *  Snapdragon
 *
 *  Created by Paul Hansen on 5/17/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Test3D

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <vector>
#include "geometry.h"
#include "RLEOperations.h"

using namespace std;
using namespace RLE;


BOOST_AUTO_TEST_CASE( TestDownsample2x2x2 )
{
    Rect3i allCells(0,0,0,1,1,1);
    Rect3i outCells(0,0,0,0,0,0);
    
    DynamicRLE3<int> rle;
    rle.mark(0,0,0,0,0,0,1);    // rle(0,0,0) = 1
    rle.mark(1,1,0,1,1,0,2);    // rle(1,1,0) = 2
    rle.mark(0,1,1,0,1,1,3);    // rle(0,1,1) = 3
    
    DynamicRLE3<int> out;
    DynamicRLE3<int> outSlow;
    
    // Get (0,0,0)
    out = downsample2(rle, allCells, outCells);
    outSlow = downsample2_slow(rle, allCells, outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK_EQUAL(out(0,0,0), 1);
    BOOST_CHECK_EQUAL(outSlow(0,0,0), 1);
    
    // Get (1,0,0)
    out = downsample2(rle, allCells+Vector3i(1,0,0), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(1,0,0), outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK(!out.isMarked(0,0,0));
    BOOST_CHECK(!outSlow.isMarked(0,0,0));
    
    // Get (0,1,0)
    out = downsample2(rle, allCells+Vector3i(0,1,0), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(0,1,0), outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK(!out.isMarked(0,0,0));
    BOOST_CHECK(!outSlow.isMarked(0,0,0));
    
    // Get (1,1,0)
    out = downsample2(rle, allCells+Vector3i(1,1,0), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(1,1,0), outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK_EQUAL(out(0,0,0), 2);
    BOOST_CHECK_EQUAL(outSlow(0,0,0), 2);
    
    // Get (0,0,1)
    out = downsample2(rle, allCells+Vector3i(0,0,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(0,0,1), outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK(!out.isMarked(0,0,0));
    BOOST_CHECK(!outSlow.isMarked(0,0,0));
    
    // Get (1,0,1)
    out = downsample2(rle, allCells+Vector3i(1,0,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(1,0,1), outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK(!out.isMarked(0,0,0));
    BOOST_CHECK(!outSlow.isMarked(0,0,0));
    
    // Get (0,1,1)
    out = downsample2(rle, allCells+Vector3i(0,1,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(0,1,1), outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK_EQUAL(out(0,0,0), 3);
    BOOST_CHECK_EQUAL(outSlow(0,0,0), 3);
    
    // Get (1,1,1)
    out = downsample2(rle, allCells+Vector3i(1,1,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(1,1,1), outCells);
    BOOST_CHECK(out == outSlow);
    BOOST_CHECK(!out.isMarked(0,0,0));
    BOOST_CHECK(!outSlow.isMarked(0,0,0));
}

BOOST_AUTO_TEST_CASE( TestDownsampleX )
{
    int outY = -3, outZ = 10;
    Rect3i allCells(-10,0,0,9,1,1);
    Rect3i outCells(-5,outY,outZ,4,outY,outZ);
    
    DynamicRLE3<int> rle;
    
    rle.mark(-10,0,0,9,0,0,1);
    rle.mark(-10,1,0,9,1,0,2);
    rle.mark(-10,0,1,9,0,1,3);
    rle.mark(-10,1,1,9,1,1,4);
    
    DynamicRLE3<int> out, outSlow;
    
    // Get row (:,0,0)
    out = downsample2(rle, allCells, outCells);
    outSlow = downsample2_slow(rle, allCells, outCells);
    BOOST_CHECK(out == outSlow);
    for (int xx = -5; xx <= 4; xx++)
    {
        BOOST_CHECK_EQUAL(out(xx,outY,outZ), 1);
        BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 1);
    }
        
    // Get row (:,1,0)
    out = downsample2(rle, allCells+Vector3i(0,1,0), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(0,1,0), outCells);
    BOOST_CHECK(out == outSlow);
    for (int xx = -5; xx <= 4; xx++)
    {
        BOOST_CHECK_EQUAL(out(xx,outY,outZ), 2);
        BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 2);
    }
    
    // Get row (:,0,1)
    out = downsample2(rle, allCells+Vector3i(0,0,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(0,0,1), outCells);
    BOOST_CHECK(out == outSlow);
    for (int xx = -5; xx <= 4; xx++)
    {
        BOOST_CHECK_EQUAL(out(xx,outY,outZ), 3);
        BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 3);
    }
    
    // Get row (:,1,1)
    out = downsample2(rle, allCells+Vector3i(0,1,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(0,1,1), outCells);
    BOOST_CHECK(out == outSlow);
    for (int xx = -5; xx <= 4; xx++)
    {
        BOOST_CHECK_EQUAL(out(xx,outY,outZ), 4);
        BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 4);
    }
}

BOOST_AUTO_TEST_CASE( TestDownsampleY )
{
    int outX = -100, outZ = -9;
    Rect3i allCells(0,-10,0,1,9,1);
    Rect3i outCells(outX,-5,outZ,outX,4,outZ);
    
    DynamicRLE3<int> rle;
    
    rle.mark(0,-10,0,0,9,0,1);
    rle.mark(1,-10,0,1,9,0,2);
    rle.mark(0,-10,1,0,9,1,3);
    rle.mark(1,-10,1,1,9,1,4);
    
    DynamicRLE3<int> out, outSlow;
    
    // Get row (0,:,0)
    out = downsample2(rle, allCells, outCells);
    outSlow = downsample2_slow(rle, allCells, outCells);
    BOOST_CHECK(out == outSlow);
    for (int yy = -5; yy <= 4; yy++)
    {
        BOOST_CHECK_EQUAL(out(outX,yy,outZ), 1);
        BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 1);
    }
        
    // Get row (1,:,0)
    out = downsample2(rle, allCells+Vector3i(1,0,0), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(1,0,0), outCells);
    BOOST_CHECK(out == outSlow);
    for (int yy = -5; yy <= 4; yy++)
    {
        BOOST_CHECK_EQUAL(out(outX,yy,outZ), 2);
        BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 2);
    }
    
    // Get row (0,:,1)
    out = downsample2(rle, allCells+Vector3i(0,0,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(0,0,1), outCells);
    BOOST_CHECK(out == outSlow);
    for (int yy = -5; yy <= 4; yy++)
    {
        BOOST_CHECK_EQUAL(out(outX,yy,outZ), 3);
        BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 3);
    }
    
    // Get row (1,:,1)
    out = downsample2(rle, allCells+Vector3i(1,0,1), outCells);
    outSlow = downsample2_slow(rle, allCells+Vector3i(1,0,1), outCells);
    BOOST_CHECK(out == outSlow);
    for (int yy = -5; yy <= 4; yy++)
    {
        BOOST_CHECK_EQUAL(out(outX,yy,outZ), 4);
        BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 4);
    }
}

BOOST_AUTO_TEST_CASE( TestDownsample2Wrap2x2x2 )
{
    Rect3i allCells(0,0,0,1,1,1);
    Rect3i outCells(0,0,0,0,0,0);
    
    DynamicRLE3<int> rle;
    rle.mark(0,0,0,0,0,0,1);    // rle(0,0,0) = 1
    rle.mark(1,1,0,1,1,0,2);    // rle(1,1,0) = 2
    rle.mark(0,1,1,0,1,1,3);    // rle(0,1,1) = 3
    
    Vector3i shifts[] = { Vector3i(-2,0,0), Vector3i(2,0,0),
        Vector3i(0,-2,0), Vector3i(0,2,0),
        Vector3i(0,0,-2), Vector3i(0,0,2) };
    
    DynamicRLE3<int> out;
    DynamicRLE3<int> outSlow;
    
    for (int ss = 0; ss < 6; ss++)
    {
        BOOST_TEST_MESSAGE("Shift " << shifts[ss]);
        
        // Get (0,0,0)
        out = downsample2Wrap(rle, allCells+shifts[ss], allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss], allCells,
            outCells);
        BOOST_CHECK(out == outSlow);
        BOOST_CHECK_EQUAL(out(0,0,0), 1);
        BOOST_CHECK_EQUAL(outSlow(0,0,0), 1);
        
        // Get (1,0,0)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(1,0,0),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(1,0,0),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        BOOST_CHECK(!out.isMarked(0,0,0));
        BOOST_CHECK(!outSlow.isMarked(0,0,0));
        
        // Get (0,1,0)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(0,1,0),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(0,1,0),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        BOOST_CHECK(!out.isMarked(0,0,0));
        BOOST_CHECK(!outSlow.isMarked(0,0,0));
        
        // Get (1,1,0)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(1,1,0),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(1,1,0),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        BOOST_CHECK_EQUAL(out(0,0,0), 2);
        BOOST_CHECK_EQUAL(outSlow(0,0,0), 2);
        
        // Get (0,0,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(0,0,1),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(0,0,1),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        BOOST_CHECK(!out.isMarked(0,0,0));
        BOOST_CHECK(!outSlow.isMarked(0,0,0));
        
        // Get (1,0,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(1,0,1),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(1,0,1),
            allCells, outCells);
        BOOST_CHECK(!out.isMarked(0,0,0));
        BOOST_CHECK(!outSlow.isMarked(0,0,0));
        
        // Get (0,1,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(0,1,1),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(0,1,1),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        BOOST_CHECK_EQUAL(out(0,0,0), 3);
        BOOST_CHECK_EQUAL(outSlow(0,0,0), 3);
        
        // Get (1,1,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(1,1,1),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(1,1,1),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        BOOST_CHECK(!out.isMarked(0,0,0));
        BOOST_CHECK(!outSlow.isMarked(0,0,0));
    }
}


BOOST_AUTO_TEST_CASE( TestDownsample2WrapX )
{
    int outY = -3, outZ = 10;
    Rect3i allCells(-10,0,0,9,1,1);
    Rect3i outCells(-5,outY,outZ,4,outY,outZ);
    
    DynamicRLE3<int> rle;
    
    rle.mark(-10,0,0,9,0,0,1);
    rle.mark(-10,1,0,9,1,0,2);
    rle.mark(-10,0,1,9,0,1,3);
    rle.mark(-10,1,1,9,1,1,4);
    
    Vector3i shifts[] = { Vector3i(-2,0,0), Vector3i(2,0,0),
        Vector3i(0,-2,0), Vector3i(0,2,0),
        Vector3i(0,0,-2), Vector3i(0,0,2) };
    
    DynamicRLE3<int> out, outSlow;
    
    for (int ss = 0; ss < 6; ss++)
    {
        // Get row (:,0,0)
        out = downsample2Wrap(rle, allCells+shifts[ss], allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss], allCells,
            outCells);
        BOOST_CHECK(out == outSlow);
        for (int xx = -5; xx <= 4; xx++)
        {
            BOOST_CHECK_EQUAL(out(xx,outY,outZ), 1);
            BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 1);
        }
            
        // Get row (:,1,0)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(0,1,0),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(0,1,0),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        for (int xx = -5; xx <= 4; xx++)
        {
            BOOST_CHECK_EQUAL(out(xx,outY,outZ), 2);
            BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 2);
        }
        
        // Get row (:,0,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(0,0,1),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(0,0,1),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        for (int xx = -5; xx <= 4; xx++)
        {
            BOOST_CHECK_EQUAL(out(xx,outY,outZ), 3);
            BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 3);
        }
        
        // Get row (:,1,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(0,1,1),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(0,1,1),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        for (int xx = -5; xx <= 4; xx++)
        {
            BOOST_CHECK_EQUAL(out(xx,outY,outZ), 4);
            BOOST_CHECK_EQUAL(outSlow(xx,outY,outZ), 4);
        }
    }
}

BOOST_AUTO_TEST_CASE( TestDownsample2WrapY )
{
    int outX = -100, outZ = -9;
    Rect3i allCells(0,-10,0,1,9,1);
    Rect3i outCells(outX,-5,outZ,outX,4,outZ);
    
    DynamicRLE3<int> rle;
    
    rle.mark(0,-10,0,0,9,0,1);
    rle.mark(1,-10,0,1,9,0,2);
    rle.mark(0,-10,1,0,9,1,3);
    rle.mark(1,-10,1,1,9,1,4);
    
    Vector3i shifts[] = { Vector3i(-2,0,0), Vector3i(2,0,0),
        Vector3i(0,-2,0), Vector3i(0,2,0),
        Vector3i(0,0,-2), Vector3i(0,0,2) };
    
    DynamicRLE3<int> out, outSlow;
    
    for (int ss = 0; ss < 6; ss++)
    {
        // Get row (0,:,0)
        out = downsample2Wrap(rle, allCells+shifts[ss], allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss], allCells,
            outCells);
        BOOST_CHECK(out == outSlow);
        for (int yy = -5; yy <= 4; yy++)
        {
            BOOST_CHECK_EQUAL(out(outX,yy,outZ), 1);
            BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 1);
        }
            
        // Get row (1,:,0)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(1,0,0),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(1,0,0),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        for (int yy = -5; yy <= 4; yy++)
        {
            BOOST_CHECK_EQUAL(out(outX,yy,outZ), 2);
            BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 2);
        }
        
        // Get row (0,:,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(0,0,1),
            allCells, outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(0,0,1),
            allCells,outCells);
        BOOST_CHECK(out == outSlow);
        for (int yy = -5; yy <= 4; yy++)
        {
            BOOST_CHECK_EQUAL(out(outX,yy,outZ), 3);
            BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 3);
        }
        
        // Get row (1,:,1)
        out = downsample2Wrap(rle, allCells+shifts[ss]+Vector3i(1,0,1),
            allCells,outCells);
        outSlow = downsample2Wrap_slow(rle, allCells+shifts[ss]+Vector3i(1,0,1),
            allCells, outCells);
        BOOST_CHECK(out == outSlow);
        for (int yy = -5; yy <= 4; yy++)
        {
            BOOST_CHECK_EQUAL(out(outX,yy,outZ), 4);
            BOOST_CHECK_EQUAL(outSlow(outX,yy,outZ), 4);
        }
    }
}


BOOST_AUTO_TEST_CASE( TestDownsample2ClampX )
{
    Rect3i allCells(0,0,0,1,1,1);
    Rect3i outCells(0,0,0,0,0,0);
    
    DynamicRLE3<int> rle;
    rle.mark(0,0,0,0,0,0,1);    // rle(0,0,0) = 1
    rle.mark(1,1,0,1,1,0,2);    // rle(1,1,0) = 2
    rle.mark(0,1,1,0,1,1,3);    // rle(0,1,1) = 3
        
    DynamicRLE3<int> out;
    DynamicRLE3<int> outSlow;
    
    vector<Vector3i> shifts;
    vector<int> values;
    static const int UNMARKED = -1;
    
    shifts.push_back(Vector3i(-10,0,0));
    shifts.push_back(Vector3i(-10,-1,0));
    shifts.push_back(Vector3i(-10,1,0));
    shifts.push_back(Vector3i(-10,0,-1));
    shifts.push_back(Vector3i(-10,0,1));
    
    values.push_back(1);
    values.push_back(1);
    values.push_back(UNMARKED);
    values.push_back(1);
    values.push_back(UNMARKED);
    
    shifts.push_back(Vector3i(0,10,0));
    shifts.push_back(Vector3i(-1,10,0));
    shifts.push_back(Vector3i(1,10,0));
    shifts.push_back(Vector3i(0,10,-1));
    shifts.push_back(Vector3i(0,10,1));
    
    values.push_back(UNMARKED);
    values.push_back(UNMARKED);
    values.push_back(2);
    values.push_back(UNMARKED);
    values.push_back(3);
    
    shifts.push_back(Vector3i(0,0,-10));
    shifts.push_back(Vector3i(-1,0,-10));
    shifts.push_back(Vector3i(1,0,-10));
    shifts.push_back(Vector3i(0,-1,-10));
    shifts.push_back(Vector3i(0,1,-10));
    
    values.push_back(1);
    values.push_back(1);
    values.push_back(UNMARKED);
    values.push_back(1);
    values.push_back(UNMARKED);
    
    for (int ss = 0; ss < shifts.size(); ss++)
    {
        out = downsample2Clamp(rle, allCells+shifts[ss], allCells, outCells);
        outSlow = downsample2Clamp_slow(rle, allCells+shifts[ss], allCells,
            outCells);
        BOOST_CHECK(out == outSlow);
        if (values[ss] != UNMARKED)
            BOOST_CHECK_EQUAL(out(0,0,0), values[ss]);
        else
            BOOST_CHECK(!out.isMarked(0,0,0));
    }
}

BOOST_AUTO_TEST_CASE( SumDownsample2Clamp_AllFilled )
{
    Rect3i allCells(-1,-1,-1,0,0,0);
    Vector3i p(10,-10,-10);
    Rect3i outCells(p[0], p[1], p[2], p[0], p[1], p[2]); // put in a weird place
    
    DynamicRLE3<int> rle;
    
    // Test a case where all the cells are marked.
    
    int mark = 0;
    for (int zz = -1; zz <= 0; zz++)
    for (int yy = -1; yy <= 0; yy++)
    for (int xx = -1; xx <= 0; xx++)
        rle.mark(xx,yy,zz,xx,yy,zz, mark++);
    
    DynamicRLE3<int> out;
    
    // Sum all eight cells
    out = sumDownsample2Clamp(rle, allCells, allCells, outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 28);
    
    // Shift x - 1
    out = sumDownsample2Clamp(rle, allCells+Vector3i(-1,0,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 24);
    
    // Shift x - 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(-10,0,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 24);
    
    // Shift x + 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(10,0,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 32);
    
    // Shift y - 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,-10,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 20);
    
    // Shift y + 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,10,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 36);
    
    // Shift z - 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,0,-10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 12);
    
    // Shift z + 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,0,10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 44);
}

BOOST_AUTO_TEST_CASE( SumDownsample2Clamp_SomeFilled )
{
    Rect3i allCells(-1,-1,-1,0,0,0);
    Vector3i p(10,-10,-10);
    Rect3i outCells(p[0], p[1], p[2], p[0], p[1], p[2]); // put in a weird place
    
    DynamicRLE3<int> rle;
    
    // Test a case where all the cells are marked.
    
    int mark = 0;
    for (int zz = -1; zz <= 0; zz++)
    for (int yy = -1; yy <= 0; yy++)
    for (int xx = -1; xx <= 0; xx++)
        rle.mark(xx,yy,zz,xx,yy,zz, mark++);
    rle.erase(0,0,0,0,0,0);
    rle.erase(-1,-1,-1,-1,-1,-1);
    BOOST_CHECK_EQUAL(rle.numRuns(), 6);
    
    DynamicRLE3<int> out;
    
    // Sum all eight cells
    out = sumDownsample2Clamp(rle, allCells, allCells, outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 21);
    
    // Shift x - 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(-10,0,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 24);
    
    // Shift x + 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(10,0,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 18);
    
    // Shift y - 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,-10,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 20);
    
    // Shift y + 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,10,0), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 22);
    
    // Shift z - 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,0,-10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 12);
    
    // Shift z + 10
    out = sumDownsample2Clamp(rle, allCells+Vector3i(0,0,10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 30);
    
    // Shift + (-10,-10,-10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(-10,-10,-10), allCells,
        outCells);
    BOOST_CHECK(!out.isMarked(p[0], p[1], p[2]));
    
    // Shift + (10, -10, -10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(10,-10,-10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 8*1);
    
    // Shift + (-10, 10, -10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(-10,10,-10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 8*2);
    
    // Shift + (10, 10, -10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(10,10,-10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 8*3);
    
    // Shift + (-10, -10, 10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(-10,-10,10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 8*4);
    
    // Shift + (10, -10, 10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(10,-10,10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 8*5);
    
    // Shift + (-10, 10, 10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(-10,10,10), allCells,
        outCells);
    BOOST_CHECK_EQUAL(out(p[0], p[1], p[2]), 8*6);
    
    // Shift + (10, 10, 10)
    out = sumDownsample2Clamp(rle, allCells+Vector3i(10,10,10), allCells,
        outCells);
    BOOST_CHECK(!out.isMarked(p[0], p[1], p[2]));
}





