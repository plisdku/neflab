/*
 *  RLEOperations.h
 *  rleFunctions
 *
 *  Created by Paul Hansen on 2/24/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifndef _RLEOPERATIONS_
#define _RLEOPERATIONS_

#include "rle/DynamicRLE3.h"
#include "rle/SupportRegion3.h"
#include "geometry.h"

template<class T>
RLE::DynamicRLE3<T> downsample2(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i outRegion);

template<class T>
RLE::DynamicRLE3<T> downsample2_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i outRegion);

// Provide the size of inRLE in inBounds.  The iteration will be over
// inRegion, and all samples will be read from wrapped coordinates.
template<class T>
RLE::DynamicRLE3<T> downsample2Wrap(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion);

// Provide the size of inRLE in inBounds.  The iteration will be over
// inRegion, and all samples will be read from wrapped coordinates.
template<class T>
RLE::DynamicRLE3<T> downsample2Wrap_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion);

template<class T>
RLE::DynamicRLE3<T> downsample2Clamp(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion);

template<class T>
RLE::DynamicRLE3<T> downsample2Clamp_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion);

template<class T>
void extrude_slow(RLE::DynamicRLE3<T> & inRLE, Rect3i innerRect,
    Rect3i outerRect);

// Sum 2x2x2 blocks of cells to obtain one output cell.
// If input dimensions are NxMxP, output dimensions are N/2 x M/2 x P/2.
template<class T>
RLE::DynamicRLE3<T> sumDownsample2Clamp(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion);

// Average 2x2x2 blocks of cells to obtain one output cell.
// If input dimensions are NxMxP, output dimensions are N/2 x M/2 x P/2.
template<class T>
RLE::DynamicRLE3<T> averageDownsample2Clamp(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion);

template<class T>
RLE::DynamicRLE3<T> averageDownsample2Clamp_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion);

template<class T>
RLE::DynamicRLE3<T> smoothClampN(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, int nSmoothing);
    
Rect3i boundingBox(const RLE::SupportRegion3 & rle);

#include "RLEOperations-inl.h"
#endif

