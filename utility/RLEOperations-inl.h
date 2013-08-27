/*
 *  RLEOperations.cpp
 *  rleFunctions
 *
 *  Created by Paul Hansen on 2/24/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifdef _RLEOPERATIONS_
#include "RLEOperations.h"


template<class T>
RLE::DynamicRLE3<T> downsample2(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i outRegion)
{
    if (inRegion.num() != outRegion.num()*2)
        throw(std::logic_error("Downsample region does not have even extent"));
    RLE::DynamicRLE3<T> outRLE(inRLE.capacity());
    
    // Choose even-numbered samples.
    int kkIn, kkOut, jjIn, jjOut, iiIn, iiOut;
    typename RLE::DynamicRLE3<T>::ConstIterator itr;
    typename RLE::DynamicRLE3<T>::ConstIterator inEnd = inRLE.end();
    itr = inRLE.location(inRegion.p1[0], inRegion.p1[1], inRegion.p1[2]);
    
    for (kkIn = inRegion.p1[2], kkOut = outRegion.p1[2]; kkIn <= inRegion.p2[2];
        kkIn += 2, kkOut++)
    for (jjIn = inRegion.p1[1], jjOut = outRegion.p1[1]; jjIn <= inRegion.p2[1];
        jjIn += 2, jjOut++)
    {
        iiIn = inRegion.p1[0];
        itr.seek(inRLE.linearIndex(iiIn, jjIn, kkIn));
        
        while (iiIn <= inRegion.p2[0] && itr < inEnd)
        {
            iiOut = outRegion.p1[0] + (iiIn - inRegion.p1[0])/2;
            int lengthRemaining = itr.runEnd() - itr.position();
            int iiInLast = iiIn + 2*(lengthRemaining/2);
            if (iiInLast >= inRegion.p2[0])
                iiInLast = inRegion.p2[0]-1; // stay on even numbers from start
            assert((iiInLast - inRegion.p1[0])%2 == 0);
            
            int iiOutLast = outRegion.p1[0] + (iiInLast - inRegion.p1[0])/2;
            assert(iiOutLast >= iiOut);
            assert(2*(iiOutLast - outRegion.p1[0]) == (iiInLast - inRegion.p1[0]));
            
            if (itr.isMarked())
            {
                outRLE.mark(iiOut, jjOut, kkOut, iiOutLast, jjOut, kkOut,
                    itr.mark());
            }
            
            iiOut = iiOutLast + 1;
            iiIn = iiInLast + 2;
            itr.seek(inRLE.linearIndex(iiIn, jjIn, kkIn));
        }
    }
    
    return outRLE;
}

template<class T>
RLE::DynamicRLE3<T> downsample2_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i outRegion)
{
    if (inRegion.num() != outRegion.num()*2)
        throw(std::logic_error("Downsample region does not have even extent"));
    RLE::DynamicRLE3<T> outRLE(inRLE.capacity());
    
    // Choose even-numbered samples.
    int kkIn, kkOut, jjIn, jjOut, iiIn, iiOut;
    
    for (kkIn = inRegion.p1[2], kkOut = outRegion.p1[2]; kkIn <= inRegion.p2[2];
        kkIn += 2, kkOut++)
    for (jjIn = inRegion.p1[1], jjOut = outRegion.p1[1]; jjIn <= inRegion.p2[1];
        jjIn += 2, jjOut++)
    for (iiIn = inRegion.p1[0], iiOut = outRegion.p1[0]; iiIn <= inRegion.p2[0];
        iiIn += 2, iiOut++)
    {
        if (inRLE.isMarked(iiIn, jjIn, kkIn))
        {
            outRLE.mark(iiOut, jjOut, kkOut, iiOut, jjOut, kkOut,
                inRLE.markAt(iiIn, jjIn, kkIn));
        }
    }
        
    return outRLE;
}


template<class T>
RLE::DynamicRLE3<T> downsample2Wrap(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion)
{
    if (inRegion.num() != outRegion.num()*2)
        throw(std::logic_error("Downsample region does not have even extent"));
    if (inBounds.num() != outRegion.num()*2)
        throw(std::logic_error("Downsample bounds do not have even extent"));
    RLE::DynamicRLE3<T> outRLE(inRLE.hasDimension(0),
        inRLE.hasDimension(1), inRLE.hasDimension(2));
    
    // Choose even-numbered samples.
    typename RLE::DynamicRLE3<T>::ConstIterator itr;
    typename RLE::DynamicRLE3<T>::ConstIterator inEnd = inRLE.end();
    itr = inRLE.location(inRegion.p1[0], inRegion.p1[1], inRegion.p1[2]);
    
    Vector3i inSize = inBounds.num();
    Vector3i inPos, outPos, wrapPos;
    
    // kkIn and jjIn iterate over the y and z dimensions of inRegion.
    // We'll wrap them/clamp them to inBounds to produce inPos.
    int kkIn, jjIn;
    for (kkIn = inRegion.p1[2], outPos[2] = outRegion.p1[2];
        kkIn <= inRegion.p2[2];
        kkIn += 2, outPos[2]++)
    for (jjIn = inRegion.p1[1], outPos[1] = outRegion.p1[1];
        jjIn <= inRegion.p2[1];
        jjIn += 2, outPos[1]++)
    {
        inPos[2] = inBounds.p1[2] + (kkIn-inBounds.p1[2] + inSize[2])%inSize[2];
        inPos[1] = inBounds.p1[1] + (jjIn-inBounds.p1[1] + inSize[1])%inSize[1];
        inPos[0] = inBounds.p1[0] +
            (inRegion.p1[0] - inBounds.p1[0] + inSize[0])%inSize[0];
        
        // now inPos is completely within inBounds.
        itr.seek(inRLE.linearIndex(inPos[0], inPos[1], inPos[2]));
        
        outPos[0] = outRegion.p1[0];
        while (outPos[0] <= outRegion.p2[0])
        {
            // 1.  Determine how far the current run extends.  It should be
            //  an even number since we're downsampling by 2.
            long inRunLength;
            if (itr.runIsBoundedAtEnd())
                inRunLength = 2*((itr.runEnd() - itr.position())/2);
            else
                inRunLength = (outRegion.p2[0] - outPos[0])*2;
            assert(inRunLength%2 == 0);
            
            // 2.  Clip the run at the end of inRegion
            if (inPos[0] + inRunLength >= inBounds.p2[0])
                inRunLength = inBounds.p2[0] - inPos[0];
            inRunLength -= inRunLength%2; // make the run length even.
            
            // 2b. Clip the run for outRegion as well.
            if (outPos[0] + inRunLength/2 >= outRegion.p2[0])
                inRunLength = 2*(outRegion.p2[0] - outPos[0]);
            
            assert(inRunLength%2 == 0);
            assert(outPos[0] + inRunLength/2 <= outRegion.p2[0]);
            
            // 3.  Mark the run!
            if (itr.isMarked())
            {
                outRLE.mark(outPos[0], outPos[1], outPos[2],
                    outPos[0]+inRunLength/2, outPos[1], outPos[2], itr.mark());
            }
            
            // 4.  Advance.  Wrap inPos to stay within inBounds as needed.
            inPos[0] += inRunLength + 2;
            if (inPos[0] > inBounds.p2[0])
                inPos[0] -= inSize[0];
            assert(inBounds.encloses(inPos));
            itr.seek(inRLE.linearIndex(inPos[0], inPos[1], inPos[2]));
            outPos[0] += inRunLength/2 + 1;
        }
    }
    
    return outRLE;
}

template<class T>
RLE::DynamicRLE3<T> downsample2Wrap_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion)
{
    if (inRegion.num() != outRegion.num()*2)
        throw(std::logic_error("Downsample region does not have even extent"));
    RLE::DynamicRLE3<T> outRLE(inRLE.capacity());
    
    Vector3i inSize = inBounds.num();
    
    // Choose even-numbered samples.
    int kkIn, kkOut, jjIn, jjOut, iiIn, iiOut;
    
    for (kkIn = inRegion.p1[2], kkOut = outRegion.p1[2]; kkIn <= inRegion.p2[2];
        kkIn += 2, kkOut++)
    for (jjIn = inRegion.p1[1], jjOut = outRegion.p1[1]; jjIn <= inRegion.p2[1];
        jjIn += 2, jjOut++)
    for (iiIn = inRegion.p1[0], iiOut = outRegion.p1[0]; iiIn <= inRegion.p2[0];
        iiIn += 2, iiOut++)
    {
        int ii = inBounds.p1[0] + (iiIn - inBounds.p1[0] + inSize[0])%inSize[0];
        int jj = inBounds.p1[1] + (jjIn - inBounds.p1[1] + inSize[1])%inSize[1];
        int kk = inBounds.p1[2] + (kkIn - inBounds.p1[2] + inSize[2])%inSize[2];
        assert(inBounds.encloses(ii,jj,kk));
        if (inRLE.isMarked(ii,jj,kk))
        {
            outRLE.mark(iiOut, jjOut, kkOut, iiOut, jjOut, kkOut,
                inRLE.markAt(ii, jj, kk));
        }
    }
        
    return outRLE;
}

template<class T>
RLE::DynamicRLE3<T> downsample2Clamp(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion)
{
    for (int xyz = 0; xyz < 3; xyz++)
    if (inRLE.hasDimension(xyz) && inRegion.num(xyz) != outRegion.num(xyz)*2)
        throw(std::logic_error("Downsample region does not have even extent"));
    
    RLE::DynamicRLE3<T> outRLE(inRLE.hasDimension(0),
        inRLE.hasDimension(1), inRLE.hasDimension(2));
    
    // Choose even-numbered samples.
    typename RLE::DynamicRLE3<T>::ConstIterator itr;
    typename RLE::DynamicRLE3<T>::ConstIterator inEnd = inRLE.end();
    itr = inRLE.location(inRegion.p1[0], inRegion.p1[1], inRegion.p1[2]);
    
    Vector3i inSize = inBounds.num();
    Vector3i inPos, outPos, wrapPos;
    
    // kkIn and jjIn iterate over the y and z dimensions of inRegion.
    // We'll wrap them/clamp them to inBounds to produce inPos.
    int kkIn, jjIn, iiIn;
    for (kkIn = inRegion.p1[2], outPos[2] = outRegion.p1[2];
        kkIn <= inRegion.p2[2];
        kkIn += 2, outPos[2]++)
    for (jjIn = inRegion.p1[1], outPos[1] = outRegion.p1[1];
        jjIn <= inRegion.p2[1];
        jjIn += 2, outPos[1]++)
    {
        iiIn = inRegion.p1[0];
        inPos = clip(inBounds, Vector3i(iiIn, jjIn, kkIn));
        // now inPos is completely within inBounds.
        itr.seek(inRLE.linearIndex(inPos[0], inPos[1], inPos[2]));
        
        outPos[0] = outRegion.p1[0];
        while (outPos[0] <= outRegion.p2[0])
        {
            // 1.  Determine how far the current run extends.  It should be
            //  an even number since we're downsampling by 2.
            long inRunLength;
            if (itr.runIsBoundedAtEnd())
                inRunLength = 2*((itr.runEnd() - itr.position())/2);
            else
                inRunLength = (outRegion.p2[0] - outPos[0])*2;
            assert(inRunLength%2 == 0);
            
            // 2.  Clip the run at the end of inRegion
            if (inPos[0] + inRunLength >= inBounds.p2[0])
                inRunLength = inBounds.p2[0] - inPos[0];
            inRunLength -= inRunLength%2; // make the run length even.
            
            // 2b. Clip the run for outRegion as well.
            if (outPos[0] + inRunLength/2 >= outRegion.p2[0])
                inRunLength = 2*(outRegion.p2[0] - outPos[0]);
            
            assert(inRunLength%2 == 0);
            assert(outPos[0] + inRunLength/2 <= outRegion.p2[0]);
            
            // 3.  Mark the run!
            if (itr.isMarked())
            {
                outRLE.mark(outPos[0], outPos[1], outPos[2],
                    outPos[0]+inRunLength/2, outPos[1], outPos[2], itr.mark());
            }
            
            // 4.  Advance.  Wrap inPos to stay within inBounds as needed.
            iiIn += inRunLength + 2;
            inPos = clip(inBounds, Vector3i(iiIn, jjIn, kkIn));
            assert(inBounds.encloses(inPos));
            itr.seek(inRLE.linearIndex(inPos[0], inPos[1], inPos[2]));
            outPos[0] += inRunLength/2 + 1;
        }
    }
    
    return outRLE;
}

template<class T>
RLE::DynamicRLE3<T> downsample2Clamp_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion)
{
    for (int xyz = 0; xyz < 3; xyz++)
    if (inRLE.hasDimension(xyz) && inRegion.num(xyz) != outRegion.num(xyz)*2)
        throw(std::logic_error("Downsample region does not have even extent"));
    RLE::DynamicRLE3<T> outRLE(inRLE.capacity());
    
    Vector3i inSize = inBounds.num();
    
    // Choose even-numbered samples.
    int kkIn, kkOut, jjIn, jjOut, iiIn, iiOut;
    
    for (kkIn = inRegion.p1[2], kkOut = outRegion.p1[2]; kkIn <= inRegion.p2[2];
        kkIn += 2, kkOut++)
    for (jjIn = inRegion.p1[1], jjOut = outRegion.p1[1]; jjIn <= inRegion.p2[1];
        jjIn += 2, jjOut++)
    for (iiIn = inRegion.p1[0], iiOut = outRegion.p1[0]; iiIn <= inRegion.p2[0];
        iiIn += 2, iiOut++)
    {
        Vector3i pos = clip(inBounds, Vector3i(iiIn, jjIn, kkIn));
        
        assert(inBounds.encloses(pos));
        if (inRLE.isMarked(pos[0], pos[1], pos[2]))
        {
            outRLE.mark(iiOut, jjOut, kkOut, iiOut, jjOut, kkOut,
                inRLE.markAt(pos[0], pos[1], pos[2]));
        }
    }
        
    return outRLE;
}

template<class T>
void extrude_slow(RLE::DynamicRLE3<T> & rle, Rect3i innerRect,
    Rect3i outerRect)
{
    if (!outerRect.encloses(innerRect))
        throw(std::logic_error("Outer rect does not enclose inner rect"));
    
    int ii, jj, kk;
    
    // Extrude X-normal faces
    // High X
    if (innerRect.p2[0] != outerRect.p2[0])
    for (kk = innerRect.p1[2]; kk <= innerRect.p2[2]; kk++)
    for (jj = innerRect.p1[1]; jj <= innerRect.p2[1]; jj++)
    if (rle.isMarked(innerRect.p2[0], jj, kk))
    {
        rle.mark(innerRect.p2[0]+1, jj, kk, outerRect.p2[0], jj, kk,
            rle.markAt(innerRect.p2[0], jj, kk));
    }
    // Low X
    if (innerRect.p1[0] != outerRect.p1[0])
    for (kk = innerRect.p1[2]; kk <= innerRect.p2[2]; kk++)
    for (jj = innerRect.p1[1]; jj <= innerRect.p2[1]; jj++)
    if (rle.isMarked(innerRect.p1[0], jj, kk))
    {
        rle.mark(outerRect.p1[0], jj, kk, innerRect.p1[0]-1, jj, kk,
            rle.markAt(innerRect.p1[0], jj, kk));
    }
    
    // Extrude Y-normal faces, including the new bits from the X extrusion
    // High Y
    if (innerRect.p2[1] != outerRect.p2[1])
    for (kk = innerRect.p1[2]; kk <= innerRect.p2[2]; kk++)
    for (ii = outerRect.p1[0]; ii <= outerRect.p2[0]; ii++)
    if (rle.isMarked(ii, innerRect.p2[1], kk))
    {
        rle.mark(ii, innerRect.p2[1]+1, kk, ii, outerRect.p2[1], kk,
            rle.markAt(ii, innerRect.p2[1], kk));
    }
    // Low Y
    if (innerRect.p1[1] != outerRect.p1[1])
    for (kk = innerRect.p1[2]; kk <= innerRect.p2[2]; kk++)
    for (ii = outerRect.p1[0]; ii <= outerRect.p2[0]; ii++)
    if (rle.isMarked(ii, innerRect.p1[1], kk))
    {
        rle.mark(ii, outerRect.p1[1], kk, ii, innerRect.p1[1]-1, kk,
            rle.markAt(ii, innerRect.p1[1], kk));
    }
    
    // Extrude Z-normal faces, including the faces from the X and Y extrusion
    // High Z
    if (innerRect.p2[2] != outerRect.p2[2])
    for (jj = outerRect.p1[1]; jj <= outerRect.p2[1]; jj++)
    for (ii = outerRect.p1[0]; ii <= outerRect.p2[0]; ii++)
    if (rle.isMarked(ii, jj, innerRect.p2[2]))
    {
        rle.mark(ii, jj, innerRect.p2[2]+1, ii, jj, outerRect.p2[2],
            rle.markAt(ii, jj, innerRect.p2[2]));
    }
    // Low Z
    if (innerRect.p1[2] != outerRect.p1[2])
    for (jj = outerRect.p1[1]; jj <= outerRect.p2[1]; jj++)
    for (ii = outerRect.p1[0]; ii <= outerRect.p2[0]; ii++)
    if (rle.isMarked(ii, jj, innerRect.p1[2]))
    {
        rle.mark(ii, jj, outerRect.p1[2], ii, jj, innerRect.p1[2]-1,
            rle.markAt(ii, jj, innerRect.p1[2]));
    }
}

template<class T>
RLE::DynamicRLE3<T> sumDownsample2Clamp(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion)
{
    RLE::DynamicRLE3<T> sum;
    sum.dimensions(inRLE.hasDimension(0), inRLE.hasDimension(1),
        inRLE.hasDimension(2));
    
    Vector3i shift;
    for (shift[0] = 0; shift[0] < 2; shift[0]++)
    for (shift[1] = 0; shift[1] < 2; shift[1]++)
    for (shift[2] = 0; shift[2] < 2; shift[2]++)
    {
        sum = sum +
            downsample2Clamp(inRLE, inRegion+shift, inBounds, outRegion);
    }
    
    return sum;
}

template<class T>
RLE::DynamicRLE3<T> averageDownsample2Clamp(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion)
{
    RLE::DynamicRLE3<T> average;
    average.dimensions(inRLE.hasDimension(0), inRLE.hasDimension(1),
        inRLE.hasDimension(2));
    
    Vector3i shift;
    for (shift[0] = 0; shift[0] < 2; shift[0]++)
    for (shift[1] = 0; shift[1] < 2; shift[1]++)
    for (shift[2] = 0; shift[2] < 2; shift[2]++)
    {
        average = average +
            downsample2Clamp(inRLE, inRegion+shift, inBounds, outRegion);
    }
    average = 0.125*average;
    
    return average;
}

template<class T>
RLE::DynamicRLE3<T> averageDownsample2Clamp_slow(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, Rect3i outRegion)
{
    for (int xyz = 0; xyz < 3; xyz++)
    if (inRLE.hasDimension(xyz) && inRegion.num(xyz) != outRegion.num(xyz)*2)
        throw(std::logic_error("Downsample region does not have even extent"));
    RLE::DynamicRLE3<T> outRLE(inRLE.capacity());
    
    // Choose even-numbered samples.
    int kkIn, kkOut, jjIn, jjOut, iiIn, iiOut;
    
    for (kkIn = inRegion.p1[2], kkOut = outRegion.p1[2]; kkIn <= inRegion.p2[2];
        kkIn += 2, kkOut++)
    for (jjIn = inRegion.p1[1], jjOut = outRegion.p1[1]; jjIn <= inRegion.p2[1];
        jjIn += 2, jjOut++)
    for (iiIn = inRegion.p1[0], iiOut = outRegion.p1[0]; iiIn <= inRegion.p2[0];
        iiIn += 2, iiOut++)
    {
        T total = 0;
        Vector3i shift;
        for (shift[0] = 0; shift[0] < 2; shift[0]++)
        for (shift[1] = 0; shift[1] < 2; shift[1]++)
        for (shift[2] = 0; shift[2] < 2; shift[2]++)
        {
            Vector3i pos(iiIn + shift[0], jjIn + shift[1], kkIn + shift[2]);
            pos = clip(inBounds, pos);
            if (inRLE.isMarked(pos[0], pos[1], pos[2]))
                total += inRLE.markAt(pos[0], pos[1], pos[2]);
        }
        
        if (total != 0)
            outRLE.mark(iiOut, jjOut, kkOut, iiOut, jjOut, kkOut, 0.125*total);
    }
    
    return outRLE;
}

template<class T>
RLE::DynamicRLE3<T> smoothClampN(const RLE::DynamicRLE3<T> & inRLE,
    Rect3i inRegion, Rect3i inBounds, int nSmoothing)
{
    RLE::DynamicRLE3<T> outRLE(inRLE.capacity());
    
    int kkIn, kkOut, jjIn, jjOut, iiIn, iiOut;
    
    double weight = 1.0/pow(2*nSmoothing+1, 3);
    
    for (int kk = inRegion.p1[2]; kk <= inRegion.p2[2]; kk++)
    for (int jj = inRegion.p1[1]; jj <= inRegion.p2[1]; jj++)
    for (int ii = inRegion.p1[0]; ii <= inRegion.p2[0]; ii++)
    {
        T total = T();
        int numSummed = 0;
        for (int kkk = kk - nSmoothing; kkk <= kk + nSmoothing; kkk++)
        for (int jjj = jj - nSmoothing; jjj <= jj + nSmoothing; jjj++)
        for (int iii = ii - nSmoothing; iii <= ii + nSmoothing; iii++)
        {
            Vector3i pos(iii, jjj, kkk);
            pos = clip(inBounds, pos);
            if (inRLE.isMarked(pos[0], pos[1], pos[2]))
            {
                total += inRLE.markAt(pos[0], pos[1], pos[2]);
                numSummed++;
            }
        }
        
        if (numSummed != 0)
            outRLE.mark(ii, jj, kk, weight*total);
    }
    
    return outRLE;
}


inline Rect3i
boundingBox(const RLE::SupportRegion3 & rle)
{
    Rect3i box(std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min());
    
    RLE::SupportRegion3::ConstIterator itr;
    
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        Vector3i p1, p2;
        rle.cartesianCoordinates(itr.runStart(), p1[0], p1[1], p1[2]);
        rle.cartesianCoordinates(itr.runEnd(), p2[0], p2[1], p2[2]);
        box.p1 = vec_min(box.p1, p1);
        box.p2 = vec_max(box.p2, p2);
    }
    
    return box;
}


#endif