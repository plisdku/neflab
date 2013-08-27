/*
 *  CGALUtilities.h
 *  Snapdragon
 *
 *  Created by Paul Hansen on 3/9/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifndef _CGALUTILITIES_
#define _CGALUTILITIES_

#include <CGAL/number_utils.h>
#include "utility/VectorMatrix.h"
#include "utility/VectorMatrix2.h"

template<class CGALPoint_3>
CGALPoint_3 toCGAL(const Vector3d & v)
{
    return CGALPoint_3(v[0], v[1], v[2]);
}

template<class CGALPoint_3>
Vector3d toVector3d(const CGALPoint_3 & p)
{
    return Vector3d(CGAL::to_double(p[0]), CGAL::to_double(p[1]),
        CGAL::to_double(p[2]));
}


#endif
