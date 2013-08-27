/*
 *  PointFacetDistance-inl.h
 *  shells
 *
 *  Created by Paul C Hansen on 6/19/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 */

#ifndef POINT_FACET_DISTANCE_INL_H
#define POINT_FACET_DISTANCE_INL_H

#include <vector>
#include <limits>
#include <CGAL/Polygon_2_algorithms.h>
#include "utility/geometry.h"

namespace PointFacetDistance
{


// Return the square of the distance to the nearest point of the polygon,
// including the interior.
template<class Kernel>
double distanceSquared(typename Kernel::Point_3 point,
    const std::vector<typename Kernel::Point_3> polygon3,
    typename Kernel::Plane_3 supportingPlane,
    const Kernel & kernel_unused);
    
// Return the square of the distance to the nearest point on the edge of the 
// polygon.
template<class Kernel>
double distanceToEdgeSquared(typename Kernel::Point_3 point,
    const std::vector<typename Kernel::Point_3> polygon3,
    const Kernel & kernel_unused);

template<class Kernel>
bool onUnboundedSide(typename Kernel::Point_3 point,
    const std::vector<typename Kernel::Point_3> polygon3,
    typename Kernel::Plane_3 supportingPlane,
    const Kernel & kernel_unused);

template<class Kernel>
bool onUnboundedSide(typename Kernel::Point_3 point,
    const std::vector<std::vector<typename Kernel::Point_3> > contours,
    int outerContourIndex,
    typename Kernel::Plane_3 supportingPlane,
    const Kernel & kernel_unused);
    
template<class Kernel>
double distanceToPlaneSquared(typename Kernel::Point_3 point,
    typename Kernel::Plane_3 plane,
    const Kernel & kernel_unused);

//
// ----- Internals
//



namespace Local
{

template<class Point3>
Vector3d toVector3d(const Point3 & cgalPoint)
{
    return Vector3d(CGAL::to_double(cgalPoint[0]),
        CGAL::to_double(cgalPoint[1]),
        CGAL::to_double(cgalPoint[2]));
}

template<class Point3>
double distanceSquared(const Point3 & point, const Point3 & edgeBegin,
    const Point3 & edgeEnd)
{
    return pointSegmentDistanceSquared(toVector3d(point),
        toVector3d(edgeBegin), toVector3d(edgeEnd));
}


} // namespace Local

template<class Kernel>
double distanceToPlaneSquared(typename Kernel::Point_3 point,
    typename Kernel::Plane_3 plane,
    const Kernel & kernel_unused)
{
    using namespace Local;
    
    Vector3d normalVector = toVector3d(plane.orthogonal_vector());
    Vector3d pointOnPlane = toVector3d(plane.point());
    
    double dotProd = dot(normalVector, toVector3d(point) - pointOnPlane);
    double distSquared = dotProd*dotProd / sumSquares(normalVector);
    
    return distSquared;
}

// A CGAL Halffacet_cycle_const_iterator is castable to
// SHalfedge_const_handle if it represenst an shalfedge.  Figure that out
// first.
template<class Kernel>
double distanceSquared(typename Kernel::Point_3 point,
    const std::vector<typename Kernel::Point_3> polygon3,
    typename Kernel::Plane_3 supportingPlane,
    const Kernel & kernel_unused)
{
    using namespace Local;
    
    if (onUnboundedSide(point, polygon3, supportingPlane, kernel_unused))
    {
        // The point is "outside" the polygon, so find the distance to the
        // closest edge of the polygon.
        
        return distanceToEdgeSquared(point, polygon3, Kernel());
    }
    else
    {
        return distanceToPlaneSquared(point, supportingPlane, Kernel());
    }
} // distanceSquared()


template<class Kernel>
double distanceToEdgeSquared(typename Kernel::Point_3 point,
    const std::vector<typename Kernel::Point_3> polygon3,
    const Kernel & kernel_unused)
{
    using namespace Local;
    
    double minDistanceSquared = std::numeric_limits<double>::max();
    for (int ee = 0; ee < polygon3.size(); ee++)
    {
        double distToEdgeSquared = distanceSquared(point,
            polygon3[ee], polygon3[(ee+1)%polygon3.size()]);
        if (distToEdgeSquared < minDistanceSquared)
            minDistanceSquared = distToEdgeSquared;
    }
    return minDistanceSquared;
} // distanceToEdgeSquared()


template<class Kernel>
bool onUnboundedSide(typename Kernel::Point_3 point,
    const std::vector<typename Kernel::Point_3> polygon3,
    typename Kernel::Plane_3 supportingPlane,
    const Kernel & kernel_unused)
{
    using namespace Local;
    
    std::vector<typename Kernel::Point_2> polygon2(polygon3.size());
    for (int nn = 0; nn < polygon3.size(); nn++)
        polygon2[nn] = supportingPlane.to_2d(polygon3[nn]);
    
    typename Kernel::Point_2 p2(supportingPlane.to_2d(point));
    
    return (CGAL::bounded_side_2(polygon2.begin(), polygon2.end(), p2, Kernel()) ==
        CGAL::ON_UNBOUNDED_SIDE);
}

template<class Kernel>
bool onUnboundedSide(typename Kernel::Point_3 point,
    const std::vector<std::vector<typename Kernel::Point_3> > contours,
    int outerContourIndex,
    typename Kernel::Plane_3 supportingPlane,
    const Kernel & kernel_unused)
{
    if (onUnboundedSide(point, contours[outerContourIndex], supportingPlane,
        Kernel()))
    {
        return true;
    }
    
    for (int nn = 0; nn < contours.size(); nn++)
    if (nn != outerContourIndex)
    {
        if (!onUnboundedSide(point, contours[nn], supportingPlane, Kernel()))
            return true;
    }
    
    return false;
}

} // namespace PointFacetDistance




#endif