/*
 *  InteriorVolumes-inl.h
 *  shells
 *
 *  Created by Paul C Hansen on 6/18/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 */

#ifndef INTERIOR_VOLUMES_INL_H
#define INTERIOR_VOLUMES_INL_H

#include <vector>
#include <limits>
#include "PointFacetDistance-inl.h"
#include "utility/geometry.h"

namespace InteriorVolumes
{

// This is the sole function offered by InteriorVolumes to the outside.
// Helper functions are in here too because of templatization.
template<class NefPolyhedron>
std::vector<typename NefPolyhedron::Volume_const_handle>
interiorVolumes(const NefPolyhedron & poly);

namespace Local
{
    bool interiorEncloses(const CGAL::Bbox_3 & outer, const CGAL::Bbox_3 & inner)
    {
        for (int ii = 0; ii < 3; ii++)
        if (outer.min(ii) >= inner.min(ii) || outer.max(ii) <= inner.max(ii))
            return false;
        
        return true;
    }
    
    bool encloses(const CGAL::Bbox_3 & outer, const CGAL::Bbox_3 & inner)
    {
        for (int ii = 0; ii < 3; ii++)
        if (outer.min(ii) > inner.min(ii) || outer.max(ii) < inner.max(ii))
            return false;
        
        return true;
    }
    
    CGAL::Bbox_3 emptyBbox()
    {
        return CGAL::Bbox_3(std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max());
    }
    
    // This class follows the Visitor pattern for shells of Nef polyhedra.
    template<class NefPolyhedron>
    class ClosestFacetFinder
    {
        typedef typename NefPolyhedron::Kernel Kernel;
        typedef typename NefPolyhedron::Halffacet_cycle_const_iterator
            HalffacetCycleConstIterator;
    public:
        ClosestFacetFinder(typename NefPolyhedron::Point_3 testPoint) :
            mPoint(testPoint),
            mMinDistanceSquared(std::numeric_limits<double>::max())
        {
        }
        
        typename NefPolyhedron::Halffacet_const_iterator
        closestFacet() const
        {
            return mClosestFacet;
        }
        
        double distanceSquared() const { return mMinDistanceSquared; }
        
        void visit(typename NefPolyhedron::Halffacet_const_handle h)
        {
            CGAL::Bbox_3 largestContourBounds(emptyBbox());
            
            std::vector<std::vector<typename Kernel::Point_3> > contours;
            int largestContourIndex;
            
            int currentContour = 0;
            
            typename NefPolyhedron::Halffacet_cycle_const_iterator itr;
            for (itr = h->facet_cycles_begin(); itr != h->facet_cycles_end();
                itr++)
            if (itr.is_shalfedge())
            {
                std::vector<typename Kernel::Point_3> polygon3;
                typename NefPolyhedron::SHalfedge_const_handle edge0(itr);
                typename NefPolyhedron::SHalfedge_const_handle currentEdge(edge0);
                
                CGAL::Bbox_3 contourBounds(emptyBbox());
                
                // PAY ATTENTION: I need to use currentEdge->prev() to get the
                // facet oriented outwards and not inwards.  It doesn't matter
                // for just finding a distance at least.
                do {
                    currentEdge = currentEdge->prev();
                    typename Kernel::Point_3 p3 =
                        currentEdge->source()->center_vertex()->point();
                    polygon3.push_back(p3);
                    contourBounds = contourBounds + p3.bbox();
                } while (currentEdge != edge0);
                
                contours.push_back(polygon3);
                if (encloses(contourBounds, largestContourBounds))
                {
                    largestContourBounds = contourBounds;
                    largestContourIndex = currentContour;
                }
                currentContour++;
            }
            
            assert(largestContourIndex == 0); // Maybe true; I want to check.
            
            double distToFacetSquared = std::numeric_limits<double>::max();
            
            if (PointFacetDistance::onUnboundedSide(mPoint, contours,
                largestContourIndex, h->plane(), Kernel()))
            {
                for (int ff = 0; ff < contours.size(); ff++)
                {
                    double dist2 = PointFacetDistance::distanceToEdgeSquared(
                        mPoint, contours[ff], Kernel());
                    if (dist2 < distToFacetSquared)
                        distToFacetSquared = dist2;
                }
            }
            else
            {
                distToFacetSquared = PointFacetDistance::distanceToPlaneSquared(
                    mPoint, h->plane(), Kernel());
            }
                            
            if (distToFacetSquared < mMinDistanceSquared)
            {
                mMinDistanceSquared = distToFacetSquared;
                mClosestFacet = h;
            }
            
//            for (int cc = 0; cc < contours.size(); cc++)
//            {
//                std::cout << "contour " << cc << ": ";
//                for (int vv = 0; vv < contours[cc].size(); vv++)
//                {
//                    std::cout << contours[cc][vv] << " ";
//                }
//                std::cout << "\n";
//            }
        }
        
        void visit(typename NefPolyhedron::Vertex_const_handle h) {}
        void visit(typename NefPolyhedron::Halfedge_const_handle h) {}
        void visit(typename NefPolyhedron::SHalfedge_const_handle h) {}
        void visit(typename NefPolyhedron::SHalfloop_const_handle h) {}
        void visit(typename NefPolyhedron::SFace_const_handle) {}
        
    private:
        typename NefPolyhedron::Point_3 mPoint;
        double mMinDistanceSquared;
        typename NefPolyhedron::Halffacet_const_handle mClosestFacet;
    };
    
    
    template<class NefPolyhedron>
    class ShellBoundsCalculator
    {
    public:
        ShellBoundsCalculator() :
            mBoundingBox(emptyBbox())
        {
        }
        
        CGAL::Bbox_3 boundingBox() const { return mBoundingBox; }
        typename NefPolyhedron::Point_3 anyPoint() const { return mAnyPoint; }
    
        void visit(typename NefPolyhedron::Vertex_const_handle h)
        {
            CGAL::Bbox_3 vertBounds(h->point().bbox());
            mBoundingBox = mBoundingBox + vertBounds;
            mAnyPoint = h->point();
            
//            std::cout << "\t\t\tvisiting " << h->point() << "\n";
        }
        
        void visit(typename NefPolyhedron::Halffacet_const_handle h) {}
        void visit(typename NefPolyhedron::Halfedge_const_handle h) {}
        void visit(typename NefPolyhedron::SHalfedge_const_handle h) {}
        void visit(typename NefPolyhedron::SHalfloop_const_handle h) {}
        void visit(typename NefPolyhedron::SFace_const_handle) {}
        
    private:
        CGAL::Bbox_3 mBoundingBox;
        typename NefPolyhedron::Point_3 mAnyPoint;
    };
    
    template<class NefPolyhedron>
    class NestedVolume
    {
        typedef typename NefPolyhedron::Kernel Kernel;
    public:
        NestedVolume(const NefPolyhedron & poly,
            typename NefPolyhedron::Volume_const_handle vol) :
            mPolyhedron(&poly),
            mNefVolume(vol),
            mBoundingBox(emptyBbox()),
            mParentVolume(0L)
        {
            // 1. Find outermost shell and its signed volume
            // 
            // old approach: assume first shell is outermost.  (Wrong!  Sigh.)
            
//            ShellBoundsCalculator<NefPolyhedron> visitor;
//            poly.visit_shell_objects(
//                typename NefPolyhedron::SFace_const_handle(vol->shells_begin()),
//                visitor);
//            mBoundingBox = visitor.boundingBox();
//            mAnyPoint = visitor.anyPoint();
            
            // Debug step: iterate over all shells...
            typename NefPolyhedron::Shell_entry_const_iterator itr;
            for (itr = vol->shells_begin(); itr != vol->shells_end(); itr++)
            {
                ShellBoundsCalculator<NefPolyhedron> vis;
                poly.visit_shell_objects(itr, vis);
                
//                std::cout << "\tshell bounds " << vis.boundingBox() << "\n";                
                
                if (Local::encloses(vis.boundingBox(), mBoundingBox))
                {
                    mOuterShell = itr;
//                    std::cout << "\t(this is now outer shell)\n";
                }
                

                mBoundingBox = mBoundingBox + vis.boundingBox();
                mAnyPoint = vis.anyPoint();
            }
            
//            std::cout << "* Bounding box is " << mBoundingBox << "\n";
        }
        
        NestedVolume() : mPolyhedron(0L)
        {}
        
        typename NefPolyhedron::Volume_const_handle nefVolume() const
        {
            return mNefVolume;
        }
        
        typename NefPolyhedron::Point_3 anyPoint() const
        {
            return mAnyPoint;
        }
        
        // assume that we've discarded the outermost volume already.
        bool isInterior() const
        {
            if (mParentVolume == 0L)
                return true;
            else
                return !mParentVolume->isInterior();
        }
        
        CGAL::Bbox_3 boundingBox() const { return mBoundingBox; }
        
        const NestedVolume<NefPolyhedron>* parent() const
        {
            return mParentVolume;
        }
        void parent(NestedVolume<NefPolyhedron>* p) { mParentVolume = p; }
        
        bool encloses(const NestedVolume<NefPolyhedron> & p) const
        {
            // little bit of fast culling.
            if (!Local::interiorEncloses(boundingBox(), p.boundingBox()))
                return false;
            
            return encloses(p.anyPoint());
        }
                
        bool encloses(const typename NefPolyhedron::Point_3 & p) const
        {
            // 1. Find halffacet closest to p
            // 2. Check which side the point is on
            
            assert(mPolyhedron != 0L);
            
            double minDistanceSquared = std::numeric_limits<double>::max();
            
            // Get outermost shell somehow!
            // cannot assume it's the first shell (why not??)
//            typename NefPolyhedron::Shell_entry_const_iterator outerShell =
//                nefVolume()->shells_begin();
            typename NefPolyhedron::Shell_entry_const_iterator outerShell =
                mOuterShell;
            
            ClosestFacetFinder<NefPolyhedron> finder(p);
            mPolyhedron->visit_shell_objects(
                typename NefPolyhedron::SFace_const_handle(
                    outerShell), finder);
            
            // Now check which side the point is on.
            // CAREFUL: The outer shell of a volume faces INTO the volume!
            // This behavior is counterintuitive to me... the point is that
            // if the point is "below" the closest facet, it is OUTSIDE the
            // volume.
            if (finder.closestFacet()->plane().has_on_positive_side(p))
                return true;
            else
                return false;
        }
    private:
        const NefPolyhedron* mPolyhedron;
        typename NefPolyhedron::Volume_const_handle mNefVolume;
        typename NefPolyhedron::Shell_entry_const_iterator mOuterShell;
        typename NefPolyhedron::Point_3 mAnyPoint;
        CGAL::Bbox_3 mBoundingBox;
        const NestedVolume<NefPolyhedron>* mParentVolume;
    };
}

template<class NefPolyhedron>
std::vector<typename NefPolyhedron::Volume_const_handle>
interiorVolumes(const NefPolyhedron & poly)
{
    typename NefPolyhedron::Volume_const_iterator itr = poly.volumes_begin();
    itr++; // SKIP outside volume
    
    std::vector<Local::NestedVolume<NefPolyhedron> > volumes;
    
    for (; itr != poly.volumes_end(); itr++)
    {
         volumes.push_back(Local::NestedVolume<NefPolyhedron>(poly,
            typename NefPolyhedron::Volume_const_handle(itr)));
    }
    
    for (int vv = 0; vv < volumes.size(); vv++)
    {
//        std::cout << "Volume " << vv << " is bounded by "
//            << volumes[vv].boundingBox() << ".\n";
    }
    
    for (int vv = 0; vv < volumes.size(); vv++)
    for (int ww = 0; ww < volumes.size(); ww++)
    if (ww != vv)
    if (volumes[ww].encloses(volumes[vv]))
    {
//        std::cout << "Volume  " << volumes[ww].boundingBox()
//            << " encloses " << volumes[vv].boundingBox() << "\n";
         if (volumes[vv].parent() == 0L)
             volumes[vv].parent(&volumes[ww]);
         else if (volumes[vv].parent()->encloses(volumes[ww])) 
             volumes[vv].parent(&volumes[ww]);
    }
    else
    {
//        std::cout << "Volume " << volumes[ww].boundingBox()
//            << " does NOT enclose " << volumes[vv].boundingBox() << "\n";
//        std::cout << "\tNon-enclosed point " << volumes[vv].anyPoint() << "\n";
//        volumes[ww].encloses(volumes[vv]);
    }
    
    std::vector<typename NefPolyhedron::Volume_const_handle> out;
    for (int vv = 0; vv < volumes.size(); vv++)
    if (volumes[vv].isInterior())
    {
        out.push_back(volumes[vv].nefVolume());
//        std::cout << "Include volume bounded by "
//            << volumes[vv].boundingBox() << "\n";
    }
    else
    {
//        std::cout << "Exclude volume bounded by "
//            << volumes[vv].boundingBox() << "\n";
    }
    return out;
}

}





#endif
