/*
 *  ShellTriangulator.h
 *  SnapdragonTopLevel
 *
 *  Created by Paul C Hansen on 6/24/11.
 *  Copyright 2011 Stanford University. All rights reserved.
 *
 */

#ifndef SHELLVISITOR_H
#define SHELLVISITOR_H

//#include "TrackedPolyhedron.h"
//#include "Triangulate/Triangulator.h"
#include "utility/geometry.h"
#include "CGALUtilities.h"
#include <limits>
#include <set>

class ShellVisitor
{
public:
//    typedef SolidModel::TrackedPolyhedron::NefPolyhedron NefPolyhedron;
    typedef NefPolyhedron::Halffacet_cycle_const_iterator
        HalffacetCycleConstIterator;
    typedef NefPolyhedron::Nef_polyhedron_S2 NefPolyhedronS2;
    typedef NefPolyhedronS2::SHalfedge_const_handle SHalfedgeConstHandle;
    
    ShellVisitor()
    {
    }
    
    /*ShellTriangulator(std::vector<SimpleMesh::Triangle> & outTriangles,
        std::vector<NefPolyhedron::Plane_3> & outTrianglePlanes) :
        mTriangleVertices(3*MAXPOLYTRIANGLES),
        mTriangulator(),
        mTriangles(outTriangles),
        mTrianglePlanes(outTrianglePlanes)
    {
    }*/
    
    void visit(NefPolyhedron::Halffacet_const_handle facetItr)
    {
        std::vector<Rect3d> boundingBoxes;
        std::vector<std::vector<Vector3d> > contours;
        NefPolyhedron::Plane_3 plane(facetItr->plane().opposite());
        // I flip the plane because shells point into volumes, not out of them.
        // This seems ridiculous.
        
        HalffacetCycleConstIterator itr;
        for (itr = facetItr->facet_cycles_begin();
            itr != facetItr->facet_cycles_end();
            itr++)
        {
            if (itr.is_shalfedge())
            {
                SHalfedgeConstHandle handle(itr);
                SHalfedgeConstHandle edgeHandle = handle;
                
                std::vector<Vector3d> contour;
                Rect3d boundingBox(std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max());
                
//                cout << "contour\n";
                do {
                    edgeHandle = edgeHandle->prev(); // to make it face outwards
                    NefPolyhedron::Point_3 p3(
                        edgeHandle->source()->center_vertex()->point());
                    Vector3d vert(toVector3d(p3));
//                    cout << "\tcgal " << p3 << ", me " << vert << "\n";
                    boundingBox.p1 = vec_min(boundingBox.p1, vert);
                    boundingBox.p2 = vec_max(boundingBox.p2, vert);
                    contour.push_back(vert);
                } while (edgeHandle != handle);
                
                boundingBoxes.push_back(boundingBox);
                contours.push_back(contour);
            }
            else if (itr.is_shalfloop())
            {
                //cerr << "Found isolated vertex; ignoring it!\n";
            }
            else
            {
                //cerr << "\tFound cycle, it's a mys it's a mys it's a mystery\n";
            }
        }
        
        // Find the largest contour.  Its bounding box is biggest...
        int biggestBoxIndex = 0;
        for (int bb = 1; bb < boundingBoxes.size(); bb++)
        {
            if (sumSquares(boundingBoxes[bb].size()) >
                sumSquares(boundingBoxes[biggestBoxIndex].size()))
            {
                biggestBoxIndex = bb;
            }
        }
        swap(contours[0], contours[biggestBoxIndex]);
        
        // Now the contours are ordered with the outside box first, and all
        // the "hole" contours are after that.
        
        // The vertices are in contours
        // The normal vector is normalVector
        
        std::vector<int> contourSizes;
        std::vector<Vector3d> vertices;
        for (int cc = 0; cc < contours.size(); cc++)
        {
            contourSizes.push_back(contours[cc].size());
            copy(contours[cc].begin(), contours[cc].end(),
                back_inserter(vertices));
        }
        
        mContourSizes.push_back(contourSizes);
        mContourVertices.push_back(vertices);
        
        std::set<Vector3d> facetVertexSet;
        for (int vv = 0; vv < vertices.size(); vv++)
        {
//            std::cerr << "\t\tvert " << vertices[vv] << "\n";
            mVertexSet.insert(vertices[vv]);
            facetVertexSet.insert(vertices[vv]);
        }
        if (facetVertexSet.size() != vertices.size())
            std::cerr << "Warning: got a repeated vertex!\n";
        
        /*
        int numTriangles = mTriangulator.triangulate(contourSizes, vertices,
            toVector3d(plane.orthogonal_vector()), mTriangleVertices);
        
        for (int tt = 0; tt < numTriangles; tt++)
        {
            mTriangles.push_back(SimpleMesh::Triangle(
                Triangle3d(*mTriangleVertices[3*tt],
                    *mTriangleVertices[3*tt+1],
                    *mTriangleVertices[3*tt+2])));
            mTrianglePlanes.push_back(plane);
            
//            sCheckTriangleNormal(mTriangles[mTriangles.size()-1].triangle(),
//                toVector3d(plane.orthogonal_vector()));
        }
        */
    }
    
    void visit(NefPolyhedron::Vertex_const_handle h) {}
    void visit(NefPolyhedron::Halfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfloop_const_handle h) {}
    void visit(NefPolyhedron::SFace_const_handle) {}
    
    const std::vector<std::vector<int> > & contourSizes() const
    {
        return mContourSizes;
    }
    
    const std::vector<std::vector<Vector3d> > & contourVertices() const
    {
        return mContourVertices;
    }
    
    const std::set<Vector3d> & vertexSet() const
    {
        return mVertexSet;
    }
    
private:
//    enum { MAXPOLYTRIANGLES = 10000 };
//    std::vector<Vector3d*> mTriangleVertices;
//    Triangulate::Triangulator mTriangulator;
//    std::vector<SimpleMesh::Triangle> & mTriangles;
//    std::vector<NefPolyhedron::Plane_3> & mTrianglePlanes;
    std::vector<std::vector<int> > mContourSizes;
    std::vector<std::vector<Vector3d> > mContourVertices;
    std::set<Vector3d> mVertexSet;
};


#endif

