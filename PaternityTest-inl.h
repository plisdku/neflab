//
//  PaternityTest.h
//  NefLab
//
//  Created by Paul Hansen on 4/9/16.
//
//

#ifndef PaternityTest_h
#define PaternityTest_h

#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
//#include <CGAL/squared_distance_3.h>

#include "PointFacetDistance-inl.h"

typedef std::vector<std::vector<unsigned int> >::const_iterator FaceIndex;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 2, FaceIndex> Box;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 3, FaceIndex> Box3;

std::vector<std::vector<unsigned int> > facetInheritance(
    const std::vector<Point_3> & v1,
    const std::vector<std::vector<unsigned int> > & f1,
    const std::vector<Point_3> & v2,
    const std::vector<std::vector<unsigned int> > & f2);

int paternityTest(const std::vector<Point_3> & vOld,
    const std::vector<unsigned int> & fOld,
    const std::vector<Point_3> & vNew,
    const std::vector<unsigned int> & fNew);


// this struct SHOULD be inside inheritFormulasFast, but I can't use a local
// type as a template parameter.  Lame.
struct CallbackTri3
{
    const std::vector<Point_3> & vOld_;
    const std::vector<std::vector<unsigned int> > & fOld_;
    const std::vector<Point_3> & vNew_;
    const std::vector<std::vector<unsigned int> > & fNew_;
    std::vector<std::vector<unsigned int> > & paternity_;
    int & numInherited_;
    
    CallbackTri3(const std::vector<Point_3> & vOld,
        const std::vector<std::vector<unsigned int> > & fOld,
        const std::vector<Point_3> & vNew,
        const std::vector<std::vector<unsigned int> > & fNew,
        std::vector<std::vector<unsigned int> > & paternity,
        int & outNumInheritedFormulae) :
        vOld_(vOld),
        fOld_(fOld),
        vNew_(vNew),
        fNew_(fNew),
        paternity_(paternity),
        numInherited_(outNumInheritedFormulae)
    {}
    
    void operator()( const Box3 & newBox, const Box3 & oldBox )
    {
        int iNew = newBox.handle() - fNew_.begin();
        int iOld = oldBox.handle() - fOld_.begin();
        
//        std::cerr << "(" << iNew << " <= " << iOld << ")\n";
//        std::cerr << "\t" << newBox.bbox() << " and " << oldBox.bbox() << "\n";
        
        int doIntersect = paternityTest(
            vOld_,
            fOld_[iOld],
            vNew_,
            fNew_[iNew]);
        
        numInherited_ += doIntersect;
        
        if (doIntersect)
            paternity_.at(iNew).push_back(iOld);
    }
};


/**
 * Return a weighted average of p0, p1, and p2.
 */
static inline Kernel::Point_2 sInset(const Kernel::Point_2 & p0,
    const Kernel::Point_2 & p1,
    const Kernel::Point_2 & p2)
{
    const Kernel::FT a = 1.0 - 1e-7;
    const Kernel::FT b = 0.5*(1.0 - a);
    
    return Kernel::Point_2(a*p0[0] + b*p1[0] + b*p2[0],
        a*p0[1] + b*p1[1] + b*p2[1]);
}

/**
 * Return a slightly inset version of the triangle.
 */
static inline CGAL::Triangle_2<Kernel> sInsetTri(const std::vector<Kernel::Point_2> & pts)
{
    return CGAL::Triangle_2<Kernel>(
        sInset(pts[0], pts[1], pts[2]),
        sInset(pts[1], pts[2], pts[0]),
        sInset(pts[2], pts[0], pts[1]));
}


int paternityTest(const std::vector<Point_3> & vOld,
    const std::vector<unsigned int> & fOld,
    const std::vector<Point_3> & vNew,
    const std::vector<unsigned int> & fNew)
{
    int numInherited = 0;
    
    if (fOld.size() != 3 || fNew.size() != 3)
        std::cerr << "Warning: not operating on triangles\n";
    
//    std::cerr << vOld.size() << " verts and " << vNew.size() << " verts:\n";
//    std::cerr << "\told tri " << fOld.size() << " new tri " << fNew.size() << "\n";
    
//    Vector_3 nOld = CGAL::unit_normal(vOld[fOld[0]], vOld[fOld[1]], vOld[fOld[2]]);
//    Vector_3 nNew = CGAL::unit_normal(vNew[fNew[0]], vNew[fNew[1]], vNew[fNew[2]]);
    Vector_3 nOld = CGAL::normal(vOld[fOld[0]], vOld[fOld[1]], vOld[fOld[2]]);
    Vector_3 nNew = CGAL::normal(vNew[fNew[0]], vNew[fNew[1]], vNew[fNew[2]]);
    
    Kernel::FT dotProd = nOld * nNew;
    Kernel::FT dotProdSquared = dotProd * dotProd;
    Kernel::FT maxDotSquared = nOld.squared_length() * nNew.squared_length();
    
    
    // TEST 1: are the facets coplanar?  Let's just do a numerical test.
    // CGAL is doing exact arithmetic on the inside but it began with double-
    // precision coordinates from the input file, so I can't ask it if the
    // two facets are EXACTLY coplanar.  Gotta do all this silly stuff.
//    if ( dotProd > 0.999999 || dotProd < -0.999999)
    if (dotProdSquared > 0.999999*maxDotSquared)
    {
//        std::cerr << "Facets appear to be parallel\n";
//        std::cerr << "\tNew normal " << nNew << " and old " << nOld << " dot2 " << dotProdSquared << " max2 " << maxDotSquared << "\n";
        
        Kernel::Plane_3 oldPlane(vOld[fOld[0]], vOld[fOld[1]], vOld[fOld[2]]);
        
//        double distSquared = PointFacetDistance::distanceToPlaneSquared(
//            vNew[fNew[0]], oldPlane, Kernel());
        Kernel::FT distSquared = CGAL::squared_distance(vNew[fNew[0]], oldPlane);
//        std::cerr << "\tdist " << distSquared << "\n";
        
        if (distSquared < 1e-9)
        {
//            std::cerr << "Facets appear to be coplanar\n";
            
            // Project both faces—presumed triangles—to the old face's
            // supporting plane.
            std::vector<Kernel::Point_2> oldFace(3), newFace(3);
            for (int ii = 0; ii < 3; ii++)
            {
                oldFace[ii] = oldPlane.to_2d(vOld[fOld[ii]]);
                newFace[ii] = oldPlane.to_2d(vNew[fNew[ii]]);
            }
            
            // Shrink the new face by a tiny fraction to exclude very very
            // small intersections.  This is NOT a CGAL-style way to do this!
            // But in practice my input geometry frequently has triangles that
            // should only share an edge being marked as intersecting.  I
            // think CGAL's do_intersect looks for intersection in the interior,
            // not on the edge, so it's an issue of having floating-point
            // coordinates passed in I guess.
            
            CGAL::Triangle_2<Kernel> tri1 = sInsetTri(oldFace);
            CGAL::Triangle_2<Kernel> tri2(newFace[0], newFace[1], newFace[2]);
            
//            CGAL::Triangle_2<Kernel> tri1(oldFace[0], oldFace[1], oldFace[2]);
//            CGAL::Triangle_2<Kernel> tri2(newFace[0], newFace[1], newFace[2]);
            
            if (CGAL::do_intersect(tri1, tri2))
            {
//                std::cerr << "Facets appear to intersect\n";
//                std::cerr << "*hit*\n";
                numInherited++;
            }
        }
        
    }
    
    /*
    // Check if the triangles are coplanar.  The ancestor triangle may face the
    // opposite direction from the descendent triangle.
    if (cgalTrianglePlanes()[newTri] == ancestor.cgalTrianglePlanes()[oldTri] ||
        cgalTrianglePlanes()[newTri] ==
            ancestor.cgalTrianglePlanes()[oldTri].opposite())
    {
        Vector3d normal = toVector3d(
            cgalTrianglePlanes()[newTri].orthogonal_vector());
        
        int projectionDir = 0;
        if (fabs(normal[1]) > fabs(normal[0]))
            projectionDir = 1;
        if (fabs(normal[2]) > fabs(normal[projectionDir]))
            projectionDir = 2;
        
        int u2 = projectionDir;
        int u0 = (u2+1)%3;
        int u1 = (u2+2)%3;
        
        Triangle2d tNew = sProject(triangles()[newTri].triangle(), u0, u1);
        Triangle2d tOld = sProject(ancestor.triangles()[oldTri].triangle(), u0, u1);
        Vector2d centerNew = (tNew[0] + tNew[1] + tNew[2])/3;
        
        // First inherit the control vertices that define the plane.  Pick a
        // triangle that does the job.  (It's possibly that several triangles
        // may all be ok here.)
        if (tOld.encloses(centerNew))
        {
            numInherited++;
            mTriangles[newTri].controlVertices(ancestor.triangles()[oldTri].
                controlVertices());
        }
        
        // Next inherit the control vertices attached to each edge of each
        // triangle.  Edges that are collinear with edges of the input
        // polyhedra will retain the control vertices of the input edges.
        // However, new edges arising from Boolean operations on polyhedra
        // will not be attached to any control vertices directly.
        //
        // The collinear edge test checks whether endpoints lie within some
        // distance of a given line.
        for (int oldEdge = 0; oldEdge < 3; oldEdge++)
        for (int newEdge = 0; newEdge < 3; newEdge++)
        if (sCollinearEdges(triangles()[newTri].triangle(), newEdge,
            ancestor.triangles()[oldTri].triangle(), oldEdge, 1e-6))
        {
            mTriangles[newTri].edgeControlVertices(newEdge,
                ancestor.triangles()[oldTri].edgeControlVertices(oldEdge));
        }
    }
    */
    
    return numInherited;
}

/**
 * Calculate bounding box of a facet.
 */
CGAL::Bbox_3 faceBounds(const std::vector<Point_3> & vertices,
    const std::vector<unsigned int> & face)
{
    std::vector<Point_3> facePoints(face.size());
    for (int ii = 0; ii < face.size(); ii++)
    {
        facePoints[ii] = vertices[face[ii]];
    }
    
    CGAL::Bbox_3 box = CGAL::bbox_3(facePoints.begin(), facePoints.end());
//    std::cerr << box << "\n";
    return box;
}

std::vector<std::vector<unsigned int> > facetInheritance(
    const std::vector<Point_3> & vertices1,
    const std::vector<std::vector<unsigned int> > & faces1,
    const std::vector<Point_3> & vertices2,
    const std::vector<std::vector<unsigned int> > & faces2)
{
    std::vector<Box3> ancestorBoxes(faces1.size());
    std::vector<Box3> newBoxes(faces2.size());
    
    
    for (int nn = 0; nn < faces1.size(); nn++)
    {
        ancestorBoxes[nn] = Box3(faceBounds(vertices1, faces1[nn]), faces1.begin() + nn);
//        std::cerr << "Old " << nn << ": " << ancestorBoxes[nn].bbox() << "\n";
    }
    
    for (int nn = 0; nn < faces2.size(); nn++)
    {
        newBoxes[nn] = Box3(faceBounds(vertices2, faces2[nn]), faces2.begin() + nn);
//        std::cerr << "New " << nn << ": " << newBoxes[nn].bbox() << "\n";
    }
    
    std::vector<std::vector<unsigned int> > paternity(faces2.size());
    
    int outNumInherited = 0;
    CallbackTri3 callback(vertices1, faces1, vertices2, faces2, paternity,
        outNumInherited);
    
    CGAL::box_intersection_d(newBoxes.begin(), newBoxes.end(),
        ancestorBoxes.begin(), ancestorBoxes.end(), callback);
//    cout << myBoxes.size() << " triangles inherited " << outNumInherited
//        << " formulas from " << ancestorBoxes.size() << " triangles.\n";
    
    return paternity;
}




#endif /* PaternityTest_h */
