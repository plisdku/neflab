//
//  FacetIntersectionTest-inl.h
//  NefLab
//
//  Created by Paul Hansen on 4/24/17.
//
//

#ifndef FacetIntersectionTest_inl_h
#define FacetIntersectionTest_inl_h

#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>

#include "PointFacetDistance-inl.h"

namespace FacetIntersectionTest
{

typedef std::vector<std::vector<unsigned int> >::const_iterator FaceIndex;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 2, FaceIndex> Box;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 3, FaceIndex> Box3;

int numIntersections(
    const std::vector<Point_3> & v1,
    const std::vector<std::vector<unsigned int> > & f1,
    const std::vector<Point_3> & v2,
    const std::vector<std::vector<unsigned int> > & f2);

int triIntersectionTest(const std::vector<Point_3> & vOld,
    const std::vector<unsigned int> & fOld,
    const std::vector<Point_3> & vNew,
    const std::vector<unsigned int> & fNew);


// this struct SHOULD be inside inheritFormulasFast, but I can't use a local
// type as a template parameter.  Lame.
struct CallbackTri3
{
    const std::vector<Point_3> & v1_;
    const std::vector<std::vector<unsigned int> > & f1_;
    const std::vector<Point_3> & v2_;
    const std::vector<std::vector<unsigned int> > & f2_;
    int & numIntersections_;
//    int & numInherited_;
    
    CallbackTri3(const std::vector<Point_3> & v1,
        const std::vector<std::vector<unsigned int> > & f1,
        const std::vector<Point_3> & v2,
        const std::vector<std::vector<unsigned int> > & f2,
        int & numIntersections) :
        v1_(v1),
        f1_(f1),
        v2_(v2),
        f2_(f2),
        numIntersections_(numIntersections)
    {}
    
    void operator()( const Box3 & b1, const Box3 & b2 )
    {
        int i1 = b1.handle() - f1_.begin();
        int i2 = b2.handle() - f2_.begin();
        
//        std::cerr << "(" << i1 << " <= " << i2 << ")\n";
//        std::cerr << "\t" << b1.bbox() << " and " << b2.bbox() << "\n";
        
        int doIntersect = triIntersectionTest(
            v1_,
            f1_[i1],
            v2_,
            f2_[i2]);
        
        numIntersections_ += doIntersect;
    }
};


///**
// * Return a weighted average of p0, p1, and p2.
// */
//static inline Kernel::Point_2 sInset(const Kernel::Point_2 & p0,
//    const Kernel::Point_2 & p1,
//    const Kernel::Point_2 & p2)
//{
//    const Kernel::FT a = 1.0 - 1e-7;
//    const Kernel::FT b = 0.5*(1.0 - a);
//    
//    return Kernel::Point_2(a*p0[0] + b*p1[0] + b*p2[0],
//        a*p0[1] + b*p1[1] + b*p2[1]);
//}

///**
// * Return a slightly inset version of the triangle.
// */
//static inline CGAL::Triangle_2<Kernel> sInsetTri(const std::vector<Kernel::Point_2> & pts)
//{
//    return CGAL::Triangle_2<Kernel>(
//        sInset(pts[0], pts[1], pts[2]),
//        sInset(pts[1], pts[2], pts[0]),
//        sInset(pts[2], pts[0], pts[1]));
//}


int triIntersectionTest(const std::vector<Point_3> & v1,
    const std::vector<unsigned int> & f1,
    const std::vector<Point_3> & v2,
    const std::vector<unsigned int> & f2)
{
    int numIntersections = 0;
    
    if (f1.size() != 3 || f2.size() != 3)
    {
        std::cerr << "Warning: not operating on triangles\n";
    }
    
    CGAL::Triangle_3<Kernel> tri1(v1[f1[0]], v1[f1[1]], v1[f1[2]]);
    CGAL::Triangle_3<Kernel> tri2(v2[f2[0]], v2[f2[1]], v2[f2[2]]);
    
//    std::cerr << "t1 " << tri1 << "\n";
//    std::cerr << "t2 " << tri2 << "\n";
    
    if (CGAL::do_intersect(tri1, tri2))
    {
        numIntersections++;
    }
    
    return numIntersections;
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

int numIntersections(
    const std::vector<Point_3> & vertices1,
    const std::vector<std::vector<unsigned int> > & faces1,
    const std::vector<Point_3> & vertices2,
    const std::vector<std::vector<unsigned int> > & faces2)
{
    std::vector<Box3> boxes1(faces1.size());
    std::vector<Box3> boxes2(faces2.size());
    
    for (int nn = 0; nn < faces1.size(); nn++)
    {
        boxes1[nn] = Box3(faceBounds(vertices1, faces1[nn]), faces1.begin() + nn);
//        std::cerr << "Old " << nn << ": " << boxes1[nn].bbox() << "\n";
    }
    
    for (int nn = 0; nn < faces2.size(); nn++)
    {
        boxes2[nn] = Box3(faceBounds(vertices2, faces2[nn]), faces2.begin() + nn);
//        std::cerr << "New " << nn << ": " << boxes2[nn].bbox() << "\n";
    }
    
    int outNumIntersections = 0;
    CallbackTri3 callback(vertices1, faces1, vertices2, faces2,
        outNumIntersections);
    
    CGAL::box_intersection_d(boxes1.begin(), boxes1.end(),
        boxes2.begin(), boxes2.end(), callback);
//    cout << myBoxes.size() << " triangles inherited " << outNumInherited
//        << " formulas from " << ancestorBoxes.size() << " triangles.\n";
    
    return outNumIntersections;
}


}; // namespace

#endif /* FacetIntersectionTest_inl_h */
