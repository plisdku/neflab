//
//  BuildMesh-inl.h
//  NefLab
//
//  Created by Paul Hansen on 4/9/16.
//
//

#ifndef BuildMesh_inl_h
#define BuildMesh_inl_h

template<class HalfedgeDS>
class BuildMesh: public CGAL::Modifier_base<HalfedgeDS>
{
    typedef typename HalfedgeDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
public:
    BuildMesh(const std::vector<Point> & vertices,
        const std::vector<std::vector<unsigned int> > & faces) :
        mVertices(vertices),
        mFaces(faces)
    {
//        std::cerr << "The faces?\n";
//        for (int ff = 0; ff < faces.size(); ff++)
        {
//            std::cout << "face " << ff << ": ";
//            for (int nn = 0; nn < faces.at(ff).size(); nn++)
//                std::cout << faces.at(ff).at(nn) << " ";
//            std::cout << "\n";
        }
    }
    
    void operator() (HalfedgeDS & halfedges)
    {
        bool verbose = false;
        CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> builder(halfedges,
            verbose);
        
        builder.begin_surface(mVertices.size(), mFaces.size());
        for (int vv = 0; vv < mVertices.size(); vv++)
        {
            builder.add_vertex(mVertices.at(vv));
//            std::cout << "Vertex " << vv << ": " << mVertices.at(vv) << "\n";
        }
        for (int ff = 0; ff < mFaces.size(); ff++)
        {
//            std::cout << "Face " << ff << ": ";
            builder.begin_facet();
            for (int vv = 0; vv < mFaces.at(ff).size(); vv++)
            {
                builder.add_vertex_to_facet(mFaces.at(ff).at(vv));
//                std::cout << mFaces.at(ff).at(vv) << " ";
            }
//            std::cout << "\n";
            builder.end_facet();
        }
        builder.end_surface();
    }
private:
    const std::vector<Point> & mVertices;
    const std::vector<std::vector<unsigned int> > & mFaces;
};


#endif /* BuildMesh_inl_h */
