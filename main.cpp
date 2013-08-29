#include <iostream>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/OFF_to_nef_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> NefPolyhedron;
typedef NefPolyhedron::Point_3 Point_3;
typedef NefPolyhedron::Traits Traits;
typedef CGAL::Polyhedron_3<Traits> Polyhedron;

#include "ShellVisitor.h"
#include "InteriorVolumes-inl.h"

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


using namespace std;

void printHelp();
void doIntersection();
void doUnion();
void doDifference();
bool testIntersection();
void printError();

NefPolyhedron readMultiOFF();
Polyhedron readPolyhedron();
void writeNefPolyhedron(const NefPolyhedron & poly);

int main(int argc, char const** argv)
{
    for (int nn = 0; nn < argc; nn++)
    {
        cerr << "Arg " << nn << " = " << argv[nn] << "\n";
    }
    
    if (argc == 1)
        printHelp();
    else
    {
        string token(argv[1]);
        
        if (token == "intersection")
            doIntersection();
        else if (token == "union")
            doUnion();
        else if (token == "difference")
            doDifference();
        else if (token == "testIntersection")
            return testIntersection();
        else
        {
            printError();
            return 1;
        }
    }
    
    
    return 0;
}


void printHelp()
{
    cout << "Usage:\n";
    cout << "\tNefLab intersection\n";
    cout << "\tNefLab union\n";
    cout << "\tNefLab difference\n";
    cout << "\tNefLab testIntersection\n";
}

void doIntersection()
{
    cerr << "intersection\n";
    
//    Polyhedron p1, p2;
//    cerr << "Reading first polyhedron.\n";
//    p1 = readPolyhedron();
//    cerr << "Reading second polyhedron.\n";
//    p2 = readPolyhedron();
//    NefPolyhedron n1(p1);
//    NefPolyhedron n2(p2);

    NefPolyhedron n1, n2;
    n1 = readMultiOFF();
    n2 = readMultiOFF();
    
    cerr << "Poly1:\n" << n1 << "\n";
    cerr << "Poly2:\n" << n2 << "\n";
    
    NefPolyhedron nefIntersection = (n1*n2).regularization();
    writeNefPolyhedron(nefIntersection);
    
    cerr << nefIntersection;
}

void doUnion()
{
    cerr << "union\n";
    
//    Polyhedron p1, p2;
//    cerr << "Reading first polyhedron.\n";
//    p1 = readPolyhedron();
//    cerr << "Reading second polyhedron.\n";
//    p2 = readPolyhedron();
//    NefPolyhedron n1(p1);
//    NefPolyhedron n2(p2);

    NefPolyhedron n1, n2;
    n1 = readMultiOFF();
    n2 = readMultiOFF();
    
    NefPolyhedron nefUnion = (n1+n2).regularization();
    writeNefPolyhedron(nefUnion);
}

void doDifference()
{
    cerr << "difference\n";
    
//    Polyhedron p1, p2;
//    cerr << "Reading first polyhedron.\n";
//    p1 = readPolyhedron();
//    cerr << "Reading second polyhedron.\n";
//    p2 = readPolyhedron();
//    NefPolyhedron n1(p1);
//    NefPolyhedron n2(p2);

    NefPolyhedron n1, n2;
    n1 = readMultiOFF();
    n2 = readMultiOFF();
    
    NefPolyhedron nefDifference = (n1 - n2).regularization();
    writeNefPolyhedron(nefDifference);
}

bool testIntersection()
{
    cerr << "testIntersection\n";
    
//    Polyhedron p1, p2;
//    cerr << "Reading first polyhedron.\n";
//    p1 = readPolyhedron();
//    cerr << "Reading second polyhedron.\n";
//    p2 = readPolyhedron();
//    NefPolyhedron n1(p1);
//    NefPolyhedron n2(p2);

    NefPolyhedron n1, n2;
    n1 = readMultiOFF();
    n2 = readMultiOFF();
    
    cerr << "Doing nothing!\n";
    return 0;
}

void printError()
{
    cout << "I don't understand what you want.\n";
    printHelp();
}

NefPolyhedron readMultiOFF()
{
    int numPositiveShells, numNegativeShells;
    int numIgnoredFacets;
    
    cin >> numPositiveShells >> numNegativeShells;
    
    NefPolyhedron nef;
    
    for (int nn = 0; nn < numPositiveShells; nn++)
    {
        NefPolyhedron addend;
        
        numIgnoredFacets = CGAL::OFF_to_nef_3(cin, addend);
        if (numIgnoredFacets > 0)
            cerr << "Ignored " << numIgnoredFacets << " facets!\n";
        nef = nef + addend;
    }
    
    for (int nn = 0; nn < numNegativeShells; nn++)
    {
        NefPolyhedron subtrahend;
        
        numIgnoredFacets = CGAL::OFF_to_nef_3(cin, subtrahend);
        if (numIgnoredFacets > 0)
            cerr << "Ignored " << numIgnoredFacets << " facets!\n";
        nef = nef - subtrahend;
    }
    
//    NefPolyhedron nef;
//    int numIgnoredFacets;
//    
//    numIgnoredFacets = CGAL::OFF_to_nef_3(cin, nef);
//    
//    if (numIgnoredFacets > 0)
//        cerr << "Ignoring " << numIgnoredFacets << " facets!!\n";
    
    return nef;
}

Polyhedron readPolyhedron()
{
    int numVertices, numFaces;
    cin >> numVertices >> numFaces;
    cerr << "Reading " << numVertices << " vertices and " << numFaces << " faces.\n";
    
    vector<Point_3> vertices(numVertices);
    vector<vector<unsigned int> > faces(numFaces);
    
    for (int vv = 0; vv < numVertices; vv++)
    {
        double x, y, z;
        cin >> x >> y >> z;
        vertices[vv] = Point_3(x, y, z);
    }
    
    for (int ff = 0; ff < numFaces; ff++)
    {
        unsigned int v1, v2, v3;
        cin >> v1 >> v2 >> v3;
        vector<unsigned int> tri(3);
        tri[0] = v1;
        tri[1] = v2;
        tri[2] = v3;
        
        faces[ff] = tri;
    }
    
    Polyhedron poly;
    BuildMesh<Polyhedron::HalfedgeDS> builder(vertices, faces);
    poly.delegate(builder);
    
    return poly;
}

void writeNefPolyhedron(const NefPolyhedron & poly)
{
    typedef NefPolyhedron::Halffacet_cycle_const_iterator
        HalffacetCycleConstIterator;
    typedef NefPolyhedron::Volume_const_handle VolumeConstHandle;
    typedef NefPolyhedron::Shell_entry_const_iterator ShellEntryConstIterator;
    typedef NefPolyhedron::SFace_const_handle SFaceConstHandle;
    
    ShellVisitor shellVisitor;
    
    vector<VolumeConstHandle> insideVolumes =
        InteriorVolumes::interiorVolumes(poly);
    
    for (int loopVertex = 0; loopVertex < insideVolumes.size(); loopVertex++)
    {
        ShellEntryConstIterator itr;
        for (itr = insideVolumes[loopVertex]->shells_begin();
            itr != insideVolumes[loopVertex]->shells_end(); itr++)
        {
            poly.visit_shell_objects(SFaceConstHandle(itr),
                shellVisitor);
        }
    }
    
    map<Vector3d, unsigned int> vertexIds;
    int vertId = 0;
    
    const std::vector<std::vector<int> > & facets =
        shellVisitor.contourSizes();
    const std::vector<std::vector<Vector3d> > & contourVertices =
        shellVisitor.contourVertices();
    const std::set<Vector3d> & vertexSet = shellVisitor.vertexSet();
    
    cout << "numVertices " << vertexSet.size() << "\n";
    
    set<Vector3d>::const_iterator itr;
    for (itr = vertexSet.begin(); itr != vertexSet.end(); itr++)
    {
        vertexIds[*itr] = vertId++;
        Vector3d v = *itr;
        cout << vertId-1 << " " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }
    
    cout << "numFacets " << facets.size() << "\n";
    
    for (int cc = 0; cc < facets.size(); cc++)
    {
        cout << "\tnumContours " << facets[cc].size() << "\n";
        
        int vertexNum = 0;
        
        // Iterate over facet loops
        for (int ll = 0; ll < facets[cc].size(); ll++)
        {
            cout << "\tnumVertices " << facets[cc][ll] << "\n";
            cout << "\t";
            
            for (int loopVertex = 0; loopVertex < facets[cc][ll]; loopVertex++)
            {
                cout << vertexIds[contourVertices[cc][vertexNum++]] << " ";
            }
            cout << "\n";
        }
    }
    
    
    
    cerr << "I have now written a NefPolyhedron.\n";
}
