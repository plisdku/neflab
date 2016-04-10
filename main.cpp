#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <set>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/OFF_to_nef_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> NefPolyhedron;
typedef NefPolyhedron::Point_3 Point_3;
typedef NefPolyhedron::Vector_3 Vector_3;
typedef NefPolyhedron::Traits Traits;
typedef CGAL::Polyhedron_3<Traits> Polyhedron;

#include "BuildMesh-inl.h"
#include "ShellVisitor.h"
#include "InteriorVolumes-inl.h"
#include "Boolean-inl.h"
#include "PaternityTest-inl.h"

using namespace std;

void printHelp();
void printError();

/**
 * Read a single polyhedral shell from OFF format into a face-vertex structure.
 */
void readPolyhedron(istream & instr,
    vector<Point_3> & outVertices,
    vector<vector<unsigned int> > & outFaces);

/**
 * Helper for readPolyhedron.
 */
void readOFF(istream & instr,
    vector<Point_3> & outVertices,
    vector<vector<unsigned int> > & outFaces);

/**
 * Determine which facets of polyhedron 2 overlap facets of polyhedron 1.
 *
 * For each facet of polyhedron 2, return a list of facet numbers in polyhedron
 * 1 which are coplanar and overlap with it.
 */
vector<vector<unsigned int> > findIntersectingFacets(
    const vector<Point_3> & vertices1,
    const vector<vector<unsigned int> > & faces1,
    const vector<Point_3> & vertices2,
    const vector<vector<unsigned int> > & faces2);

/**
 * Main function for intersection, union and difference.
 */
int handleBooleans(int argc, char const** argv);

/**
 * Main function for facet inheritance (the only operation offered here that
 * doesn't use Nef polyhedra).
 */
int handleFacetInheritance(int argc, char const** argv);

/**
 * Load a Nef polyhedron from a single file specifying all its shells as
 * separate OFF files.  A brief header gives the number of shells, and the
 * remainder of the file contains concatenated OFF files.
 */
NefPolyhedron readMultiOFF(istream & instr);

/**
 * Write a Nef polyhedron in a face-vertex format.
 */
void writeNefPolyhedron(const NefPolyhedron & poly);




int main(int argc, char const** argv)
{
    cout << setprecision(15);
    for (int nn = 0; nn < argc; nn++)
    {
        cerr << "Arg " << nn << " = " << argv[nn] << "\n";
    }
    
    if (argc == 1)
    {
        printHelp();
        return 0;
    }
    
    int returnVal;
    string algorithmName(argv[1]);
    
    if (algorithmName == "inherit")
    {
        returnVal = handleFacetInheritance(argc, argv);
    }
    else
    {
        returnVal = handleBooleans(argc, argv);
    }
    
    return returnVal;
}



int handleBooleans(int argc, char const **argv)
{
    string algorithmName(argv[1]);
    
    // ------ Read two polyhedra, either from a file or from standard input
    
    NefPolyhedron p1, p2;
    if (argc == 3)
    {
        string fname(argv[2]);
        cerr << "Opening " << fname << "\n";
        
        ifstream in(fname.c_str());
        
        if (!in)
        {
            cerr << "Cannot open file.\n";
            return 1;
        }

        p1 = readMultiOFF(in);
        p2 = readMultiOFF(in);
    }
    else
    {
        p1 = readMultiOFF(cin);
        p2 = readMultiOFF(cin);
    }
    
    // ------ Go do the calculation
    
    NefPolyhedron result;
    if (algorithmName == "intersection")
    {
        result = doIntersection(p1, p2);
    }
    else if (algorithmName == "union")
    {
        result = doUnion(p1, p2);
    }
    else if (algorithmName == "difference")
    {
        result = doDifference(p1, p2);
    }
    else
    {
        printError();
        return 1;
    }
    
    writeNefPolyhedron(result);
    
    return 0;
}



int handleFacetInheritance(int argc, char const **argv)
{
    // ------ Read two polyhedra, either from a file or from standard input
    
    // Vertices and faces of the two polyhedra
    vector<Point_3> vertices1, vertices2;
    vector<vector<unsigned int> > faces1, faces2;
    
    // Indices into polyhedron 1, for each face of polyhedron 2
    vector<vector<unsigned int> > ancestorFaces;
    
    // For the facet inheritance algorithm, I just want to deal with
    // collections of facets.  NefPolyhedron will actually make this job
    // impossible because it automatically merges coplanar input facets.
    if (argc == 3)
    {
        string fname(argv[2]);
        cerr << "Opening " << fname << "\n";
        
        ifstream in(fname.c_str());
        
        if (!in)
        {
            cerr << "Cannot open file.\n";
            return 1;
        }
        
        readPolyhedron(in, vertices1, faces1);
        readPolyhedron(in, vertices2, faces2);
    }
    else
    {
        readPolyhedron(cin, vertices1, faces1);
        readPolyhedron(cin, vertices2, faces2);
    }
    
//    std::cerr << "From " << vertices1.size() << " verts and "
//        << faces1.size() << " faces to "
//        << vertices2.size() << " verts and "
//        << faces2.size() << " faces.\n";
    ancestorFaces = facetInheritance(vertices1, faces1, vertices2, faces2);
    
    for (int aa = 0; aa < ancestorFaces.size(); aa++)
    {
        std::cout << aa << ":\t";
//        std::cerr << "New face " << aa << " paternity: [ ";
        for (int bb = 0; bb < ancestorFaces[aa].size(); bb++)
        {
//            std::cerr << ancestorFaces[aa][bb] << " ";
            std::cout << ancestorFaces[aa][bb] << " ";
        }
//        std::cerr << "]\n";
        std::cout << "\n";
    }
    
    return 0;
}


void printHelp()
{
    cout << "Usage:\n";
    cout << "\tNefLab intersection\n";
    cout << "\tNefLab union\n";
    cout << "\tNefLab difference\n";
}

void printError()
{
    cout << "I don't understand what you want.\n";
    printHelp();
}

NefPolyhedron readMultiOFF(istream & instr)
{
    int numPositiveShells, numNegativeShells;
    int numIgnoredFacets;
    
//    cerr << "Bad() = " << instr.bad() << "\n";
//    cerr << "Fail() = " << instr.fail() << "\n";
//    cerr << "EOF() = " << instr.eof() << "\n";
//    cerr << "Good() = " << instr.good() << "\n";
    
    instr >> numPositiveShells;
    
//    cerr << "Bad() = " << instr.bad() << "\n";
//    cerr << "Fail() = " << instr.fail() << "\n";
//    cerr << "EOF() = " << instr.eof() << "\n";
//    cerr << "Good() = " << instr.good() << "\n";
    
    instr >> numNegativeShells;
    
    NefPolyhedron nef;
    
    for (int nn = 0; nn < numPositiveShells; nn++)
    {
        NefPolyhedron addend;
        
        numIgnoredFacets = CGAL::OFF_to_nef_3(instr, addend);
        if (numIgnoredFacets > 0)
            cerr << "Ignored " << numIgnoredFacets << " facets!\n";
        nef = nef + addend;
    }
    
    for (int nn = 0; nn < numNegativeShells; nn++)
    {
        NefPolyhedron subtrahend;
        
        numIgnoredFacets = CGAL::OFF_to_nef_3(instr, subtrahend);
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

void readPolyhedron(istream & instr,
    vector<Point_3> & outVertices,
    vector<vector<unsigned int> > & outFaces)
{
    int numPositiveShells, numNegativeShells;
    instr >> numPositiveShells >> numNegativeShells;
    
    int numShells = numPositiveShells + numNegativeShells;
    
    // Empty out the vertex and face vectors to begin.
    outVertices = vector<Point_3>();
    outFaces = vector<vector<unsigned int> >();
    
    int nextVertexNumber = 0;
    for (int ss = 0; ss < numShells; ss++)
    {
        std::vector<Point_3> shellVerts;
        std::vector<std::vector<unsigned int> > shellFaces;
        
        readOFF(instr, shellVerts, shellFaces);
        
        // increment vertex index into total number of vertices so we can
        // concatenate all the shells into one vertex list and one face list.
        for (int ff = 0; ff < shellFaces.size(); ff++)
        for (int ii = 0; ii < shellFaces[ff].size(); ii++)
        {
            shellFaces[ff][ii] += nextVertexNumber;
        }
        
        nextVertexNumber += shellVerts.size();
        
        outVertices.insert(outVertices.end(), shellVerts.begin(), shellVerts.end());
        outFaces.insert(outFaces.end(), shellFaces.begin(), shellFaces.end());
//        std::copy(shellVerts.begin(), shellVerts.end(), std::back_inserter(outVertices));
//        std::copy(shellFaces.begin(), shellFaces.end(), std::back_inserter(outFaces));
    }
}

void readOFF(istream & instr,
    vector<Point_3> & outVertices,
    vector<vector<unsigned int> > & outFaces)
{
    std::string shouldBeOFF;
    instr >> shouldBeOFF;
    assert(shouldBeOFF == "OFF");
    
    int numVertices, numFaces, numEdges_unused;
    instr >> numVertices >> numFaces >> numEdges_unused;
//    cerr << "Reading " << numVertices << " vertices and " << numFaces << " faces.\n";
    
    outVertices = vector<Point_3>(numVertices);
    outFaces = vector<vector<unsigned int> >(numFaces);
    
    for (int vv = 0; vv < numVertices; vv++)
    {
        double x, y, z;
        instr >> x >> y >> z;
        outVertices[vv] = Point_3(x, y, z);
    }
    
    for (int ff = 0; ff < numFaces; ff++)
    {
        unsigned int numVerts, v1, v2, v3;
        instr >> numVerts >> v1 >> v2 >> v3;
        vector<unsigned int> tri(3);
        tri[0] = v1;
        tri[1] = v2;
        tri[2] = v3;
        
        outFaces[ff] = tri;
    }
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
//        cerr << "Volume " << loopVertex << "\n";
        
        int shell = 0;
        ShellEntryConstIterator itr;
        for (itr = insideVolumes[loopVertex]->shells_begin();
            itr != insideVolumes[loopVertex]->shells_end(); itr++)
        {
//            cout << "\tShell " << shell++ << "\n";
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
    
    cout << setprecision(15);
    set<Vector3d>::const_iterator itr;
    for (itr = vertexSet.begin(); itr != vertexSet.end(); itr++)
    {
        if (vertexIds.count(*itr))
        {
            cerr << "Warning: duplicate vertex detected!";
        }
        vertexIds[*itr] = vertId++;
        Vector3d v = *itr;
        cout << vertId-1 << " " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }
    
    cout << "numFacets " << facets.size() << "\n";
    
    for (int cc = 0; cc < facets.size(); cc++)
    {
        cout << "\tnumContours " << facets[cc].size() << "\n";
        
        int vertexNum = 0;
        
        std::set<int> vertIdSet;
        
        // Iterate over facet loops
        for (int ll = 0; ll < facets[cc].size(); ll++)
        {
            cout << "\tnumVertices " << facets[cc][ll] << "\n";
            cout << "\t";
            
            
            for (int loopVertex = 0; loopVertex < facets[cc][ll]; loopVertex++)
            {
                int vId = vertexIds[contourVertices[cc][vertexNum++]];
                cout << vId << " ";
                
                if (vertIdSet.count(vId))
                {
                    cerr << "Warning: duplicate vertex in facet!";
                }
                vertIdSet.insert(vId);
            }
            cout << "\n";
        }
    }
}

