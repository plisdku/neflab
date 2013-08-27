#include "Triangulator.h"
#include <cstdlib>
#include <algorithm>
#include <stdint.h>

using namespace std;

namespace Triangulate
{

const double Triangulator::kINFINITY = 1e20;

const size_t VERTEX_INDEX_STACK_SIZE = 100000;

static ptrdiff_t myRandom(ptrdiff_t ii)
{
    return lrand48()%ii;
}

Triangulator::
Triangulator():
    mNumSegments(0),
    mNumTrapezoids(0),
    mNumQueryNodes(0),
    mNumQueryNodeLinks(0),
    mSegmentStore(100000),
    mTrapezoidStore(100000),
    mQueryNodeStore(100000),
    mQueryNodeLinkStore(1000000),
    mNextSegment(100000),
    mStartingTrapezoid(100000),
    mVertexChain(100000),
    mVertexIndexStack(),
    mTrapezoidUsed(100000, 0),
    mVertexInTrapezoidation(100000, 0),
    mSegmentInTriangulation(100000,0)
{
    mVertexIndexStack.reserve(VERTEX_INDEX_STACK_SIZE);
    
    // The randomness of this algorithm is part of how it can achieve high 
    // performance.  Unfortunately, randomness also means that I won't get the
    // exact same triangulation every time I triangulate a given polygon.  This
    // ends up causing problems for unit tests down the road (when I perturb
    // polygons and really want the triangulation to have proceeded in the same
    // way for two very similar polygons).  So, using the same random seed here
    // is necessary and, alas, will screw up my runtimes.
    //
    // A better solution might expose this re-seeding as an option that the
    // unit tests can tweak.
    srand48(0L);
    //cerr << "Warning: resetting random seed for triangulator.\n";
}

struct Projection
{
    Projection(Vector3d surfaceNormal)
    {
        int normalDirection = 0;
        if (fabs(surfaceNormal[1]) > fabs(surfaceNormal[0]))
            normalDirection = 1;
        if (fabs(surfaceNormal[2]) > fabs(surfaceNormal[normalDirection]))
            normalDirection = 2;
        
        // Pick a projection that preserves handedness.  The right projection
        // depends on which way the surface normal points.
        if (surfaceNormal[normalDirection] > 0)
        {
            //cout << "Up!\n";
            // three unit vectors
            m_u0 = (normalDirection+1)%3;
            m_u1 = (normalDirection+2)%3;
            m_u2 = normalDirection;
        }
        else
        {
            //cout << "Down!\n";
            m_u0 = (normalDirection+2)%3;
            m_u1 = (normalDirection+1)%3;
            m_u2 = normalDirection;
        }
    }
    
    Vector2d operator()(const Vector3d & v) const
    {
        return Vector2d(v[m_u0], v[m_u1]);
    }
    
    int m_u0;
    int m_u1;
    int m_u2;
};

int Triangulator::
triangulate(const vector<int> & contourSizes,
    const vector<Vector3d> & vertices,
    Vector3d surfaceNormal,
    vector<Vector3d*> & outTriangleIndices)
{
    assert(vertices.size() > 0);
    
    int numVertices, numTriangles;
    countPrimitives(contourSizes, numVertices, numTriangles);
    
    vector<Vector2d> vertices2d;
    vector<Vector2d*> triangleIndices(outTriangleIndices.size());
    
    // Project the polygon onto the XY, YZ or ZX plane.  Choose based on the
    // smallest dimension of the bounding box over all the vertices.
    
    transform(vertices.begin(), vertices.end(), back_inserter(vertices2d),
        Projection(surfaceNormal));
    
    triangulate(contourSizes, vertices2d, triangleIndices);
    
    for (int nn = 0; nn < 3*numTriangles; nn++)
    {
        ptrdiff_t offset = triangleIndices[nn] - &vertices2d[0];
        outTriangleIndices.at(nn) = const_cast<Vector3d*>(&vertices[0]) + offset;
    }
    
    return numTriangles;
}

int Triangulator::
triangulate(const vector<int> & contourSizes,
    const vector<Vector2d> & vertices,
    vector<Vector2d*> & outTriangleIndices)
{
    int numVertices, numTriangles;
    countPrimitives(contourSizes, numVertices, numTriangles);
    
    if (signedArea(vertices.begin(), vertices.begin()+contourSizes[0]) < 0)
    {
        signedArea(vertices.begin(), vertices.begin()+contourSizes[0]);
    }
    assert(signedArea(vertices.begin(), vertices.begin()+contourSizes[0]) > 0);
    
    initSegments(contourSizes, vertices);
    //cout << "Shuffling segments.\n";
    random_shuffle(segmentsBegin(), segmentsEnd(), myRandom);
    //reverse(segmentsBegin(), segmentsEnd());
    
    for (vector<Segment>::iterator itr = segmentsBegin();
        itr != segmentsEnd(); itr++)
    {
        mNextSegment[&itr->v0() - mFirstUserVertex] = &(*itr);
    }
    
    mNumTrapezoids = 0;
    mNumQueryNodes = 0;
    mNumQueryNodeLinks = 0;
    
    /*
    for (int mm = 0; mm < mVertexInTrapezoidation.size(); mm++)
        assert(mVertexInTrapezoidation[mm] == 0);
    for (int mm = 0; mm < mTrapezoidUsed.size(); mm++)
        assert(mTrapezoidUsed[mm] == 0);
    for (int mm = 0; mm < mSegmentInTriangulation.size(); mm++)
        assert(mSegmentInTriangulation[mm] == 0);
    */
    
    // Initialize the segment array.  Walk each contour and put the segments
    // in order.
    //
    // Anyway, the segments all go from one vertex to the next...
    
    // Now trapezoidalize!
    
    Trapezoid* firstTrapezoid = newTrapezoid(); // Represents the WHOLE PLANE.
    QueryNode* firstQueryNode = newQueryNode();
    firstQueryNode->trapezoid(firstTrapezoid);
    firstTrapezoid->queryNode(firstQueryNode);
    
    vector<Segment>::const_iterator itr;
    
//    for (itr = segmentsBegin(); itr != segmentsEnd(); itr++)
//    {
//        int num = itr - segmentsBegin();
//        
//        cout << "Segment " << num << " [" << itr->v0() << ", "
//            << itr->v1() << "]\n";
//    }
    
    
    int nSeg = 0;
    for (itr = segmentsBegin(); itr != segmentsEnd(); itr++)
    {
        //cout << "Add segment " << nSeg++ << endl;
        addSegmentToTrapezoidation(*itr);
    }
    
    //cout << "Done with trapezoidation.\n";
    //printTrapezoids();
    //checkPointers();
    
    //plotSegments();
    //plotTrapezoids();
    
    //cout << "Traverse!\n";
    traverseChains(outTriangleIndices);
    
    // Erase the flags we used
    fill(mVertexInTrapezoidation.begin(),
        mVertexInTrapezoidation.begin() + vertices.size(), 0);
    fill(mTrapezoidUsed.begin(),
        mTrapezoidUsed.begin() + mNumTrapezoids, 0);
    fill(mSegmentInTriangulation.begin(),
        mSegmentInTriangulation.begin() + mNumSegments, 0);
    
    return numTriangles;
}

void Triangulator::
countPrimitives(const vector<int> & contourSizes, int & outNumVertices,
    int & outNumTriangles) const
{
    int numVerts = 0, numHoles = contourSizes.size()-1;
    for (int mm = 0; mm < contourSizes.size(); mm++)
        numVerts += contourSizes[mm];
    
    outNumVertices = numVerts;
    outNumTriangles = numVerts - 2 + 2*numHoles;
}

void Triangulator::
initSegments(const vector<int> & contourSizes,
    const vector<Vector2d> & vertices)
{
    mNumSegments = 0;
    mFirstUserVertex = &(vertices[0]);
    
    int vv = 0;
    for (uint32_t cc = 0; cc < contourSizes.size(); cc++)
    {
        int vFirst = vv;
        int vLast = vv + contourSizes[cc] - 1;
        
        for (vv = vFirst; vv <= vLast; vv++)
        {
            const Vector2d* v0 = &(vertices[vv]);
            const Vector2d* v1;
            
            if (vv+1 <= vLast)
                v1 = &(vertices[vv+1]);
            else
                v1 = &(vertices[vFirst]);
            
            Segment* segment = newSegment();
            segment->v0(v0);
            segment->v1(v1);
        }
    }
}

void Triangulator::
addSegmentToTrapezoidation(const Segment & segment)
{
    vector<Trapezoid>::const_iterator itr;
    
//    cout << "QNL " << mNumQueryNodeLinks << "\n";
    
    assert(validateQueryTree());
//    cout << "Add " << segment.v0() << " to " << segment.v1() << endl;
    
    // If first vertex is new, split the trapezoid it hits.
    if (false == isVertexInTrapezoidation(&segment.v0()))
    {
        Trapezoid* trap = enclosingTrapezoid(segment.v0());
        assert(trap != 0L);
//        assert(cout << hex << enclosingTrapezoid(segment.v0()) << dec << "\n");
//        assert(cout << hex << enclosingTrapezoid_safe(segment.v0()) << dec
//            << "\n");
        assert(trap == enclosingTrapezoid_safe(segment.v0()));
        assert(trap->pointersAreValid());
        
        Trapezoid* bottomHalf = splitVertically(trap, segment.v0());
        assert(validateQueryTree());
        assert(bottomHalf->pointersAreValid());
        assert(trap->pointersAreValid());
        
//        cout << "Vertex " << segment.v0() << " induced trapezoids:\n";
//        print(*trap);
//        print(*bottomHalf);
        
        flagVertexInTrapezoidation(&segment.v0());
    }
    
    // If second vertex is new, split the trapezoid it hits.
    if (false == isVertexInTrapezoidation(&segment.v1()))
    {
        Trapezoid* trap = enclosingTrapezoid(segment.v1());
        assert(trap != 0L);
        assert(trap == enclosingTrapezoid_safe(segment.v1()));
        assert(trap->pointersAreValid());
        
        Trapezoid* bottomHalf = splitVertically(trap, segment.v1());
        assert(validateQueryTree());
        assert(bottomHalf->pointersAreValid());
        assert(trap->pointersAreValid());
        
//        cout << "Vertex " << segment.v1() << " induced trapezoids:\n";
//        print(*trap);
//        print(*bottomHalf);
        
        flagVertexInTrapezoidation(&segment.v1());
    }
    
    assert(checkPointers());
    
    // Split the trapezoids we're crossing.  But which one are we crossing?
    Trapezoid* lowestTrapToSplit = lowestIntersectingTrapezoid(segment);
    assert(lowestTrapToSplit);
    assert(lowestTrapToSplit == lowestIntersectingTrapezoid_safe(segment));
    Trapezoid* newRightHalf;
    
    if (lowestTrapToSplit->bottomRight())
    {
        newRightHalf = splitHorizontally(
            lowestTrapToSplit, lowestTrapToSplit->bottomLeft(),
            lowestTrapToSplit->bottomRight(), segment);
        assert(validateQueryTree());
    }
    else
    {
        // use redundant bottomLeft neighbor to indicate that there's only one
        // neighbor underneath
        newRightHalf = splitHorizontally(
            lowestTrapToSplit, lowestTrapToSplit->bottomLeft(),
            lowestTrapToSplit->bottomLeft(), segment);
        assert(validateQueryTree());
    }
    
//    cout << "Added!\n";
//    printTrapezoids();
    assert(checkPointers());
    
    // Now re-merge trapezoids.
    Trapezoid* lowL = lowestTrapToSplit;
    Trapezoid* previous = lowL;
    Trapezoid* tt;
    
    // Merging on the left of the segment.  The current trapezoid is tt, and
    // it will absorb trapezoids below it.
    for (tt = lowL->upperRightNeighbor();
        tt != 0L && tt->isBelow(segment.topVertex());
        tt = tt->upperRightNeighbor())
    {
        if (tt->leftSegment() == previous->leftSegment())
        {
            // tt will absorb previous.
            assert(tt->pointersAreValid());
            
            assert(validateQueryTree());
            // Step 0: connect the query nodes
            QueryNode* q = tt->queryNode();
            QueryNode* prevNode = previous->queryNode();
            prevNode->replaceWith(q);
            for (QueryNodeLink* link = prevNode->parentListHead();
                link != 0L; link = link->next)
            {
                q->addParentNode(link->node, newQueryNodeLink());
            }
            assert(validateQueryTree());
            
            // Step 1: connect the trapezoid topology
            tt->bottomVertex(previous->bottomVertex());
            tt->bottomLeft(previous->bottomLeft());
            tt->bottomRight(previous->bottomRight());
            if (previous->bottomLeft())
            {
                if (previous == previous->bottomLeft()->topLeft())
                    previous->bottomLeft()->topLeft(tt);
                else
                {
                    assert(previous == previous->bottomLeft()->topRight());
                    previous->bottomLeft()->topRight(tt);
                }
            }
            if (previous->bottomRight())
            {
                if (previous == previous->bottomRight()->topLeft())
                    previous->bottomRight()->topLeft(tt);
                else
                {
                    assert(previous == previous->bottomRight()->topRight());
                    previous->bottomRight()->topRight(tt);
                }
            }
            assert(tt->pointersAreValid());
            assert(validateQueryTree());
            deleteTrapezoid(previous);
            assert(validateQueryTree());
            assert(tt->pointersAreValid());
        }
        previous = tt;
    }
    
    // Merging on the right of the segment.
    previous = newRightHalf;
    
    for (tt = newRightHalf->topLeft();
        tt != 0L && tt->isBelow(segment.topVertex());
        tt = tt->topLeft())
    {
        if (tt->rightSegment() == previous->rightSegment())
        {
            // tt will absorb previous
            assert(tt->pointersAreValid());
            
            // Step 0: connect the query nodes
            assert(validateQueryTree());
            QueryNode* q = tt->queryNode();
            QueryNode* prevNode = previous->queryNode();
            prevNode->replaceWith(q);
            for (QueryNodeLink* link = prevNode->parentListHead();
                link != 0L; link = link->next)
            {
                q->addParentNode(link->node, newQueryNodeLink());
            }
            assert(validateQueryTree());
            
            // Step 1: connect the trapezoid topology
            tt->bottomVertex(previous->bottomVertex());
            tt->bottomLeft(previous->bottomLeft());
            tt->bottomRight(previous->bottomRight());
            if (previous->bottomLeft())
            {
                if (previous == previous->bottomLeft()->topLeft())
                    previous->bottomLeft()->topLeft(tt);
                else
                {
                    assert(previous == previous->bottomLeft()->topRight());
                    previous->bottomLeft()->topRight(tt);
                }
            }
            if (tt->bottomRight())
            {
                if (previous == previous->bottomRight()->topLeft())
                    previous->bottomRight()->topLeft(tt);
                else
                {
                    assert(previous == previous->bottomRight()->topRight());
                    previous->bottomRight()->topRight(tt);
                }
            }
            assert(tt->pointersAreValid());
//            cout << "tt = " << tt << endl;
            assert(validateQueryTree());
            deleteTrapezoid(previous);
            assert(validateQueryTree());
            assert(tt->pointersAreValid());
        }
        previous = tt;
    }
    
//    cout << "Merged!\n";
//    printTrapezoids();
    assert(checkPointers());
}

bool Triangulator::
isVertexInTrapezoidation(const Vector2d* vertex) const
{
    return mVertexInTrapezoidation.at(vertex - mFirstUserVertex);
}

void Triangulator::
flagVertexInTrapezoidation(const Vector2d * vertex)
{
    mVertexInTrapezoidation.at(vertex - mFirstUserVertex) = true;
}

void Triangulator::
traverseChains(vector<Vector2d*> & outTriangleIndices)
{
    /*
        The rules: at the beginning of the do-loop,
        
        currentVertex       has been added to a vertex chain already
        currentSegment      begins at currentVertex (so, not traversed yet!)
        currentTrapezoid    is bounded by currentSegment (left or right) and
                            currentVertex (top or bottom)
        
        Which vertex do we move to?
        - If there is a diagonal that we haven't gone across, go across it.
          Obtain new currentSegment as mNextSegment[newVertex].
          Obtain new currentTrapezoid by graph traversal (complicated-ish).
        - Otherwise, proceed along boundary of polygon.
          Obtain new currentSegment as mNextSegment[currentSegment->v1()].
          Obtain new currentVertex as (new) currentSegment->v0().
          Obtain new currentTrapezoid by graph traversal (complicated-ish).
    */
    
    const Segment* segBegin = &(*segmentsBegin());
    
    const Vector2d* currentVertex;
    const Vector2d* lastVertex;
    const Segment* currentSegment;
    const Segment* lastSegment;
    const Trapezoid* currentTrapezoid;
    
    int triangleNum = 0;
    
    // Set up a segment-to-trapezoid mapping.
    vector<Trapezoid>::const_iterator itr;
    for (itr = trapezoidsBegin(); itr != trapezoidsEnd(); itr++)
    if (itr->isInside())
    {
        if (itr->topVertex() == &itr->leftSegment()->v0())
            mStartingTrapezoid[itr->leftSegment() - segBegin] = &(*itr);
        if (itr->bottomVertex() == &itr->rightSegment()->v0())
            mStartingTrapezoid[itr->rightSegment() - segBegin] = &(*itr);
    }
    
    // Pick the first segment, trapezoid and vertex
    vector<Segment>::const_iterator startingSegmentItr = segmentsBegin();
    currentSegment = &(*startingSegmentItr++);
    currentTrapezoid = mStartingTrapezoid[currentSegment - segBegin];
    currentVertex = &(currentSegment->v0());
    lastSegment = 0L;
    lastVertex = 0L;
    
    int numVertexChainsDone = 0;
    while ("bad design")
    {
//        cout << "******** NEW POLYGON from " << *currentVertex << "\n";
        
        int vertexChainLength = 0;
        const Vector2d* firstVertex = currentVertex;
        do {
//            cout << "----------------------\n";
//            cout << "Current trapezoid:\n";
//            print(*currentTrapezoid);
            
            // if (the trapezoid is non-diagonal, OR
            //     this is the first segment of this vertex chain, OR
            //     there's a diagonal but we just used it)
            //      advance along current segment
            //  else
            //      advance along the diagonal
            
            mVertexChain[vertexChainLength++] = currentVertex;
            
            if ( (false == currentTrapezoid->willBeCutDiagonally()) ||
                (lastVertex == 0L) ||
                (lastVertex == currentTrapezoid->oppositeVertex(currentVertex)))
            {
//                cout << "Move along edge." << endl;
                // either
                //      we were moving along the boundary and there is no diagonal
                // or
                //      we just crossed a diagonal
                // so
                //      continue forward along the boundary (don't cross diagonally)
                lastVertex = currentVertex;
                lastSegment = currentSegment;
                currentVertex = &currentSegment->v1();
                currentSegment = mNextSegment[currentVertex - mFirstUserVertex];
                
//                cout << "Marking ";
//                print(*lastSegment);
//                cout << "\n";
                mSegmentInTriangulation[lastSegment - segBegin] = true;
                
                if (isBelow(*lastVertex, *currentVertex)) // if ascending
                {
//                    cout << "\tascending:" << endl;
                    while (currentTrapezoid->topVertex() != currentVertex)
                    {
                        currentTrapezoid = currentTrapezoid->upperRightNeighbor();
                        assert(currentTrapezoid);
                    }
                    if ( (false == currentTrapezoid->willBeCutDiagonally()) &&
                        currentTrapezoid->numTopNeighbors() != 0)
                    {
                        currentTrapezoid = currentTrapezoid->upperRightNeighbor();
                        assert(currentTrapezoid);
                    }
                }
                else // if descending
                {
//                    cout << "\tdescending:" << endl;
                    while (currentTrapezoid->bottomVertex() != currentVertex)
                    {
                        currentTrapezoid = currentTrapezoid->bottomLeft();
                        assert(currentTrapezoid);
                    }
                    if ( (false == currentTrapezoid->willBeCutDiagonally()) &&
                        currentTrapezoid->numBottomNeighbors() != 0)
                    {
                        currentTrapezoid = currentTrapezoid->bottomLeft();
                        assert(currentTrapezoid);
                    }
                }
            }
            else
            {
                // we were moving along the polygon boundary.
                // we have reached a diagonal "cut" across the middle.
                // move diagonally across the middle.
                
//                cout << "Move across diagonal." << endl;
                assert(currentTrapezoid->willBeCutDiagonally());
                
                lastVertex = currentVertex;
                lastSegment = currentSegment;
                currentVertex = currentTrapezoid->oppositeVertex(currentVertex);
                currentSegment = mNextSegment[currentVertex - mFirstUserVertex];
                
                if (currentVertex == &currentTrapezoid->rightSegment()->v1())
                {
//                    cout << "\tascending right:";
                    if (currentTrapezoid->topLeft())
                    {
//                        cout << " top has neighbors" << endl;
                        currentTrapezoid = currentTrapezoid->topLeft();
                    }
                    else
                    {
//                        cout << " top is triangular" << endl;
                    }
                }
                else if (currentVertex ==&currentTrapezoid->leftSegment()->v1())
                {
//                    cout << "\tdescending left:";
                    if (currentTrapezoid->bottomRightNeighbor())
                    {
//                        cout << " bottom has neighbors" << endl;
                        currentTrapezoid =
                            currentTrapezoid->bottomRightNeighbor();
                    }
                    else
                    {
//                        cout << " bottom is triangular" << endl;
                    }
                }
                else if (currentVertex == &currentTrapezoid->rightSegment()->v0())
                {
//                    cout << "\tdescending right:" << endl;
                    // stay in same trapezoid
                }
                else if (currentVertex == &currentTrapezoid->leftSegment()->v0())
                {
//                    cout << "\tascending left:" << endl;
                    // stay in same trapezoid
                }
                else if (isBelow(*lastVertex, *currentVertex))
                {
//                    cout << "\tascending center:" << endl;
                    currentTrapezoid = currentTrapezoid->topLeft();
                }
                else
                {
//                    cout << "\tdescending center:" << endl;
                    currentTrapezoid = currentTrapezoid->bottomRightNeighbor();
                }
                
                assert(currentTrapezoid);
            }
//            cout << "Vertex: " << *currentVertex << "\n";
        } while (currentVertex != firstVertex);
        
        numVertexChainsDone++;
        
        triangleNum = triangulateChain(vertexChainLength, triangleNum,
            outTriangleIndices);
        
        while ( (startingSegmentItr != segmentsEnd()) &&
            (mSegmentInTriangulation[startingSegmentItr-segmentsBegin()] == 1))
        {
            startingSegmentItr++;
        }
        
        if (startingSegmentItr == segmentsEnd())
        {
//            cout << "We're done: " << numVertexChainsDone << " chains.\n";
            return;
        }
        
        currentSegment = &(*startingSegmentItr++);
        currentTrapezoid = mStartingTrapezoid[currentSegment - segBegin];
        currentVertex = &(currentSegment->v0());
        lastSegment = 0L;
        lastVertex = 0L;
        
//        cout << "New current segment ";
//        print(*currentSegment);
//        cout << "\n";
//        cout << "Marked? " << mSegmentInTriangulation[currentSegment - segBegin]
//            << "\n";
    } // while 1
}

static bool vectorLessThan(const Vector2d* lhs, const Vector2d* rhs)
{
    return Triangulator::isBelow(*lhs, *rhs);
}

int Triangulator::
triangulateChain(int chainLength, int triangleNum,
    vector<Vector2d*> & outTriangleIndices)
{
    //mVertexIndexStack.clear();
    assert(mVertexIndexStack.size() == 0);
    assert(mVertexIndexStack.capacity() == VERTEX_INDEX_STACK_SIZE);
//    cout << "CHAIN:\n";
//    for (int nn = 0; nn < chainLength; nn++)
//    {
//        cout << *mVertexChain[nn] << "\n";
//    }
    
    // Sort the chain
    sort(mVertexChain.begin(), mVertexChain.begin() + chainLength, vectorLessThan);
//    cout << "SORTED CHAIN:\n";
//    for (int nn = 0; nn < chainLength; nn++)
//    {
//        cout << *mVertexChain[nn] << "\n";
//    }
    
    // Is it flat on the right or flat on the left?
    if (cross(*mVertexChain[1]-*mVertexChain[0],
        *mVertexChain[chainLength-1] - *mVertexChain[0]) > 0)
    {
        // Second vertex is to the right of top vertex.  Flat on the left!
//        cout << "Flat on the left.\n";
        
        int nextVertIndex = 0;
//        cout << "Push " << nextVertIndex << endl;
        mVertexIndexStack.push_back(nextVertIndex++);
        while (nextVertIndex < chainLength-1)
        {
//            cout << "Push " << nextVertIndex << endl;
            mVertexIndexStack.push_back(nextVertIndex++);
            
            if (mVertexIndexStack.size() > 1)
            {
                const Vector2d* prevVert = mVertexChain[mVertexIndexStack.at(
                    mVertexIndexStack.size()-2)];
                const Vector2d* topVert = mVertexChain[mVertexIndexStack.at(
                    mVertexIndexStack.size()-1)];
                const Vector2d* nextVert = mVertexChain[nextVertIndex];
                
                while (mVertexIndexStack.size() > 1 &&
                    cross(*nextVert - *topVert, *topVert - *prevVert) < 0)
                {
//                    cout << mVertexIndexStack[mVertexIndexStack.size()-1] << " is ear" << endl;
                    
                    int ii = 3*triangleNum;
                    outTriangleIndices[ii] = const_cast<Vector2d*>(prevVert);
                    outTriangleIndices[ii+1] = const_cast<Vector2d*>(topVert);
                    outTriangleIndices[ii+2] = const_cast<Vector2d*>(nextVert);
                    triangleNum++;
                    
                    mVertexIndexStack.pop_back();
                    
                    if (mVertexIndexStack.size() > 1)
                    {
                        prevVert = mVertexChain[mVertexIndexStack.at(
                            mVertexIndexStack.size()-2)];
                        topVert = mVertexChain[mVertexIndexStack.at(
                            mVertexIndexStack.size()-1)];
                        nextVert = mVertexChain[nextVertIndex];
                    }
                }
            }
        }
    }
    else
    {
        // Second vertex is to the left of top vertex.  Flat on the right!
//        cout << "Flat on the right.\n";
        
        int nextVertIndex = 0;
//        cout << "Push " << nextVertIndex << endl;
        mVertexIndexStack.push_back(nextVertIndex++);
        while (nextVertIndex < chainLength-1)
        {
//            cout << "Push " << nextVertIndex << endl;
            mVertexIndexStack.push_back(nextVertIndex++);
            
            if (mVertexIndexStack.size() > 1)
            {
                const Vector2d* prevVert = mVertexChain[mVertexIndexStack.at(mVertexIndexStack.size()-2)];
                const Vector2d* topVert = mVertexChain[mVertexIndexStack.at(mVertexIndexStack.size()-1)];
                const Vector2d* nextVert = mVertexChain[nextVertIndex];
                
//                cout << "try (" << mVertexIndexStack.at(mVertexIndexStack.size()-2)
//                    << ", " << mVertexIndexStack.at(mVertexIndexStack.size()-1) << ", "
//                    << nextVertIndex << ")" << endl;
                
                while (mVertexIndexStack.size() > 1 &&
                    cross(*nextVert - *topVert, *topVert - *prevVert) > 0)
                {
//                    cout << mVertexIndexStack[mVertexIndexStack.size()-1] << " is ear" << endl;
                    
                    int ii = 3*triangleNum;
                    outTriangleIndices[ii] = const_cast<Vector2d*>(prevVert);
                    outTriangleIndices[ii+1] = const_cast<Vector2d*>(nextVert);
                    outTriangleIndices[ii+2] = const_cast<Vector2d*>(topVert);
                    triangleNum++;
                    
                    mVertexIndexStack.pop_back();
                    
                    if (mVertexIndexStack.size() > 1)
                    {
                        prevVert = mVertexChain[mVertexIndexStack.at(mVertexIndexStack.size()-2)];
                        topVert = mVertexChain[mVertexIndexStack.at(mVertexIndexStack.size()-1)];
                        nextVert = mVertexChain[nextVertIndex];
                    }
                }
            }
        }
    }
    
    assert(mVertexIndexStack.size() == 1);
    mVertexIndexStack.pop_back();
    
    // Chain is unimonotone, flat on the left?
    //  sort vertices bottom-to-top
    // Chain is unimonotone, flat on the right?
    //  sort vertices top-to-bottom
    
    return triangleNum;
}


void Triangulator::
plotSegments()
{
    vector<Segment>::const_iterator itr;
    
    cout << "xx = [";
    for (itr = segmentsBegin(); itr != segmentsEnd(); itr++)
        cout << itr->v0()[0] << " " << itr->v1()[0] << "; ";
    cout << "];\n";
    cout << "yy = [";
    for (itr = segmentsBegin(); itr != segmentsEnd(); itr++)
        cout << itr->v0()[1] << " " << itr->v1()[1] << "; ";
    cout << "];\n";
    cout << "plot(xx',yy');\n";
}

void Triangulator::
plotTrapezoids()
{
    vector<Trapezoid>::const_iterator itr;
    
    cout << "tx = [";
    for (itr = trapezoidsBegin(); itr != trapezoidsEnd(); itr++)
    if (mTrapezoidUsed[&(*itr) - &(mTrapezoidStore[0])])
    if (itr->isInside())
    {
        if (itr->leftSegment()->height() != 0)
        {
            cout << itr->leftSegment()->x((*itr->bottomVertex())[1]) << " ";
            cout << itr->leftSegment()->x((*itr->topVertex())[1]) << " ";
        }
        else
        {
            cout << itr->leftSegment()->v0()[0] << " "
                << itr->leftSegment()->v1()[0] << " ";
        }
        
        if (itr->rightSegment()->height() != 0)
        {
            cout << itr->rightSegment()->x((*itr->topVertex())[1]) << " ";
            cout << itr->rightSegment()->x((*itr->bottomVertex())[1]) << " ";
        }
        else
        {
            cout << itr->rightSegment()->v0()[0] << " "
                << itr->rightSegment()->v1()[0] << " ";
        }
        
        if (itr->leftSegment()->height() != 0)
            cout << itr->leftSegment()->x((*itr->bottomVertex())[1]);
        else
            cout << itr->leftSegment()->v0()[0];
        
        cout << "; ";
    }
    cout << "];\n";
    
    cout << "ty = [";
    for (itr = trapezoidsBegin(); itr != trapezoidsEnd(); itr++)
    if (mTrapezoidUsed[&(*itr) - &(mTrapezoidStore[0])])
    if (itr->isInside())
    {
        cout << (*itr->bottomVertex())[1] << " "
            << (*itr->topVertex())[1] << " "
            << (*itr->topVertex())[1] << " "
            << (*itr->bottomVertex())[1] << " "
            << (*itr->bottomVertex())[1] << "; ";
    }
    cout << "];\n";
    cout << "plot(tx', ty');\n";
}

Trapezoid* Triangulator::
enclosingTrapezoid_safe(const Vector2d & vertex)
{
    Trapezoid* found = 0L;
    
    vector<Trapezoid>::iterator itr;
    for (itr = trapezoidsBegin(); itr != trapezoidsEnd(); itr++)
    if (mTrapezoidUsed[&(*itr) - &(mTrapezoidStore[0])])
    {
        if (itr->encloses(vertex))
        {
            if (found == 0L)
                found = &(*itr);
            else
                cerr << "Duplicate enclosing trapezoid!\n";
        }
    }
    if (found == 0L)
        cerr << "No enclosing trapezoid found!\n";
    return found;
}

Trapezoid* Triangulator::
enclosingTrapezoid(const Vector2d & vertex)
{
    assert(mQueryNodeStore.size() > 0);
    QueryNode* currentNode = &(mQueryNodeStore[0]);
    
    while (currentNode->type() != QueryNode::kLeafNode)
    {
        if (currentNode->type() == QueryNode::kYNode)
        {
            assert(currentNode->keyY());
            assert(vertex != *currentNode->keyY());
            if (isBelow(vertex, *currentNode->keyY()))
                currentNode = currentNode->lesserNode();
            else
                currentNode = currentNode->greaterNode();
        }
        else if (currentNode->type() == QueryNode::kXNode)
        {
            assert(currentNode->keyX());
            if (isToLeft(vertex, *currentNode->keyX()))
                currentNode = currentNode->lesserNode();
            else
                currentNode = currentNode->greaterNode();
        }
    }

    assert(mTrapezoidUsed[currentNode->trapezoid() - &(mTrapezoidStore[0])]);
    assert(currentNode->trapezoid()->encloses(vertex));
    assert(currentNode->trapezoid());
    return currentNode->trapezoid();
}

Trapezoid* Triangulator::
lowestIntersectingTrapezoid_safe(const Segment & segment)
{
    Trapezoid* found = 0L;
    const Vector2d* bottomVertex = &(segment.bottomVertex());
    
    Vector2d segmentVector = segment.topVertex() - segment.bottomVertex();
    
    vector<Trapezoid>::iterator itr;
    for (itr = trapezoidsBegin(); itr != trapezoidsEnd(); itr++)
    if (mTrapezoidUsed[&(*itr) - &(mTrapezoidStore[0])])
    {
        bool thisIsIt = 0;
        if (itr->bottomVertex() == bottomVertex)
        {
            // No more than two trapezoids should share this vertex, I think.
            if (itr->isBoundedLeft() &&
                itr->bottomVertex() == &(itr->leftSegment()->bottomVertex()))
            {
                // segment shares bottom vertex with itr->leftSegment(), so
                // segment must be to the RIGHT of itr->leftSegment().
                Vector2d leftVector = itr->leftSegment()->topVertex() -
                    itr->leftSegment()->bottomVertex();
                if (cross(segmentVector, leftVector) > 0)
                    thisIsIt = 1;
            }
            else if (itr->isBoundedRight() && 
                itr->bottomVertex() == &(itr->rightSegment()->bottomVertex()))
            {
                // segment shares bottom vertex with itr->rightSegment(), so
                // segment must be to the LEFT if itr->rightSegment().
                Vector2d rightVector = itr->rightSegment()->topVertex() -
                    itr->rightSegment()->bottomVertex();
                if (cross(rightVector, segmentVector) > 0)
                    thisIsIt = 1;
            }
            else
            {
                thisIsIt = 1;
            }
        }
        if (thisIsIt)
        {
            if (found == 0L)
                found = &(*itr);
            else
                cerr << "Found more than one lowest intersecting trapezoids!\n";
        }
    }
    
    return found;
}

Trapezoid* Triangulator::
lowestIntersectingTrapezoid(const Segment & segment)
{
    const Vector2d* bottomVertex = &(segment.bottomVertex());
    Vector2d segmentVector = segment.topVertex() - segment.bottomVertex();
    
//    cout << "Testing ";
//    print(segment);
//    cout << "\n";
    
    assert(mQueryNodeStore.size() > 0);
    QueryNode* currentNode = &(mQueryNodeStore[0]);
    
    while (currentNode->type() != QueryNode::kLeafNode)
    {
        if (currentNode->type() == QueryNode::kXNode)
        {
            assert(currentNode->keyX());
            
//            cout << "X node: ";
//            print(*currentNode->keyX());
//            cout << ": ";
            
            // A precondition for isToLeft(vertex, segment) is that the vertex
            // does not coincide with the segment.
            if (&currentNode->keyX()->bottomVertex() == bottomVertex)
            {
                Vector2d keyVector = currentNode->keyX()->topVertex()
                    - currentNode->keyX()->bottomVertex();
                if (cross(keyVector, segmentVector) < 0)
                {
                    //cout << "To the right, shared bottom vertex.\n";
                    currentNode = currentNode->greaterNode();
                }
                else
                {
                    //cout << "To the left, shared bottom vertex.\n";
                    currentNode = currentNode->lesserNode();
                }
            }
            else if (isToLeft(*bottomVertex, *currentNode->keyX()))
            {
                //cout << "To the left.\n";
                currentNode = currentNode->lesserNode();
            }
            else
            {
                //cout << "To the right.\n";
                currentNode = currentNode->greaterNode();
            }
        }
        else if (currentNode->type() == QueryNode::kYNode)
        {
            assert(currentNode->keyY());
            //cout << "Y node: " << *currentNode->keyY() << ": ";
            if (isBelow(*bottomVertex, *currentNode->keyY()))
            {
                //cout << "Below.\n";
                currentNode = currentNode->lesserNode();
            }
            else
            {
                //cout << "Above.\n";
                currentNode = currentNode->greaterNode();
            }
        }
    }
    assert(currentNode->trapezoid());
    assert(currentNode->trapezoid()->bottomVertex() == bottomVertex);
    
//    cout << "Found: ";
//    print(*currentNode->trapezoid());
//    cout << endl;
    
    return currentNode->trapezoid();
}

Trapezoid* Triangulator::
splitVertically(Trapezoid* t, const Vector2d & v)
{
    Trapezoid* bottomHalf = newTrapezoid();
    
    // Step 0: split the query node
    QueryNode* q = t->queryNode();
    assert(q);
    assert(q->type() == QueryNode::kLeafNode);
    q->type(QueryNode::kYNode);
    q->keyY(&v);
    QueryNode* bottomNode = newQueryNode();
    QueryNode* topNode = newQueryNode();
    q->lesserNode(bottomNode);
    q->greaterNode(topNode);
    bottomNode->addParentNode(q, newQueryNodeLink());
    //bottomNode->parentNode(q);
    topNode->addParentNode(q, newQueryNodeLink());
    //topNode->parentNode(q);
    bottomNode->trapezoid(bottomHalf);
    bottomHalf->queryNode(bottomNode);
    topNode->trapezoid(t);
    t->queryNode(topNode);
    assert(bottomNode->type() == QueryNode::kLeafNode);
    assert(topNode->type() == QueryNode::kLeafNode);
    
    // Step 1: copy segments and top/bottom from old trap to new bottom half
    bottomHalf->leftSegment(t->leftSegment());
    bottomHalf->rightSegment(t->rightSegment());
    bottomHalf->topVertex(&v);
    bottomHalf->bottomVertex(t->bottomVertex());
    t->bottomVertex(&v);
    
    // Step 2: re-link bottom neighbors to the new bottom half trapezoid
    if (t->bottomLeft())
    {
        if (t == t->bottomLeft()->topLeft())
            t->bottomLeft()->topLeft(bottomHalf);
        else if (t == t->bottomLeft()->topRight())
            t->bottomLeft()->topRight(bottomHalf);
    }
    if (t->bottomRight())
    {
        if (t == t->bottomRight()->topLeft())
            t->bottomRight()->topLeft(bottomHalf);
        if (t == t->bottomRight()->topRight())
            t->bottomRight()->topRight(bottomHalf);
    }
    
    // Step 3: link new bottom half trapezoid to bottom neighbors; link the
    // old trap and the new trap
    bottomHalf->bottomLeft(t->bottomLeft());
    bottomHalf->bottomRight(t->bottomRight());
    bottomHalf->topLeft(t);
    t->bottomLeft(bottomHalf);
    t->bottomRight(0L);
    
    return bottomHalf;
}

Trapezoid* Triangulator::
splitHorizontally(Trapezoid* t, Trapezoid* oldBottomLeft,
    Trapezoid* newBottomRight, const Segment & segment)
{
    assert(t->pointersAreValid());
    Trapezoid* rightHalf = newTrapezoid();
    const Segment* originalRightSegment = t->rightSegment();
    const Segment* originalLeftSegment = t->leftSegment();
    
    // Step 0.  split the query node
    QueryNode* q = t->queryNode();
    assert(q);
    assert(q->type() == QueryNode::kLeafNode);
    q->type(QueryNode::kXNode);
    q->keyX(&segment);
    QueryNode* leftNode = newQueryNode();
    QueryNode* rightNode = newQueryNode();
    q->lesserNode(leftNode);
    q->greaterNode(rightNode);
    //leftNode->parentNode(q);
    leftNode->addParentNode(q, newQueryNodeLink());
    //rightNode->parentNode(q);
    rightNode->addParentNode(q, newQueryNodeLink());
    leftNode->trapezoid(t);
    t->queryNode(leftNode);
    rightNode->trapezoid(rightHalf);
    rightHalf->queryNode(rightNode);
    assert(leftNode->type() == QueryNode::kLeafNode);
    assert(rightNode->type() == QueryNode::kLeafNode);
    
    // Step 1.  copy segments and top/bottom from old trap to new right half
    rightHalf->rightSegment(t->rightSegment());
    rightHalf->leftSegment(&segment);
    rightHalf->topVertex(t->topVertex());
    rightHalf->bottomVertex(t->bottomVertex());
    t->rightSegment(&segment);
    
    // Step 2.  pass pointers along to new right-hand trapezoid
    rightHalf->topLeft(t->topLeft());
    rightHalf->topRight(t->topRight());
    rightHalf->bottomLeft(t->bottomLeft());
    rightHalf->bottomRight(t->bottomRight());
    
    // Step 3.  specific new links.  There are several cases.
    
    // Switch: (1 of 2)
    // t has one neighbor below (implies that t is the lowest trap to split)
    //  is t, the new leftHalf, a triangle?
    //  is the new rightHalf a triangle?
    //  is neither one a triangle?
    // t has two neighbors below (may or may not be the bottom trap to split)
    //  t's bottom left neighbor split already (t might be the bottom trap)
    //  t's bottom right neighbor split already (t is not the bottom trap)
    // else t has one neighbor below that already split
    
    // Switch: (2 of 2)
    // t is not the top trapezoid to split
    //  split t's top neighbor (recursion)
    // else t is the top trapezoid to split
    //  is t, the new leftHalf, a triangle?
    //  is the new rightHalf a triangle?
    //  is neither one a triangle?
    
    if (oldBottomLeft == newBottomRight) // implies that t is the lowest trap!
    {
        // Test to see if rightHalf or t (the old left half) is a triangle.
        if (segment.sharesBottomEndpointWith(*originalLeftSegment))
        {
            // t is a triangle, and I think t == oldBottomLeft->topRight()
            assert(t == oldBottomLeft->topRight());
            t->bottomLeft(0L);
            oldBottomLeft->topRight(rightHalf);
        }
        else if (segment.sharesBottomEndpointWith(*originalRightSegment))
        {
            // rightHalf is a triangle, and so t == oldBottomLeft->topLeft()
            assert(t == oldBottomLeft->topLeft());
            rightHalf->bottomLeft(0L);
        }
        else // neither one is a triangle
        {
            oldBottomLeft->topLeft(t);
            oldBottomLeft->topRight(rightHalf);
        }
    }
    else if (t->bottomRight()) // i.e. old trapezoid had two bottom neighbors
    {
        assert(t->bottomLeft());
        if (oldBottomLeft == t->bottomLeft())
        {
            if (newBottomRight == t->bottomRight())
            {
                t->bottomRight(0L);
                rightHalf->bottomLeft(newBottomRight);
                rightHalf->bottomRight(0L);
                newBottomRight->topLeft(rightHalf);
            }
            else
            {
                rightHalf->bottomLeft(newBottomRight);
                rightHalf->bottomRight(t->bottomRight());
                t->bottomRight()->topLeft(rightHalf);
                t->bottomRight(0L);
                newBottomRight->topLeft(rightHalf);
            }
        }
        else if (oldBottomLeft == t->bottomRight())
        {
            // this case cannot possibly occur when splitting the lowest trap.
            // when t already had two bottom neighbors, and t->bottomRight
            // splits horizontally, t will appear to have three lower neighbors:
            // t->bottomLeft(), t->bottomRight() == oldBottomLeft, and
            // newBottomRight.  Connect rightHalf to newBottomRight to fix this.
            rightHalf->bottomLeft(newBottomRight);
            newBottomRight->topLeft(rightHalf);
            rightHalf->bottomRight(0L);
        }
    }
    else // old trapezoid had one neighbor below, on the left, but it split
    {
        assert(t->bottomRight() == 0L);
        rightHalf->bottomLeft(newBottomRight);
        newBottomRight->topLeft(rightHalf);
    }
    
    // Step 4.  Possible recursion.  Check if the top of the splitting segment
    // stops here.  Use the isBelow predicate, which is a lexicographic
    // comparison (not just a y-coordinate comparison).  This guarantees that
    // for any two distinct points, one is "below" the other.
    if (isBelow( (*t->topVertex()), segment.topVertex() ))
    {
        if (t->topRight() != 0L) // Two top neighbors
        {
            // Pick the correct upper neighbor to split
            if (t->topRight()->bottomCrossesSegment(segment))
            {
                rightHalf->topLeft(t->topRight());
                rightHalf->topRight(0L);
                splitHorizontally(t->topRight(), t, rightHalf, segment);
            }
            else
            {
                assert(t->topLeft()->bottomCrossesSegment(segment));
                t->topRight()->bottomLeft(rightHalf);
                t->topRight(0L);
                splitHorizontally(t->topLeft(), t, rightHalf, segment);
            }
        }
        else
        {
            splitHorizontally(t->topLeft(), t, rightHalf, segment);
        }
    }
    else // signifies that t is the top trapezoid to split.
    {
        if (segment.sharesTopEndpointWith(*originalLeftSegment))
        {
            // t is a triangle
            assert(t->topRight() == 0L);
            t->topLeft()->bottomRight(rightHalf);
            t->topLeft(0L);
        }
        else if (segment.sharesTopEndpointWith(*originalRightSegment))
        {
            // rightHalf is a triangle
            assert(rightHalf->topRight() == 0L);
            rightHalf->topLeft(0L);
        }
        else if (t->topRight() == 0L)
        {
            // neither one is a triangle
            t->topLeft()->bottomRight(rightHalf);
        }
        else
        {
            // neither is a triangle, and there were two trapezoids above.
            // this will happen when segment shares a vertex with the segment
            // that splits the upper trapezoids.  (yuck.)
            t->topRight(0L);
            rightHalf->topLeft(rightHalf->topRight());
            rightHalf->topRight(0L);
            rightHalf->topLeft()->bottomLeft(rightHalf);
        }
    }
    
    assert(t->pointersAreValid());
    assert(rightHalf->pointersAreValid());
    assert(oldBottomLeft->pointersAreValid());
    assert(newBottomRight->pointersAreValid());
    
    return rightHalf;
}

double Triangulator::
signedArea(vector<Vector2d>::const_iterator iBegin,
    vector<Vector2d>::const_iterator iEnd) const
{
    double twiceArea = 0.0;
    vector<Vector2d>::const_iterator itr(iBegin);
    Vector2d vCurrent;
    Vector2d vLast = *iBegin;
    while (itr != iEnd)
    {
        vCurrent = *itr;
        twiceArea -= (vCurrent[0]-vLast[0])*(vCurrent[1]+vLast[1]);
        vLast = vCurrent;
        itr++;
    }
    vCurrent = *iBegin;
    twiceArea -= (vCurrent[0]-vLast[0])*(vCurrent[1]+vLast[1]);
    
    return 0.5*twiceArea;
}

void Triangulator::
printTrapezoids()
{
    vector<Trapezoid>::const_iterator itr;
    for (itr = trapezoidsBegin(); itr != trapezoidsEnd(); itr++)
    if (mTrapezoidUsed[&(*itr) - &(mTrapezoidStore[0])])
        print(*itr);
}

bool Triangulator::
checkPointers()
{
    vector<Trapezoid>::const_iterator itr;
    for (itr = trapezoidsBegin(); itr != trapezoidsEnd(); itr++)
    if (mTrapezoidUsed[&(*itr) - &(mTrapezoidStore[0])])
        assert(itr->pointersAreValid());
    return true;
}

void Triangulator::
print(const Trapezoid & trap) const
{
    int num = &trap - &(mTrapezoidStore[0]);
    cout << "Trapezoid " << num << " [";
    if (trap.leftSegment())
        cout << int(trap.leftSegment() - &(mSegmentStore[0]));
    else
        cout << "unbounded";
    
    cout << ", ";
    if (trap.rightSegment())
        cout << int(trap.rightSegment() - &(mSegmentStore[0]));
    else
        cout << "unbounded";
    
    cout << "] ";
    if (trap.bottomVertex())
        cout << *(trap.bottomVertex()) << " to ";
    else
        cout << "(unbounded) to ";
    if (trap.topVertex())
        cout << *(trap.topVertex());
    else
        cout << "(unbounded)";
    
    cout << " bottom (";
    if (trap.bottomLeft())
        cout << int(trap.bottomLeft() - &(mTrapezoidStore[0]));
    else
        cout << "NULL";
    cout << ", ";
    if (trap.bottomRight())
        cout << int(trap.bottomRight() - &(mTrapezoidStore[0]));
    else
        cout << "NULL";
    cout << ") top (";
    if (trap.topLeft())
        cout << int(trap.topLeft() - &(mTrapezoidStore[0]));
    else
        cout << "NULL";
    cout << ", ";
    if (trap.topRight())
        cout << int(trap.topRight() - &(mTrapezoidStore[0]));
    else
        cout << "NULL";
    cout << ")\n";
}


void Triangulator::
print(const Segment & seg) const
{
    std::cout << "[" << seg.v0() << ", " << seg.v1() << "]";
}

bool Triangulator::
validateQueryTree() const
{
    validateQueryNode(&(mQueryNodeStore[0]));
    return true;
}

void Triangulator::
validateQueryNode(const QueryNode* q) const
{
    assert(q);
    if (q->type() == QueryNode::kLeafNode)
    {
        int trapNum = q->trapezoid() - &(mTrapezoidStore[0]);
        if (mTrapezoidUsed.at(trapNum) == false)
            print(*q->trapezoid());
        assert(q->trapezoid()->queryNode() == q);
        assert(mTrapezoidUsed.at(trapNum));
    }
    else //if (q->type() == QueryNode::kYNode)
    {
        assert(q->lesserNode());
        assert(q->greaterNode());
        assert(q->lesserNode() != q->greaterNode());
        validateQueryNode(q->lesserNode());
        validateQueryNode(q->greaterNode());
    }
//    else if (q->type() == QueryNode::kXNode)
//    {
//        validateQueryNode(q->lesserNode());
//        validateQueryNode(q->greaterNode());
//    }
}


#pragma mark *** Trapezoid ***


bool Trapezoid::
pointersAreValid() const
{
    if (topLeft())
    {
        if (this != topLeft()->bottomLeft() &&
            this != topLeft()->bottomRight())
        {
            return 0;
        }
        if (topLeft() == topRight())
            return 0;
    }
    if (topRight())
    {
        if (this != topRight()->bottomLeft() &&
            this != topRight()->bottomLeft())
        {
            return 0;
        }
        if (topRight() == topLeft())
            return 0;
    }
    if (bottomLeft())
    {
        if (this != bottomLeft()->topLeft() &&
            this != bottomLeft()->topRight())
        {
            return 0;
        }
        if (bottomLeft() == bottomRight())
            return 0;
    }
    if (bottomRight())
    {
        if (this != bottomRight()->topLeft() &&
            this != bottomRight()->topRight())
        {
            return 0;
        }
        if (bottomRight() == bottomLeft())
            return 0;
    }
    
    return 1;
}


// Definition of inclusion:
//
// on the bottom boundary: inside
// on the top boundary: outside
// on the left boundary: inside
// on the right boundary: outside
bool Trapezoid::
encloses(const Vector2d & v) const
{
    /*
        What we've got: the bounding segments, the top vertex, the bottom
        vertex.
        
        First cull by top/bottom.
    */
    if (mBottomVertex && Triangulator::isBelow(v, *mBottomVertex))
        return false;
    if (mTopVertex && !Triangulator::isBelow(v, *mTopVertex))
        return false;
    
    if (isBoundedLeft())
    {
        assert(mTopVertex);
        assert(mBottomVertex);
        if ( (*mTopVertex)[1] == (*mBottomVertex)[1])
        {
            // If it's a degenerate trapezoid, then we know that the vertex is
            // in the half-open segment (mBottomVertex, mTopVertex].  If the
            // top vertex is on the left hand segment, then v is on the left
            // boundary and is included.  Otherwise, v is in the degenerate
            // interior region, and is ALSO included.  (-:
            return true;
        }
        else if (v[0] < mLeft->x(v[1])) // nondegenerate trap: safe to call x() 
            return false;
        else if (isBoundedRight() && v[0] >= mRight->x(v[1]))
            return false;
    }
    else if (isBoundedRight())
    {
        assert(mTopVertex);
        assert(mBottomVertex);
        if ( (*mTopVertex)[1] == (*mBottomVertex)[1])
        {
            // This is a degenerate trapezoid; we know that the vertex is in the
            // half-open segment (mBottomVertex, mTopVertex].  The top vertex
            // is a point on mRight.  If the bottom vertex is also on mRight,
            // then v lies on the right boundary and is excluded.  Otherwise,
            // v is somewhere in the degenerate "interior" region.
            
            return &(rightSegment()->bottomVertex()) != mBottomVertex;
        }
        else if (v[0] >= mRight->x(v[1]))
            return false;
    }
    
    // If the trapezoid is unbounded (no mRight or mLeft) it encompasses the
    // entire plane, so the point is contained in the trapezoid.
    
    return true;
}

}; // namespace Triangulate




