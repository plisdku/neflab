#ifndef _TRIANGULATOR_
#define _TRIANGULATOR_

#include <vector>
#include <stack>
#include <stdexcept>
#include "geometry2.h"
#include "geometry.h"

namespace Triangulate
{

class Trapezoid;
class QueryNode;
class QueryNode;
struct QueryNodeLink;
class Segment;

/*
    Known bug:
    In the event that the triangulation routine throws an exception, it does
    not clean up the flags (mTrapezoidUsed, mVertexInTrapezoidation,
    and mSegmentInTriangulation).  This will very likely cause problems if
    the caller recovers after the exception and tries to triangulate another
    polygon.
*/

class Triangulator
{
public:
    Triangulator();
    
    // returns number of triangles
    int triangulate(
        const std::vector<int> & contourSizes,
        const std::vector<Vector3d> & vertices,
        Vector3d surfaceNormal,
        std::vector<Vector3d*> & outTriangleIndices);
    
    // returns number of triangles
    int triangulate(
        const std::vector<int> & contourSizes,
        const std::vector<Vector2d> & vertices,
        std::vector<Vector2d*> & outTriangleIndices);
    
private:
    static const double kINFINITY;
    
    void countPrimitives(const std::vector<int> & contourSizes,
        int & outNumVertices, int & outNumTriangles) const;
    
    void initSegments(const std::vector<int> & contourSizes,
        const std::vector<Vector2d> & vertices);
    
    void addSegmentToTrapezoidation(const Segment & segment);
    bool isVertexInTrapezoidation(const Vector2d* vertex) const;
    void flagVertexInTrapezoidation(const Vector2d* vertex);
    
    //void addDiagonals();
    void traverseChains(std::vector<Vector2d*> & outTriangleIndices);
    int triangulateChain(int chainLength, int triangleNum,
        std::vector<Vector2d*> & outTriangleIndices); // returns new triangleNum
    
    void plotSegments();
    void plotTrapezoids();
    
    // Find the trapezoid that encloses the given vertex.  If more than one
    // trapezoid is found, an assertion fails (so all trapezoids are tested
    // against the vertex).  This method operates by VERY safe brute force:
    // it searches every trapezoid, and will also report to cerr if more than
    // one trapezoid is found that encloses the vertex (this should not happen).
    Trapezoid* enclosingTrapezoid_safe(const Vector2d & vertex);
    
    // Quick search using the query tree.
    Trapezoid* enclosingTrapezoid(const Vector2d & vertex);
    
    // Find the lowest trapezoid that is crossed by this segment.
    Trapezoid* lowestIntersectingTrapezoid_safe(const Segment & segment);
    Trapezoid* lowestIntersectingTrapezoid(const Segment & segment);
public:
    // lexicographic comparison
    inline static bool isBelow(const Vector2d & lhs, const Vector2d & rhs);
        
    // This comparison is worth thinking about.  In the case of a horizontal
    // segment, there is NO LEGAL INPUT that permits the point lhs to lie in
    // the relative interior of segment rhs.  This function assumes that the
    // input is legal in this sense!
    inline static bool isToLeft(const Vector2d & lhs, const Segment & rhs);
    
private:
    
    Trapezoid* splitVertically(Trapezoid* t, const Vector2d & v);
    Trapezoid* splitHorizontally(Trapezoid* t, Trapezoid* oldBottomLeft,
        Trapezoid* newBottomRight, const Segment & segment);
    
    double signedArea(std::vector<Vector2d>::const_iterator iBegin,
        std::vector<Vector2d>::const_iterator iEnd) const;
    
    void printTrapezoids();
    bool checkPointers();
    void print(const Trapezoid & trap) const;
    void print(const Segment & seg) const;
    
    inline Segment* newSegment();
    inline QueryNode* newQueryNode();
    inline QueryNodeLink* newQueryNodeLink();
    inline Trapezoid* newTrapezoid();
    inline void deleteTrapezoid(Trapezoid* t);
    bool validateQueryTree() const;
    void validateQueryNode(const QueryNode* q) const;
    
    std::vector<Segment>::const_iterator segmentsBegin() const
        { return mSegmentStore.begin(); }
    std::vector<Segment>::const_iterator segmentsEnd() const
        { return mSegmentStore.begin() + mNumSegments; }
        
    std::vector<Segment>::iterator segmentsBegin()
        { return mSegmentStore.begin(); }
    std::vector<Segment>::iterator segmentsEnd()
        { return mSegmentStore.begin() + mNumSegments; }
    
    std::vector<Trapezoid>::const_iterator trapezoidsBegin() const
        { return mTrapezoidStore.begin(); }
    std::vector<Trapezoid>::const_iterator trapezoidsEnd() const
        { return mTrapezoidStore.begin() + mNumTrapezoids; }
    
    std::vector<Trapezoid>::iterator trapezoidsBegin()
        { return mTrapezoidStore.begin(); }
    std::vector<Trapezoid>::iterator trapezoidsEnd()
        { return mTrapezoidStore.begin() + mNumTrapezoids; }
    
    const Vector2d* mFirstUserVertex;
    
    // mNumSegments is the number of segments that are used.  It will be less
    // than mSegmentStore.size().
    int mNumSegments;
    int mNumTrapezoids;
    int mNumQueryNodes;
    int mNumQueryNodeLinks;
    std::vector<Segment> mSegmentStore;
    std::vector<Trapezoid> mTrapezoidStore;
    std::vector<QueryNode> mQueryNodeStore;
    std::vector<QueryNodeLink> mQueryNodeLinkStore;
    
    std::vector<Segment*> mNextSegment; // mNextSegment[vertexIndex] = ...
    std::vector<const Trapezoid*> mStartingTrapezoid; // mStartingTrapezoid[segNum]
    std::vector<const Vector2d*> mVertexChain;
    std::vector<int> mVertexIndexStack;
    
    // Markers!
    std::vector<char> mTrapezoidUsed;
    std::vector<char> mVertexInTrapezoidation;
    std::vector<char> mSegmentInTriangulation;
};

class Segment
{
public:
    Segment() : m_v0(0L), m_v1(0L) {}
    Segment(const Vector2d* v0, const Vector2d* v1) :
        m_v0(v0), m_v1(v1)
    {
    }
    
    const Vector2d & v0() const { assert(m_v0 != 0L); return *m_v0; }
    const Vector2d & v1() const { assert(m_v1 != 0L); return *m_v1; }
    
    double height() const { return fabs( (*m_v0)[1] - (*m_v1)[1]); }
    double width() const { return fabs( (*m_v0)[0] - (*m_v1)[0]); }
    double left() const { return fmin( (*m_v0)[0], (*m_v1)[0] ); }
    double right() const { return fmax( (*m_v0)[0], (*m_v1)[0] ); }
    
    inline const Vector2d & topVertex() const;
    inline const Vector2d & bottomVertex() const;
    
    void v0(const Vector2d* v0) { m_v0 = v0; }
    void v1(const Vector2d* v1) { m_v1 = v1; }
    
    // The segment lies on a supporting line.  At what x is the supporting
    // line crossed by a line of constant y?
    inline double x(double y) const;
    inline bool sharesEndpointWith(const Segment & rhs) const;
    inline bool sharesTopEndpointWith(const Segment & rhs) const;
    inline bool sharesBottomEndpointWith(const Segment & rhs) const;
    inline bool hasEndpoint(const Vector2d* v) const;
    
private:
    const Vector2d* m_v0;
    const Vector2d* m_v1;
};

struct QueryNodeLink
{
    QueryNodeLink() : node(0L), next(0L) {}
    QueryNode* node;
    QueryNodeLink* next;
};

class QueryNode
{
public:
    QueryNode() : mType(kLeafNode), mKeyX(0L), mKeyY(0L), mParents(0L),
        mLesser(0L), mGreater(0L), mTrapezoid(0L) {}
    
    enum Type
    {
        kXNode,
        kYNode,
        kLeafNode
    };
    
    Type type() const { return mType; }
    const Segment* keyX() const { return mKeyX; }
    const Vector2d* keyY() const { return mKeyY; }
    //QueryNode* parentNode() const { return mParent; }
    QueryNodeLink* parentListHead() { return mParents; }
    const QueryNodeLink* parentListHead() const { return mParents; }
    QueryNode* lesserNode() const { return mLesser; }
    QueryNode* greaterNode() const { return mGreater; }
    
    void type(Type t) { mType = t; }
    void keyX(const Segment* s) { mKeyX = s; }
    void keyY(const Vector2d* v) { mKeyY = v; }
    //void parentNode(QueryNode* q) { mParent = q; }
    inline void addParentNode(QueryNode* q, QueryNodeLink* linkNode);
    void lesserNode(QueryNode* q) { mLesser = q; }
    void greaterNode(QueryNode* q) { mGreater = q; }
    
    // Iterate through the parents of this node.  If this == parent->mLesser,
    // set parent->mLesser = q.  Same for mGreater.  Doesn't change
    // this->mParents at all.
    inline void replaceWith(QueryNode* q);
    
    
    Trapezoid* trapezoid() { return mTrapezoid; }
    const Trapezoid* trapezoid() const { return mTrapezoid; }
    
    void trapezoid(Trapezoid* t) { mTrapezoid = t; }
    
private:
    Type mType;
    const Segment* mKeyX;
    const Vector2d* mKeyY;
    //QueryNode* mParent;
    QueryNodeLink* mParents;
    QueryNode* mLesser;
    QueryNode* mGreater;
    Trapezoid* mTrapezoid;
};

class Trapezoid
{
public:
    Trapezoid() : mLeft(0L), mRight(0L),
        mTopLeft(0L), mTopRight(0L), mBottomLeft(0L), mBottomRight(0L),
        mBottomVertex(0L), mTopVertex(0L), mQueryNode(0L) {}
    
    Trapezoid(const Segment* left, const Segment* right,
        const Vector2d* bottomVertex, const Vector2d* topVertex) :
        mLeft(left), mRight(right),
        mTopLeft(0L), mTopRight(0L), mBottomLeft(0L), mBottomRight(0L),
        mBottomVertex(bottomVertex), mTopVertex(topVertex), mQueryNode(0L)
    {
        // used for testing only
    }
    bool pointersAreValid() const;
    
    inline bool isInside() const;    
    Trapezoid* topLeft() const { return mTopLeft; }
    Trapezoid* topRight() const { return mTopRight; }
    Trapezoid* bottomLeft() const { return mBottomLeft; }
    Trapezoid* bottomRight() const { return mBottomRight; }
    
    inline Trapezoid* upperRightNeighbor() const;
    inline Trapezoid* bottomRightNeighbor() const;
    
    void topLeft(Trapezoid* t) { mTopLeft = t; }
    void topRight(Trapezoid* t) { mTopRight = t; }
    void bottomLeft(Trapezoid* t) { mBottomLeft = t; }
    void bottomRight(Trapezoid* t) { mBottomRight = t; }
    
    inline int numTopNeighbors() const;    
    inline int numBottomNeighbors() const;
    const Segment* leftSegment() const { return mLeft; }
    const Segment* rightSegment() const { return mRight; }
    void leftSegment(const Segment* left) { mLeft = left; }
    void rightSegment(const Segment* right) { mRight = right; }
    QueryNode* queryNode() { return mQueryNode; }
    const QueryNode* queryNode() const { return mQueryNode; }
    void queryNode(QueryNode* queryNode) { mQueryNode = queryNode; }
    
    bool isBoundedLeft() const { return mLeft != 0L; }
    bool isBoundedRight() const { return mRight != 0L; }
    
    void topVertex(const Vector2d* newTop) { mTopVertex = newTop; }
    void bottomVertex(const Vector2d* newBottom) { mBottomVertex = newBottom; }
    const Vector2d* topVertex() const { return mTopVertex; }
    const Vector2d* bottomVertex() const { return mBottomVertex; }
    
    inline const Vector2d* oppositeVertex(const Vector2d* v) const;
    inline bool willBeCutDiagonally() const;
    bool encloses(const Vector2d & v) const;
    inline bool bottomCrossesSegment(const Segment & seg) const;
    inline bool isBelow(const Vector2d & vertex) const;
    
private:
    const Segment* mLeft;
    const Segment* mRight;
    
    // If the trapezoid has only one neighbor above (below), it is 
    // stored in mTopLeft (mBottomLeft), and mTopRight (mBottomRight) is
    // null.
    Trapezoid* mTopLeft;
    Trapezoid* mTopRight;
    Trapezoid* mBottomLeft;
    Trapezoid* mBottomRight;
    
    const Vector2d* mBottomVertex;
    const Vector2d* mTopVertex;
    
    QueryNode* mQueryNode;
};

}; // namespace Triangulate

#include "Triangulator-inl.h"

#endif
