/*
 *  Triangulator-inl.h
 *  triangulate
 *
 *  Created by Paul Hansen on 2/13/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifdef _TRIANGULATOR_

namespace Triangulate
{

#pragma mark *** Triangulator ***

bool Triangulator::
isBelow(const Vector2d & lhs, const Vector2d & rhs)
{
    if (lhs[1] < rhs[1])
        return 1;
    else if (lhs[1] > rhs[1])
        return 0;
    return lhs[0] < rhs[0];
}

bool Triangulator::
isToLeft(const Vector2d & lhs, const Segment & rhs)
{
    if (rhs.height() > 0)
        return lhs[0] < rhs.x(lhs[1]);
    return lhs[0] < rhs.left();
}

Segment* Triangulator::
newSegment()
{
    assert(mNumSegments+1 < mSegmentStore.size());
    mSegmentStore[mNumSegments] = Segment();
    return &(mSegmentStore[mNumSegments++]);
}

QueryNode* Triangulator::
newQueryNode()
{
//    if (mNumQueryNodes+1 >= mQueryNodeStore.size())
//        throw(std::runtime_error("Query node store too small"));
    assert(mNumQueryNodes+1 < mQueryNodeStore.size());
    mQueryNodeStore[mNumQueryNodes] = QueryNode();
    return &(mQueryNodeStore[mNumQueryNodes++]);
}

QueryNodeLink* Triangulator::
newQueryNodeLink()
{
    //std::cout << "Link " << mNumQueryNodeLinks << "\n";
    assert(mNumQueryNodeLinks+1 < mQueryNodeLinkStore.size());
    mQueryNodeLinkStore[mNumQueryNodeLinks] = QueryNodeLink();
    return &(mQueryNodeLinkStore[mNumQueryNodeLinks++]);
}

Trapezoid* Triangulator::
newTrapezoid()
{
    assert(mNumTrapezoids+1 < mTrapezoidStore.size());
    mTrapezoidStore[mNumTrapezoids] = Trapezoid();
    mTrapezoidUsed[mNumTrapezoids] = true;
    return &(mTrapezoidStore[mNumTrapezoids++]);
}

void Triangulator::
deleteTrapezoid(Trapezoid* t)
{
    mTrapezoidUsed[int(t - &(mTrapezoidStore[0]))] = false;
}


#pragma mark *** Segment ***

const Vector2d & Segment::
topVertex() const
{
    assert(m_v0 != 0L);
    assert(m_v1 != 0L);
    if (Triangulator::isBelow(*m_v0, *m_v1))
        return *m_v1;
    return *m_v0;
}

const Vector2d & Segment::
bottomVertex() const
{
    assert(m_v0 != 0L);
    assert(m_v1 != 0L);
    if (Triangulator::isBelow(*m_v0, *m_v1))
        return *m_v0;
    return *m_v1;
}

double Segment::
x(double y) const
{
    assert((*m_v1)[1] != (*m_v0)[1]);
    return (*m_v0)[0] + (y-(*m_v0)[1]) * ((*m_v1)[0]-(*m_v0)[0])
        / ((*m_v1)[1] - (*m_v0)[1]);
}

bool Segment::
sharesEndpointWith(const Segment & rhs) const
{
    if (&rhs != 0L)
        return (m_v0 == rhs.m_v1 || m_v1 == rhs.m_v0);
    return 0;
}

bool Segment::
sharesTopEndpointWith(const Segment & rhs) const
{
    if (&rhs != 0L)
        return (topVertex() == rhs.topVertex());
    return 0;
}

bool Segment::
sharesBottomEndpointWith(const Segment & rhs) const
{
    if (&rhs != 0L)
        return (bottomVertex() == rhs.bottomVertex());
    return 0;
}

bool Segment::
hasEndpoint(const Vector2d* v) const
{
    return (m_v0 == v || m_v1 == v);
}


#pragma mark *** QueryNode ***

void QueryNode::
addParentNode(QueryNode* q, QueryNodeLink* linkNode)
{
    linkNode->node = q;
    linkNode->next = mParents;
    mParents = linkNode;
}

void QueryNode::
replaceWith(QueryNode* q)
{
    QueryNodeLink* parent = mParents;
    while(parent)
    {
        if (this == parent->node->lesserNode())
        {
            parent->node->lesserNode(q);
        }
        else
        {
            assert(this == parent->node->greaterNode());
            parent->node->greaterNode(q);
        }
        parent = parent->next;
    }
}

#pragma mark *** Trapezoid ***

bool Trapezoid::
isInside() const
{
    if (false == isBoundedLeft() || false == isBoundedRight())
        return false;
    
    const Vector2d& v0 = rightSegment()->v0();
    const Vector2d& v1 = rightSegment()->v1();
    return Triangulator::isBelow(v0, v1);
}

Trapezoid* Trapezoid::
upperRightNeighbor() const
{
    if (mTopRight)
        return mTopRight;
    return mTopLeft;
}

Trapezoid* Trapezoid::
bottomRightNeighbor() const
{
    if (mBottomRight)
        return mBottomRight;
    return mBottomLeft;
}

int Trapezoid::
numTopNeighbors() const
{
    if (mTopRight)
    {
        assert(mTopLeft);
        return 2;
    }
    else if (mTopLeft)
        return 1;
    return 0;
}

int Trapezoid::
numBottomNeighbors() const
{
    if (mBottomRight)
    {
        assert(mBottomLeft);
        return 2;
    }
    else if (mBottomLeft)
        return 1;
    return 0;
}


const Vector2d* Trapezoid::
oppositeVertex(const Vector2d* v) const
{
    if (v == mTopVertex)
        return mBottomVertex;
    else if (v == mBottomVertex)
        return mTopVertex;
    throw(std::logic_error("Vertex not in trapezoid"));
}

bool Trapezoid::
willBeCutDiagonally() const
{
//    if (numBottomNeighbors() == 0 || numTopNeighbors() == 0)
//        return false;
    
    if (&leftSegment()->v0() == topVertex() &&
        &leftSegment()->v1() == bottomVertex())
    {
        return false;
    }
    
    if (&rightSegment()->v0() == bottomVertex() &&
        &rightSegment()->v1() == topVertex())
    {
        return false;
    }
    
    return true;
}

bool Trapezoid::
bottomCrossesSegment(const Segment & seg) const
{
    if (mBottomVertex == 0L)
        return false;
    
    const Vector2d & vBot = seg.bottomVertex();
    const Vector2d & vTop = seg.topVertex();
    double height = vTop[1] - vBot[1];
    double width = vTop[0] - vBot[0];
    double x;
    if (height != 0)
        x = vBot[0] + ( (*mBottomVertex)[1] -vBot[1] )*width/height;
    else
        x = vBot[0]; // think this'll work?
    
    if (isBoundedLeft())
    {
        if (leftSegment()->height() != 0 &&
            leftSegment()->x((*mBottomVertex)[1]) > x)
            return 0;
        else if (leftSegment()->left() > x)
            return 0;
    }
    
    if (isBoundedRight())
    {
        if (rightSegment()->height() != 0 &&
            rightSegment()->x((*mBottomVertex)[1]) < x)
            return 0;
        else if (rightSegment()->right() < x)
            return 0;
    }
    return 1;
}

bool Trapezoid::
isBelow(const Vector2d & vertex) const
{
    assert(&vertex != 0L);
    if (mTopVertex)
        return !Triangulator::isBelow(vertex, *mTopVertex); // use "<="
    return 0;
}



}; // namespace Triangulate

#endif