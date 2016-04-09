//
//  Boolean-inl.h
//  NefLab
//
//  Created by Paul Hansen on 4/9/16.
//
//

#ifndef Boolean_inl_h
#define Boolean_inl_h

NefPolyhedron doIntersection(const NefPolyhedron & p1, const NefPolyhedron & p2);
NefPolyhedron doUnion(const NefPolyhedron & p1, const NefPolyhedron & p2);
NefPolyhedron doDifference(const NefPolyhedron & p1, const NefPolyhedron & p2);
//bool testIntersection(const NefPolyhedron & p1, const NefPolyhedron & p2);


NefPolyhedron doIntersection(const NefPolyhedron & p1, const NefPolyhedron & p2)
{
    std::cerr << "intersection\n";
    
    NefPolyhedron nefIntersection = (p1*p2).regularization();
    return nefIntersection;
}

NefPolyhedron doUnion(const NefPolyhedron & p1, const NefPolyhedron & p2)
{
    std::cerr << "union\n";
    
//    cout << "P1:\n";
//    writeNefPolyhedron(p1);
//    cout << "P1:\n";
//    writeNefPolyhedron(p1);
//    cout << "Union:\n";
    
    NefPolyhedron nefUnion = (p1+p2).regularization();
    return nefUnion;
}

NefPolyhedron doDifference(const NefPolyhedron & p1, const NefPolyhedron & p2)
{
    std::cerr << "difference\n";
    
    NefPolyhedron nefDifference = (p1 - p2).regularization();
    return nefDifference;
}




#endif /* Boolean_inl_h */
