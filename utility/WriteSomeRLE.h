/*
 *  WriteSomeRLE.h
 *  Snapdragon
 *
 *  Created by Paul Hansen on 4/26/10.
 *  Copyright 2010 Stanford University. All rights reserved.
 *
 */

#ifndef _WRITESOMERLE_
#define _WRITESOMERLE_

#include <string>
#include <fstream>
#include <sstream>
#include "geometry.h"
#include "rle/DynamicRLE3.h"

inline void writeRLE(const RLE::DynamicRLE3<double> & rle, Rect3i voxels,
    std::string filePrefix)
{
    std::ostringstream str;
    str << "_" << voxels.num(0) << "_" << voxels.num(1) << "_" << voxels.num(2);
    std::string fileName = filePrefix + str.str();
    std::ofstream file(fileName.c_str(), std::ofstream::binary);
    
    for (int kk = voxels.p1[2]; kk <= voxels.p2[2]; kk++)
    for (int jj = voxels.p1[1]; jj <= voxels.p2[1]; jj++)
    for (int ii = voxels.p1[0]; ii <= voxels.p2[0]; ii++)
    {
        double val = rle.at(ii,jj,kk);
        file.write((char*)&val, sizeof(double));
    }
    
    file.close();
}

inline void writeRLE(const RLE::DynamicRLE3<float> & rle, Rect3i voxels,
    std::string filePrefix)
{
    std::ostringstream str;
    str << "_" << voxels.num(0) << "_" << voxels.num(1) << "_" << voxels.num(2);
    std::string fileName = filePrefix + str.str();
    std::ofstream file(fileName.c_str(), std::ofstream::binary);
    
    for (int kk = voxels.p1[2]; kk <= voxels.p2[2]; kk++)
    for (int jj = voxels.p1[1]; jj <= voxels.p2[1]; jj++)
    for (int ii = voxels.p1[0]; ii <= voxels.p2[0]; ii++)
    {
        float val = rle.at(ii,jj,kk);
        file.write((char*)&val, sizeof(float));
    }
    
    file.close();
}

inline void writeRLE(const RLE::DynamicRLE3<Vector3d> & rle, Rect3i voxels,
    std::string filePrefix)
{
    std::ostringstream str;
    str << "_" << voxels.num(0) << "_" << voxels.num(1) << "_" << voxels.num(2);
    std::string fileName = filePrefix + str.str();
    std::ofstream file(fileName.c_str(), std::ofstream::binary);
    
    for (int kk = voxels.p1[2]; kk <= voxels.p2[2]; kk++)
    for (int jj = voxels.p1[1]; jj <= voxels.p2[1]; jj++)
    for (int ii = voxels.p1[0]; ii <= voxels.p2[0]; ii++)
    {
        Vector3d val = rle.at(ii,jj,kk);
        double val2[] = { val[0], val[1], val[2] };
        file.write((char*)val2, 3*sizeof(double));
    }
    
    file.close();
}

inline void writeRLE(const RLE::DynamicRLE3<Matrix3d> & rle, Rect3i voxels,
    std::string filePrefix)
{
    std::ostringstream str;
    str << "_" << voxels.num(0) << "_" << voxels.num(1) << "_" << voxels.num(2);
    std::string fileName = filePrefix + str.str();
    std::ofstream file(fileName.c_str(), std::ofstream::binary);
    
    for (int kk = voxels.p1[2]; kk <= voxels.p2[2]; kk++)
    for (int jj = voxels.p1[1]; jj <= voxels.p2[1]; jj++)
    for (int ii = voxels.p1[0]; ii <= voxels.p2[0]; ii++)
    {
        Matrix3d val = rle.at(ii,jj,kk);
        double val2[] = { val[0], val[1], val[2], val[3], val[4], val[5],
            val[6], val[7], val[8] };
        file.write((char*)val2, 9*sizeof(double));
    }
    
    file.close();
}


template<class T>
void writeOne(const T & val, std::ostream & stream)
{
    stream << val;
}

inline void writeOne(const Vector3d & val, std::ostream & stream)
{
    stream << "[" << val[0] << " " << val[1] << " " << val[2] << "]";
}

inline void writeOne(const Matrix3d & val, std::ostream & stream)
{
    writeOne(diagonalPartOf(val), stream);
}

template<class T>
inline void matlabRLE(const RLE::DynamicRLE3<T> & rle, Rect3i voxels,
    std::string variableName, std::ostream & outStream)
{
    outStream << variableName << ".bounds = ["
        << voxels.p1[0] << " " << voxels.p1[1] << " " << voxels.p1[2] << " "
        << voxels.p2[0] << " " << voxels.p2[1] << " " << voxels.p2[2] << "];\n";
    
    outStream << variableName << ".indices = [...\n";
    
    typename RLE::DynamicRLE3<T>::ConstIterator itr;
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        long ii,jj,kk;
        rle.cartesianCoordinates(itr.runStart(), ii, jj, kk);
        
        outStream << "\t" << ii-voxels.p1[0]
            << " " << jj-voxels.p1[1]
            << " " << kk-voxels.p1[2]
            << " ";
        rle.cartesianCoordinates(itr.runEnd(), ii, jj, kk);
        outStream << ii << " " << jj << " " << kk << ";...\n";
    }
    outStream << "];\n";
    
    outStream << variableName << ".values = [...\n";
    for (itr = rle.begin(); itr != rle.end(); itr.nextMarkedRun())
    {
        outStream << "\t";
        writeOne(itr.mark(), outStream);
        outStream << ";\n";
    }
    outStream << "];\n";
}





#endif
