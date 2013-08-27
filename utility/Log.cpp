/*
 *  Log.cpp
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 12/10/07.
 *  Copyright 2007 Stanford University. All rights reserved.
 *
 *  $Rev:: 185                           $:  Revision of last commit
 *  $Author:: pch                        $:  Author of last commit
 *
 *  $Date: 2009-10-20 10:49:03 -0700 (Tue, 20 Oct 2009) $:
 *  $Id: Log.cpp 185 2009-10-20 17:49:03Z pch $:
 *
 */

#include "Log.h"

std::ofstream TrogLog::sLogfile("log.out.txt");
StreamTee TrogLog::sTee(std::cout, TrogLog::sLogfile);

std::string TrogLog::
stripArgs(const std::string & str)
{
    int c2 = str.find_first_of("(");  // find the start of the arguments block
    std::string noArgs = str.substr(0,c2);
    int c1 = noArgs.find_last_of(" ");  // find the beginning of the invocation
    if (c1 == std::string::npos)
        c1 = -1;
    return str.substr(c1+1, c2-c1) + ")";
}
