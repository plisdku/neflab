/*
 *  TimeWrapper.h
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 7/12/05.
 *  Copyright 2005 Stanford University. All rights reserved.
 *
 *  $Rev:: 185                           $:  Revision of last commit
 *  $Author:: pch                        $:  Author of last commit
 *
 *  $Date: 2009-10-20 10:49:03 -0700 (Tue, 20 Oct 2009) $:
 *  $Id: TimeWrapper.h 185 2009-10-20 17:49:03Z pch $:
 *
 */

#ifndef _TIMEWRAPPER_
#define _TIMEWRAPPER_

#include <string>

//#include <sys/time.h>

double timeInMicroseconds();
//{
//    timeval tv;
//    struct timezone tz;
//    
//    gettimeofday(&tv, &tz);
//    
//    return ((double)tv.tv_usec + 1e6*tv.tv_sec);
//}

std::string timestamp();

#endif
