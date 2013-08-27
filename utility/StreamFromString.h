/*
 *  StreamFromString.h
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 12/11/07.
 *  Copyright 2007 Stanford University. All rights reserved.
 *
 *  $Rev:: 35                            $:  Revision of last commit
 *  $Author:: pch                        $:  Author of last commit
 *
 *  $Date: 2009-05-13 15:11:45 -0700 (Wed, 13 May 2009) $:
 *  $Id: StreamFromString.h 35 2009-05-13 22:11:45Z pch $:
 *
 */

#ifndef _STREAMFROMSTRING_
#define _STREAMFROMSTRING_

#include <iostream>
#include <sstream>

template<typename T>
std::string operator >> (std::string inString, T & value)
{
	std::string remainder;
	std::istringstream istr(inString);
	istr >> value;
	getline(istr, remainder, '\0');
	return remainder;
}

#endif
