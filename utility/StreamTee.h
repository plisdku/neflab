/*
 *  StreamTee.h
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 12/5/07.
 *  Copyright 2007 Stanford University. All rights reserved.
 *
 *  $Rev:: 35                            $:  Revision of last commit
 *  $Author:: pch                        $:  Author of last commit
 *
 *  $Date: 2009-05-13 15:11:45 -0700 (Wed, 13 May 2009) $:
 *  $Id: StreamTee.h 35 2009-05-13 22:11:45Z pch $:
 *
 */

#ifndef _STREAMTEE_
#define _STREAMTEE_

#include <iostream>

class StreamTee
{
public:
	StreamTee();
	StreamTee(std::ostream & str);
	StreamTee(std::ostream & str1, std::ostream & str2);
	
	void setStream(std::ostream & str);
	void setStreams(std::ostream & str1, std::ostream & str2);
	
	// Normal output
	template <class T>
	StreamTee & operator<<(T arg);
	
	// Special case: catch iomanips like endl and flush
	StreamTee & operator<<(std::ostream& (*x)(std::ostream&));
	
private:
	int m_numStreams;
	std::ostream * m_str1;
	std::ostream * m_str2;
};

template <class T>
StreamTee&
StreamTee::operator<<(T arg)
{
	if (m_numStreams > 0)
		*m_str1 << arg;
	if (m_numStreams > 1)
		*m_str2 << arg;
	
	return *this;
}


#endif
