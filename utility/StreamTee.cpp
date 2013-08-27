/*
 *  StreamTee.cpp
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 12/5/07.
 *  Copyright 2007 Stanford University. All rights reserved.
 *
 *  $Rev:: 35                            $:  Revision of last commit
 *  $Author:: pch                        $:  Author of last commit
 *
 *  $Date: 2009-05-13 15:11:45 -0700 (Wed, 13 May 2009) $:
 *  $Id: StreamTee.cpp 35 2009-05-13 22:11:45Z pch $:
 *
 */

#include "StreamTee.h"

StreamTee::StreamTee() :
	m_numStreams(0),
	m_str1(0L),
	m_str2(0L)
{
}

StreamTee::StreamTee(std::ostream & str) :
	m_numStreams(1),
	m_str1(&str),
	m_str2(&str)
{
}

StreamTee::StreamTee(std::ostream & str1, std::ostream & str2) :
	m_numStreams(2),
	m_str1(&str1),
	m_str2(&str2)
{
}


void
StreamTee::setStream(std::ostream & str)
{
	m_str1 = &str;
	m_numStreams = 1;
}

void
StreamTee::setStreams(std::ostream & str1, std::ostream & str2)
{
	m_str1 = &str1;
	m_str2 = &str2;
	m_numStreams = 2;
}

StreamTee&
StreamTee::operator<<(std::ostream& (*x)(std::ostream&))
{
	if (m_numStreams > 0)
		*m_str1 << x;
	if (m_numStreams > 1)
		*m_str2 << x;
	
	return *this;
}