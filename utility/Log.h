/*
 *  Log.h
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 12/5/07.
 *
 */

#ifndef _LOGGER_
#define _LOGGER_

#include "StreamTee.h"

#include <boost/current_function.hpp>
#include "TimeWrapper.h"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <libgen.h>

//std::ofstream sLogfile("log.out.txt");
//StreamTee sLogTee(std::cout, sLogfile);

// Due to the static variables in the TrogLog class, it's possible to get into
// static initialization purgatory if any static members in the project get
// initialized before TrogLog's members AND then make use of the logfile
// macros.

class TrogLog
{
public:
	static StreamTee & tee() { return sTee; }
	static std::ofstream & logf() { return sLogfile; }
    static std::string stripArgs(const std::string & str);
private:
	static std::ofstream sLogfile;
	static StreamTee sTee;
	TrogLog() {}
	~TrogLog() {}
};


#define LOGF (TrogLog::logf() << timestamp() << " " << basename(const_cast<char*>(__FILE__)) << " " << std::setw(4) << __LINE__ << ": ")
#define LOGFMORE TrogLog::logf()
#define LOG (TrogLog::tee() << timestamp() << " " << basename(const_cast<char*>(__FILE__)) << " " << std::setw(4) << __LINE__ << ": ")
#define LOGMORE TrogLog::tee()

//#define LOGF (TrogLog::logf() << timestamp() << " [" << TrogLog::stripArgs((BOOST_CURRENT_FUNCTION)) << ", " << __LINE__ << "]: ")
//#define LOGFMORE TrogLog::logf()
//#define LOG (TrogLog::tee() << timestamp() << " [" << TrogLog::stripArgs((BOOST_CURRENT_FUNCTION)) << ", " << __LINE__ << "]: ")
//#define LOGMORE TrogLog::tee()

#endif

