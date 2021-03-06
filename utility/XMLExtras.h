/*
 *  XMLExtras.h
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 5/12/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 */

#ifndef _XMLEXTRAS_
#define _XMLEXTRAS_

#include "tinyxml.h"
#include <string>
#include <sstream>
#include <stdexcept>
#include "Map.h"


#pragma mark *** Implementations ***


static std::string
sErr(const std::string & str, const TiXmlElement* elem)
{
	std::ostringstream outStr;
	outStr << str << " (row " << elem->Row() << ", col " << elem->Column()
		<< ")";
	return outStr.str();
}


static Map<std::string, std::string>
sGetAttributes(const TiXmlElement* elem)
{
    Map<std::string, std::string> attribs;
    if (elem)
    {
        const TiXmlAttribute* attrib = elem->FirstAttribute();
        while (attrib)
        {
            attribs[attrib->Name()] = attrib->Value();
            attrib = attrib->Next();
        }
    }
    return attribs;
}

static void sSetAttributes(TiXmlElement* elem,
    const Map<std::string, std::string> & attribs) throw(std::logic_error)
{
    if (elem)
    {
        std::map<std::string, std::string>::const_iterator itr;
        for (itr = attribs.begin(); itr != attribs.end(); itr++)
            elem->SetAttribute(itr->first, itr->second);\
    }
}

template <class T>
void sGetMandatoryAttribute(const TiXmlElement* elem,
	const std::string & attribute, T & val) throw(std::logic_error)
{
	assert(elem);
	std::ostringstream err;
	
	if (elem->Attribute(attribute.c_str()) == 0L)
	{
		err << elem->Value() << " needs attribute \"" << attribute << "\"";
		throw(std::logic_error(sErr(err.str(), elem)));
	}
	
	std::istringstream istr(elem->Attribute(attribute.c_str()));
	
	if (!(istr >> val))
	{
//        assert(!"Trouble.");
		err << elem->Value() << " has invalid \"" << attribute
			<<"\" attribute";
		throw(std::logic_error(sErr(err.str(), elem)));
	}
}

template <class T>
bool sTryGetAttribute(const TiXmlElement* elem,
	const std::string & attribute, T & val)
{
	assert(elem);
	std::ostringstream err;
	
	if (elem->Attribute(attribute.c_str()) == 0L)
		return 0;
	else
	{
		std::istringstream istr(elem->Attribute(attribute.c_str()));
		
		if (!(istr >> val))
		{
            assert(!"Trouble.");
			err << elem->Value() << " has invalid \"" << attribute
				<< "\" attribute";
			throw(std::logic_error(sErr(err.str(), elem)));
		}
	}
	return 1;
}

template <class T>
void sGetOptionalAttribute(const TiXmlElement* elem,
	const std::string & attribute, T & val, const T & defaultVal)
{
	assert(elem);
	std::ostringstream err;
	
	if (elem->Attribute(attribute.c_str()) == 0L)
	{
		val = defaultVal;
	}
	else
	{
		std::istringstream istr(elem->Attribute(attribute.c_str()));
		
		if (!(istr >> val))
		{
            assert(!"Trouble.");
			err << elem->Value() << " has invalid \"" << attribute
				<< "\" attribute";
			throw(std::logic_error(sErr(err.str(), elem)));
		}
	}
}


template<class T>
void sSetAttribute(TiXmlElement* elem,
    const std::string & attribute, T & val) throw(std::logic_error)
{
    std::istringstream istr;
    istr << val;
    if (elem)
        elem->SetAttribute(attribute, istr.str().c_str());
}

//template <>
//void sGetMandatoryAttribute<std::string>(const TiXmlElement* elem,
//	const std::string & attribute, std::string & val) throw(std::logic_error)
static void sGetMandatoryAttribute(const TiXmlElement* elem,
	const std::string & attribute, std::string & val) throw(std::logic_error)
{
	assert(elem);
	std::ostringstream err;
	
	if (elem->Attribute(attribute.c_str()) == 0L)
	{
		err << elem->Value() << " needs attribute \"" << attribute << "\"";
		throw(std::logic_error(sErr(err.str(), elem)));
	}
	val = elem->Attribute(attribute.c_str());
}

/*
static bool sTryGetAttribute(const TiXmlElement* elem, 
    const std::string & attribute, std::string & val)
{
	assert(elem);
	std::ostringstream err;
	
	if (elem->Attribute(attribute.c_str()) == 0L)
		return 0;
	else
	{
		val = elem->Attribute(attribute.c_str());
	}
	return 1;
}
*/

/*
static void sGetOptionalAttribute(const TiXmlElement* elem,
    const std::string & attribute, std::string & val,
    const std::string & defaultVal)
{
	assert(elem);
	std::ostringstream err;
	
	if (elem->Attribute(attribute.c_str()) == 0L)
	{
		val = defaultVal;
	}
	else
	{
		val = elem->Attribute(attribute.c_str());
	}
}
*/



#endif
