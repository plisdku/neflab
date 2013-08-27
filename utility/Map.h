/*
 *  Map.h
 *  SNAPDRAGON
 *
 *  Created by Paul Hansen on 6/14/06.
 *  Copyright 2007 Stanford University. All rights reserved.
 *
 *  $Rev:: 35                            $:  Revision of last commit
 *  $Author:: pch                        $:  Author of last commit
 *
 *  $Date: 2009-05-13 15:11:45 -0700 (Wed, 13 May 2009) $:
 *  $Id: Map.h 35 2009-05-13 22:11:45Z pch $:
 *
 *  This wraps the std::map class to provide a convenient operator:
 *
 *      const operator[] (Key) const
 *
 *  which bypasses the need for using find() and all that.  It crashes
 *  on failure to find the key, so caveat cryptor!
 *
 */

#ifndef _MAP_
#define _MAP_

#include <map>
#include <cassert>

template <typename Key, typename Value>
class Map : public std::map<Key, Value>
{
public:
    Map();
    virtual ~Map();
    
    Map(const std::map<Key, Value> & rhs);
    Map<Key, Value> & operator=(const std::map<Key, Value> & rhs);
    
    template<typename Key2>
    const Value& operator[] (const Key2& key) const;
    
    template<typename Key2>
    Value& operator[] (const Key2& key);
};


template <typename Key, typename Value>
Map<Key, Value> ::
Map()
{
}

template <typename Key, typename Value>
Map<Key, Value> ::
~Map()
{
}

template <typename Key, typename Value>
Map<Key, Value> ::
Map(const std::map<Key, Value> & rhs) :
    std::map<Key, Value>(rhs)
{
}

template <typename Key, typename Value>
Map<Key, Value> & Map<Key, Value>::
operator=(const std::map<Key, Value> & rhs)
{
    *this = rhs;
    return *this;
}

template <typename Key, typename Value>
template <typename Key2>
const Value& Map<Key, Value>::
operator[] (const Key2& key) const
{
    assert( count(Key(key)) != 0 );
    return (* find(Key(key))).second;
}

template <typename Key, typename Value>
template <typename Key2>
Value& Map<Key, Value>::
operator[] (const Key2& key)
{
    return std::map<Key, Value>::operator[](Key(key));
}





#endif

