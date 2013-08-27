

add_library( utility
calc.hh
geometry.h
geometry-inl.h
geometry2.h
geometry2-inl.h
Log.cpp
Log.h
Pointer.h
Polynomial.h
Polynomial-inl.h
RationalFunction.h
RationalFunction-inl.h
RLEOperations.h
RLEOperations-inl.h
STLOutput.h
StreamFromString.h
StreamTee.cpp
StreamTee.h
TimeWrapper.cpp
TimeWrapper.h
VectorMatrix.h
VectorMatrix-inl.h
VectorMatrix2.h
VectorMatrix2-inl.h
XMLExtras.h
WriteSomeRLE.h
)

find_package( Boost COMPONENTS date_time )

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message("\nBoost includes are at ${Boost_INCLUDE_DIRS}")
    #message("\nBoost libraries are ${Boost_LIBRARIES}")
    target_link_libraries(utility ${Boost_LIBRARIES})
else (Boost_FOUND)
    message(FATAL_ERROR "\nCan't find Boost!")
endif (Boost_FOUND)

add_subdirectory(tests)