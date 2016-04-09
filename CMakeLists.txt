project(NefLab)

add_subdirectory( utility )
#add_subdirectory( Triangulate )

# ------------------ Boilerplate

cmake_minimum_required(VERSION 2.4.5)
if (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

# ---- done with the sillies

#set(Boost_USE_MULTITHREADED ON)

find_package(CGAL QUIET COMPONENTS Core )
if (NOT CGAL_FOUND)
    message(FATAL_ERROR "Cannot find CGAL.")
endif (NOT CGAL_FOUND)

include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )


set(CGAL_3RD_PARTY_LIBRARIES "/usr/local/lib/libboost_thread-mt.dylib;/usr/local/lib/libboost_system-mt.dylib" )

message(ÒCGAL libraries are  ${CGAL_LIBRARIES}Ó)
message(ÒCGAL use file is ${CGAL_USE_FILE}Ó)
message(ÒCGAL 3rd party is ${CGAL_3RD_PARTY_LIBRARIES}Ó)

include_directories(
#    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

add_executable( NefLab
    Boolean-inl.h
    BuildMesh-inl.h
    CGALUtilities.h
    InteriorVolumes-inl.h
    main.cpp
    PaternityTest-inl.h
    PointFacetDistance-inl.h
    ShellVisitor.h
)

target_link_libraries( NefLab
    ${CGAL_LIBRARIES}
)

#add_subdirectory( tests )

