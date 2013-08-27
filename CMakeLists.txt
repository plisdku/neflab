project(NefLab)

add_subdirectory( utility )
#add_subdirectory( Triangulate )

# ------------------ Boilerplate

cmake_minimum_required(VERSION 2.4.5)
if (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

# ---- done with the sillies

find_package(CGAL QUIET COMPONENTS Core )
if (NOT CGAL_FOUND)
    message(FATAL_ERROR "Cannot find CGAL.")
endif (NOT CGAL_FOUND)

include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

include_directories(
#    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

add_executable( NefLab
    main.cpp
    InteriorVolumes-inl.h
    PointFacetDistance-inl.h
    ShellVisitor.h
    CGALUtilities.h
)

target_link_libraries( NefLab
    ${CGAL_LIBRARIES}
)

#add_subdirectory( tests )

