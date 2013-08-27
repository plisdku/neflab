
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

set(Boost_ADDITIONAL_VERSIONS "1.41" "1.41.0")
find_package( Boost COMPONENTS unit_test_framework REQUIRED  )
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(
    /usr/local/include
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

add_executable(testGeometry
    testGeometry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../geometry.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../geometry-inl.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../geometry2.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../geometry2-inl.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../VectorMatrix.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../VectorMatrix-inl.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../VectorMatrix2.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../VectorMatrix2-inl.h
)

add_executable(testPolynomial
    testPolynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Polynomial.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../Polynomial-inl.h
)

add_executable(testRationalFunction
    testRationalFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Polynomial.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../Polynomial-inl.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../RationalFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../RationalFunction-inl.h
)

add_executable(testRLEOperations
    testRLEOperations.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../RLEOperations.h
    ${CMAKE_CURRENT_SOURCE_DIR}/../RLEOperations-inl.h
)
target_link_libraries(testRLEOperations utility)

macro(setupUnitTest executableName)
    target_link_libraries(${executableName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    get_target_property( path ${executableName} LOCATION)
    add_custom_command( TARGET ${executableName} POST_BUILD
        COMMAND ${path} --log_level=message)
endmacro()

setupUnitTest(testGeometry)
setupUnitTest(testPolynomial)
setupUnitTest(testRationalFunction)
setupUnitTest(testRLEOperations)
