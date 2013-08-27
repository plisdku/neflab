add_library( Triangulate
	Triangulator.cpp
	Triangulator.h
	Triangulator-inl.h
)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/../utility
)

target_link_libraries(
	utility
)


