cmake_minimum_required(VERSION 2.4 FATAL_ERROR)

project(DT_Reader)

set(CMAKE_BUILD_TYPE debug)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{DATATREE_HOME})

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(DataTree REQUIRED)

set(CMAKE_CXX_FLAGS "-Wall")

include(${ROOT_USE_FILE})

include_directories(
    $ENV{HADDIR}/include
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${DataTree_INCLUDE_DIR}
    )
		
link_directories(
    ${CMAKE_SOURCE_DIR}
    ${DataTree_LIBRARY_DIR}
    )


set (HEADERS
		src/Centrality.h
		src/Selector.h
		src/Qvector.h
		src/Qvector3SE.h
		src/TrackQA.h
		src/EventQA.h
		src/Flow.h
		src/Flow3SE.h
		src/Reader.h
	)

set (SOURCES 
		src/Centrality.cxx
		src/Selector.cxx
		src/Qvector.cxx
		src/Qvector3SE.cxx
		src/TrackQA.cxx
		src/EventQA.cxx
		src/Flow.cxx
		src/Flow3SE.cxx
		src/Reader.cxx
		src/main.cxx
	)

add_definitions(${ROOT_CXX_FLAGS})

add_executable(DT_Reader ${HEADERS} ${SOURCES})
target_link_libraries(DT_Reader ${ROOT_LIBRARIES} DataTree)
