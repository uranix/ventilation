cmake_minimum_required(VERSION 3.0)

project(Ventilation)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(GetGitRevisionDescription)
git_describe(GIT_VERSION --tags)
message(STATUS "Configured version is ${GIT_VERSION}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GitVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/GitVersion.h @ONLY)

set(vent_SOURCES src/flux.cpp src/contact.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=5 -std=c++11 -Wall -Werror -Wextra")

if(NOT DEFINED SECOND_ORDER)
	set(SECOND_ORDER FALSE)
endif()

if(NOT DEFINED RECONSTRUCT)
	set(RECONSTRUCT "${SECOND_ORDER}")
endif()

message(STATUS "Use output reconstruction: ${RECONSTRUCT}")
message(STATUS "Use second order method: ${SECOND_ORDER}")

if (${RECONSTRUCT})
	add_definitions("-DRECONSTRUCTION=1")
else()
	add_definitions("-DRECONSTRUCTION=0")
endif()

if (${SECOND_ORDER})
	add_definitions("-DSECOND_ORDER=1")
else()
	add_definitions("-DSECOND_ORDER=0")
endif()

find_package(SWIG 3.0 REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include_directories(include/)
add_library(vent SHARED ${vent_SOURCES})

SET_SOURCE_FILES_PROPERTIES(vent.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(vent.i PROPERTIES SWIG_FLAGS "-Wall")
SWIG_ADD_MODULE(vent python vent.i)
SWIG_LINK_LIBRARIES(vent ${PYTHON_LIBRARIES} vent)

add_executable(main main.cpp)
target_link_libraries(main vent)

# add_executable(col col.cpp)
# target_link_libraries(col vent)
