cmake_minimum_required(VERSION 3.0)

project(Ventilation)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(GetGitRevisionDescription)
git_describe(GIT_VERSION --tags)
message(STATUS "Configured version is ${GIT_VERSION}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GitVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/GitVersion.h @ONLY)

set(vent_SOURCES
	src/contact.cpp
	src/vtk.cpp
	src/object.cpp
	src/pipe.cpp
	src/solver.cpp
	)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -Wextra -g")
set(max_errors 5)
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=${max_errors}")
endif()
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ferror-limit=${max_errors}")
endif()

if(NOT DEFINED PRECISE_RIEMANN)
	set(PRECISE_RIEMANN FALSE)
endif()

if(NOT DEFINED FP_TRAP)
	set(FP_TRAP TRUE)
endif()

if(NOT DEFINED NC)
	set(NC 2)
endif()

message(STATUS "Configuerd with following params")
message(STATUS " ")
message(STATUS "Use precise Riemann solver (-DPRECISE_RIEMANN): ${PRECISE_RIEMANN}")
message(STATUS "Use floating point traps (-DFP_TRAP): ${FP_TRAP}")
message(STATUS "Number of components (-DNC): ${NC}")
message(STATUS " ")

if (${PRECISE_RIEMANN})
	add_definitions("-DPRECISE_RIEMANN=1")
	list(APPEND vent_SOURCES src/flux_riemann.cpp)
else()
	add_definitions("-DPRECISE_RIEMANN=0")
	list(APPEND vent_SOURCES src/flux_cir.cpp)
endif()

if (${FP_TRAP})
	add_definitions("-DFP_TRAP=1")
else()
	add_definitions("-DFP_TRAP=0")
endif()

add_definitions("-DNC=${NC}")

find_package(SWIG 3.0 REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include_directories(include/)
include_directories(include/eigen/)
add_library(vent SHARED ${vent_SOURCES})

SET_SOURCE_FILES_PROPERTIES(pyvent.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(pyvent.i PROPERTIES SWIG_FLAGS "-Wall")
SWIG_ADD_MODULE(pyvent python pyvent.i)
SWIG_LINK_LIBRARIES(pyvent ${PYTHON_LIBRARIES} vent)

#add_executable(main main.cpp)
#target_link_libraries(main vent)

# add_executable(col col.cpp)
# target_link_libraries(col vent)
