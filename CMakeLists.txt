cmake_minimum_required(VERSION 3.14)

project(Ventilation)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

string(TIMESTAMP _date "%d.%m.%Y, %H:%M:%S")
message(STATUS "Build timestamp is ${_date}")

include(GetGitRevisionDescription)
git_describe(GIT_VERSION --tags)
message(STATUS "Configured version is ${GIT_VERSION}")
file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/GitVersion.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GitVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/GitVersion.h @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(vent_SOURCES
	src/state.cpp
	src/gasinfo.cpp
	src/slope.cpp
	src/contact.cpp
	src/vtk.cpp
	src/object.cpp
	src/pipe.cpp
	src/solver.cpp
	src/atm.cpp
	src/fan.cpp
	src/room.cpp
	src/pulse.cpp
	src/corrector.cpp
	)

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -Wextra -g -funroll-loops")
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -Wextra -g -funroll-loops")
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Werror -Wextra -g")
endif()

#find_package(OpenMP)
#if (${OPENMP_FOUND})
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#endif()

set(max_errors 5)
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=${max_errors}")
endif()
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ferror-limit=${max_errors}")
endif()

if(NOT DEFINED ENTROPY_FIX)
	set(ENTROPY_FIX FALSE)
endif()

if(NOT DEFINED PRECISE_RIEMANN)
	set(PRECISE_RIEMANN FALSE)
endif()

if(NOT DEFINED TURBULENCE)
	set(TURBULENCE FALSE)
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
message(STATUS "Use entropy fix (-DENTROPY_FIX): ${ENTROPY_FIX}")
message(STATUS "Enable turbulence model (-DTURBULENCE): ${TURBULENCE}")
message(STATUS "Use floating point traps (-DFP_TRAP): ${FP_TRAP}")
message(STATUS "Number of components (-DNC): ${NC}")
message(STATUS " ")

if (${ENTROPY_FIX})
	add_definitions("-DENTROPY_FIX=1")
else()
	add_definitions("-DENTROPY_FIX=0")
endif()

if (${PRECISE_RIEMANN})
	add_definitions("-DPRECISE_RIEMANN=1")
	list(APPEND vent_SOURCES src/predictor_godunov.cpp)
else()
	add_definitions("-DPRECISE_RIEMANN=0")
	list(APPEND vent_SOURCES src/predictor_cir.cpp)
endif()

if (${TURBULENCE})
	add_definitions("-DTURBULENCE=1")
else()
	add_definitions("-DTURBULENCE=0")
endif()

if (${FP_TRAP})
	add_definitions("-DFP_TRAP=1")
else()
	add_definitions("-DFP_TRAP=0")
endif()

add_definitions("-DNC=${NC}")

find_package(SWIG 3.0 REQUIRED)
include(UseSWIG)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Python3 COMPONENTS Development)

include_directories(include/)
add_library(vent SHARED ${vent_SOURCES})

set_source_files_properties(pyvent.i PROPERTIES CPLUSPLUS ON)
set_source_files_properties(pyvent.i PROPERTIES SWIG_FLAGS "-Wall")
swig_add_library(pyvent LANGUAGE python SOURCES pyvent.i)
swig_link_libraries(pyvent Python3::Python vent)

#add_executable(main main.cpp)
#target_link_libraries(main vent)

# add_executable(col col.cpp)
# target_link_libraries(col vent)
