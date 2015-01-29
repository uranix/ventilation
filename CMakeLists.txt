cmake_minimum_required(VERSION 3.0)

project(Ventilation)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(GetGitRevisionDescription)
git_describe(GIT_VERSION --tags)
message(STATUS "Configured version is ${GIT_VERSION}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GitVersion.h.in ${CMAKE_CURRENT_BINARY_DIR}/GitVersion.h @ONLY)

set(vent_SOURCES src/flux.cpp src/contact.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=5 -std=c++11 -Wall -Werror -Wextra -g")

if(NOT DEFINED SECOND_ORDER)
	set(SECOND_ORDER FALSE)
endif()

if(NOT DEFINED RECONSTRUCT)
	set(RECONSTRUCT "${SECOND_ORDER}")
endif()

if(NOT DEFINED PRECISE_RIEMANN)
	set(PRECISE_RIEMANN FALSE)
endif()

if(NOT DEFINED FP_TRAP)
	set(FP_TRAP TRUE)
endif()

if(NOT DEFINED LIMITERS)
	set(LIMITERS "${SECOND_ORDER}")
endif()

message(STATUS "Configuerd with following params")
message(STATUS "Use output reconstruction (-DRECONSTRUCT): ${RECONSTRUCT}")
message(STATUS "Use second order method (-DSECOND_ORDER): ${SECOND_ORDER}")
message(STATUS "Use slope limiters (-DLIMITERS): ${LIMITERS}")
message(STATUS "Use precise Riemann solver (-DPRECISE_RIEMANN): ${PRECISE_RIEMANN}")
message(STATUS "Use floating point traps (-DFP_TRAP): ${FP_TRAP}")

if (${RECONSTRUCT})
	add_definitions("-DRECONSTRUCTED_OUTPUT=1")
else()
	add_definitions("-DRECONSTRUCTED_OUTPUT=0")
endif()

if (${SECOND_ORDER})
	add_definitions("-DSECOND_ORDER=1")
else()
	add_definitions("-DSECOND_ORDER=0")
endif()

if (${PRECISE_RIEMANN})
	add_definitions("-DPRECISE_RIEMANN=1")
else()
	add_definitions("-DPRECISE_RIEMANN=0")
endif()

if (${FP_TRAP})
	add_definitions("-DFP_TRAP=1")
else()
	add_definitions("-DFP_TRAP=0")
endif()

if (${LIMITERS})
	add_definitions("-DUSE_LIMITERS=1")
else()
	add_definitions("-DUSE_LIMITERS=0")
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

# add_executable(main main.cpp)
# target_link_libraries(main vent)

# add_executable(col col.cpp)
# target_link_libraries(col vent)
