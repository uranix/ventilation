cmake_minimum_required(VERSION 2.8)

project(Ventilation)

set(vent_SOURCES src/flux.cpp src/contact.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ferror-limit=5 -std=c++11 -Wall -Wno-deprecated-register")

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(include/)
add_library(vent SHARED ${vent_SOURCES})

SET_SOURCE_FILES_PROPERTIES(vent.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(vent.i PROPERTIES SWIG_FLAGS "-Wall")
SWIG_ADD_MODULE(vent python vent.i)
SWIG_LINK_LIBRARIES(vent ${PYTHON_LIBRARIES} vent)

add_executable(main main.cpp)
target_link_libraries(main vent)

add_executable(col col.cpp)
target_link_libraries(col vent)
