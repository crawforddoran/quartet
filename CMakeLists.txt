cmake_minimum_required(VERSION 3.0)
project(quartet)

file(GLOB SOURCES src/*.cpp)

#message("${SOURCES}")

list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/gluvi.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/view_tet.cpp")

message("${SOURCES}")

add_executable(quartet ${SOURCES})
