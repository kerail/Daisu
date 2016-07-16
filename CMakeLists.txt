cmake_minimum_required(VERSION 3.3)
project(CameraCalibration)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

add_definitions(-DBOOST_TEST_DYN_LINK)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(.)

add_subdirectory(geom)