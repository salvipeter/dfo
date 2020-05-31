cmake_minimum_required(VERSION 3.1.0)

project(dfo)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(CMAKE_COMPILER_IS_GNUCXX)

add_library(dfo STATIC
  direct.hh
  direct.cc
  mads.hh
  mads.cc
  nelder-mead.hh
  nelder-mead.cc
  powell.hh
  powell.cc
  cross-entropy.hh
  cross-entropy.cc
  )

add_executable(dfo-test
  test.cc
  )

target_link_libraries(dfo-test dfo)

if(WIN32)
  set(CMAKE_CXX_FLAGS "/std:c++17")
endif()
