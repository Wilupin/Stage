cmake_minimum_required(VERSION 2.8.9)
set(CMAKE_BUILD_TYPE Release)
project(BMW)
set(SOURCE_DIR src)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

file(GLOB_RECURSE SOURCE_FILES
  ${SOURCE_DIR}/*
  )

add_executable(
  exe
  ${SOURCE_FILES}
  )
