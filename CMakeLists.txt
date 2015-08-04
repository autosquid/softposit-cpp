cmake_minimum_required(VERSION 3.0)

project(SoftPosit)

IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR CMAKE_VERSION VERSION_GREATER "3.0.0")
  CMAKE_POLICY(SET CMP0045 OLD)  # suppress set_target_properties empty target warnning
  CMAKE_POLICY(SET CMP0026 OLD)  # suppress set_target_properties read LOCATION warnning
ENDIF()

# set (CMAKE_VERBOSE_MAKEFILE 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=c++1y")

# use libc++ for osx
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(BUILD_STATIC_LIBS ON)
# default out of source build
if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( FATAL_ERROR "Please select another Build Directory ! (and give it a clever name, like bin_Visual2012_64bits/)" )
endif()

if( CMAKE_SOURCE_DIR MATCHES " " )
  message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()

if( CMAKE_BINARY_DIR MATCHES " " )
  message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()


# required thirdy party libraries
INCLUDE_DIRECTORIES (/System/Library/Frameworks)
find_library(Accelerate_LIBRARIES Accelerate)
set(ALL_LIBS ${Accelerate_LIBRARIES})


# find boost
# set(Boost_USE_STATIC_LIBS ON)
# set(Boost_USE_MULTITHREADED ON)
# set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost REQUIRED COMPONENTS filesystem system unit_test_framework random date_time log_setup log regex thread program_options)


# armadillo
find_package(Armadillo REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
  ${Boost_INCLUDE_DIRS}
  ${ARMADILLO_INCLUDE_DIRS}
  )

link_directories(${Boost_LIBRARY_DIRS})

set(ALL_LIBS ${ALL_LIBS} ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES})

add_definitions("-D_USE_MATH_DEFINES")

# include(flycheck-clang)
# write_dir_locals_el(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(src)
