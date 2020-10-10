cmake_minimum_required(VERSION 3.0)

project(SoftPosit)

IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR CMAKE_VERSION VERSION_GREATER "3.0.0")
  CMAKE_POLICY(SET CMP0045 OLD)  # suppress set_target_properties empty target warnning
  CMAKE_POLICY(SET CMP0026 OLD)  # suppress set_target_properties read LOCATION warnning
ENDIF()

# use libc++ for osx
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=c++1y")

if (APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(ARCHIVE_OUTPUT_PATH ${PROJECT_BINARY_DIR})


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
if (APPLE)
  INCLUDE_DIRECTORIES (/System/Library/Frameworks)
  find_library(Accelerate_LIBRARIES Accelerate)
  set(ALL_LIBS ${Accelerate_LIBRARIES})
endif()


# find boost
find_package(Boost)
link_directories(${Boost_LIBRARY_DIRS})

# armadillo
find_package(Armadillo REQUIRED)
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
  ${Boost_INCLUDE_DIRS}
  ${ARMADILLO_INCLUDE_DIRS}
  "/usr/include/eigen3/"
  )

add_definitions("-D_USE_MATH_DEFINES")

set(ALL_LIBS ${ALL_LIBS} ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES})


add_subdirectory(src)
