cmake_minimum_required(VERSION 3.0)

project(SoftPosit)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)

macro(get_gcc_compile_flags target flags)
  string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" name)
  set(flags "${${name}} ${CMAKE_CXX_COMPILER_ARG1}")
  get_target_property(value ${target} COMPILE_FLAGS)
  if (value)
    list(APPEND flags ${value})
  endif()
  get_target_property(value ${target} TYPE)
  get_target_property(value ${target} COMPILE_DEFINITIONS)
  if (value)
    foreach(item ${value})
      list(APPEND flags "-D${item}")
    endforeach()
  endif()
  STRING(TOUPPER "COMPILE_DEFINITIONS_${CMAKE_BUILD_TYPE}" name)
  get_target_property(value ${target} ${name})

  if (value)
    foreach(item ${value})
      list(APPEND flags "-D${item}")
    endforeach()
  endif()
  get_directory_property(value DEFINITIONS)
  if (value)
    set( value_list ${value} )
    separate_arguments(value_list)
    foreach (item ${value_list})
      list(APPEND flags ${item})
    endforeach()
  endif()
  get_directory_property(value INCLUDE_DIRECTORIES)
  if (value)
    foreach(item ${value})
      list(APPEND flags "-I${item}")
    endforeach()
  endif()
  #separate_arguments(flags)
endmacro()


IF(CMAKE_VERSION VERSION_EQUAL "3.0.0" OR CMAKE_VERSION VERSION_GREATER "3.0.0")
  CMAKE_POLICY(SET CMP0045 OLD)  # suppress set_target_properties empty target warnning
  CMAKE_POLICY(SET CMP0026 OLD)  # suppress set_target_properties read LOCATION warnning
ENDIF()

# set (CMAKE_VERBOSE_MAKEFILE 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -std=c++1y")

# use libc++ for osx
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")


set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)

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

# install in install directory by default
SET(MY_OWN_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install CACHE PATH "Prefix prepended to
install directories")
SET(CMAKE_INSTALL_PREFIX "${MY_OWN_INSTALL_PREFIX}" CACHE INTERNAL "Prefix
prepended to install directories" FORCE)

# set output place for executable and lib and static lib
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(ARCHIVE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/static)
# MyWorkingDirectory is working directory
set(MyWorkingDirectory ${PROJECT_BINARY_DIR})

# required thirdy party libraries

INCLUDE_DIRECTORIES (/System/Library/Frameworks)
find_library(Accelerate_LIBRARIES Accelerate)
set(ALL_LIBS ${Accelerate_LIBRARIES})
set(OpenCV_DIR "/usr/local/opt/opencv3/share/OpenCV" CACHE PATH "opencv path" FORCE)

set(CML_DIR "${3rd_root}/cml")
# set(LAPACK_DIR "/usr/local/opt/lapack")

#add qt support
SET(QT_SEARCH_PATH "" CACHE PATH "Prefix prepended to install directories")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${QT_SEARCH_PATH})
# find boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(TBB REQUIRED)
# find_package(lapack REQUIRED)

find_package(VTK REQUIRED)

find_package(OpenCL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OpenCV REQUIRED)

find_package(Armadillo REQUIRED)
find_package(OpenMesh REQUIRED)

find_package(Boost REQUIRED COMPONENTS filesystem system unit_test_framework random date_time log_setup log regex thread program_options)

set(3rd_bin_dir "${_3rd_dir}/bin")
list(APPEND 3rd_bin_dir "${3rd_root}/common/bin")
list(APPEND 3rd_bin_dir "${Boost_LIBRARY_DIRS}")
list(APPEND 3rd_bin_dir "${OpenCV_LIB_DIR}/../bin")
list(APPEND 3rd_bin_dir "${CMAKE_PREFIX_PATH}/bin")

include(${VTK_USE_FILE})

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
  ${Boost_INCLUDE_DIRS}
  ${NLOPT_INCLUDE_DIRS}
  # ${LAPACK_INCLUDE_DIRS}
  ${ARMADILLO_INCLUDE_DIRS}
  ${OPENCL_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${OPENMESH_INCLUDE_DIR}
  ${CML_DIR}
  ${ICP_INCLUDE_DIRS}
  ${LEVMAR_INCLUDE_DIRS}
  ${OPENGL_INCLUDE_DIR}
  ${TBB_INCLUDE_DIRS}
  )


link_directories(${Boost_LIBRARY_DIRS})
get_filename_component(tbb_lib_dir "${TBB_INCLUDE_DIRS}/../lib" ABSOLUTE)
link_directories(${tbb_lib_dir})

# armadillo
# add_definitions("-DARMA_DONT_USE_LAPACK")
# add_definitions("-DARMA_DONT_USE_BLAS")
add_definitions("-DWITH_OPENGL=ON")
set(ALL_LIBS
  ${ALL_LIBS}
  ${Boost_LIBRARIES}
  ${ARMADILLO_LIBRARIES}
  ${OPENCL_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${OpenCV_LIBS}
  ${ICP_LIBRARIES}
  ${TBB_LIBRARIES}
  ${VTK_LIBRARIES}
  )

add_definitions("-D_USE_MATH_DEFINES")

set (proj_include)
# include(flycheck-clang)
# write_dir_locals_el(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(sposit demo.cpp softposit.cpp)
target_link_libraries(sposit ${ALL_LIBS})


get_gcc_compile_flags(sposit flags)

if (proj_include)
  foreach(item ${proj_include})
    list(APPEND flags "-I${item}")
  endforeach()
endif()

list(REMOVE_DUPLICATES flags)
list(REMOVE_ITEM flags "-I${CMAKE_CURRENT_SOURCE_DIR}" "-I${CMAKE_CURRENT_BINARY_DIR}")

file(WRITE "${CMAKE_CURRENT_SOURCE_DIR}/.clang_complete" "")
foreach (aflag ${flags})
  if (aflag)
    file(APPEND "${CMAKE_CURRENT_SOURCE_DIR}/.clang_complete" ${aflag} "\n")
  endif()
endforeach(aflag)

############ test
