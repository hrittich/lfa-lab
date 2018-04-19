# LFA Lab - Library to simplify local Fourier analysis.
# Copyright (C) 2018  Hannah Rittich
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

cmake_minimum_required(VERSION 3.0)
project("LFA Lab"
  VERSION 0.4.0
  LANGUAGES CXX)
# project status is one of: alpha, beta, rc, stable
set(PROJECT_STATUS "alpha")
# project(LFA CXX C Fortran)

if(NOT PROJECT_STATUS STREQUAL "stable")
  set(PROJECT_VERSION "${PROJECT_VERSION}-${PROJECT_STATUS}")
endif()

# Do not include RPATH into executable. Otherwise build fails on MacOS
set(CMAKE_SKIP_RPATH ON)


# make a Release Build by default
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  # enable warnings
  set(CMAKE_CXX_FLAGS "-Wall ${CMAKE_CXX_FLAGS}")
endif()


# ===== Check C++ 11 Features ======
include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES("
#include <memory>
std::shared_ptr<int> smart_ptr;
int main() { return 0; }
" HAVE_STD_SHARED_PTR)
if(HAVE_STD_SHARED_PTR)
    set(HAVE_STD_SHARED_PTR TRUE)
endif()
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES("
#include <cstddef>
void* p = nullptr;
int main() { return 0; }
" HAVE_NULLPTR)


# ===== Check System Features =====

include(CheckSymbolExists)
check_symbol_exists(feenableexcept "fenv.h" HAVE_FEENABLEEXCEPT)

# ===== PkgConfig =====
find_package(PkgConfig)

# ===== Doxygen =====
find_package(Doxygen)

# ===== Eigen3 =====
pkg_check_modules(EIGEN REQUIRED eigen3)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})

# ====== LAPACK ======
find_package(LAPACK)
# ${LAPACK_FOUND}
option(WITH_LAPACK "Use LAPACK Library" ${LAPACK_FOUND})
if(WITH_LAPACK)
    list(APPEND LIBS ${LAPACK_LIBRARIES})
endif()

# ====== ARPACK ======
find_library(ARPACK_LIBRARIES arpack)
if(ARPACK_LIBRARIES)
    set(ARPACK_FOUND ON)
else()
    set(ARPACK_FOUND OFF)
endif()

option(WITH_ARPACK "Use ARPACK Library" OFF)
if(WITH_ARPACK)
    list(APPEND LIBS ${ARPACK_LIBRARIES})
endif()

# ====== Boost ======

if(NOT HAVE_STD_SHARED_PTR)
  # If we do not have C++11 shared_ptr we need Boost
  message(WARNING "Could not find a C++11 compatible compiler. Trying to find boost.")

  set(Boost_USE_STATIC_LIBS        OFF)
  set(Boost_USE_MULTITHREADED      ON)
  find_package(Boost 1.40.0 REQUIRED)

  include_directories(${Boost_INCLUDE_DIR})
  list(APPEND LIBS ${Boost_LIBRARIES})
endif()


# ====== Google Test ======
find_package(Threads)

# try to find gtest sources
find_path(GTEST_INCLUDE_DIR NAMES "gtest/gtest.h")
find_path(GTEST_SRC_DIR NAMES "src/gtest-all.cc"
        PATH_SUFFIXES "src/gtest")
if(GTEST_INCLUDE_DIR AND GTEST_SRC_DIR)
    include_directories(${GTEST_INCLUDE_DIR} ${GTEST_SRC_DIR})
    set(GTEST_FOUND TRUE)

    # build gtest
    add_library(gtest
        ${GTEST_SRC_DIR}/src/gtest-all.cc
        ${GTEST_SRC_DIR}/src/gtest_main.cc)

    set(GTEST_BOTH_LIBRARIES gtest)
endif()

# try finding it by gtest module
if(NOT GTEST_FOUND)
    find_package(GTest)
    # defines the GTEST_FOUND variable
endif()

option(WITH_TESTS "Build unit tests" ${GTEST_FOUND})

# ====== OpenMP ======
find_package(OpenMP)
#option(WITH_OPENMP "Use OpenMP" ${OPENMP_FOUND})
option(WITH_OPENMP "Use OpenMP" OFF)
if(WITH_OPENMP)
  # OpenMP
  set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

  # GCC specific...
  #if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_LINK_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_LINK_FLAGS}")
  #endif()

endif()

# ====== SWIG Python Wrapper ======
find_program(SWIG_EXECUTABLE swig)
if(NOT SWIG_EXECUTABLE)
    message(SEND_ERROR "SWIG not found.")
endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ====== Python + NumPy ======
find_package(PythonInterp)
include(python.cmake)

# use python to generate a configuration file
execute_process(COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/python_config.py
    RESULT_VARIABLE PYTHON_CONFIG
    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/python_config.cmake)
# include the configuration
include(${CMAKE_CURRENT_BINARY_DIR}/python_config.cmake)

if(NOT PYTHON_CONFIG EQUAL 0)
    message(SEND_ERROR "Python configuration failed.")
endif()

option(USER_INSTALL "Install into per user directory")
if(NOT USER_INSTALL STREQUAL USER_INSTALL_OLD)
  message("USER INSTALL SETTING CHANGED")
  if(USER_INSTALL)
    set(NEW_PYTHON_INSTALL_DIR ${PYTHON_USER_INSTALL_DIR})
  else()
    set(NEW_PYTHON_INSTALL_DIR ${PYTHON_GLOBAL_INSTALL_DIR})
  endif()
  set(PYTHON_INSTALL_DIR "${NEW_PYTHON_INSTALL_DIR}" CACHE PATH
    "Installation directory of python packages" FORCE)
  set(USER_INSTALL_OLD ${USER_INSTALL} CACHE INTERNAL "Internal")
endif()

include_directories(${NUMPY_INCLUDES})
include_directories(${PYTHON_INCLUDE_DIRS})

if(NOT NUMPY_FOUND)
    message(SEND_ERROR "NumPy not found")
endif()

# Configure GCC
option(BOUND_CHECKS "Enable bound checks" OFF)
if(BOUND_CHECKS)
    if(CMAKE_COMPILER_IS_GNUCXX)
        # enable STL Debug mode
        set(GCC_BOUND_CHECKS ON)
    else()
        message(FATAL_ERROR "Bound checks not supported")
    endif()
endif()


# Display Summary
message("
----------------------------------------
 Feature    Status
----------------------------------------
 Version    ${PROJECT_NAME} ${PROJECT_VERSION}
 Build      ${CMAKE_BUILD_TYPE}
 Eigen      ${EIGEN_VERSION}
 LAPACK     ${WITH_LAPACK}
 ARPACK     ${WITH_ARPACK}
 OpenMP     ${WITH_OPENMP}
 Python     ${PYTHON_VERSION}
 Tests      ${WITH_TESTS}
----------------------------------------

Install to ${PYTHON_INSTALL_DIR}
")


# Write configuration
configure_file(Config.h.in Config.h)
include_directories(BEFORE ${CMAKE_BINARY_DIR})

# Environment source file
configure_file("setup-env.sh.in" "setup-env.sh")

# Tags
add_custom_target(tags COMMAND ctags *.cpp *.h)

# Create distribution directory.
set(DIST_DIR "lfa-lab-${PROJECT_VERSION}")
add_custom_target(dist
  COMMAND git clone ${CMAKE_SOURCE_DIR} ${DIST_DIR}
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${DIST_DIR}/.git
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_BINARY_DIR}/doc/html ${DIST_DIR}/doc/html
  COMMAND ${CMAKE_COMMAND} -E remove_directory
    ${DIST_DIR}/doc/html/publications
  COMMAND ${CMAKE_COMMAND} -E remove_directory
    ${DIST_DIR}/doc/publications
  COMMAND tar -zcvf "${DIST_DIR}.tar.gz" ${DIST_DIR}
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${DIST_DIR})
add_dependencies(dist sphinx-doc)


# Target to run all checks
add_custom_target(check)

# Subdirectories
add_subdirectory(lfa_lab)
add_subdirectory(doc)
add_subdirectory(test)
add_subdirectory(demo)

