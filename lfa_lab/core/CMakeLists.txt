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


add_library(lfa STATIC
  ../../Config.h
  Common.cpp Common.h
  MathUtil.h MathUtil.cpp
  EigenSolver.cpp EigenSolver.h
  DenseStencil.cpp DenseStencil.h
  Stencil2d.cpp Stencil2d.h
  Stencil3d.h
  SystemStencil.cpp SystemStencil.h
  StencilGallery.cpp StencilGallery.h
  ExEigenSolver.cpp ExEigenSolver.h
  VectorizedIndex.h
  MultiArray.h
  BlockStencil.cpp BlockStencil.h
  NMatrixView.h
  BlockStencil.cpp BlockStencil.h
  FoContext.cpp FoContext.h
  FoStencil.cpp FoStencil.h
  SystemSymbol.cpp SystemSymbol.h
  SystemSymbolProperties.cpp SystemSymbolProperties.h
  SystemClusterSymbol.cpp SystemClusterSymbol.h
  FrequencyIndices.cpp FrequencyIndices.h
  HarmonicClusters.cpp HarmonicClusters.h
  HarmonicIndices.cpp HarmonicIndices.h
  NdRange.h
  BdMatrix.cpp BdMatrix.h
  ClusterSymbol.cpp ClusterSymbol.h
  Symbol.cpp Symbol.h
  Grid.cpp Grid.h
  SplitFrequencyDomain.cpp SplitFrequencyDomain.h
  CartIterator.h
  SamplingProperties.cpp SamplingProperties.h
  FoProperties.cpp FoProperties.h
  DiscreteDomain.cpp DiscreteDomain.h
  SymbolBuilder.cpp SymbolBuilder.h
  BlockSb.cpp BlockSb.h
  NdArray.cpp NdArray.h
  SparseStencil.cpp SparseStencil.h
  ConstantSb.cpp ConstantSb.h
  HpFilterSb.cpp HpFilterSb.h
)
set_property(TARGET lfa PROPERTY POSITION_INDEPENDENT_CODE ON)

if(WITH_TESTS)
    add_executable(lfa_test
        test_Stencil.cpp
        test_General.cpp
        test_Harmonics.cpp
        test_BdMatrix.cpp
        test_Symbol.cpp
        test_CartIterator.cpp
        test_SystemSymbol.cpp
        test_NdArray.cpp
        test_SparseStencil.cpp
        test_FoProperties.cpp)
    target_link_libraries(lfa_test lfa ${LIBS} ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif()

if(PYTHON_VERSION_MAJOR EQUAL 3)
  set(SWIG_PYTHON_VERSION_SWITCH "-py3")
else()
  set(SWIG_PYTHON_VERSION_SWITCH "")
endif()

# Create wrapper using Swig
file(GLOB SWIG_INPUT_FILES *.i)
add_custom_command(
    OUTPUT
        ${CMAKE_CURRENT_BINARY_DIR}/extension.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/extension.py
    COMMAND
        ${SWIG_EXECUTABLE} ARGS -c++ -python ${SWIG_PYTHON_VERSION_SWITCH}
        -o ${CMAKE_CURRENT_BINARY_DIR}/extension.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/lfa.i"
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/lfa.i
    DEPENDS ${SWIG_INPUT_FILES})

add_library(extension SHARED extension.cpp)
target_link_libraries(extension lfa ${LIBS})
# set required output name
set_target_properties(extension
  PROPERTIES LIBRARY_OUTPUT_NAME "_extension"
  PREFIX "")

if(APPLE)
  set_target_properties(extension
    PROPERTIES LINK_FLAGS "-undefined dynamic_lookup"
               SUFFIX ".so")
endif()

install(TARGETS extension DESTINATION ${PYTHON_INSTALL_DIR}/lfa_lab/core)
add_python(lfa_lab_core_py lfa_lab/core
  ${CMAKE_CURRENT_BINARY_DIR}/extension.py
  __init__.py)

