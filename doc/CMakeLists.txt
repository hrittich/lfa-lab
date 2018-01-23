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

# Target to build all documentation.
add_custom_target(doc)

configure_file("conf.py.in" "conf.py")
add_custom_target(sphinx-doc
  ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}
  ${PYTHON_EXECUTABLE} -m sphinx
  -c ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/html)
add_dependencies(sphinx-doc lfa_lab_py extension)
add_dependencies(doc sphinx-doc)

configure_file("Doxyfile.in" "Doxyfile")
if(DOXYGEN_FOUND)
  add_custom_target(doxygen-doc
    #COMMAND
    #  ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/html"
    COMMAND
      ${DOXYGEN_EXECUTABLE} Doxyfile)

  add_dependencies(doc doxygen-doc)
endif()


