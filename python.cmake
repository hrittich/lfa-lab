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


function(ADD_PYTHON TARGET MODULE_PATH)

  foreach(_PYTHON_FILE ${ARGN})
    get_filename_component(_SOURCE_NAME ${_PYTHON_FILE} NAME)
    get_filename_component(_SOURCE_STRIPPED ${_SOURCE_NAME} NAME_WE)

    if(IS_ABSOLUTE ${_PYTHON_FILE})
      set(_ORIGINAL_SOURCE_FILE "${_PYTHON_FILE}")
    else()
      set(_ORIGINAL_SOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_SOURCE_NAME}")
    endif()
    set(_SOURCE_FILE "${CMAKE_CURRENT_BINARY_DIR}/${_SOURCE_NAME}")
    set(_TARGET_FILE "${CMAKE_CURRENT_BINARY_DIR}/${_SOURCE_STRIPPED}.pyc")

    # Copy file to binary dir
    # configure_file("${_ORIGINAL_SOURCE_FILE}" ${_SOURCE_FILE} COPYONLY)

    if(NOT (_SOURCE_FILE STREQUAL _ORIGINAL_SOURCE_FILE))
      add_custom_command(
        OUTPUT ${_SOURCE_FILE}
        COMMAND ${CMAKE_COMMAND} -E copy
                ${_ORIGINAL_SOURCE_FILE}
                ${_SOURCE_FILE}
        MAIN_DEPENDENCY ${_ORIGINAL_SOURCE_FILE})
    endif()

    list(APPEND PYTHON_OBJECTS "${_TARGET_FILE}")

    add_custom_command(
      OUTPUT ${_TARGET_FILE}
      COMMAND
        ${PYTHON_EXECUTABLE} ARGS -m py_compile ${_SOURCE_FILE}
      DEPENDS ${_SOURCE_FILE})

    install(FILES "${_SOURCE_FILE}" "${_TARGET_FILE}"
            DESTINATION ${PYTHON_INSTALL_DIR}/${MODULE_PATH})
  endforeach()

  add_custom_target(${TARGET} ALL DEPENDS ${PYTHON_OBJECTS})
endfunction(ADD_PYTHON)


