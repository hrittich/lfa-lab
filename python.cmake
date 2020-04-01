# LFA Lab - Library to simplify local Fourier analysis.
# Copyright (C) 2018-2020  Hannah Rittich
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
  set(PYTHON_SOURCES "")

  foreach(_PYTHON_FILE ${ARGN})
    get_filename_component(_SOURCE_NAME ${_PYTHON_FILE} NAME)

    if(IS_ABSOLUTE ${_PYTHON_FILE})
      set(_ORIGINAL_SOURCE_FILE "${_PYTHON_FILE}")
    else()
      set(_ORIGINAL_SOURCE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_SOURCE_NAME}")
    endif()
    set(_SOURCE_FILE "${CMAKE_CURRENT_BINARY_DIR}/${_SOURCE_NAME}")
    set(_INSTALL_FILE "${PYTHON_INSTALL_DIR}/${MODULE_PATH}/${_SOURCE_NAME}")

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

    list(APPEND PYTHON_SOURCES "${_SOURCE_FILE}")

    install(FILES "${_SOURCE_FILE}"
            DESTINATION ${PYTHON_INSTALL_DIR}/${MODULE_PATH})

    # We create the Python cache in the installation directory. This is
    # easier, because we do not have to figure out the names of the generated
    # files.
    # The problem is, that on Python 2 uninstalling does not remove the
    # cache. Since Python 2, however, is not supported anymore. I will not
    # fix this. Python 2 support will be dropped in the future.
    install(CODE "message(\"Compiling ${_INSTALL_FILE}\")")
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m py_compile ${_INSTALL_FILE})")
  endforeach()

  add_custom_target(${TARGET} ALL DEPENDS ${PYTHON_SOURCES})
endfunction(ADD_PYTHON)


