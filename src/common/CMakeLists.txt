# Copyright (C) 2025 Andrew D Smith
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

file(GLOB cpp_files "*.cpp")

set(LIBRARY_OBJECTS "")
foreach(cpp_file ${cpp_files})
  get_filename_component(BASE_NAME ${cpp_file} NAME_WE)
  add_library(${BASE_NAME} OBJECT ${cpp_file})
  target_link_libraries(${BASE_NAME} PRIVATE
    bamxx
    smithlab_cpp
  )
  target_include_directories(${BASE_NAME} PUBLIC
    ${CMAKE_BINARY_DIR}
  )
  list(APPEND LIBRARY_OBJECTS ${BASE_NAME})
endforeach()

# Create static library linking the individual objects
add_library(dnmtools_objs STATIC)
target_include_directories(dnmtools_objs PUBLIC
  ${CMAKE_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(dnmtools_objs PUBLIC
  ${LIBRARY_OBJECTS}
  smithlab_cpp
  bamxx
)
