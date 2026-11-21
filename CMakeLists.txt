cmake_minimum_required(VERSION 3.16)

project(ntmesh_nearest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -D wins, then the environment variable, then a sibling checkout.
function(nt_resolve_root var marker sibling hint)
  if(NOT ${var} AND DEFINED ENV{${var}} AND NOT "$ENV{${var}}" STREQUAL "")
    set(${var} "$ENV{${var}}" PARENT_SCOPE)
    return()
  endif()
  if(NOT ${var} AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${sibling}/${marker}")
    get_filename_component(found "${CMAKE_CURRENT_SOURCE_DIR}/${sibling}" ABSOLUTE)
    set(${var} "${found}" PARENT_SCOPE)
    return()
  endif()
  if(NOT ${var})
    message(FATAL_ERROR "${var} is not set. Pass -D${var}= or export ${var}. ${hint}")
  endif()
endfunction()

nt_resolve_root(NTV_ROOT "src/nvec3.h" "../ntv/ntv" "NTV_ROOT is the ntv repo that contains src/nvec3.h.")
nt_resolve_root(NTMESH_ROOT "src/ntmesh.h" "../ntmesh" "NTMESH_ROOT is the ntmesh repo that contains src/ntmesh.h.")

if(NOT IS_ABSOLUTE "${NTV_ROOT}")
  get_filename_component(NTV_ROOT "${NTV_ROOT}" ABSOLUTE)
endif()
if(NOT IS_ABSOLUTE "${NTMESH_ROOT}")
  get_filename_component(NTMESH_ROOT "${NTMESH_ROOT}" ABSOLUTE)
endif()
if(NOT EXISTS "${NTV_ROOT}/src/nvec3.h")
  message(FATAL_ERROR "NTV_ROOT does not contain src/nvec3.h: ${NTV_ROOT}")
endif()
if(NOT EXISTS "${NTMESH_ROOT}/src/ntmesh.h")
  message(FATAL_ERROR "NTMESH_ROOT does not contain src/ntmesh.h: ${NTMESH_ROOT}")
endif()

set(LIGHT_BUNNY "${CMAKE_CURRENT_SOURCE_DIR}/samplegeo/light_bunny.geo")
if(NOT EXISTS "${LIGHT_BUNNY}")
  message(FATAL_ERROR "missing samplegeo/light_bunny.geo")
endif()

find_package(Threads REQUIRED)

function(nt_near_target name)
  add_executable(${name} ${ARGN})
  target_include_directories(${name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/test
    ${NTMESH_ROOT}/src
    ${NTV_ROOT}/src
  )
  target_compile_definitions(${name} PRIVATE LIGHT_BUNNY="${LIGHT_BUNNY}")
  target_compile_options(${name} PRIVATE -Wall -Wextra -Wpedantic)
  target_link_libraries(${name} PRIVATE Threads::Threads)
endfunction()

nt_near_target(ntmesh_nearest_tests test/main.C test/test_near.C)
nt_near_target(nearpc_scene samples/nearpc/src/main.C)
target_compile_definitions(nearpc_scene PRIVATE OUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/samples/nearpc")

message(STATUS "NTV_ROOT: ${NTV_ROOT}")
message(STATUS "NTMESH_ROOT: ${NTMESH_ROOT}")
