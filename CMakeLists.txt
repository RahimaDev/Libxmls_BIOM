# XMls root CMakeLists
# Defines all variables global to the project
# author: Bruno Vallet
# date: 11/2016

cmake_minimum_required(VERSION 2.8)
PROJECT(XMls)

set (CMAKE_CXX_STANDARD 11)



# enable C++11
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
#set (CMAKE_CXX_STANDARD 14)

####################################
### DEFAULT BUILD CONFIGURATION  ###
####################################
# project variables to export in the ${PROJECT_NAME}Config.cmake
string(TOUPPER ${PROJECT_NAME} UPPER_PROJECT_NAME)
set(${UPPER_PROJECT_NAME}_ROOT_DIR "${PROJECT_SOURCE_DIR}")
set(${UPPER_PROJECT_NAME}_DEFINITIONS "")
set(${UPPER_PROJECT_NAME}_INCLUDE_DIRS "")
set(${UPPER_PROJECT_NAME}_LIBRARIES "")
set(CMAKE_VERBOSITY 1 CACHE INT "0: none, 1: normal, 2: verbose")

# build type
if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING
	" Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# global definitions
if (UNIX AND NOT APPLE)
	add_definitions(-D__LINUX__)
endif (UNIX AND NOT APPLE)
if (WIN32 AND MSVC)
	add_definitions(-D__WIN32__)
endif (WIN32 AND MSVC)

# output directories
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin CACHE PATH
	" Where to put all the RUNTIME targets when built" FORCE)
endif (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)

if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib CACHE PATH 
	" Output directory in which to build LIBRARY target files" FORCE)
endif (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)

if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
	set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib CACHE PATH
	" Where to put all the ARCHIVE targets when built" FORCE)
endif (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

if (NOT EXECUTABLE_OUTPUT_PATH)
	set(EXECUTABLE_OUTPUT_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} CACHE PATH
        " old variable not used with cmake 2.6+" FORCE)
endif (NOT EXECUTABLE_OUTPUT_PATH)

if (NOT LIBRARY_OUTPUT_PATH)
	set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} CACHE PATH
        " old variable not used with cmake 2.6+" FORCE)
endif (NOT LIBRARY_OUTPUT_PATH)

mark_as_advanced(
CMAKE_ARCHIVE_OUTPUT_DIRECTORY
EXECUTABLE_OUTPUT_PATH
LIBRARY_OUTPUT_PATH)

#include directory
include_directories(include)

######################################
### DEFAULT INSTALL CONFIGURATION  ###
######################################
# install path
if(UNIX)
	set(INSTALL_PREFIX "/usr/local" CACHE PATH " install path" )
endif(UNIX)
if(WIN32)
	set(INSTALL_PREFIX "C:/Program Files/MATIS" CACHE PATH " install path")
endif(WIN32)
set(CMAKE_INSTALL_PREFIX "${INSTALL_PREFIX}" CACHE INTERNAL " real install path" FORCE)

# documentation install default path
if(UNIX)
	set(MATIS_DOC_DIR "/usr/local/doc" CACHE PATH " install path" )
endif(UNIX)
if(WIN32)
	set(MATIS_DOC_DIR "C:/Program Files/MATIS/doc" CACHE PATH " install path")
endif(WIN32)
SET(MATIS_DOC_COMMON_DIR ${MATIS_DOC_DIR}/matis_doc_template )
SET(MATIS_DOC_TRAC_DIR ${MATIS_DOC_DIR}/trac )

# prefix path (mainly used to tell find_package where to look)
set(PREFIX_PATH "${INSTALL_PREFIX}" CACHE PATHS " where find_package looks")
set(CMAKE_PREFIX_PATH "${PREFIX_PATH}" CACHE INTERNAL " where find_package looks" FORCE)

#######################
### GENERAL OUTPUT  ###
#######################
if(CMAKE_VERBOSITY GREATER 0)
	message(STATUS "---------------- GENERAL --------------------")
	message(STATUS "- Project ${PROJECT_NAME} under ${CMAKE_SYSTEM} in ${CMAKE_BUILD_TYPE} mode")
	message(STATUS "- Project root: ${${UPPER_PROJECT_NAME}_ROOT_DIR}")
	message(STATUS "- Makefile/solutions are in: ${CMAKE_BINARY_DIR}")
	message(STATUS "- Executables are built in: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
	message(STATUS "- Libraries are built in: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
	message(STATUS "- Project is installed in: ${CMAKE_INSTALL_PREFIX}")
	if(CMAKE_VERBOSITY GREATER 1)
		string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)
		message(STATUS "- Default cxx flags are: ${CMAKE_CXX_FLAGS_${UPPER_CMAKE_BUILD_TYPE}}")
	endif(CMAKE_VERBOSITY GREATER 1)
	message(STATUS "---------------------------------------------")
	message(STATUS "Looking for packages...")
endif(CMAKE_VERBOSITY GREATER 0)

##############################
### Find required packages ###
##############################
# most commonly used packages, remove the ones you do not need
# the packages found here will be accessible project wide.
# If you need a package for an exe or lib only, put the find_package in the appropriate CMakeList

set(XMls_LIBRARIES "")
set(XMls_DEFINITIONS "")

find_package(PkgConfig REQUIRED)
pkg_check_modules(TINYXML tinyxml)
if(TINYXML_FOUND)
  message(STATUS "Tinyxml found")
  include_directories(${TINYXML_INCLUDE_DIRS})
  set(XMls_DEFINITIONS ${XMls_DEFINITIONS} ${TINYXML_DEFINITIONS})
  set(XMls_LIBRARIES ${XMls_LIBRARIES} ${TINYXML_LIBRARIES})
else(TINYXML_FOUND)
  message(ERROR "Tinyxml not found")
endif(TINYXML_FOUND)

option(USE_Ori "Use libOri" 0)
if(USE_Ori)
    find_package(Ori REQUIRED)
endif()

find_package( Boost REQUIRED COMPONENTS date_time system filesystem )
include_directories(Boost_INCLUDE_DIRS)
set(XMls_DEFINITIONS ${XMls_DEFINITIONS} ${Boost_DEFINITIONS})
set(XMls_LIBRARIES ${XMls_LIBRARIES} ${Boost_LIBRARIES})

option(USE_CGAL "Use CGAL" 1)
if(USE_CGAL)
#find_package(CGAL REQUIRED Components Qt5 Core)
    find_package(CGAL REQUIRED Qt5 Core)
    message(STATUS "Found CGAL: CGAL_INCLUDE_DIRS=" ${CGAL_INCLUDE_DIRS} " and CGAL_LIBRARY=" ${CGAL_LIBRARY} " and CGAL_USE_FILE=" ${CGAL_USE_FILE})
    include(${CGAL_USE_FILE})
endif()


find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt5Gui)
find_package(Qt5OpenGL)
find_package(Qt5Xml)

get_target_property(QtCore_location Qt5::Core LOCATION)
find_package(PCL 1.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})






find_path(PROJ4_INCLUDE_DIRS NAMES proj_api.h PATHS ${PROJ4_DIR}/include)
if(WIN32)
    find_library(PROJ4_LIB NAMES proj_4_8 proj4 proj PATHS ${PROJ4_DIR}/${MSVC_PREFIX}_lib)
    find_library(PROJ4_LIB_DBG NAMES proj_4_8 proj4 projd PATHS ${PROJ4_DIR}/${MSVC_PREFIX}_lib)
    if(PROJ4_LIB_DBG)
       set(PROJ4_LIBRARY optimized ${PROJ4_LIB} debug ${PROJ4_LIB_DBG})
    else(PROJ4_LIB_DBG)
        set(PROJ4_LIBRARY ${PROJ4_LIB})
    endif(PROJ4_LIB_DBG)
else(WIN32)
    find_library(PROJ4_LIBRARY NAMES proj4 proj)
    #set(PROJ4_LIBRARY  proj4)
endif(WIN32)
set(XMls_LIBRARIES ${XMls_LIBRARIES} ${PROJ4_LIBRARY})
set(XMls_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)

add_subdirectory(src)

install (EXPORT XMls-targets DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)
set(XMls_LIBRARIES XBase XMls ${XMls_LIBRARIES})
configure_file(XMlsConfig.cmake.in ${CMAKE_BINARY_DIR}/XMlsConfig.cmake.in2 @ONLY)
configure_file(${CMAKE_BINARY_DIR}/XMlsConfig.cmake.in2 ${CMAKE_BINARY_DIR}/XMlsConfig.cmake @ONLY)
install (FILES ${CMAKE_BINARY_DIR}/XMlsConfig.cmake DESTINATION ${CMAKE_INSTALL_PREFIX}/cmake)
