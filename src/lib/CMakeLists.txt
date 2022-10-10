# XMls/src/lib CMakeLists
# adds all executables to your project
# you should only modify this file when you add a new library
# you can also put here cmake code common to all libraries
# author: Bruno Vallet
# date: 11/2016

set(LIB_TYPE "STATIC")
option(BUILD_SHARED_LIBS "Build ${PROJECT_NAME} libraries shared." ON)
if(BUILD_SHARED_LIBS)
	set(LIB_TYPE "SHARED")
endif(BUILD_SHARED_LIBS)

add_subdirectory(libXBase)
add_subdirectory(libXMls)
