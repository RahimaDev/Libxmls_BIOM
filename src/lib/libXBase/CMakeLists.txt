# XMls/src/lib/XBase CMakeLists
# adds library XBase to the project
# you should only modify this file when you add new files or dependencies to XBase
# author: Bruno Vallet
# date: 11/2016

set(XBase_SRCS
    XArchiGeoref.cpp
    XMat3D.cpp
    XPt2D.cpp
    XPt3D.cpp
    XQuaternion.cpp
    XRect.cpp
    XStringTools.cpp
    XArchiXMLBaseTools.cpp
    XArchiXMLTools.cpp
    XArchiGeorefXML.cpp
    XPath.cpp
    XSystemInfo.cpp
    XAbsoluteTime.cpp
 )

file(GLOB XBase_HEADERS ${PROJECT_SOURCE_DIR}/include/libXBase/*.h)

if(WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32)

add_definitions(-DBOOST_NO_CXX11_SCOPED_ENUMS)
include_directories(${Boost_INCLUDE_DIRS})

add_library(XBase ${LIB_TYPE} ${XBase_SRCS} ${XBase_HEADERS})
target_link_libraries(XBase ${XMls_LIBRARIES})

# copy raf98.txt into the bin folder
# can't use $<TARGET_FILE_DIR:XGeodesy> as it is the lib folder
add_custom_command(TARGET XBase PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                        ${CMAKE_CURRENT_SOURCE_DIR}/raf98.txt ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIGURATION>)

if(CMAKE_VERBOSITY GREATER 0)
message(STATUS "Added target XBase (lib)")
endif(CMAKE_VERBOSITY GREATER 0)

install(TARGETS XBase
  EXPORT XMls-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)

install(FILES ${XBase_HEADERS} DESTINATION "${XMls_INCLUDE_DIR}/libXBase")
