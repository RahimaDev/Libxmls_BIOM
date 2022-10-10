# XMls/src/lib/XMls CMakeLists
# adds library XMls to the project
# you should only modify this file when you add new files or dependencies to XMls
# author: Bruno Vallet
# date: 11/2016

set(XMls_SOURCES XAttribute.cpp XEchoPulseTables.cpp XMls.cpp XTrajecto.cpp XCam.cpp XMCam.cpp)

file(GLOB XMls_HEADERS ${PROJECT_SOURCE_DIR}/include/libXMls/*.h PARENT_SCOPE)
set(XMLS_HEADERS ${XMls_HEADERS} PARENT_SCOPE)

add_library(XMls ${LIB_TYPE} ${XMls_SOURCES} ${XMls_HEADERS})
target_link_libraries(XMls XBase ${PROJ4_LIBRARY})
	
if(CMAKE_VERBOSITY GREATER 0)
message(STATUS "Added target XMls (lib)")
endif(CMAKE_VERBOSITY GREATER 0)

##############################################################################
# Install
##############################################################################
install(TARGETS XMls
  EXPORT XMls-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)

install(FILES ${XMls_HEADERS} DESTINATION "${XMls_INCLUDE_DIR}/libXMls")

