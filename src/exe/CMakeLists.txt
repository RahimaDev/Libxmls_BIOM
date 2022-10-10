# XMls/src/exe CMakeLists
# adds all executables to your project
# you should only modify this file when you add a new executable
# you can also put here cmake code common to all executables
# author: Bruno Vallet
# date: 11/2016

add_executable(EptReader EptReader.cpp)
target_link_libraries(EptReader XMls)

add_executable(RieglFlexport RieglFlexport.cpp)
target_link_libraries(RieglFlexport XMls)

add_executable(XMlsInfo XMlsInfo.cpp)
target_link_libraries(XMlsInfo XMls)

add_executable(XMlsAttrib XMlsAttrib.cpp)
target_link_libraries(XMlsAttrib XMls)

add_executable(XTrajectoSample XTrajectoSample.cpp)
target_link_libraries(XTrajectoSample XMls)

add_executable(GenerateTimePivot GenerateTimePivot.cpp)
target_link_libraries(GenerateTimePivot XMls)

add_executable(SensorMesh SensorMesh.cpp)
target_link_libraries(SensorMesh XMls)

add_executable(SensorMesh2 SensorMesh2.cpp)
target_link_libraries(SensorMesh2 XMls)

add_executable(SensorMesh3 SensorMesh3.cpp)
target_link_libraries(SensorMesh3 XMls)

add_executable(SensorEdges SensorEdges.cpp)
target_link_libraries(SensorEdges XMls)

add_executable(ListIntersect ListIntersect.cpp)
target_link_libraries(ListIntersect XMls)

add_executable(SelfRegistration SelfRegistration.cpp)
target_link_libraries(SelfRegistration XMls)

if(USE_CGAL)
include_directories(${CGAL_INCLUDE_DIRS})
add_definitions(-DCGAL_USE_BASIC_VIEWER)
#add_executable(Decimate Decimate.cpp)
#target_link_libraries(Decimate ${CGAL_LIBRARY})

add_executable(Polygon Polygon.cpp)
#target_link_libraries(Polygon XMls Qt5::Widgets Qt5::Xml Qt5::Gui Qt5::OpenGL CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_Core)
target_link_libraries(Polygon XMls Qt5::Widgets Qt5::Xml Qt5::Gui Qt5::OpenGL CGAL::CGAL  CGAL::CGAL_Core)
target_link_libraries (Polygon  ${PCL_LIBRARIES} )
add_subdirectory(KinWa)
endif()
