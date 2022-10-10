# XMls/src/exe/Kinwa CMakeLists
# adds Kinwa executable to your project
# author: Bruno Vallet
# date: 02/2019

add_executable(KinWa2d KinWa2d.cpp PlyExport.cpp PriorityVertex.cpp PriorityEdge.cpp)
target_link_libraries(KinWa2d ${CGAL_LIBRARY} XMls)
