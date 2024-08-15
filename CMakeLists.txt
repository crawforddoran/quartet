cmake_minimum_required(VERSION 3.8)

project(quartet LANGUAGES CXX)

option(QUARTET_ENABLE_VIEWER "Build the mesh viewer" ON)

add_library(quartetcore
    src/predicates.cpp
    src/geometry_queries.cpp
    src/sdf.cpp
    src/trimesh.cpp
    src/tet_mesh.cpp
    src/feature.cpp
    src/read_obj.cpp
    src/tet_quality.cpp
    src/match_features.cpp
    src/optimize_tet_mesh.cpp
    src/make_signed_distance.cpp
    src/make_tet_mesh.cpp
)

add_executable(quartet src/main.cpp)
target_link_libraries(quartet quartetcore)

if (QUARTET_ENABLE_VIEWER)

    set(Open_GL_PREFERENCE "GLVND")
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)

    add_library(quartetview src/gluvi.cpp)
    target_link_libraries(quartetview PUBLIC OpenGL::OpenGL OpenGL::GLU GLUT::GLUT)

    add_executable(view-tet src/view_tet.cpp)
    target_link_libraries(view-tet quartetcore quartetview)
endif()
