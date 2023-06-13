cmake_minimum_required(VERSION 3.1)
# get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(PROJECT_NAME "Viewer")
project(${PROJECT_NAME})
set(CMAKE_BUILD_TYPE release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
polyscope
GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
GIT_TAG master
GIT_SHALLOW TRUE
GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(polyscope)

FetchContent_Declare(
libigl
GIT_REPOSITORY https://github.com/libigl/libigl.git
GIT_TAG main
GIT_SHALLOW TRUE
GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(libigl)

set(LIBIGL_GLFW   OFF CACHE BOOL "Build libigl with GLFW")
set(LIBIGL_IMGUI  OFF CACHE BOOL "Build libigl with ImGui")
set(LIBIGL_OPENGL OFF CACHE BOOL "Build libigl with OpenGL")
# add_subdirectory("./libigl")
# add_subdirectory("./polyscope")

find_package(OpenMP REQUIRED)

# Add your project files
include_directories("." "./PoissonCurveRecon/" "./Misha/" "./Include/")
file(GLOB SOURCEFILES  ./Include/*.cpp) #        ./ExteriorPoissonRecon/*.cpp)
file(GLOB HEADERS  ./Misha/*.h ./Include/*.h ./ExteriorPoissonRecon/*.h)
file(GLOB INLINE_HEADERS  ./Misha/*.inl ./Include/*.inl)

add_executable(${PROJECT_NAME} ./Viewer/Viewer.cpp ${SOURCEFILES} ${HEADERS} ${INLINE_HEADERS})
target_link_libraries(${PROJECT_NAME} igl::core polyscope OpenMP::OpenMP_CXX)
