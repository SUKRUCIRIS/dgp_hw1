cmake_minimum_required(VERSION 3.12)
project(sukru_hw)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}$<0:>)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(Easy3D_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(Easy3D_BUILD_TUTORIALS OFF CACHE BOOL "" FORCE)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/third_party/Easy3D)
set(Easy3D_RESOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/Easy3D/resources)
add_subdirectory(third_party/Easy3D)
add_executable(sukru_hw ./src/main.cpp ./src/hw1.cpp ./src/hw2.cpp)
target_compile_definitions(sukru_hw PRIVATE 
    Easy3D_RESOURCE_DIR="${Easy3D_RESOURCE_DIR}"
	EASY3D_STATIC_DEFINE
)
set(Easy3D_DIR "${CMAKE_BINARY_DIR}/third_party/Easy3D")
target_link_libraries(sukru_hw
	easy3d::viewer
    easy3d::gui 
    easy3d::renderer 
    easy3d::algo 
    easy3d::core
	3rd_imgui
)
