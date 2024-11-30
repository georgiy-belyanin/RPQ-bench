cmake_minimum_required(VERSION 3.17)
project(rpq-bench)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vendor/LAGraph)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/vendor/GraphBLAS)

enable_language(C)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_executable(rpq-bench
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.c)

set_target_properties ( PROPERTIES
    C_STANDARD_REQUIRED ON
    C_STANDARD 11 )

target_link_libraries(rpq-bench GraphBLAS LAGraph LAGraphX)
target_include_directories(rpq-bench PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

