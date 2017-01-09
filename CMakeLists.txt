cmake_minimum_required(VERSION 2.8)
enable_testing()
project(mjolnir)

include_directories(${PROJECT_SOURCE_DIR})

option(ENABLE_PARALLEL_COMPUTING "enable multi thread computing" ON)
if(ENABLE_PARALLEL_COMPUTING)
    find_package(Threads REQUIRED)
    add_definitions(-DMJOLNIR_PARALLEL_THREAD)
endif(ENABLE_PARALLEL_COMPUTING)

add_subdirectory(src)
add_subdirectory(test)
