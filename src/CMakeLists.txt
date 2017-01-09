set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

add_definitions("-std=c++11 -O3")
add_executable(mjolnir main.cpp)
if(ENABLE_PARALLEL_COMPUTING)
    target_link_libraries(mjolnir ${CMAKE_THREAD_LIBS_INIT})
endif(ENABLE_PARALLEL_COMPUTING)


