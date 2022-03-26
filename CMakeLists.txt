cmake_minimum_required(VERSION 3.12)
project(mainIsoElpsPackTestEntropy VERSION 1.0)

add_subdirectory(src)

add_executable("${PROJECT_NAME}" "${PROJECT_NAME}.cpp")

target_link_libraries("${PROJECT_NAME}" PRIVATE xDEM)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries("${PROJECT_NAME}" PUBLIC OpenMP::OpenMP_CXX)
endif()

# find_package(OpenMP REQUIRED)
# if(OpenMP_FOUND)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif()

# cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release